# plumber.R
library(plumber)
require(tidyr)
require(dplyr)
require(sf)
require(sp)
require(ggVennDiagram)
require(RColorBrewer)
require(ggplot2)
require(gbifdb)
require(lubridate)
require(jsonlite)

Sys.setenv("AWS_ACCESS_KEY_ID" = "e9190baae65b40a38bf43ade883b04a6",
           "AWS_SECRET_ACCESS_KEY" = "7aa129a7d84744efa76183cc9cf4b0a5")

#* Echo the parameter that was sent in
#* @param ProjectID:string Project ID
#* @param First_Date:string YYYY-MM-DD
#* @param Last_Date:string YYYY-MM-DD
#* @param Marker:string Target marker name
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @param Geographic_Scale:string Local, State, or Nation
#* @get /venn
venn <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,Geographic_Scale){
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
  Primer <- Marker
  #Read in parsed Tronko database.
  TronkoInput <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/",Primer,"/Taxa_Parsed.tsv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  TronkoInput <- read.table(text = TronkoInput,header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  
  #Coerce data type.
  TronkoDB <- as.data.frame(TronkoInput)
  
  #Filter by maximum number of read mismatches.
  TronkoDB <- TronkoDB[TronkoDB$Mismatch <= Num_Mismatch & !is.na(TronkoDB$Mismatch),]
  
  #Read in Metadata
  Metadata <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/Metadata.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  Metadata <- read.table(text = Metadata,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  #Set date type
  Metadata$sample_date <- as.Date(lubridate::ymd(Metadata$sample_date))
  
  #Filter metadata by date range.
  First_Date <- as.Date(lubridate::ymd(First_Date))
  Last_Date <- as.Date(lubridate::ymd(Last_Date))
  Metadata <- Metadata[Metadata$sample_date >= First_Date & Metadata$sample_date <= Last_Date,]
  
  #Filter Tronko-assign data by remaining samples.
  TronkoDB <- TronkoDB[TronkoDB$SampleID %in% Metadata$sample_id,]
  
  #Filter Tronko-assign by read counts per sample
  TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID) %>% dplyr::filter(n() > CountThreshold)
  
  #Filter by relative abundance per taxon per sample.
  TronkoDB <- TronkoDB[!is.na(TronkoDB[,TaxonomicRank]),]
  TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID,!!sym(TaxonomicRank)) %>% 
    dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n/sum(n)) %>% 
    dplyr::ungroup() %>% dplyr::filter(freq > FilterThreshold) %>% select(-n,-freq)
  TronkoDB <- as.data.frame(TronkoDB)
  
  #Get unique taxa list from Tronko-assign
  Tronko_Taxa <- na.omit(unique(TronkoDB[,TaxonomicRank]))
  
  #Read in GBIF occurrences.
  gbif <- gbif_local()
  
  #Get unique states and nations in project.
  country_list <- na.omit(unique(Metadata$Nation))
  state_province_list <- na.omit(unique(Metadata$State))
  
  if(Geographic_Scale=="Local"){
    #Get local bounds for sample locations, add 0.5 degree buffer.
    Local_East <- max(na.omit(Metadata$longitude))+0.5
    Local_West <- min(na.omit(Metadata$longitude))-0.5
    Local_South <- min(na.omit(Metadata$latitude))-0.5
    Local_North <- max(na.omit(Metadata$latitude))+0.5
    
    #Clip GBIF occurrence locations by local boundaries.
    Taxa_Local <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                  coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                  occurrencestatus=="PRESENT",
                                  decimallongitude >= Local_West & decimallongitude <= Local_East & decimallatitude >= Local_South & decimallatitude <= Local_North) %>% 
      select(!!sym(TaxonomicRank))
    Taxa_Local <- as.data.frame(Taxa_Local)
    Taxa_Local <- na.omit(unique(Taxa_Local[,TaxonomicRank]))
    venn_list <- toJSON(list(GBIF=Taxa_Local,eDNA=Tronko_Taxa))
  }
  if(Geographic_Scale=="State"){
    #Clip GBIF occurrence locations by state/province boundaries.
    Taxa_State <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                  coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                  occurrencestatus=="PRESENT", stateprovince %in% state_province_list) %>% select(!!sym(TaxonomicRank))
    Taxa_State <- as.data.frame(Taxa_State)
    Taxa_State <- na.omit(unique(Taxa_State[,TaxonomicRank]))
    venn_list <- toJSON(list(GBIF=Taxa_State,eDNA=Tronko_Taxa))
  }
  if(Geographic_Scale=="Nation"){
    #Clip GBIF occurrence locations by national boundaries.
    Taxa_Nation <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                   coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                   occurrencestatus=="PRESENT", countrycode %in% country_list) %>% select(!!sym(TaxonomicRank))
    Taxa_Nation <- as.data.frame(Taxa_Nation)
    Taxa_Nation <- na.omit(unique(Taxa_Nation[,TaxonomicRank]))
    venn_list <- toJSON(list(GBIF=Taxa_Nation,eDNA=Tronko_Taxa))
  }
  return(venn_list)
}
