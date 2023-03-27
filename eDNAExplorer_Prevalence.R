# plumber.R
library(plumber)
require(tidyr)
require(dplyr)
require(lubridate)
require(jsonlite)

#* Echo the parameter that was sent in
#* @param ProjectID:string
#* @param First_Date:string YYYY-MM-DD
#* @param Last_Date:string YYYY-MM-DD
#* @param Marker:string Target marker name
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @get /prevalence

prevalence <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold){
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
  First_Date <<- as.Date(lubridate::ymd(First_Date))
  Last_Date <<- as.Date(lubridate::ymd(Last_Date))
  Metadata <- Metadata[Metadata$sample_date >= First_Date & Metadata$sample_date <= Last_Date,]
  
  #Filter Tronko-assign data by remaining samples.
  TronkoDB <- TronkoDB[TronkoDB$SampleID %in% Metadata$sample_id,]
  
  #Filter Tronko-assign by read counts per sample
  TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID) %>% dplyr::filter(n() > CountThreshold)
  
  #Filter by relative abundance per taxon per sample.
  TronkoDB <- TronkoDB[!is.na(TronkoDB[,TaxonomicRank]),]
  KingdomMatch <- TronkoDB[,c("kingdom",TaxonomicRank)]
  KingdomMatch <- KingdomMatch[!duplicated(KingdomMatch),]
  TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID,!!sym(TaxonomicRank)) %>% 
    dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n/sum(n)) %>% 
    dplyr::ungroup() %>% dplyr::filter(freq > FilterThreshold) %>% select(-n,-freq)
  TronkoDB <- TronkoDB %>% dplyr::group_by(!!sym(TaxonomicRank)) %>% dplyr::summarise(per=n()/length(unique(TronkoDB$SampleID)))
  TronkoDB <- dplyr::left_join(TronkoDB,KingdomMatch)
  TronkoDB <- toJSON(TronkoDB)
  return(TronkoDB)  
}
