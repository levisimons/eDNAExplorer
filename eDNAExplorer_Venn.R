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
require(data.table)
require(DBI)
require(RPostgreSQL)
require(digest)

Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")

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
#* @param SpeciesList:string Name of csv file containing selected species list.
#* @get /venn
venn <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,Geographic_Scale,SpeciesList){
  CategoricalVariables <- c("grtgroup","biome_type","IUCN_CAT","ECO_NAME","HYBAS_ID")
  ContinuousVariables <- c("bio01","bio12","gHM","elevation","NDVI","Average_Radiance")
  FieldVars <- c("FastqID","Sample Date","Latitude","Longitude","Spatial Uncertainty")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  First_Date <- lubridate::ymd(First_Date)
  Last_Date <- lubridate::ymd(Last_Date)
  Num_Mismatch <- as.numeric(Num_Mismatch)
  CountThreshold <- as.numeric(CountThreshold)
  FilterThreshold <- as.numeric(FilterThreshold)
  SelectedSpeciesList <- as.character(paste(SpeciesList,".csv",sep=""))
           
  #Read in species list
  if(SelectedSpeciesList != "None.csv"){
    SpeciesList_df <- system(paste("aws s3 cp s3://ednaexplorer/specieslists/",SelectedSpeciesList," - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    SpeciesList_df <- read.table(text = paste(SpeciesList_df,sep = ","),header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  }
  
  #Establish sql connection
  Database_Driver <- dbDriver("PostgreSQL")
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Read in metadata and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  Metadata <- tbl(con,"TronkoMetadata")
  #Keep_Vars <- c(CategoricalVariables,ContinuousVariables,FieldVars)[c(CategoricalVariables,ContinuousVariables,FieldVars) %in% dbListFields(con,"TronkoMetadata")]
  Metadata <- Metadata %>% filter(sample_date >= First_Date & sample_date <= Last_Date) %>%
    filter(ProjectID == ProjectID) %>% filter(!is.na(latitude) & !is.na(longitude))
  Metadata <- as.data.frame(Metadata)
  dbDisconnect(con)
  
  #Create sample metadata matrix
  Sample <- Metadata[!is.na(Metadata$fastqid),]
  rownames(Sample) <- Sample$fastqid
  Sample$fastqid <- NULL
  Sample <- sample_data(Sample)
  remaining_Samples <- rownames(Sample)
  
  #Read in Tronko output and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  TronkoInput <- tbl(con,"TronkoOutput")
  if(sample_TaxonomicRank != "species"){
    TronkoInput <- TronkoInput %>% filter(ProjectID == sample_ProjectID) %>% filter(Primer == sample_Primer) %>% 
      filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(sample_TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > sample_CountThreshold) %>% 
      select(SampleID,species,sample_TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    if(SelectedSpeciesList != "None.csv"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample) & TronkoDB$species %in% SpeciesList_df$Species,]}
    if(SelectedSpeciesList == "None.csv"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample),]}
    TronkoDB$species <- NULL
  } else{
    TronkoInput <- TronkoInput %>% filter(ProjectID == sample_ProjectID) %>% filter(Primer == sample_Primer) %>% 
      filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(sample_TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > sample_CountThreshold) %>% 
      select(SampleID,sample_TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    if(SelectedSpeciesList != "None.csv"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample) & TronkoDB$species %in% SpeciesList_df$Species,]}
    if(SelectedSpeciesList == "None.csv"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample),]}
  }
           
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Filter by relative abundance per taxon per sample.
  if(nrow(TronkoDB) > 1){
             TronkoDB <- TronkoDB[!is.na(TronkoDB[,TaxonomicRank]),]
             TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID,!!sym(TaxonomicRank)) %>% 
               dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n/sum(n)) %>% 
               dplyr::ungroup() %>% dplyr::filter(freq > FilterThreshold) %>% select(-n,-freq)
             TronkoDB <- TronkoDB %>% dplyr::group_by(!!sym(TaxonomicRank)) %>% dplyr::summarise(per=n()/length(unique(TronkoDB$SampleID)))
             TronkoDB <- as.data.frame(TronkoDB)

             #Get unique taxa list from Tronko-assign
             Tronko_Taxa <- na.omit(unique(TronkoDB[,TaxonomicRank]))
  } else {
           Tronko_Taxa <- c()
  }
  
  #Read in GBIF occurrences.
  gbif <- gbif_local()
  
  #Get unique states and nations in project.
  country_list <- na.omit(unique(Metadata$nation))
  state_province_list <- na.omit(unique(Metadata$state))
  
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
