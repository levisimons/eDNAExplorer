# plumber.R
library(plumber)
require(aws.s3)
require(tidyr)
require(dplyr)
require(ggplot2)
require(rgbif)
require(gbifdb)
require(RPostgreSQL)
require(lubridate)
require(plotly)
require(jsonlite)

Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")

#* Echo the parameter that was sent in
#* @param Marker:string Target marker name
#* @param Taxon_name:string Scientific taxon name
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @get /timeline

timeline <- function(Marker,Taxon_name,TaxonomicRank,Num_Mismatch,CountThreshold,FilterThreshold){
  
  #Establish sql connection
  Database_Driver <- dbDriver("PostgreSQL")
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Select taxon to map.
  #User input
  Taxon <- Taxon_name
  #Get GBIF taxonomy key for taxon.
  Taxon_GBIF <- name_backbone(name=Taxon,rank=TaxonomicRank)$usageKey
  #Ensure numeric values.
  Num_Mismatch <- as.numeric(Num_Mismatch)
  CountThreshold <- as.numeric(CountThreshold)
  FilterThreshold <- as.numeric(FilterThreshold)
  
  #Read in GBIF occurrences.
  gbif <- gbif_local()
  
  #Filter GBIF occurrences to a particular taxon.
  GBIFDB <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                occurrencestatus=="PRESENT",taxonkey==Taxon_GBIF) %>% select(year)
  GBIFDB <- as.data.frame(GBIFDB)
  
  if(nrow(GBIFDB) < 1){
    Taxa_Time <- data.frame(matrix(nrow=1,ncol=2))
    colnames(Taxa_Time) <- c("year","gbif_occurrences")
  }
  if(nrow(GBIFDB) >=1){
    #Get GBIF taxon occurrences over time.
    Taxa_Time <- as.data.frame(table(GBIFDB))
    colnames(Taxa_Time) <- c("year","gbif_occurrences") 
  }
  
  #Read in Tronko output.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  TronkoInput <- tbl(con,"TronkoOutput")
  
  #Get samples where taxon occurs in general.
  TaxonDB <- TronkoInput %>% filter(!!sym(TaxonomicRank) == Taxon) %>% 
    select(ProjectID,SampleID) %>% distinct_all()
  TaxonDB <- as.data.frame(TaxonDB)
  
  #Filter Tronko output by mismatches, sample count, and relative abundance
  TronkoDB <- TronkoInput %>% filter(Primer == Marker) %>% 
    filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>% 
    group_by(SampleID) %>% filter(n() > CountThreshold) %>% 
    summarise(n=n()) %>% mutate(freq=n/sum(n)) %>% 
    ungroup() %>% filter(freq > FilterThreshold) %>% select(-n,-freq)
  TronkoDB <- as.data.frame(TronkoDB)
  
  #Get samples where taxon occurs and meets Tronko filters.
  TaxonDB <- TaxonDB[TaxonDB$SampleID %in% TronkoDB$SampleID,]
  taxon_samples <- unique(TaxonDB$SampleID)
  taxon_projects <- unique(TaxonDB$ProjectID)
  
  #Read in metadata and filter it.
  Metadata <- tbl(con,"TronkoMetadata")
  Metadata_Filtered <- Metadata %>% filter(!is.na(latitude) & !is.na(longitude)) %>%
    filter(projectid %in% taxon_projects) %>% filter(fastqid %in% taxon_samples) %>%
    select(projectid,sample_date)
  Metadata_Filtered <- as.data.frame(Metadata_Filtered)
  Metadata_Filtered$year <- lubridate::year(lubridate::ymd(Metadata_Filtered$sample_date))
  
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Get projects and years where taxon was detected via eDNA.
  eDNA <- Metadata_Filtered[,c("year","projectid")]
  colnames(eDNA) <- c("year","eDNA_projects")
  eDNA <- eDNA[complete.cases(eDNA),]
  eDNA <- eDNA[!duplicated(eDNA),]
  
  #Merge eDNA and GBIF timeline data.
  Timeline <- merge(Taxa_Time,eDNA,by="year",all=T)
  Timeline <- Timeline[rowSums(is.na(Timeline)) != ncol(Timeline), ]
  
  Timeline <- toJSON(Timeline)
  return(Timeline)
}
