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
require(data.table)

readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")

#* Echo the parameter that was sent in
#* @param ProjectID:string
#* @param Marker:string Target marker name
#* @param Taxon_name:string Scientific taxon name
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @get /map


map <- function(ProjectID,Marker,Taxon_name,TaxonomicRank,Num_Mismatch,CountThreshold,FilterThreshold){
  #Establish sql connection
  Database_Driver <- dbDriver("PostgreSQL")
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Select taxon to map.
  #User input
  Project_ID <- as.character(ProjectID)
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
                            occurrencestatus=="PRESENT",taxonkey==Taxon_GBIF) %>% select(decimallongitude,decimallatitude)
  GBIFDB <- as.data.frame(GBIFDB)
  colnames(GBIFDB) <- c("longitude","latitude")
  
  #Get unique taxon locations
  if(nrow(GBIFDB) < 1){
    TaxonMap <- data.frame(matrix(nrow=1,ncol=3))
    colnames(TaxonMap) <- c("source","longitude","latitude")
    TaxonMap$source <- "GBIF"
  }
  if(nrow(GBIFDB) >=1){
    TaxonMap <- GBIFDB[,c("longitude","latitude")]
    TaxonMap <- TaxonMap[complete.cases(TaxonMap),]
    TaxonMap <- TaxonMap[!duplicated(TaxonMap),]
    TaxonMap$source <- "GBIF"
    TaxonMap <- TaxonMap[,c("source","longitude","latitude")]
  }
  
  #Read in Tronko output and filter it.
  TronkoFile <- paste(Marker,".csv",sep="")
  system(paste("aws s3 cp s3://ednaexplorer/tronko_output/",Project_ID,"/",TronkoFile," ",TronkoFile," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
  system(paste("cut -d ',' -f 2,6,7,8,9,10,11,12,13,14,16 ",TronkoFile," > subset.csv",sep=""))
  TronkoInput <- fread(file="subset.csv",header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
  if(TaxonomicRank != "species"){
    TronkoInput <- TronkoInput %>% filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > CountThreshold) %>% 
      select(ProjectID,SampleID,species,TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    TronkoDB$species <- NULL
  } else{
    TronkoInput <- TronkoInput %>% filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > CountThreshold) %>% 
      select(ProjectID,SampleID,TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
  }
  TaxonDB <- TronkoDB[!duplicated(TronkoDB),]
  TaxonDB$SampleID <- gsub("-","_",TaxonDB$SampleID)
  
  #Get samples where taxon occurs and meets Tronko filters.
  taxon_samples <- unique(TaxonDB$SampleID)
  taxon_projects <- unique(TaxonDB$ProjectID)
  
  #Read in metadata and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  Metadata <- tbl(con,"TronkoMetadata")
  Metadata_Filtered <- Metadata %>% filter(!is.na(latitude) & !is.na(longitude)) %>%
    filter(projectid %in% taxon_projects) %>% filter(fastqid %in% taxon_samples) %>%
    select(longitude,latitude)
  Metadata_Filtered <- as.data.frame(Metadata_Filtered)
  
  Metadata_Filtered$source <- "eDNA"
  Metadata_Filtered <- Metadata_Filtered[,c("source","longitude","latitude")]
  
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  Taxon_Map_Data <- rbind(TaxonMap,Metadata_Filtered)
  #Return results
  Taxon_Map_Data <- jsonlite::toJSON(Taxon_Map_Data)
  filename <- paste("Map_Metabarcoding_Marker_",Marker,"_Taxon_",Taxon,"_Rank_",TaxonomicRank,"_Mismatch_",Num_Mismatch,"_CountThreshold_",CountThreshold,"_AbundanceThreshold_",format(FilterThreshold,scientific=F),".json",sep="")
  filename <- tolower(filename)
  filename <- gsub(" ","_",filename)
  write(Taxon_Map_Data,filename)
  system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  system(paste("rm ",filename,sep=""))
}
