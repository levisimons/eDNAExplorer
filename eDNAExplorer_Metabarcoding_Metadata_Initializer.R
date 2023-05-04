rm(list=ls())
require(tidyr)
require(sf)
require(sp)
require(lubridate)
require(httr)
require(curl)
httr::set_config(httr::config(http_version = 2))
curl::handle_setopt(new_handle(),http_version=2)
require(gbifdb)
require(jsonlite)
require(rgbif)
require(data.table)
require(dplyr)
require(DBI)
require(RPostgreSQL)
require(digest)

#Establish database credentials.
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
Database_Driver <- dbDriver("PostgreSQL")
#Force close any possible postgreSQL connections.
sapply(dbListConnections(Database_Driver), dbDisconnect)

ProjectID <- "LARiverRound1" #This is hard-coded for now.

#Read in initial metadata.
Metadata_Initial <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/InputMetadata.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Initial <- read.table(text = Metadata_Initial,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
Metadata_Initial$`Sample Date` <- lubridate::mdy(Metadata_Initial$`Sample Date`)
#Get field variables from initial metadata.
Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% c("Sample ID","Longitude","Latitude","Sample Date","Spatial Uncertainty"))]
#Read in extracted metadata.
Metadata_Extracted <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/MetadataOutput.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Extracted <- read.table(text = Metadata_Extracted,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
Metadata_Extracted$Sample_Date <- lubridate::ymd_hms(Metadata_Extracted$Sample_Date)

#Merge metadata
Metadata <- dplyr::left_join(Metadata_Initial[,c("Sample ID","Sample Date","Latitude","Longitude","Spatial Uncertainty",Field_Variables)],Metadata_Extracted,by=c("Sample ID"="name","Sample Date"="Sample_Date","Latitude","Longitude","Spatial Uncertainty"="Spatial_Uncertainty"))

#Add project ID
Metadata$ProjectID <- ProjectID

#Add Fastq ID to make sure metadata and Tronko-assign output lines up.
Metadata$FastqID <- gsub("_R.*","",Metadata$`Fastq Forward Reads Filename`)

#Read in state/province boundaries.
#Boundaries are from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
sf_use_s2(FALSE)
SpatialBucket <- system("aws s3 ls s3://ednaexplorer/spatial --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",intern=TRUE)
SpatialBucket <- read.table(text = paste(SpatialBucket,sep = ""),header = FALSE)
colnames(SpatialBucket) <- c("Date", "Time", "Size","Filename")
SpatialFiles <- unique(SpatialBucket$Filename)
for(SpatialFile in SpatialFiles){
  system(paste("aws s3 cp s3://ednaexplorer/",SpatialFile," . --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
}
GADM_1_Boundaries <- sf::st_read("ne_10m_admin_1_states_provinces.shp")
#Determine the unique list of national and state/proving boundaries sample locations cover.
GADM_Boundaries <- st_join(st_as_sf(Metadata[!is.na(Metadata$Latitude) & !is.na(Metadata$Latitude),], coords = c("Longitude", "Latitude"), crs = 4326), GADM_1_Boundaries[c('iso_a2','woe_name')],join = st_intersects)
GADM_Boundaries <- GADM_Boundaries %>% st_drop_geometry()
GADM_Boundaries <- as.data.frame(GADM_Boundaries[,c("Sample ID","Sample Date","iso_a2","woe_name")])
names(GADM_Boundaries)[names(GADM_Boundaries) == "woe_name"] <- "State"
names(GADM_Boundaries)[names(GADM_Boundaries) == "iso_a2"] <- "Nation"
Metadata <- dplyr::left_join(Metadata,GADM_Boundaries,by=c("Sample ID","Sample Date"))
country_list <- na.omit(unique(GADM_Boundaries$Nation))
state_province_list <- na.omit(unique(GADM_Boundaries$State))

#Remove rows without an associated sequence file.
Metadata <- Metadata[!is.na(Metadata$`Fastq Forward Reads Filename`) & !is.na(Metadata$`Fastq Reverse Reads Filename`),]

#Generate unique code for each Tronko-assign output
Metadata$UniqueID <- sapply(paste(Metadata$ProjectID,Metadata$FastqID,Metadata$`Sample Date`,Metadata$Latitude,Metadata$Longitude,Metadata$`Spatial Uncertainty`),digest,algo="md5")

#Match metadata column names to format in SQL database.
colnames(Metadata) <- gsub(" ","_",tolower(colnames(Metadata)))

#Create Metadata database.
con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)

#Check for redundant data.
#Add new metadata.
if(dbExistsTable(con,"TronkoMetadata")){
  Metadata_Check <-  tbl(con,"TronkoMetadata")
  Metadata_IDs <- Metadata$uniqueid
  Metadata_Check <- Metadata_Check %>% filter(uniqueid %in% Metadata_IDs)
  Metadata_Check <- as.data.frame(Metadata_Check)
  Metadata_Check_IDs <- Metadata_Check$uniqueid
  Metadata_Append <- Metadata[!(Metadata_IDs %in% Metadata_Check_IDs),]
  dbWriteTable(con,"TronkoMetadata",Metadata_Append,row.names=FALSE,append=TRUE)
} else{
  dbWriteTable(con,"TronkoMetadata",Metadata,row.names=FALSE,append=TRUE)
}
RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
system("rm ne_10m_admin_1_states_provinces.*")
