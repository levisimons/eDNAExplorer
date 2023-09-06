#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
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
require(anytime)

# Fetch project ID early so we can use it for error output when possible.
ProjectID <- args[1]

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, error = error_message))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  
  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://ednaexplorer/errors/metadata/", new_filename, sep = "")
  } else {
    paste("s3://ednaexplorer/tronko_output/", ProjectID, "/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system("rm ne_10m_admin_1_states_provinces.*")
  system(paste("rm ",filename,sep=""))
  RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  stop(error_message)
}

#Establish database credentials.
readRenviron(".env")
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

tryCatch(
  {
    #Get project ID.
    #Rscript --vanilla eDNAExplorer_Metabarcoding_Metadata_Initializer.R "project ID string"
    if (length(args)==0) {
      stop("Need a project ID", call.=FALSE)
    } else if (length(args)==1) {
      ProjectID <- args[1]
    }
  },
  error = function(e) {
    process_error(e)
  }
)

# Update taxonomy tables and cache processed Tronko-assign output.
tryCatch(
  {
    #Find metabarcoding project data file and read it into a dataframe.
    Project_Data <- system(paste("aws s3 cp s3://ednaexplorer/projects",ProjectID,"METABARCODING.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep="/"),intern=TRUE)
    Project_Data <- gsub("[\r\n]", "", Project_Data)
    if(length(Project_Data)==0) {
      stop("Error: No initial metadata present.")
    }
    Project_Data <- read.table(text = Project_Data,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
    names(Project_Data) <- gsub(x = names(Project_Data), pattern = "ForwardPS", replacement = "Forward PS")
    names(Project_Data) <- gsub(x = names(Project_Data), pattern = "ReversePS", replacement = "Reverse PS")
    addFormats(c("%m/%d/%y","%m-%d-%y","%d/%m/%y","%y/%m/%d"))
    Project_Data$`Sample Date` <- as.Date(as.character(parse_date_time(Project_Data$`Sample Date`, orders = c("ymd", "dmy", "mdy"))))
    Project_Data$`Data type` <- NULL
    Project_Data$`Additional environmental metadata....` <- NULL
    #Remove zero length variable names
    Project_Data <- Project_Data[,nchar(colnames(Project_Data))>0]
    Project_Data <- Project_Data %>% dplyr::mutate_at(c("Latitude","Longitude","Spatial Uncertainty"),as.numeric)
    Project_Data <- as.data.frame(Project_Data)
    Metadata_Initial <- Project_Data
    
    Required_Variables <- c("Site","Sample ID","Sample Type","Longitude","Latitude","Sample Date","Sequencing Platform","Sequence Length","Adapter type","Fastq Forward Reads Filename","Fastq Reverse Reads Filename",grep("^Marker [[:digit:]]$",colnames(Metadata_Initial),value=T),grep("^Marker [[:digit:]] Forward PS$",colnames(Metadata_Initial),value=T),grep("^Marker [[:digit:]] Reverse PS$",colnames(Metadata_Initial),value=T))
    #Get field variables from initial metadata.  These are generally project-specific non-required variables.
    Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% Required_Variables)]
    #Read in extracted metadata.
    Metadata_Extracted <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/MetadataOutput_Metabarcoding.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    if(length(Metadata_Extracted)==0) {
      stop("Error: No extracted metadata present.")
    }
    Metadata_Extracted <- read.table(text = Metadata_Extracted,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
    Metadata_Extracted$Sample_Date <- lubridate::ymd_hms(Metadata_Extracted$Sample_Date)
    Metadata_Extracted$Sample_Date <- as.Date(as.character(as.POSIXct(Metadata_Extracted$Sample_Date)))
    #Set no data results.
    Metadata_Extracted[Metadata_Extracted==-999999] <- NA
    Metadata_Extracted[Metadata_Extracted==-32768] <- NA
    
    #Merge metadata
    Metadata <- dplyr::left_join(Metadata_Initial[,Required_Variables],Metadata_Extracted,by=c("Sample ID"="name","Sample Date"="Sample_Date","Latitude","Longitude"),na_matches="never")
    
    #Add project ID
    Metadata$ProjectID <- ProjectID
    
    #Add Fastq ID to make sure metadata and Tronko-assign output lines up.
    Metadata$FastqID <- gsub("R1_001.fastq.gz","",Metadata$`Fastq Forward Reads Filename`)
    Metadata$FastqID <- gsub("_$","",Metadata$FastqID)
    
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
    system("rm ne_10m_admin_1_states_provinces.*")
    
    #Remove rows without an associated sequence file.
    Metadata <- Metadata[!is.na(Metadata$`Fastq Forward Reads Filename`) & !is.na(Metadata$`Fastq Reverse Reads Filename`),]
    
    #Generate unique code for each Tronko-assign output
    Metadata$UniqueID <- sapply(paste(Metadata$ProjectID,Metadata$FastqID,Metadata$`Sample Date`,Metadata$Latitude,Metadata$Longitude,Metadata$`Spatial Uncertainty`),digest,algo="md5")
    
    #Match metadata column names to format in SQL database.
    colnames(Metadata) <- gsub(" ","_",tolower(colnames(Metadata)))
    
    #Set character and numeric columns.
    col_non_numeric <- c("adapter_type","name","biome_type","eco_name","fastq_forward_reads_filename","fastqid","fastq_reverse_reads_filename",
                         "grtgroup","hybas_id","marker_1","nation","projectid","realm","sample_date","sample_id","sample_type","sequencing_platform",
                         "site","state","desig_eng","gov_type","iucn_cat","uniqueid","marker_1_forward_ps","marker_1_reverse_ps","landform",
                         "wdpa_pid","marker_10","marker_10_forward_ps","marker_10_reverse_ps","marker_2","marker_2_forward_ps","marker_2_reverse_ps",
                         "marker_3","marker_3_forward_ps","marker_3_reverse_ps","marker_4","marker_4_forward_ps","marker_4_reverse_ps","marker_5",
                         "marker_5_forward_ps","marker_5_reverse_ps","marker_6","marker_6_forward_ps","marker_6_reverse_ps","marker_7",
                         "marker_7_forward_ps","marker_7_reverse_ps","marker_8","marker_8_forward_ps","marker_8_reverse_ps","marker_9",
                         "marker_9_forward_ps","marker_9_reverse_ps")
    col_numeric <- colnames(Metadata)[!(colnames(Metadata) %in% col_non_numeric)]
    Metadata[, col_numeric] <- lapply(Metadata[, col_numeric], as.numeric)
    
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
  },
  error = function(e) {
    process_error(e, filename)
  }
)
