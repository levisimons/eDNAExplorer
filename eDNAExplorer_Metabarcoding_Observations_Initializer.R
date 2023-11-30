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
require(uuid)

# Fetch project ID early so we can use it for error output when possible.
ProjectID <- args[1]

#Establish database credentials.
readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
gbif_dir <- Sys.getenv("GBIF_HOME")
bucket <- Sys.getenv("S3_BUCKET")
taxonomy_home <- Sys.getenv("taxonomy_home")
Database_Driver <- dbDriver("PostgreSQL")
ENDPOINT_URL <- Sys.getenv("ENDPOINT_URL")

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, lastRanAt = Sys.time(), error = error_message))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  dest_filename <- sub("\\.json$", ".build", filename)
  
  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://",bucket,"/errors/observations/", new_filename, sep = "")
  } else {
    paste("s3://",bucket,"/projects/", ProjectID, "/plots/error_", dest_filename, " --endpoint-url ",ENDPOINT_URL, sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}

#Force close any possible postgreSQL connections.
sapply(dbListConnections(Database_Driver), dbDisconnect)

#Set taxonomic rank column names.
TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")

tryCatch(
  {
    #Get project ID.
    #Rscript --vanilla eDNAExplorer_Metabarcoding_Observations_Initializer.R "ProjectID" "TaxonomicRank"
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

tryCatch(
  {
    #Read in project metadata.
    con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)
    Metadata <- tbl(con,"TronkoMetadata")
    Metadata <- Metadata %>% filter(projectid == ProjectID)
    Metadata <- as.data.frame(Metadata)
    
    #Get primers
    Markers <- grep("^marker_[[:digit:]]$",colnames(Metadata),value=T)
    Primers <- na.omit(unique(unlist(Metadata[,Markers])))
    
    #Read in Tronko-assign output files.  Standardize sample IDs within them.
    TronkoBucket <- system(paste("aws s3 ls s3://",bucket,"/tronko_output/",ProjectID," --recursive --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    TronkoBucket <- read.table(text = paste(TronkoBucket,sep = ""),header = FALSE)
    colnames(TronkoBucket) <- c("Date", "Time", "Size","Filename")
    TronkoFiles <- unique(TronkoBucket$Filename)
    TronkoFiles <- TronkoFiles[grepl("*.csv$",TronkoFiles)]
    if(length(TronkoFiles) <= 0){
      stop("Tronko-assign output files needed", call.=FALSE)
    }
    
    #Loop over Tronko-assign data to database of observations
    TronkoObservations <- c()
    i=1
    for(Primer in Primers){
      TronkoFile <- paste("s3://",bucket,"/tronko_output/",ProjectID,"/",Primer,".csv",sep="")
      TronkoFile_tmp <- gsub(".csv",paste("_",UUIDgenerate(),".csv",sep=""),basename(TronkoFile))
      system(paste("aws s3 cp ",TronkoFile," ",TronkoFile_tmp," --endpoint-url ",ENDPOINT_URL,sep=""))
      
      # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
      SubsetFile <- gsub(".csv","_subset.csv",TronkoFile_tmp)
      awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"ProjectID\"], $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",TronkoFile_tmp, SubsetFile)
      system(awk_command, intern = TRUE)
      
      TronkoInput <- fread(file=SubsetFile,header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
      
      for(TaxonomicRank in TaxonomicRanks){
        #Filter Tronko-assign output by columns, mismatches, and taxa at a selected rank.
        TronkoFiltered <- TronkoInput %>% select("ProjectID","SampleID",!!sym(TaxonomicRank),"Mismatch") %>% 
          filter(Mismatch <= 5) %>% filter(!is.na(!!sym(TaxonomicRank))) %>% select(-Mismatch) %>%
          distinct(.keep_all = TRUE)
        
        if(nrow(TronkoFiltered)>0){
          TronkoDB <- as.data.frame(TronkoFiltered)
          names(TronkoDB)[names(TronkoDB) == TaxonomicRank] <- "Taxon"
          TronkoDB$rank <- TaxonomicRank
          TronkoObservations[[i]] <- TronkoDB
          i=i+1
        }
      }
      
      system(paste("rm ",TronkoFile_tmp,sep=""))
      system(paste("rm ",SubsetFile,sep=""))
    }
    
    #Aggregate unique taxa detected from Tronko-assign outputs
    TronkoObservations <- rbindlist(TronkoObservations, use.names=TRUE, fill=TRUE)
    TronkoObservations <- TronkoObservations[!duplicated(TronkoObservations),]
    
    #Get the year and coordinates for each sample in a project.
    MetadataObservations <- Metadata %>% select(projectid,fastqid,sample_id,sample_date,latitude,longitude)
    #Convert the date string to year.
    MetadataObservations$year <- as.integer(format(as.Date(MetadataObservations$sample_date,format("%Y-%m-%d")), "%Y"))
    #Make a sample ID column to match the format with Tronko-assign outputs.
    MetadataObservations$SampleID <- gsub("_","-",MetadataObservations$fastqid)
    
    #Create merged observation data set.
    ObservationsExport <- dplyr::right_join(TronkoObservations,MetadataObservations,by=c("ProjectID"="projectid","SampleID"="SampleID"))
    #Remove entries without taxa
    ObservationsExport <- ObservationsExport[!is.na(ObservationsExport$Taxon),]
    #Rename columns for export
    names(ObservationsExport)[names(ObservationsExport) == "ProjectID"] <- "projectId"
    names(ObservationsExport)[names(ObservationsExport) == "sample_id"] <- "sampleId"
    
    #Designate a unique id
    ObservationsExport$id <- sapply(paste(ObservationsExport$Taxon,ObservationsExport$rank,ObservationsExport$projectId,ObservationsExport$sampleId,ObservationsExport$latitude,ObservationsExport$longitude,ObservationsExport$year),digest,algo="md5")
    
    #Get common names for each taxon observation.
    if (nrow(ObservationsExport) > 0) {
      TaxaList <- na.omit(unique(ObservationsExport$Taxon))
    }
    if (nrow(ObservationsExport) == 0) {
      TaxaList <- c()
    }
    TaxonomicRanks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
    TaxonomicKeyRanks <- c("kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")
    TaxonomyInput <- tbl(con, "Taxonomy")
    TaxonomyInput <- TaxonomyInput %>%
      filter(Taxon %in% TaxaList) %>%
      select(TaxonomicRanks[2:length(TaxonomicRanks)],TaxonomicKeyRanks,Taxon,rank,Common_Name)
    TaxonomyDB <- as.data.frame(TaxonomyInput)
    #Figure out which taxonomy version is more complete.
    TaxonomyDB$rankCount <- rowSums(!is.na(TaxonomyDB[,colnames(TaxonomyDB) %in% TaxonomicRanks]))
    TaxonomyDB$rankKeyCount <- rowSums(!is.na(TaxonomyDB[,colnames(TaxonomyDB) %in% TaxonomicKeyRanks]))
    TaxonomyDB <- TaxonomyDB %>%
      group_by(Taxon) %>%
      slice_max(order_by = rankCount, n = 1) %>%
      ungroup()
    TaxonomyDB <- TaxonomyDB %>%
      group_by(Taxon) %>%
      slice_max(order_by = rankKeyCount, n = 1) %>%
      ungroup()
    #Figure out which common_name is most common per taxon.
    TaxonomyDB <- TaxonomyDB %>%
      group_by(Taxon) %>%
      mutate(Most_Common_Name = ifelse(all(is.na(Common_Name)), NA, names(which.max(table(Common_Name[!is.na(Common_Name)]))))) %>%
      ungroup()
    TaxonomyDB$Common_Name <- TaxonomyDB$Most_Common_Name
    TaxonomyDB$Most_Common_Name <- NULL
    TaxonomyDB <- as.data.frame(TaxonomyDB)
    TaxonomyDB$rankCount <- NULL
    TaxonomyDB <- TaxonomyDB[!duplicated(TaxonomyDB),]
    #Merge in common names to occurrence taxonomies.
    ObservationsExport <- dplyr::left_join(ObservationsExport,TaxonomyDB[,c("Taxon","rank","Common_Name")],by=c("Taxon","rank"))
    
    #Retain key columns
    ObservationsExport <- ObservationsExport[,c("id","Taxon","rank","Common_Name","projectId","sampleId","latitude","longitude","year")]
    
    #Check for redundant data.
    #Add new metadata.
    if(dbExistsTable(con,"Occurence")){
      Observations_Check <-  tbl(con,"Occurence")
      Observations_IDs <- ObservationsExport$id
      Observations_Check <- Observations_Check %>% filter(id %in% Observations_IDs)
      Observations_Check <- as.data.frame(Observations_Check)
      Observations_Check_IDs <- Observations_Check$id
      Observations_Append <- ObservationsExport[!(Observations_IDs %in% Observations_Check_IDs),]
      dbWriteTable(con,"Occurence",Observations_Append,row.names=FALSE,append=TRUE)
    } else{
      dbWriteTable(con,"Occurence",ObservationsExport,row.names=FALSE,append=TRUE)
    }
    
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  },
  error = function(e) {
    process_error(e, filename)
  }
)