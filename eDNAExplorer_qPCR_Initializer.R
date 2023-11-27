#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
require(tidyr)
require(sf)
#require(sp)
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

#Establish database credentials.
readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"),
           "PROJ_LIB" = Sys.getenv("PROJ_LIB"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
bucket <- Sys.getenv("S3_BUCKET")
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
    paste("s3://",bucket,"/errors/qpcr/", new_filename, sep = "")
  } else {
    paste("s3://",bucket,"/projects/", ProjectID, "/plots/", dest_filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}

#Force close any possible postgreSQL connections.
sapply(dbListConnections(Database_Driver), dbDisconnect)

tryCatch(
  {
    #Get project ID.
    #Rscript --vanilla eDNAExplorer_qPCR_Initializer.R "project ID string"
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
    #Read in qPCR project data.
    Project_Data <- system(paste("aws s3 cp s3://",bucket,"/projects/",ProjectID,"/QPCR.csv - --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    if(length(Project_Data)==0) {
      stop("Error: No initial metadata present.")
    }
    Project_Data <- read.table(text = Project_Data,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
    Project_Data$`Sample Date` <- as.Date(as.character(parse_date_time(Project_Data$`Sample Date`, orders = c("ymd","mdy","dmy"))))
    Project_Data$`Data type` <- NULL
    Project_Data$`Additional environmental metadata....` <- NULL
    #Remove zero length variable names
    Project_Data <- Project_Data[,nchar(colnames(Project_Data))>0]
    #Define numerical variables
    Project_Data <- Project_Data %>% dplyr::mutate_at(c("Latitude","Longitude","Spatial Uncertainty"),as.numeric)
    Project_Data <- as.data.frame(Project_Data)
    colnames(Project_Data) <- gsub('qPCR Probe Fluorophore \\(dye\\)','qPCR Probe Fluorophore',colnames(Project_Data))
    colnames(Project_Data) <- gsub('Cycle Threshold \\(ct\\)','Cycle Threshold',colnames(Project_Data))
    Metadata_Initial <- Project_Data
    
    Required_Variables <- c("Site","Sample ID","Sample Replicate Number","Sample Type","Longitude","Latitude","Sample Date",grep("^Target [[:digit:]] Organism$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Taxonomic Rank of Organism$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] qPCR Probe Fluorophore$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Primer Probe Set Name$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] ForwardPS$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] ReversePS$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Probe$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Cycle Threshold$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Min RFU threshold$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Final RFU$",colnames(Metadata_Initial),value=T),grep("^Target [[:digit:]] Was organism detected$",colnames(Metadata_Initial),value=T))
    #Get field variables from initial metadata.  These are generally project-specific non-required variables.
    Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% Required_Variables)]
    #Get target organism variables from initial metadata.
    Target_Variables <- colnames(Metadata_Initial)[grep("^Target [[:digit:]]",colnames(Metadata_Initial))]
    #Get non-target variables from initial metadata.
    Non_Target_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% Target_Variables)]
    #Get number of target organisms.
    Target_Numbers <- unique(as.numeric(gsub("\\D", "", Target_Variables)))
    
    #Restructure initial metadata to keep all target organism results within a single set of columns
    Metadata_Initial_List <- c()
    i=1
    for(Target_Number in Target_Numbers){
      Target_Num_Variables <- colnames(Metadata_Initial)[grep(paste("^Target ",Target_Number,sep=""),colnames(Metadata_Initial))]
      Metadata_Subset <- Metadata_Initial[,c(Non_Target_Variables,Target_Num_Variables)]
      colnames(Metadata_Subset) <- c(Non_Target_Variables,gsub(paste("^Target ",Target_Number," ",sep=""),"Target ",Target_Num_Variables))
      Metadata_Initial_List[[i]] <- Metadata_Subset
      i=i+1
    }
    Metadata_Initial <- rbindlist(Metadata_Initial_List, use.names=TRUE, fill=TRUE)
    Metadata_Initial <- as.data.frame(Metadata_Initial)
    
    Metadata_Extracted <- system(paste("aws s3 cp s3://",bucket,"/projects/",ProjectID,"/MetadataOutput_qPCR.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    if(length(Metadata_Extracted)==0) {
      stop("Error: No extracted metadata present.")
    }
    Metadata_Extracted <- read.table(text = Metadata_Extracted,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A","n/a","na"))
    Metadata_Extracted$Sample_Date <- lubridate::ymd_hms(Metadata_Extracted$Sample_Date)
    Metadata_Extracted$Sample_Date <- as.Date(as.character(as.POSIXct(Metadata_Extracted$Sample_Date)))
    #Set no data results.
    Metadata_Extracted[Metadata_Extracted==-999999] <- NA
    Metadata_Extracted[Metadata_Extracted==-32768] <- NA
    #Deduplicate
    Metadata_Extracted <- Metadata_Extracted[!duplicated(Metadata_Extracted),]
    
    #Merge metadata
    Required_Variables <- c("Site","Sample ID","Sample Replicate Number","Sample Type","Longitude","Latitude","Sample Date","Target Organism","Target Taxonomic Rank of Organism","Target qPCR Probe Fluorophore","Target Primer Probe Set Name","Target ForwardPS","Target Probe","Target ReversePS","Target Cycle Threshold","Target Min RFU threshold","Target Final RFU","Target Was organism detected")
    Metadata <- dplyr::right_join(Metadata_Initial[,Required_Variables],Metadata_Extracted,by=c("Sample ID"="name","Sample Date"="Sample_Date","Latitude","Longitude"),na_matches="never",multiple="all",relationship="many-to-one")
    Metadata <- Metadata[!duplicated(Metadata),]
    
    #Add project ID
    MergedData <- Metadata
    MergedData$ProjectID <- ProjectID
    
    #Read in state/province boundaries.
    #Boundaries are from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
    sf_use_s2(FALSE)
    SpatialBucket <- system(paste("aws s3 ls s3://",bucket,"/spatial --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    SpatialBucket <- read.table(text = paste(SpatialBucket,sep = ""),header = FALSE)
    colnames(SpatialBucket) <- c("Date", "Time", "Size","Filename")
    SpatialFiles <- unique(SpatialBucket$Filename)
    for(SpatialFile in SpatialFiles){
      system(paste("aws s3 cp s3://",bucket,"/",SpatialFile," . --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
    }
    GADM_1_Boundaries <- sf::st_read("ne_10m_admin_1_states_provinces.shp")
    #Determine the unique list of national and state/proving boundaries sample locations cover.
    GADM_Boundaries <- st_join(st_as_sf(MergedData[!is.na(MergedData$Latitude) & !is.na(MergedData$Latitude),], coords = c("Longitude", "Latitude"), crs = 4326), GADM_1_Boundaries[c('iso_a2','woe_name')],join = st_intersects)
    GADM_Boundaries <- GADM_Boundaries %>% st_drop_geometry()
    GADM_Boundaries <- as.data.frame(GADM_Boundaries[,c("Sample ID","Sample Date","iso_a2","woe_name")])
    names(GADM_Boundaries)[names(GADM_Boundaries) == "woe_name"] <- "State"
    names(GADM_Boundaries)[names(GADM_Boundaries) == "iso_a2"] <- "Nation"
    GADM_Boundaries <- GADM_Boundaries[!duplicated(GADM_Boundaries),]
    MergedData <- dplyr::left_join(MergedData,GADM_Boundaries,by=c("Sample ID","Sample Date"),na_matches = "never")
    country_list <- na.omit(unique(GADM_Boundaries$Nation))
    state_province_list <- na.omit(unique(GADM_Boundaries$State))
    #Remove spatial files
    system("rm ne_10m_admin_1_states_provinces.*")
    
    #Generate unique code for each sample
    MergedData <- MergedData[!duplicated(MergedData),]
    MergedData$UniqueID <- sapply(apply(MergedData[,Required_Variables],1,paste,collapse = "" ),digest,algo="md5")
    
    #Read in GBIF occurrences.
    gbif <- gbif_local()
    
    #Get local bounds for sample locations, add 0.5 degree buffer.
    Local_East <- max(na.omit(MergedData$Longitude))+0.5
    Local_West <- min(na.omit(MergedData$Longitude))-0.5
    Local_South <- min(na.omit(MergedData$Latitude))-0.5
    Local_North <- max(na.omit(MergedData$Latitude))+0.5
    
    #Clip GBIF occurrence locations by local boundaries.
    Taxa_Local <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                  coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                  occurrencestatus=="PRESENT",
                                  decimallongitude >= Local_West & decimallongitude <= Local_East & decimallatitude >= Local_South & decimallatitude <= Local_North) %>% 
      select(decimallongitude,decimallatitude,taxonkey,taxonrank,species,genus,family,order,class,phylum,kingdom)
    Taxa_Local <- as.data.frame(Taxa_Local)
    #Assign weight value for species, genera, and families.
    Taxa_Local <- Taxa_Local %>% dplyr:: mutate(Local_GBIFWeight = dplyr::case_when(taxonrank=="SPECIES" ~ 4, taxonrank=="SUBSPECIES" ~ 4, taxonrank=="GENUS" ~ 2, taxonrank=="FAMILY" ~ 1, !(taxonrank %in% c("SPECIES","GENUS","FAMILY"))~0))
    
    #Clip GBIF occurrence locations by state/province boundaries.
    Taxa_State <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                  coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                  occurrencestatus=="PRESENT", stateprovince %in% state_province_list) %>% select(decimallongitude,decimallatitude,taxonkey,taxonrank,species,genus,family,order,class,phylum,kingdom)
    Taxa_State <- as.data.frame(Taxa_State)
    #Assign weight value for species, genera, and families.
    Taxa_State <- Taxa_State %>% dplyr:: mutate(State_GBIFWeight = dplyr::case_when(taxonrank=="SPECIES" ~ 4, taxonrank=="SUBSPECIES" ~ 4, taxonrank=="GENUS" ~ 2, taxonrank=="FAMILY" ~ 1, !(taxonrank %in% c("SPECIES","GENUS","FAMILY"))~0))
    
    #Clip GBIF occurrence locations by national boundaries.
    Taxa_Nation <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                   coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                   occurrencestatus=="PRESENT", countrycode %in% country_list) %>% select(decimallongitude,decimallatitude,taxonkey,taxonrank,species,genus,family,order,class,phylum,kingdom)
    Taxa_Nation <- as.data.frame(Taxa_Nation)
    #Assign weight value for species, genera, and families.
    Taxa_Nation <- Taxa_Nation %>% dplyr:: mutate(Ecoregion_GBIFWeight = dplyr::case_when(taxonrank=="SPECIES" ~ 4, taxonrank=="SUBSPECIES" ~ 4, taxonrank=="GENUS" ~ 2, taxonrank=="FAMILY" ~ 1, !(taxonrank %in% c("SPECIES","GENUS","FAMILY"))~0))
    
    #Get unique taxa with full taxonomy
    TaxaList <- unique(MergedData[,"Target Organism"])
    Taxa <- list()
    i=1
    TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
    for(Taxon in TaxaList){
      Taxon_Rank <- unique(MergedData[MergedData$`Target Organism`==Taxon,"Target Taxonomic Rank of Organism"])
      Taxon_GBIF <- name_backbone(name=Taxon,rank=Taxon_Rank,verbose=T,strict=F,curlopts=list(http_version=2))
      Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",colnames(Taxon_GBIF) %in% c(TaxonomicRanks)]
      if(nrow(Taxon_GBIF)>0){
        Taxon_GBIF <- Taxon_GBIF[1,]
        Taxon_GBIF$Taxon <- Taxon
        Taxa[[i]] <- Taxon_GBIF
        i=i+1
      }
    }
    TaxaDB <- rbindlist(Taxa, use.names=TRUE, fill=TRUE)
    TaxaDB <- as.data.frame(TaxaDB[!duplicated(TaxaDB),])
    #Count eDNA taxonomic resolution and weigh them.
    #Species = 4, Genus = 2, Family = 1.  Everything else = 0.
    tmp <- TaxaDB[,c("Taxon","family","genus","species")]
    tmp <- tmp[!duplicated(tmp),]
    tmp$eDNAWeight <- 1*as.numeric(!is.na(tmp$family))+2*as.numeric(!is.na(tmp$genus))+4*as.numeric(!is.na(tmp$species))
    #Check if any of the eDNA reads show up in the local set of GBIF family observations.
    tmp$LocalFamilyPresentGBIF <- as.numeric(lapply(tmp$family,is.element,unique(na.omit(Taxa_Local$family))))
    #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
    tmp$StateFamilyPresentGBIF <- as.numeric(lapply(tmp$family,is.element,unique(na.omit(Taxa_State$family))))
    #Check if any of the eDNA reads show up in the realm set of GBIF family observations.
    tmp$NationFamilyPresentGBIF <- as.numeric(lapply(tmp$family,is.element,unique(na.omit(Taxa_Nation$family))))
    #Check if any of the eDNA reads show up in the local set of GBIF family observations.
    tmp$LocalGenusPresentGBIF <- as.numeric(lapply(tmp$genus,is.element,unique(na.omit(Taxa_Local$genus))))
    #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
    tmp$StateGenusPresentGBIF <- as.numeric(lapply(tmp$genus,is.element,unique(na.omit(Taxa_State$genus))))
    #Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
    tmp$NationGenusPresentGBIF <- as.numeric(lapply(tmp$genus,is.element,unique(na.omit(Taxa_Nation$genus))))
    #Check if any of the eDNA reads show up in the local set of GBIF family observations.
    tmp$LocalSpeciesPresentGBIF <- as.numeric(lapply(tmp$species,is.element,unique(na.omit(Taxa_Local$species))))
    #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
    tmp$StateSpeciesPresentGBIF <- as.numeric(lapply(tmp$species,is.element,unique(na.omit(Taxa_State$species))))
    #Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
    tmp$NationSpeciesPresentGBIF <- as.numeric(lapply(tmp$species,is.element,unique(na.omit(Taxa_Nation$species))))
    
    #Assign TOS scores for GBIF results.
    tmp$TOS_Local <- (1*tmp$LocalFamilyPresentGBIF+2*tmp$LocalGenusPresentGBIF+4*tmp$LocalSpeciesPresentGBIF)/tmp$eDNAWeight
    tmp$TOS_State <- (1*tmp$StateFamilyPresentGBIF+2*tmp$StateGenusPresentGBIF+4*tmp$StateSpeciesPresentGBIF)/tmp$eDNAWeight
    tmp$TOS_Nation <- (1*tmp$NationFamilyPresentGBIF+2*tmp$NationGenusPresentGBIF+4*tmp$NationSpeciesPresentGBIF)/tmp$eDNAWeight
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    tmp[is.nan(tmp)] <- 0
    #Calculate a geographically weighted TOS.
    tmp$gw_TOS <- (4*tmp$TOS_Local+2*tmp$TOS_State+1*tmp$TOS_Nation)/7
    
    #Merge TOS results back into merged data set.
    MergedData <- dplyr::left_join(MergedData,tmp[,c("Taxon","TOS_Local","TOS_State","TOS_Nation","gw_TOS")],by=c("Target Organism"="Taxon"))
    
    #Match metadata column names to format in SQL database.
    colnames(MergedData) <- gsub(" ","_",tolower(colnames(MergedData)))
    
    #Create database connection.
    con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)
    
    #Clear old metadata entries.
    dbExecute(con,paste('DELETE FROM "QPCRSample_test" WHERE "projectid" = \'',ProjectID,'\'',sep=""))
    
    #Check for redundant data.
    #Add new qPCR data.
    if(dbExistsTable(con,"QPCRSample_test")){
      Metadata_Check <-  tbl(con,"QPCRSample_test")
      Metadata_IDs <- MergedData$uniqueid
      Metadata_Check <- Metadata_Check %>% filter(uniqueid %in% Metadata_IDs)
      Metadata_Check <- as.data.frame(Metadata_Check)
      Metadata_Check_IDs <- Metadata_Check$uniqueid
      Metadata_Append <- MergedData[!(Metadata_IDs %in% Metadata_Check_IDs),]
      dbWriteTable(con,"QPCRSample_test",Metadata_Append,row.names=FALSE,append=TRUE)
    } else{
      dbWriteTable(con,"QPCRSample_test",MergedData,row.names=FALSE,append=TRUE)
    }
    
    #Prepare data for export to occurence table.
    ObservationsExport <- MergedData[MergedData$target_was_organism_detected==1,c("target_organism","target_taxonomic_rank_of_organism","projectid","sample_id","latitude","longitude","sample_date")]
    #Get sample years
    ObservationsExport$year <- as.numeric(format(ObservationsExport$sample_date,'%Y'))
    #Rename columns for export
    names(ObservationsExport)[names(ObservationsExport) == "target_organism"] <- "Taxon"
    names(ObservationsExport)[names(ObservationsExport) == "target_taxonomic_rank_of_organism"] <- "rank"
    names(ObservationsExport)[names(ObservationsExport) == "projectid"] <- "projectId"
    names(ObservationsExport)[names(ObservationsExport) == "sample_id"] <- "sampleId"
    ObservationsExport$sample_date <- NULL
    #Designate a unique id
    ObservationsExport$id <- sapply(paste(ObservationsExport$Taxon,ObservationsExport$rank,ObservationsExport$projectId,ObservationsExport$sampleId,ObservationsExport$latitude,ObservationsExport$longitude,ObservationsExport$year),digest,algo="md5")
    #Change case on taxonomic rank entries
    ObservationsExport$rank <- tolower(ObservationsExport$rank)
    
    #Check for redundant qPCR data.
    #Add new qPCR data.
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
