#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
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
gbif_dir <- Sys.getenv("GBIF_HOME")
bucket <- Sys.getenv("S3_BUCKET")

Taxon_name <- args[1]

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, lastRanAt = Sys.time(), error = error_message))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  dest_filename <- sub("\\.json$", ".build", filename)
  
  s3_path <- if (is.null(Taxon_name) || Taxon_name == "") {
    paste("s3://",bucket,"/errors/timelines/", new_filename, sep = "")
  } else {
    paste("s3://",bucket,"/timelines/", dest_filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}

# Get filtering parameters.
# Taxon_name:string Scientific taxon name
# TaxonomicRank:string Taxonomic level to aggregate results to
# Rscript --vanilla eDNAExplorer_Timeline_Metabarcoding.R "Taxon_name" "TaxonomicRank"

tryCatch(
  {
    if (length(args) != 2) {
      stop("Need the following inputs: Taxon_name, TaxonomicRank.", call. = FALSE)
    } else if (length(args) == 2) {
      Taxon_name <- args[1]
      TaxonomicRank <- args[2]
      #Define variables.
      Taxon_name <- as.character(Taxon_name)
      #Get GBIF taxonomy key for taxon.
      Taxon_GBIF <- name_backbone(name=Taxon_name,rank=TaxonomicRank)$usageKey
      
      #Define output filename.
      filename <- paste("Timeline_Metabarcoding_Taxon_",Taxon_name,"_Rank_",TaxonomicRank,".json",sep="")
      filename <- tolower(filename)
      filename <- gsub(" ","_",filename)
    }
  },
  error = function(e) {
    process_error(e)
  }
)

# Generate the output filename for cached plots.
tryCatch(
  {
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    dest_filename <- sub("\\.json$", ".build", filename) # Write to a temporary file first as .build
    system(paste("aws s3 cp ", filename, " s3://",bucket,"/timelines/", dest_filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""), intern = TRUE)
    system(paste("rm ", filename, sep = ""))
    
    #Establish sql connection
    Database_Driver <- dbDriver("PostgreSQL")
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
    
    #Read in GBIF occurrences.
    gbif <- gbif_local()
    
    #Filter GBIF occurrences to a particular taxon.
    GBIFDB <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                              coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                              occurrencestatus=="PRESENT",taxonkey==Taxon_GBIF) %>% select(year)
    GBIFDB <- as.data.frame(GBIFDB)
    
    if(nrow(GBIFDB) < 1){
      Taxa_Time <- data.frame(matrix(nrow=1,ncol=2))
      colnames(Taxa_Time) <- c("year","occurrences")
      Taxa_Time$year <- " "
      Taxa_Time$occurrences <- " "
    }
    if(nrow(GBIFDB) >=1){
      #Get GBIF taxon occurrences over time.
      Taxa_Time <- as.data.frame(table(GBIFDB))
      colnames(Taxa_Time) <- c("year","occurrences") 
    }
    Taxa_Time$year <- as.numeric(as.character(Taxa_Time$year))
    
    # Read in Tronko output and filter it.
    TronkoInput <- tbl(con, "Occurence")
    TronkoFiltered <- TronkoInput %>% filter(Taxon == Taxon_name) %>%
      filter(rank == TaxonomicRank) %>% select(c(year))
    TronkoDB <- as.data.frame(TronkoFiltered)
    TronkoDB <- as.data.frame(TronkoDB[complete.cases(TronkoDB),])
    colnames(TronkoDB) <- c("year")
    if(nrow(TronkoDB) < 1){
      Tronko_Time <- data.frame(matrix(nrow=1,ncol=2))
      colnames(Tronko_Time) <- c("year","occurrences")
      Tronko_Time$year <- " "
      Tronko_Time$occurrences <- " "
    }
    if(nrow(TronkoDB) >=1){
      #Get GBIF taxon occurrences over time.
      Tronko_Time <- as.data.frame(table(TronkoDB))
      colnames(Tronko_Time) <- c("year","occurrences") 
    }
    
    # Generate JSON object for export and mapping.
    datasets <- list(datasets = list(eDNA=Tronko_Time,GBIF=Taxa_Time))
    
    # Export file for mapping
    filename <- paste("Timeline_Metabarcoding_Taxon_",Taxon_name,"_Rank_",TaxonomicRank,".json",sep="")
    filename <- tolower(filename)
    filename <- gsub(" ","_",filename)
    write(toJSON(datasets,auto_unbox = TRUE),filename)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/timelines/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  },
  error = function(e) {
    process_error(e, filename)
  }
)