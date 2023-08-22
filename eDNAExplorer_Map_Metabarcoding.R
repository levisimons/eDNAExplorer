#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
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
    paste("s3://ednaexplorer/errors/map/", new_filename, sep = "")
  } else {
    paste("s3://ednaexplorer/projects/", ProjectID, "/plots/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  stop(error_message)
}

readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
gbif_dir <- Sys.getenv("GBIF_HOME")

# Get filtering parameters.
# ProjectID:string
# Marker:string Target marker name
# Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
# TaxonomicRank:string Taxonomic level to aggregate results to
# CountThreshold:numeric Read count threshold for retaining samples
# FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
# Taxon_name:string Scientific taxon name
# Rscript --vanilla eDNAExplorer_Map_Metabarcoding.R "ProjectID" "Marker" "Num_Mismatch" "TaxonomicRank" "CountThreshold" "FilterThreshold" "Taxon_name"

tryCatch(
  {
    if (length(args) != 7) {
      stop("Need the following inputs: ProjectID, Marker, Num_Mismatch, TaxonomicRank, CountThreshold, FilterThreshold, Taxon_name.", call. = FALSE)
    } else if (length(args) == 7) {
      ProjectID <- args[1]
      Marker <- args[2]
      Num_Mismatch <- args[3]
      TaxonomicRank <- args[4]
      CountThreshold <- args[5]
      FilterThreshold <- args[6]
      Taxon_name <- args[7]
    }
    Project_ID <- as.character(ProjectID)
    Taxon <- as.character(Taxon_name)
    #Get GBIF taxonomy key for taxon.
    Taxon_GBIF <- name_backbone(name=Taxon,rank=TaxonomicRank)$usageKey
    #Ensure numeric values.
    Num_Mismatch <- as.numeric(Num_Mismatch)
    CountThreshold <- as.numeric(CountThreshold)
    FilterThreshold <- as.numeric(FilterThreshold)
  },
  error = function(e) {
    process_error(e)
  }
)

# Generate the output filename for cached plots.
tryCatch(
  {
    #Establish sql connection
    Database_Driver <- dbDriver("PostgreSQL")
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    
    #Read in GBIF occurrences.
    gbif <- gbif_local(dir=gbif_dir)
    
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
    
    # Read in Tronko output and filter it.
    TronkoFile <- paste(Marker, ".csv", sep = "")
    system(paste("aws s3 cp s3://ednaexplorer/tronko_output/", Project_ID, "/", TronkoFile, " ", TronkoFile, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""))
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- "subset.csv"
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",TronkoFile, SubsetFile)
    system(awk_command, intern = TRUE)
    # Filter on the number of mismatches.  Remove entries with NA for mismatches and for the selected taxonomic rank.
    TronkoInput <- fread(file = "subset.csv", header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    TronkoInput$Mismatch <- as.numeric(as.character(TronkoInput$Mismatch))
    TronkoInput <- TronkoInput %>%
      filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>%
      filter(!is.na(!!sym(TaxonomicRank))) %>%
      group_by(SampleID) %>%
      filter(n() > CountThreshold) %>%
      select(SampleID, kingdom, phylum, class, order, family, genus, species)
    TronkoDB <- as.data.frame(TronkoInput)
    TronkoDB <- TronkoDB[TronkoDB[,TaxonomicRank]==Taxon,]
    
    #Get samples where taxon occurs and meets Tronko filters.
    TaxonDB <- TronkoDB[!duplicated(TronkoDB),]
    TaxonDB$SampleID <- gsub("-","_",TaxonDB$SampleID)
    system(paste("rm",TronkoFile,sep=" "))
    system("rm subset.csv")
    taxon_samples <- unique(TaxonDB$SampleID)
    taxon_projects <- ProjectID
    
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
  },
  error = function(e) {
    process_error(e, filename)
  }
)
