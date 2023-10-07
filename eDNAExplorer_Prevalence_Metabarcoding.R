#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(tidyr)
require(dplyr)
require(lubridate)
require(jsonlite)
require(data.table)
require(DBI)
require(RPostgreSQL)
require(digest)
require(uuid)

# Fetch project ID early so we can use it for error output when possible.
ProjectID <- args[1]

# Establish database credentials.
readRenviron(".env")
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
  "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY")
)
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
bucket <- Sys.getenv("S3_BUCKET")
home_dir <- Sys.getenv("home_dir")

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
    paste("s3://",bucket,"/errors/prevalence/", new_filename, sep = "")
  } else {
    paste("s3://",bucket,"/projects/", ProjectID, "/plots/", dest_filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}
# Get filtering parameters.
# ProjectID:string
# First_Date:string YYYY-MM-DD
# Last_Date:string YYYY-MM-DD
# Marker:string Target marker name
# Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
# TaxonomicRank:string Taxonomic level to aggregate results to
# CountThreshold:numeric Read count threshold for retaining samples
# FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
# SpeciesList:string Name of csv file containing selected species list.
# Rscript --vanilla eDNAExplorer_Prevalence_Metabarcoding.R "ProjectID" "First_Date" "Last_Date" "Marker" "Num_Mismatch" "TaxonomicRank" "CountThreshold" "FilterThreshold" "SpeciesList"

tryCatch(
  {
    if (length(args) != 9) {
      stop("Need the following inputs: ProjectID, First_Date, Last_Date, Marker, Num_Mismatch, TaxonomicRank, CountThreshold, FilterThreshold, SpeciesList.", call. = FALSE)
    } else if (length(args) == 9) {
      ProjectID <- args[1]
      First_Date <- args[2]
      Last_Date <- args[3]
      Marker <- args[4]
      Num_Mismatch <- args[5]
      TaxonomicRank <- args[6]
      CountThreshold <- args[7]
      FilterThreshold <- args[8]
      SpeciesList <- args[9]
      
      CategoricalVariables <- c("site","grtgroup", "biome_type", "iucn_cat", "eco_name", "hybas_id")
      ContinuousVariables <- c("bio01", "bio12", "ghm", "elevation", "ndvi", "average_radiance")
      FieldVars <- c("fastqid", "sample_date", "latitude", "longitude", "spatial_uncertainty")
      TaxonomicRanks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
      TaxonomicKeyRanks <- c("kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")
      TaxonomicNum <- which(TaxonomicRanks == TaxonomicRank)
      First_Date <- lubridate::ymd(First_Date)
      Last_Date <- lubridate::ymd(Last_Date)
      Num_Mismatch <- as.numeric(Num_Mismatch)
      CountThreshold <- as.numeric(CountThreshold)
      FilterThreshold <- as.numeric(FilterThreshold)
      SelectedSpeciesList <- as.character(SpeciesList)
      Project_ID <- as.character(ProjectID)
      
      filename <- paste("Prevalence_Metabarcoding_FirstDate", First_Date, "LastDate", Last_Date, "Marker", Marker, "Rank", TaxonomicRank, "Mismatch", Num_Mismatch, "CountThreshold", CountThreshold, "AbundanceThreshold", format(FilterThreshold, scientific = F), "SpeciesList", SelectedSpeciesList, sep = "_")
      filename <- paste(filename, ".json", sep = "")
      filename <- tolower(filename)
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
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/plots/",dest_filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Establish sql connection
    Database_Driver <- dbDriver("PostgreSQL")
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    con <- dbConnect(Database_Driver, host = db_host, port = db_port, dbname = db_name, user = db_user, password = db_pass)
    
    # Read in species list
    if (SelectedSpeciesList != "None") {
      SpeciesList_df <- tbl(con, "SpeciesListItem")
      SpeciesList_df <- SpeciesList_df %>% filter(species_list == SelectedSpeciesList)
      SpeciesList_df <- as.data.frame(SpeciesList_df)
    }
    
    # Read in metadata and filter it.
    Metadata <- tbl(con, "TronkoMetadata")
    Keep_Vars <- c(CategoricalVariables, ContinuousVariables, FieldVars)[c(CategoricalVariables, ContinuousVariables, FieldVars) %in% dbListFields(con, "TronkoMetadata")]
    # Get the number of samples in a project before filtering.
    Metadata_Unfiltered <- Metadata %>% filter(projectid == Project_ID)
    Metadata_Unfiltered <- as.data.frame(Metadata_Unfiltered)
    total_Samples <- nrow(Metadata_Unfiltered)
    Metadata <- Metadata %>%
      filter(projectid == Project_ID) %>%
      filter(!is.na(latitude) & !is.na(longitude))
    Metadata <- as.data.frame(Metadata)
    Metadata$sample_date <- lubridate::ymd(Metadata$sample_date)
    Metadata <- Metadata %>% filter(sample_date >= First_Date & sample_date <= Last_Date)
    if(nrow(Metadata) == 0 || ncol(Metadata) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    Metadata$fastqid <- gsub("_", "-", Metadata$fastqid)
    
    # Read in Tronko output and filter it.
    TronkoFile <- paste(Marker, ".csv", sep = "")
    TronkoFile_tmp <- paste(Marker,"_prevalence_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://",bucket,"/tronko_output/", Project_ID, "/", TronkoFile, " ", TronkoFile_tmp, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""))
    #Check if file exists.
    if(file.info(TronkoFile_tmp)$size== 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- paste("subset_prevalence_",UUIDgenerate(),".csv",sep="")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",TronkoFile_tmp, SubsetFile)
    system(awk_command, intern = TRUE)
    TronkoInput <- fread(file=SubsetFile, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    TronkoInput$Mismatch <- as.numeric(as.character(TronkoInput$Mismatch))
    #Remove samples with missing coordinates, and which are outside of the date filters.
    TronkoInput <- TronkoInput <- TronkoInput[TronkoInput$SampleID %in% unique(na.omit(Metadata$fastqid)), ]
    #Store the unfiltered reads.
    Tronko_Unfiltered <- TronkoInput
    # Calculate relative abundance of taxa with a given rank in the unfiltered reads.
    Tronko_Unfiltered <- Tronko_Unfiltered %>%
      dplyr::group_by(SampleID, !!sym(TaxonomicRank)) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>%
      dplyr::ungroup() %>% select(-n)
    #Filter on the number of reads per sample, then a mismatch threshold.
    TronkoInput <- TronkoInput %>% group_by(SampleID) %>%
      filter(n() > CountThreshold) %>% filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>% select(-Mismatch)
    #Merege relative abudance results into data filtered by reads per sample and mismatches.
    #Then filtered on relative abundances.
    TronkoInput <- TronkoInput %>%
      left_join(Tronko_Unfiltered,na_matches="never") %>%
      filter(freq >= FilterThreshold) %>% select(-freq)
    TronkoDB <- as.data.frame(TronkoInput)
    #Remove taxa which are unknown at a given rank.
    TronkoDB <- TronkoDB[,c("SampleID",TaxonomicRanks[1:TaxonomicNum])]
    TronkoDB <- TronkoDB[!is.na(TronkoDB[, TaxonomicRank]), ]
    #Filter results by species list.
    if (SelectedSpeciesList != "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)) & TronkoDB$species %in% SpeciesList_df$name, ]
    }
    if (SelectedSpeciesList == "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)), ]
    }
    system(paste("rm ",TronkoFile_tmp,sep=""))
    system(paste("rm ",SubsetFile,sep=""))
    
    # Read in Taxonomy output and filter it.
    TaxonomyInput <- tbl(con, "Taxonomy")
    if (nrow(TronkoDB) > 0) {
      TaxaList <- na.omit(unique(TronkoDB[, TaxonomicRank]))
    }
    if (nrow(TronkoDB) == 0) {
      TaxaList <- c()
    }
    if (TaxonomicRank != "kingdom") {
      TaxonomyInput <- TaxonomyInput %>%
        filter(rank==TaxonomicRank) %>%
        filter(!!sym(TaxonomicRank) %in% TaxaList) %>%
        select(TaxonomicRanks[2:TaxonomicNum],TaxonomicKeyRanks,Common_Name,Image_URL)
      TaxonomyDB <- as.data.frame(TaxonomyInput)
      #Figure out which taxonomy version is more complete.
      TaxonomyDB$rankCount <- rowSums(!is.na(TaxonomyDB[,colnames(TaxonomyDB) %in% TaxonomicRanks]))
      TaxonomyDB$rankKeyCount <- rowSums(!is.na(TaxonomyDB[,colnames(TaxonomyDB) %in% TaxonomicKeyRanks]))
      TaxonomyDB <- TaxonomyDB %>%
        group_by(!!sym(TaxonomicRank)) %>%
        slice_max(order_by = rankCount, n = 1) %>%
        ungroup()
      TaxonomyDB <- TaxonomyDB %>%
        group_by(!!sym(TaxonomicRank)) %>%
        slice_max(order_by = rankKeyCount, n = 1) %>%
        ungroup()
      #Figure out which common_name is most common per taxon.
      TaxonomyDB <- TaxonomyDB %>%
        group_by(!!sym(TaxonomicRank)) %>%
        mutate(Most_Common_Name = ifelse(all(is.na(Common_Name)), NA, names(which.max(table(Common_Name[!is.na(Common_Name)]))))) %>%
        ungroup()
      TaxonomyDB$Common_Name <- TaxonomyDB$Most_Common_Name
      TaxonomyDB$Most_Common_Name <- NULL
      TaxonomyDB <- as.data.frame(TaxonomyDB)
      TaxonomyDB <- subset(TaxonomyDB, select = -grep("Key", colnames(TaxonomyDB)))
      TaxonomyDB$rankCount <- NULL
      TaxonomyDB <- TaxonomyDB[!duplicated(TaxonomyDB),]
    }
    if (TaxonomicRank == "kingdom") {
      TaxonomyDB <- data.frame(
        kingdom = c("Fungi", "Plantae", "Animalia", "Bacteria", "Archaea", "Protista", "Monera", "Chromista"),
        Image_URL = c(
          "https://images.phylopic.org/images/7ebbf05d-2084-4204-ad4c-2c0d6cbcdde1/raster/958x1536.png",
          "https://images.phylopic.org/images/573bc422-3b14-4ac7-9df0-27d7814c099d/raster/1052x1536.png",
          "https://images.phylopic.org/images/0313dc90-c1e2-467e-aacf-0f7508c92940/raster/681x1536.png",
          "https://images.phylopic.org/images/d8c9f603-8930-4973-9a37-e9d0bc913a6b/raster/1536x1128.png",
          "https://images.phylopic.org/images/7ccfe198-154b-4a2f-a7bf-60390cfe6135/raster/1177x1536.png",
          "https://images.phylopic.org/images/4641171f-e9a6-4696-bdda-e29bc4508538/raster/336x1536.png",
          "https://images.phylopic.org/images/018ee72f-fde6-4bc3-9b2e-087d060ee62d/raster/872x872.png",
          "https://images.phylopic.org/images/1fd55f6f-553c-4838-94b4-259c16f90c31/raster/1054x1536.png"
        )
      )
    }
    
    #Calculate prevalence of taxa per sample and output results for plotting.
    if(nrow(TronkoDB) > 0){
      num_filteredSamples <- length(unique(TronkoDB$SampleID))
      TronkoDB <- TronkoDB %>%
        dplyr::group_by(!!sym(TaxonomicRank)) %>%
        dplyr::summarise(per = n_distinct(SampleID)/n_distinct(TronkoDB$SampleID))
      TronkoDB <- as.data.frame(TronkoDB)
      #Merge in taxonomy data.
      TronkoDB <- dplyr::left_join(TronkoDB, TaxonomyDB,na_matches="never")
      #TronkoDB <- dplyr::left_join(TronkoDB,TaxonomyDB,by=c("
      TronkoDB$Image_URL <- ifelse(is.na(TronkoDB$Image_URL), 'https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png', TronkoDB$Image_URL)
      if (TaxonomicRank != "kingdom") {
        colnames(TronkoDB)[which(names(TronkoDB) == TaxonomicRank)] <- "Latin_Name"
      }
      if (TaxonomicRank == "kingdom") {
        TronkoDB$Latin_Name <- TronkoDB$kingdom
      }
      TronkoDB <- TronkoDB %>% mutate_all(~ifelse(is.na(.), " ", .))
      #De-duplicated dataframe for tronko table.
      TronkoDB <- TronkoDB[!duplicated(TronkoDB),]
      
      # Insert the number of samples and number of samples post-filtering as a return object.
      SampleDB <- data.frame(matrix(ncol = 2, nrow = 1))
      colnames(SampleDB) <- c("totalSamples", "filteredSamples")
      SampleDB$totalSamples <- total_Samples
      SampleDB$filteredSamples <- num_filteredSamples
      
      datasets <- list(datasets = list(results = TronkoDB, metadata = SampleDB))
      write(toJSON(datasets), filename)
      system(paste("aws s3 cp ", filename, " s3://",bucket,"/projects/", Project_ID, "/plots/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""), intern = TRUE)
      system(paste("rm ", filename, sep = ""))
    }
    if(nrow(TronkoDB) == 0){
      stop("Error: Filters are too stringent. Cannot proceed.")
    }
    sapply(dbListConnections(Database_Driver), dbDisconnect)
  },
  error = function(e) {
    process_error(e)
  }
)
