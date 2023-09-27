#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
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
require(uuid)

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
    paste("s3://ednaexplorer_staging/errors/venn/", new_filename, sep = "")
  } else {
    paste("s3://ednaexplorer_staging/projects/", ProjectID, "/plots/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
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
# First_Date:string YYYY-MM-DD
# Last_Date:string YYYY-MM-DD
# Marker:string Target marker name
# Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
# TaxonomicRank:string Taxonomic level to aggregate results to
# CountThreshold:numeric Read count threshold for retaining samples
# FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
# SpeciesList:string Name of csv file containing selected species list.
# Geographic_Scale:string Local, State, or Nation
# Rscript --vanilla ednaexplorer_staging_Venn_Metabarcoding.R "ProjectID" "First_Date" "Last_Date" "Marker" "Num_Mismatch" "TaxonomicRank" "CountThreshold" "FilterThreshold" "SpeciesList" "Geographic_Scale"

tryCatch(
  {
    if (length(args) != 10) {
      stop("Need the following inputs: ProjectID, First_Date, Last_Date, Marker, Num_Mismatch, TaxonomicRank, CountThreshold, FilterThreshold, SpeciesList, Geographic_Scale.", call. = FALSE)
    } else if (length(args) == 10) {
      ProjectID <- args[1]
      First_Date <- args[2]
      Last_Date <- args[3]
      Marker <- args[4]
      Num_Mismatch <- args[5]
      TaxonomicRank <- args[6]
      CountThreshold <- args[7]
      FilterThreshold <- args[8]
      SpeciesList <- args[9]
      Geographic_Scale <- args[10]
      
      CategoricalVariables <- c("site","grtgroup", "biome_type", "iucn_cat", "eco_name", "hybas_id")
      ContinuousVariables <- c("bio01", "bio12", "ghm", "elevation", "ndvi", "average_radiance")
      FieldVars <- c("fastqid", "sample_date", "latitude", "longitude", "spatial_uncertainty")
      TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
      TaxonomicKeyRanks <- c("kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")
      TaxonomicNum <- which(TaxonomicRanks == TaxonomicRank)
      Project_ID <- as.character(ProjectID)
      First_Date <- lubridate::ymd(First_Date)
      Last_Date <- lubridate::ymd(Last_Date)
      Marker <- as.character(Marker)
      Num_Mismatch <- as.numeric(Num_Mismatch)
      CountThreshold <- as.numeric(CountThreshold)
      FilterThreshold <- as.numeric(FilterThreshold)
      SelectedSpeciesList <- as.character(SpeciesList)
      
      filename <- paste("Venn_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,"GeographicScale",Geographic_Scale,".json",sep="_")
      filename <- gsub("_.json",".json",filename)
      filename <- tolower(filename)
    }
  },
  error = function(e) {
    process_error(e)
  }
)

tryCatch(
  {
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    system(paste("aws s3 cp ", filename, " s3://ednaexplorer_staging/projects/", Project_ID, "/plots/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""), intern = TRUE)
    system(paste("rm ", filename, sep = ""))
    
    #Establish sql connection
    Database_Driver <- dbDriver("PostgreSQL")
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
    
    #Read in species list
    if(SelectedSpeciesList != "None"){
      SpeciesList_df <- tbl(con,"SpeciesListItem")
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
    TronkoFile_tmp <- paste(Marker,"_venn_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://ednaexplorer_staging/tronko_output/", Project_ID, "/", TronkoFile, " ", TronkoFile_tmp, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""))
    #Check if file exists.
    if(file.info(TronkoFile_tmp)$size== 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- paste("subset_venn_",UUIDgenerate(),".csv",sep="")
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
    
    #Get eDNA results summarized.
    if(nrow(TronkoDB) >= 1){
      num_filteredSamples <- length(unique(TronkoDB$SampleID))
      #Get unique taxa list from Tronko-assign
      Tronko_Taxa <- na.omit(unique(TronkoDB[,TaxonomicRank]))
      Tronko_Taxa <- as.data.frame(Tronko_Taxa)
      colnames(Tronko_Taxa) <- c("eDNA")
    } else {
      Tronko_Taxa <- data.frame(matrix(ncol=1,nrow=1))
      colnames(Tronko_Taxa) <- c("eDNA")
      Tronko_Taxa$eDNA <- NA
      num_filteredSamples <- 0
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
      Taxa_GBIF <- as.data.frame(Taxa_Local)
      Taxa_GBIF <- na.omit(unique(Taxa_GBIF[,TaxonomicRank]))
      Taxa_GBIF <- as.data.frame(Taxa_GBIF)
      colnames(Taxa_GBIF) <- c("Taxa_Local")
      
      #Define results shared between eDNA and GBIF.
      Both <- as.data.frame(intersect(Tronko_Taxa$eDNA,Taxa_GBIF$Taxa_Local))
      colnames(Both) <- c("Both")
      #Define eDNA results as those unique from GBIF
      eDNA_only <- as.data.frame(Tronko_Taxa$eDNA[!(Tronko_Taxa$eDNA %in% Taxa_GBIF$Taxa_Local)])
      colnames(eDNA_only) <- c("eDNA")
      #Define GBIF results those unique from eDNA
      GBIF_only <- as.data.frame(Taxa_GBIF$Taxa_Local[!(Taxa_GBIF$Taxa_Local %in% Tronko_Taxa$eDNA)])
      colnames(GBIF_only) <- c("Taxa_Local")
       
      rm(Taxa_GBIF)
      Taxa_GBIF <- GBIF_only
      rm(Tronko_Taxa)
      Tronko_Taxa <- eDNA_only
    }
    if(Geographic_Scale=="State"){
      #Clip GBIF occurrence locations by state/province boundaries.
      Taxa_State <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                    coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                    occurrencestatus=="PRESENT", stateprovince %in% state_province_list) %>% select(!!sym(TaxonomicRank))
      Taxa_GBIF <- as.data.frame(Taxa_State)
      Taxa_GBIF <- na.omit(unique(Taxa_GBIF[,TaxonomicRank]))
      Taxa_GBIF <- as.data.frame(Taxa_GBIF)
      colnames(Taxa_GBIF) <- c("Taxa_State")
      
      #Define results shared between eDNA and GBIF.
      Both <- as.data.frame(intersect(Tronko_Taxa$eDNA,Taxa_GBIF$Taxa_State))
      colnames(Both) <- c("Both")
      #Define eDNA results as those unique from GBIF
      eDNA_only <- as.data.frame(Tronko_Taxa$eDNA[!(Tronko_Taxa$eDNA %in% Taxa_GBIF$Taxa_State)])
      colnames(eDNA_only) <- c("eDNA")
      #Define GBIF results those unique from eDNA
      GBIF_only <- as.data.frame(Taxa_GBIF$Taxa_State[!(Taxa_GBIF$Taxa_State %in% Tronko_Taxa$eDNA)])
      colnames(GBIF_only) <- c("Taxa_State")
      
      rm(Taxa_GBIF)
      Taxa_GBIF <- GBIF_only
      rm(Tronko_Taxa)
      Tronko_Taxa <- eDNA_only
    }
    if(Geographic_Scale=="Nation"){
      #Clip GBIF occurrence locations by national boundaries.
      Taxa_Nation <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                     coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                     occurrencestatus=="PRESENT", countrycode %in% country_list) %>% select(!!sym(TaxonomicRank))
      Taxa_GBIF <- as.data.frame(Taxa_Nation)
      Taxa_GBIF <- na.omit(unique(Taxa_GBIF[,TaxonomicRank]))
      Taxa_GBIF <- as.data.frame(Taxa_GBIF)
      colnames(Taxa_GBIF) <- c("Taxa_Nation")
      
      #Define results shared between eDNA and GBIF.
      Both <- as.data.frame(intersect(Tronko_Taxa$eDNA,Taxa_GBIF$Taxa_Nation))
      colnames(Both) <- c("Both")
      #Define eDNA results as those unique from GBIF
      eDNA_only <- as.data.frame(Tronko_Taxa$eDNA[!(Tronko_Taxa$eDNA %in% Taxa_GBIF$Taxa_Nation)])
      colnames(eDNA_only) <- c("eDNA")
      #Define GBIF results those unique from eDNA
      GBIF_only <- as.data.frame(Taxa_GBIF$Taxa_Nation[!(Taxa_GBIF$Taxa_Nation %in% Tronko_Taxa$eDNA)])
      colnames(GBIF_only) <- c("Taxa_Nation")
      
      rm(Taxa_GBIF)
      Taxa_GBIF <- GBIF_only
      rm(Tronko_Taxa)
      Tronko_Taxa <- eDNA_only
    }
    
    #Insert the number of samples and number of samples post-filtering as a return object.
    SampleDB <- data.frame(matrix(ncol=2,nrow=1))
    colnames(SampleDB) <- c("totalSamples","filteredSamples")
    SampleDB$totalSamples <- total_Samples
    SampleDB$filteredSamples <- num_filteredSamples
    datasets <- list(datasets = list(eDNA=Tronko_Taxa[,1],Both=Both[,1],GBIF=Taxa_GBIF[,1],metadata=SampleDB))
    filename <- paste("Venn_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,"GeographicScale",Geographic_Scale,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    write(toJSON(datasets),filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer_staging/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
  },
  error = function(e) {
    process_error(e, filename)
  }
)
