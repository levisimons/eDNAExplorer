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
require(zoo)
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
ENDPOINT_URL <- Sys.getenv("ENDPOINT_URL")

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
}

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
    paste("s3://",bucket,"/errors/spacetime/", new_filename," --endpoint-url ",ENDPOINT_URL, sep = "")
  } else {
    paste("s3://",bucket,"/projects/", ProjectID, "/plots/error_", dest_filename, " --endpoint-url ",ENDPOINT_URL, sep = "")
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
# Rscript --vanilla eDNAExplorer_Spacetime_Metabarcoding.R "ProjectID" "First_Date" "Last_Date" "Marker" "Num_Mismatch" "TaxonomicRank" "CountThreshold" "FilterThreshold" "SpeciesList"

# Generate the output filename for cached plots.
tryCatch(
  {
    #Generate the output filename for cached plots.
    filename <- paste("PresenceByTime_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    dest_filename <- sub("\\.json$", ".build", filename) # Write to a temporary file first as .build
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/plots/",dest_filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    filename <- paste("PresenceBySite_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    dest_filename <- sub("\\.json$", ".build", filename) # Write to a temporary file first as .build
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/plots/",dest_filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    filename <- paste("PresenceBySample_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    dest_filename <- sub("\\.json$", ".build", filename) # Write to a temporary file first as .build
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/plots/",dest_filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    filename <- paste("PresenceBySiteAndTime_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    dest_filename <- sub("\\.json$", ".build", filename) # Write to a temporary file first as .build
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/plots/",dest_filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank filtered taxonomy table as a default.  This gets overwritten is actual material exists.
    filename <- paste("FilteredTaxonomy_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".csv",sep="_")
    filename <- gsub("_.csv",".csv",filename)
    filename <- tolower(filename)
    write.table(data.frame(error=c("No results"),message=c("Filter too stringent")),filename,quote=FALSE,sep=",",row.names = FALSE)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",Project_ID,"/tables/",filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
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
    TronkoFile_tmp <- paste(Marker,"_spacetime_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://",bucket,"/tronko_output/", Project_ID, "/", TronkoFile, " ", TronkoFile_tmp, " --endpoint-url ",ENDPOINT_URL, sep = ""))
    #Check if file exists.
    if(file.info(TronkoFile_tmp)$size== 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- paste("subset_spacetime_",UUIDgenerate(),".csv",sep="")
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
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    
    #Determine the number of samples passing the filter.
    #Merge in taxonomic data.
    if(nrow(TronkoDB)>0){
      num_filteredSamples <- length(unique(TronkoDB$SampleID))
      #Merge in taxonomy data.
      TronkoDB <- dplyr::left_join(TronkoDB, TaxonomyDB[,c(TaxonomicRank,"Common_Name","Image_URL")],na_matches="never")
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
    }
    if(nrow(TronkoDB)==0){
      stop("Error: Filters are too stringent. Cannot proceed.")
    }
    
    #Merge Tronko output with sample metadata
    ProjectDB <- dplyr::left_join(TronkoDB,Metadata[,c(CategoricalVariables,ContinuousVariables,FieldVars,"sample_id")],by=c("SampleID"="fastqid"))
    print(paste("ProjectDB",nrow(ProjectDB),ncol(ProjectDB)))
    
    # Generate the number of samples and number of samples post-filtering as a return object,
    # to append to output json objects.
    SampleDB <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(SampleDB) <- c("totalSamples", "filteredSamples")
    SampleDB$totalSamples <- total_Samples
    SampleDB$filteredSamples <- num_filteredSamples
    
    #Aggregate merged data by the appropriate time interval to find taxa presence by time.
    ProjectDB_byTime <- ProjectDB %>% dplyr::distinct(site,SampleID,sample_date,Latin_Name,.keep_all=T) %>% 
      dplyr::group_by(sample_date,Latin_Name) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    
    #Export taxa presence by time.
    ProjectDB_byTime <- as.data.frame(ProjectDB_byTime)
    #Store date range information
    unique_dates <- unique(as.character(ProjectDB_byTime$sample_date))
    #Merge back in taxonomic information.
    tmp <- ProjectDB[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    ProjectDB_byTime <- dplyr::left_join(ProjectDB_byTime,tmp)
    # Convert to a wide data frame.
    ProjectDB_byTime <- tidyr::pivot_wider(ProjectDB_byTime[c("Latin_Name","Common_Name","Image_URL","sample_date","freq")], names_from = sample_date, values_from = freq, values_fill = 0)
    ProjectDB_byTime <- as.data.frame(ProjectDB_byTime)
    #Specify filname
    filename <- paste("PresenceByTime_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    # Convert data set to JSON object.
    json_list <- list()
    # Iterate through rows of the data frame
    df <- ProjectDB_byTime
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      # Create a list for the date ranges and Presence values
      sample_data <- list()
      for(unique_date in unique_dates){
        freq <- df[i, unique_date]
        sample_data[[unique_date]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name = Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        Sample_Date = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",Project_ID,"/plots/",filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Aggregate merged data by site to find taxa presence by site.
    ProjectDB_bySite <- ProjectDB %>% dplyr::distinct(site,SampleID,Latin_Name,.keep_all=T) %>% 
      dplyr::group_by(site,Latin_Name) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    #Export taxa presence by site.
    ProjectDB_bySite <- as.data.frame(ProjectDB_bySite)
    # Get unique sites.
    unique_sites <- unique(ProjectDB_bySite$site)
    #Merge back in taxonomic information.
    tmp <- ProjectDB[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    ProjectDB_bySite <- dplyr::left_join(ProjectDB_bySite,tmp)
    # Convert to a wide data frame.
    ProjectDB_bySite <- tidyr::pivot_wider(ProjectDB_bySite, names_from = site, values_from = freq, values_fill = 0)
    ProjectDB_bySite <- as.data.frame(ProjectDB_bySite)
    # Convert data set to JSON object.
    json_list <- list()
    # Iterate through rows of the data frame
    df <- ProjectDB_bySite
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      # Create a list for the Site and Presence values
      sample_data <- list()
      for(unique_site in unique_sites){
        freq <- df[i, unique_site]
        sample_data[[unique_site]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name= Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        Site = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    filename <- paste("PresenceBySite_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",Project_ID,"/plots/",filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Find taxa presence/absence by sample.
    ProjectDB_bySample <- ProjectDB[,c("Latin_Name","sample_id")]
    ProjectDB_bySample <- ProjectDB_bySample[!duplicated(ProjectDB_bySample),]
    ProjectDB_bySample <- ProjectDB_bySample[!is.na(ProjectDB_bySample[,"Latin_Name"]),]
    colnames(ProjectDB_bySample) <- c("taxa","SampleID")
    ProjectDB_bySample$Presence <- 1
    #Get sample names.
    unique_samples <- unique(ProjectDB$sample_id)
    # Create a reference data frame with all possible combinations of taxa and samples.
    all_combinations <- expand.grid(
      taxa = na.omit(unique(ProjectDB_bySample$taxa)),
      SampleID = unique(ProjectDB_bySample$SampleID)
    )
    # Merge the original data frame with the reference data frame
    ProjectDB_bySample <- merge(all_combinations, ProjectDB_bySample, by = c("taxa", "SampleID"), all.x = TRUE)
    #Convert NA values to 0 in Presence.
    ProjectDB_bySample["Presence"][is.na(ProjectDB_bySample["Presence"])] <- 0
    #Merge back in taxonomic information.
    tmp <- ProjectDB[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    ProjectDB_bySample <- dplyr::left_join(ProjectDB_bySample,tmp,by=c("taxa"="Latin_Name"))
    # Convert to a presence/absence data frame.
    ProjectDB_bySample <- tidyr::pivot_wider(ProjectDB_bySample, names_from = SampleID, values_from = Presence, values_fill = 0)
    ProjectDB_bySample <- as.data.frame(ProjectDB_bySample)
    #Re-insert taxonomic rank.
    colnames(ProjectDB_bySample)[which(names(ProjectDB_bySample) == "taxa")] <- "Latin_Name"
    #Specify filename
    filename <- paste("PresenceBySample_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    # Convert data set to JSON object.
    json_list <- list()
    # Iterate through rows of the data frame
    df <- ProjectDB_bySample
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      # Create a list for the SampleID and Presence values
      sample_data <- list()
      for(unique_sample in unique_samples){
        presence <- df[i, unique_sample]
        sample_data[[unique_sample]] <- presence
      }
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name = Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        SampleID = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",Project_ID,"/plots/",filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Aggregate merged data by the appropriate time interval to find taxa presence by time.
    ProjectDB_bySiteTime <- ProjectDB %>% dplyr::distinct(site,SampleID,sample_date,Latin_Name,.keep_all=T) %>% 
      dplyr::group_by(site,sample_date,Latin_Name) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    
    #Export taxa presence by time.
    ProjectDB_bySiteTime <- as.data.frame(ProjectDB_bySiteTime)
    # Get unique sites.
    unique_sites <- unique(ProjectDB_bySiteTime$site)
    # Get unique date ranges.
    unique_dates <- unique(as.character(ProjectDB_bySiteTime$sample_date))
    #Merge back in taxonomic information.
    tmp <- ProjectDB[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    ProjectDB_bySiteTime <- dplyr::left_join(ProjectDB_bySiteTime,tmp)
    # Convert to a wide data frame.
    ProjectDB_bySiteTime <- tidyr::pivot_wider(ProjectDB_bySiteTime[,c("site","Latin_Name","freq","sample_date","Common_Name","Image_URL")], names_from = site, values_from = freq, values_fill = 0)
    ProjectDB_bySiteTime <- as.data.frame(ProjectDB_bySiteTime)
    filename <- paste("PresenceBySiteAndTime_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    # Convert data set to JSON object.
    json_list <- list()
    # Iterate through rows of the data frame
    df <- ProjectDB_bySiteTime
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      Sample_Date <- df[i,"sample_date"]
      # Create a list for the site/date/frequency values
      sample_data <- list()
      for(unique_site in unique_sites){
        freq <- df[i, unique_site]
        sample_data[[unique_site]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name = Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        sample_date = Sample_Date,
        Site = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",Project_ID,"/plots/",filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Export filtered taxonomy table.
    TronkoTable <- ProjectDB[,c("SampleID","Latin_Name")]
    TronkoTable <- as.data.frame(table(TronkoTable[,c("SampleID","Latin_Name")]))
    TronkoTable <- as.data.frame(pivot_wider(TronkoTable, names_from = SampleID, values_from = Freq))
    filename <- paste("FilteredTaxonomy_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".csv",sep="_")
    filename <- gsub("_.csv",".csv",filename)
    filename <- tolower(filename)
    write.table(TronkoTable,filename,quote=FALSE,sep=",",row.names = FALSE)
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",Project_ID,"/tables/",filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
  },
  error = function(e) {
    process_error(e, filename)
  }
)
