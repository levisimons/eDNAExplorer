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

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, error = error_message))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  
  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://ednaexplorer/errors/spacetime/", new_filename, sep = "")
  } else {
    paste("s3://ednaexplorer/projects/", ProjectID, "/plots/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}

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
      First_Date <- lubridate::ymd(First_Date)
      Last_Date <- lubridate::ymd(Last_Date)
      Num_Mismatch <- as.numeric(Num_Mismatch)
      CountThreshold <- as.numeric(CountThreshold)
      FilterThreshold <- as.numeric(FilterThreshold)
      SelectedSpeciesList <- as.character(SpeciesList)
      Project_ID <- as.character(ProjectID)
    }
  },
  error = function(e) {
    process_error(e)
  }
)

# Generate the output filename for cached plots.
tryCatch(
  {
    #Generate the output filename for cached plots.
    filename <- paste("PresenceByTime_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    filename <- paste("PresenceBySite_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    filename <- paste("PresenceBySample_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    filename <- paste("PresenceBySiteAndTime_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",ProjectID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))

    # Output a blank filtered taxonomy table as a default.  This gets overwritten is actual material exists.
    filename <- paste("FilteredTaxonomy_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".csv",sep="_")
    filename <- gsub("_.csv",".csv",filename)
    filename <- tolower(filename)
    write.table(data.frame(error=c("No results"),message=c("Filter too stringent")),filename,quote=FALSE,sep=",",row.names = FALSE)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/tables/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
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
    tmp <- Metadata %>% filter(projectid == Project_ID)
    tmp <- as.data.frame(tmp)
    total_Samples <- nrow(tmp)
    Metadata <- Metadata %>%
      filter(projectid == Project_ID) %>%
      filter(!is.na(latitude) & !is.na(longitude))
    Metadata <- as.data.frame(Metadata)
    Metadata$sample_date <- lubridate::ymd(Metadata$sample_date)
    Metadata <- Metadata %>% filter(sample_date >= First_Date & sample_date <= Last_Date)
    Metadata$fastqid <- gsub("_", "-", Metadata$fastqid)
    if(nrow(Metadata) == 0 || ncol(Metadata) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    print(paste("Metadata",nrow(Metadata),ncol(Metadata)))
    
    # Read in Tronko output and filter it.
    TronkoFile <- paste(Marker, ".csv", sep = "")
    TronkoFile_tmp <- paste(Marker,"_spacetime_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://ednaexplorer/tronko_output/", Project_ID, "/", TronkoFile, " ", TronkoFile_tmp, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""))
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- paste("subset_spacetime_",UUIDgenerate(),".csv",sep="")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",TronkoFile_tmp, SubsetFile)
    system(awk_command, intern = TRUE)
    # Filter on the number of mismatches.  Remove entries with NA for mismatches and for the selected taxonomic rank.
    TronkoInput <- fread(file=SubsetFile, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    TronkoInput$Mismatch <- as.numeric(as.character(TronkoInput$Mismatch))
    TronkoInput <- TronkoInput %>%
      filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>%
      filter(!is.na(!!sym(TaxonomicRank))) %>%
      group_by(SampleID) %>%
      filter(n() > CountThreshold) %>%
      select(SampleID, kingdom, phylum, class, order, family, genus, species)
    TronkoDB <- as.data.frame(TronkoInput)
    print(paste("TronkoDB",nrow(TronkoDB),ncol(TronkoDB)))
    if (SelectedSpeciesList != "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)) & TronkoDB$species %in% SpeciesList_df$name, ]
    }
    if (SelectedSpeciesList == "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)), ]
    }
    system(paste("rm ",TronkoFile_tmp,sep=""))
    system(paste("rm ",SubsetFile,sep=""))
    print(paste("TronkoDB",nrow(TronkoDB),ncol(TronkoDB)))
    if(nrow(TronkoDB) == 0 || ncol(TronkoDB) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    print(paste("TronkoDB",nrow(TronkoDB),ncol(TronkoDB)))
    
    #Read in Taxonomy output and filter it.
    TaxonomyInput <- tbl(con,"Taxonomy")
    TaxaList <- na.omit(unique(TronkoDB[,TaxonomicRank]))
    TaxonomyInput <- TaxonomyInput %>% filter(rank == TaxonomicRank) %>% filter(Taxon %in% TaxaList) %>% select(Taxon,Common_Name,Image_URL)
    TaxonomyDB <-  as.data.frame(TaxonomyInput)
    colnames(TaxonomyDB) <- c(TaxonomicRank,"Common_Name","Image_URL")
    if(TaxonomicRank=="kingdom"){
      TaxonomyDB <- data.frame(kingdom=c("Fungi","Plantae","Animalia","Bacteria","Archaea","Protista","Monera","Chromista"),
                               Image_URL=c("https://images.phylopic.org/images/e5d32221-7ea9-46ed-8e0a-d9dbddab0b4a/raster/1536x1082.png",
                                           "https://images.phylopic.org/images/e42e270b-13ca-4871-b706-003ed1b98b8a/raster/1220x1705.png",
                                           "https://images.phylopic.org/images/9c234021-ce53-45d9-8fdd-b0ca3115a451/raster/1913x1113.png",
                                           "https://images.phylopic.org/images/891bdfc2-292b-422d-b7e5-a93823b9ce85/raster/1484x1536.png",
                                           "https://images.phylopic.org/images/a333c62f-a421-4ae0-894c-cf4dd8e1b70c/raster/1336x1536.png",
                                           "https://images.phylopic.org/images/595359bd-6638-4900-8536-82b2ed511a25/raster/535x1536.png",
                                           "https://images.phylopic.org/images/1edff864-c53b-492b-a0cf-f6fe816815a8/raster/1536x1282.png",
                                           "https://images.phylopic.org/images/4924b6bd-cfb8-4d60-a32a-442d02afbe85/raster/1460x1536.png"))
    }
    TaxonomyDB <- TaxonomyDB[!duplicated(TaxonomyDB),]
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    
    #Merge Tronko output with sample metadata
    ProjectDB <- dplyr::left_join(TronkoDB,Metadata,by=c("SampleID"="fastqid"))
    print(paste("ProjectDB",nrow(ProjectDB),ncol(ProjectDB)))
    #Get project duration
    Duration <- abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="days"))
    #Get day of year
    ProjectDB$day <- lubridate::yday(ProjectDB$sample_date)
    #Get week of year
    ProjectDB$week <- strftime(ProjectDB$sample_date,format="%V")
    #Get month of year
    ProjectDB$month <- lubridate::month(ProjectDB$sample_date)
    #Get quarter of year
    ProjectDB$quarter <- quarters(as.Date(ProjectDB$sample_date))
    #Get year
    ProjectDB$year <- lubridate::year(ProjectDB$sample_date)
    
    #Filter by relative abundance per taxon per sample.
    if(nrow(TronkoDB)>0){
      TronkoDB <- TronkoDB[!is.na(ProjectDB[,TaxonomicRank]),]
      TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n/sum(n)) %>% 
        dplyr::ungroup() %>% dplyr::filter(freq > FilterThreshold) %>% select(-n,-freq)
      TronkoDB <- as.data.frame(TronkoDB) 
      num_filteredSamples <- length(unique(TronkoDB$SampleID))
    }
    if(nrow(TronkoDB)==0){
      num_filteredSamples <- 0
    }
    
    # Generate the number of samples and number of samples post-filtering as a return object,
    # to append to output json objects.
    SampleDB <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(SampleDB) <- c("totalSamples", "filteredSamples")
    SampleDB$totalSamples <- total_Samples
    SampleDB$filteredSamples <- num_filteredSamples
    
    #Filter merged data.
    ProjectDB <- dplyr::inner_join(TronkoDB,ProjectDB,multiple="all")
    
    #Get start time of merged data.
    StartTime <- as.Date(paste(year(min(ProjectDB$sample_date)), 1, 1, sep = "-"))
    
    #Aggregate merged data by the appropriate time interval to find taxa presence by time.
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="days"))) <= 7){
      ProjectDB_byTime <- ProjectDB %>% dplyr::distinct(site,SampleID,day,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(day,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_byTime$date_range <- paste(StartTime+days(ProjectDB_byTime$day)-days(1),"to",StartTime+days(ProjectDB_byTime$day))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 1 &
       as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) <= 4){
      ProjectDB_byTime <- ProjectDB %>% dplyr::distinct(site,SampleID,week,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(week,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_byTime$date_range <- paste(StartTime+weeks(ProjectDB_byTime$week)-weeks(1),"to",StartTime+weeks(ProjectDB_byTime$week))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 4 &
       as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) <= 13){
      ProjectDB_byTime <- ProjectDB %>% dplyr::distinct(site,SampleID,month,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(month,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_byTime$date_range <- paste(StartTime+months(ProjectDB_byTime$month)-months(1),"to",StartTime+months(ProjectDB_byTime$month))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 13 &
       as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) <= 52){
      ProjectDB_byTime <- ProjectDB %>% dplyr::distinct(site,SampleID,quarter,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(quarter,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_byTime$date_range <- paste(as.Date(as.yearqtr(paste(ProjectDB_byTime$quarter,year(StartTime)), format = "Q%q %Y")),"to",as.Date(as.yearqtr(paste(ProjectDB_byTime$quarter,year(StartTime)), format = "Q%q %Y"))+months(3))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 52){
      ProjectDB_byTime <- ProjectDB %>% dplyr::distinct(site,SampleID,year,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(year,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_byTime$date_range <- paste(ProjectDB_byTime$year,"to",ProjectDB_byTime$year+1)
    }
    #Export taxa presence by time.
    ProjectDB_byTime <- as.data.frame(ProjectDB_byTime)
    #Merge in taxonomy data.
    ProjectDB_byTime <- dplyr::left_join(ProjectDB_byTime,TaxonomyDB,relationship = "many-to-many")
    if (TaxonomicRank != "kingdom") {
      colnames(ProjectDB_byTime)[which(names(ProjectDB_byTime) == TaxonomicRank)] <- "Latin_Name"
    }
    if (TaxonomicRank == "kingdom") {
      ProjectDB_byTime$Latin_Name <- ProjectDB_byTime$kingdom
    }
    #Store date range information
    unique_dates <- unique(ProjectDB_byTime$date_range)
    # Convert to a wide data frame.
    ProjectDB_byTime <- tidyr::pivot_wider(ProjectDB_byTime[c("Latin_Name","Common_Name","Image_URL","date_range","freq")], names_from = date_range, values_from = freq, values_fill = 0)
    ProjectDB_byTime <- as.data.frame(ProjectDB_byTime)
    #Convert NA values for image URLs to the default Phylopic image.
    ProjectDB_byTime["Image_URL"][is.na(ProjectDB_byTime["Image_URL"])] <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
    #Convert missing common names to blanks.
    ProjectDB_byTime["Common_Name"][is.na(ProjectDB_byTime["Common_Name"])] <- " "
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
      # Create a list for the SampleID and Presence values
      sample_data <- list()
      for(unique_date in unique_dates){
        freq <- df[i, unique_date]
        sample_data[[unique_date]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        Latin_Name = latin_name,
        Common_Name = common_name,
        Image_URL = image_url,
        Date_Range = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Aggregate merged data by site to find taxa presence by site.
    ProjectDB_bySite <- ProjectDB %>% dplyr::distinct(site,SampleID,!!sym(TaxonomicRank),.keep_all=T) %>% 
      dplyr::group_by(site,!!sym(TaxonomicRank)) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    #Export taxa presence by site.
    ProjectDB_bySite <- as.data.frame(ProjectDB_bySite)
    #Merge in taxonomy data.
    ProjectDB_bySite <- dplyr::left_join(ProjectDB_bySite,TaxonomyDB,relationship = "many-to-many")
    if (TaxonomicRank != "kingdom") {
      colnames(ProjectDB_bySite)[which(names(ProjectDB_bySite) == TaxonomicRank)] <- "Latin_Name"
    }
    if (TaxonomicRank == "kingdom") {
      ProjectDB_bySite$Latin_Name <- ProjectDB_bySite$kingdom
    }
    #Convert NA values for image URLs to the default Phylopic image.
    ProjectDB_bySite["Image_URL"][is.na(ProjectDB_bySite["Image_URL"])] <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
    #Convert missing common names to blanks.
    ProjectDB_bySite["Common_Name"][is.na(ProjectDB_bySite["Common_Name"])] <- " "
    # Get unique sites.
    unique_sites <- unique(ProjectDBby_Site$site)
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
      # Create a list for the SampleID and Presence values
      sample_data <- list()
      for(unique_site in unique_sites){
        freq <- df[i, unique_site]
        sample_data[[unique_site]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        Latin_Name = latin_name,
        Common_Name = common_name,
        Image_URL = image_url,
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
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Find taxa presence/absence by sample.
    ProjectDB_bySample <- ProjectDB[,c(TaxonomicRank,"SampleID")]
    ProjectDB_bySample <- ProjectDB_bySample[!duplicated(ProjectDB_bySample),]
    ProjectDB_bySample <- ProjectDB_bySample[!is.na(ProjectDB_bySample[,TaxonomicRank]),]
    colnames(ProjectDB_bySample) <- c("taxa","SampleID")
    ProjectDB_bySample$Presence <- 1
    # Create a reference data frame with all possible combinations of taxa and samples.
    all_combinations <- expand.grid(
      taxa = na.omit(unique(ProjectDB_bySample$taxa)),
      SampleID = unique(ProjectDB_bySample$SampleID)
    )
    # Merge the original data frame with the reference data frame
    ProjectDB_bySample <- merge(all_combinations, ProjectDB_bySample, by = c("taxa", "SampleID"), all.x = TRUE)
    #Convert NA values to 0 in Presence.
    ProjectDB_bySample["Presence"][is.na(ProjectDB_bySample["Presence"])] <- 0
    # Convert to a presence/absence data frame.
    ProjectDB_bySample <- tidyr::pivot_wider(ProjectDB_bySample, names_from = SampleID, values_from = Presence, values_fill = 0)
    ProjectDB_bySample <- as.data.frame(ProjectDB_bySample)
    #Re-insert taxonomic rank.
    colnames(ProjectDB_bySample)[which(names(ProjectDB_bySample) == "taxa")] <- TaxonomicRank
    #Merge in taxonomy data.
    ProjectDB_bySample <- dplyr::left_join(ProjectDB_bySample,TaxonomyDB,relationship = "many-to-many")
    if (TaxonomicRank != "kingdom") {
      colnames(ProjectDB_bySample)[which(names(ProjectDB_bySample) == TaxonomicRank)] <- "Latin_Name"
    }
    if (TaxonomicRank == "kingdom") {
      ProjectDB_bySample$Latin_Name <- ProjectDB_bySample$kingdom
    }
    #Convert NA values for image URLs to the default Phylopic image.
    ProjectDB_bySample["Image_URL"][is.na(ProjectDB_bySample["Image_URL"])] <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
    #Convert missing common names to blanks.
    ProjectDB_bySample["Common_Name"][is.na(ProjectDB_bySample["Common_Name"])] <- " "
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
      unique_samples <- unique(ProjectDB$SampleID)
      for(unique_sample in unique_samples){
        presence <- df[i, unique_sample]
        sample_data[[unique_sample]] <- presence
      }
      # Create the JSON object for the current row
      json_obj <- list(
        Latin_Name = latin_name,
        Common_Name = common_name,
        Image_URL = image_url,
        SampleID = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Aggregate merged data by the appropriate time interval to find taxa presence by time.
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="days"))) <= 7){
      ProjectDB_bySiteTime <- ProjectDB %>% dplyr::distinct(site,SampleID,day,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(site,day,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_bySiteTime$date_range <- paste(StartTime+days(ProjectDB_bySiteTime$day)-days(1),"to",StartTime+days(ProjectDB_bySiteTime$day))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 1 &
       as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) <= 4){
      ProjectDB_bySiteTime <- ProjectDB %>% dplyr::distinct(site,SampleID,week,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(site,week,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_bySiteTime$date_range <- paste(StartTime+weeks(ProjectDB_bySiteTime$week)-weeks(1),"to",StartTime+weeks(ProjectDB_bySiteTime$week))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 4 &
       as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) <= 13){
      ProjectDB_bySiteTime <- ProjectDB %>% dplyr::distinct(site,SampleID,month,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(site,month,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_bySiteTime$date_range <- paste(StartTime+months(ProjectDB_bySiteTime$month)-months(1),"to",StartTime+months(ProjectDB_bySiteTime$month))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 13 &
       as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) <= 52){
      ProjectDB_bySiteTime <- ProjectDB %>% dplyr::distinct(site,SampleID,quarter,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(site,quarter,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_bySiteTime$date_range <- paste(as.Date(as.yearqtr(paste(ProjectDB_bySiteTime$quarter,year(StartTime)), format = "Q%q %Y")),"to",as.Date(as.yearqtr(paste(ProjectDB_bySiteTime$quarter,year(StartTime)), format = "Q%q %Y"))+months(3))
    }
    if(as.numeric(abs(difftime(range(ProjectDB$sample_date)[1],range(ProjectDB$sample_date)[2],units="weeks"))) > 52){
      ProjectDB_bySiteTime <- ProjectDB %>% dplyr::distinct(site,SampleID,year,!!sym(TaxonomicRank),.keep_all=T) %>% 
        dplyr::group_by(site,year,!!sym(TaxonomicRank)) %>% 
        dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
      ProjectDB_bySiteTime$date_range <- paste(ProjectDB_bySiteTime$year,"to",ProjectDB_bySiteTime$year+1)
    }
    #Export taxa presence by time.
    ProjectDB_bySiteTime <- as.data.frame(ProjectDB_bySiteTime)
    #Merge in taxonomy data.
    ProjectDB_bySiteTime <- dplyr::left_join(ProjectDB_bySiteTime,TaxonomyDB,relationship = "many-to-many")
    #Merge in taxonomy data.
    if (TaxonomicRank != "kingdom") {
      colnames(ProjectDB_bySiteTime)[which(names(ProjectDB_bySiteTime) == TaxonomicRank)] <- "Latin_Name"
    }
    if (TaxonomicRank == "kingdom") {
      ProjectDB_bySiteTime$Latin_Name <- ProjectDB_bySiteTime$kingdom
    }
    #Convert NA values for image URLs to the default Phylopic image.
    ProjectDB_bySiteTime["Image_URL"][is.na(ProjectDB_bySiteTime["Image_URL"])] <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
    #Convert missing common names to blanks.
    ProjectDB_bySiteTime["Common_Name"][is.na(ProjectDB_bySiteTime["Common_Name"])] <- " "
    # Get unique sites.
    unique_sites <- unique(ProjectDB_bySiteTime$site)
    # Get unique date ranges.
    unique_dates <- unique(ProjectDB_bySiteTime$date_range)
    # Convert to a wide data frame.
    ProjectDB_bySiteTime <- tidyr::pivot_wider(ProjectDB_bySiteTime[,c("site","Latin_Name","freq","date_range","Common_Name","Image_URL")], names_from = site, values_from = freq, values_fill = 0)
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
      Date_Range <- df[i,"date_range"]
      # Create a list for the site/date/frequency values
      sample_data <- list()
      for(unique_site in unique_sites){
        freq <- df[i, unique_site]
        sample_data[[unique_site]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        Latin_Name = latin_name,
        Common_Name = common_name,
        Image_URL = image_url,
        Date_Range = Date_Range,
        Site = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    
    datasets <- list(datasets = list(results = json_list, metadata = SampleDB))
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    #Export filtered taxonomy table.
    TronkoTable <- ProjectDB[,c("SampleID",TaxonomicRank)]
    #TronkoTable$sum.taxonomy <- apply(ProjectDB[ ,TaxonomicRanks] , 1 , paste , collapse = ";" )
    TronkoTable <- as.data.frame(table(TronkoTable[,c("SampleID",TaxonomicRank)]))
    TronkoTable <- as.data.frame(pivot_wider(TronkoTable, names_from = SampleID, values_from = Freq))
    filename <- paste("FilteredTaxonomy_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,".csv",sep="_")
    filename <- gsub("_.csv",".csv",filename)
    filename <- tolower(filename)
    write.table(TronkoTable,filename,quote=FALSE,sep=",",row.names = FALSE)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/tables/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
  },
  error = function(e) {
    process_error(e, filename)
  }
)
