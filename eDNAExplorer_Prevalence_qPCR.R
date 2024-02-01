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
require(uuid)
require(ggplot2)
require(plotly)
require(jsonlite)

# Fetch project ID early so we can use it for error output when possible.
project_id <- args[1]

#Establish database credentials.
readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
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
  
  s3_path <- if (is.null(project_id) || project_id == "") {
    paste("s3://",bucket,"/errors/qpcr_prevalence/", new_filename, sep = "")
  } else {
    paste("s3://",bucket,"/projects/", project_id, "/plots/error+", dest_filename, " --endpoint-url ",ENDPOINT_URL, sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}

# Get filtering parameters.
# project_id:string
# first_date:string YYYY-MM-DD
# last_date:string YYYY-MM-DD
# species_list:string Name of csv file containing selected species list.
# Rscript --vanilla eDNAExplorer_Prevalence_qPCR.R "project_id" "first_date" "last_date" "marker" "species_list"

tryCatch(
  {
    if (length(args) != 4) {
      stop("Need the following inputs: project_id, first_date, last_date, species_list.", call. = FALSE)
    } else if (length(args) == 4) {
      project_id <- args[1]
      first_date <- args[2]
      last_date <- args[3]
      species_list <- args[4]
      
      first_date <- lubridate::ymd(first_date)
      last_date <- lubridate::ymd(last_date)
      selected_species_list <- as.character(species_list)
      project_id <- as.character(project_id)
      
      categorical_variables <- c("site","grtgroup","biome_type","iucn_cat","eco_name","hybas_id")
      continuous_variables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
      FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
      taxonomic_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
      taxonomic_key_ranks <- c("kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")
      
      #Save plot name.
      filename <- paste("Prevalence_qPCR_FirstDate",first_date,"LastDate",last_date,"species_list",selected_species_list,".json",sep="_")
      filename <- gsub("_.json",".json",filename)
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
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",project_id,"/plots/",dest_filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Establish sql connection
    Database_Driver <- dbDriver("PostgreSQL")
    sapply(dbListConnections(Database_Driver), dbDisconnect)
    con <- dbConnect(Database_Driver, host = db_host, port = db_port, dbname = db_name, user = db_user, password = db_pass)
    
    # Read in species list
    if (selected_species_list != "None") {
      SpeciesList_df <- tbl(con, "SpeciesListItem")
      SpeciesList_df <- SpeciesList_df %>% filter(species_list == selected_species_list)
      SpeciesList_df <- as.data.frame(SpeciesList_df)
    }
    
    #Read in qPCR data and filter it.
    qPCR_input <- tbl(con,"QPCRSample_test")
    Keep_Vars <- c(categorical_variables, continuous_variables, FieldVars)[c(categorical_variables, continuous_variables, FieldVars) %in% dbListFields(con, "QPCRSample_test")]
    # Get the number of samples in a project before filtering.
    qPCR_file <- paste("s3://",bucket,"/projects/",project_id,"/QPCR.csv",sep="")
    total_Samples <- as.numeric(system(paste("aws s3 cp ",qPCR_file," --endpoint-url ",ENDPOINT_URL," - | wc -l",sep=""),intern=T))-1
    #Filter samples
    qPCR_input <- qPCR_input %>%
      filter(projectid == project_id) %>%
      filter(!is.na(latitude) & !is.na(longitude))
    qPCR_input <- as.data.frame(qPCR_input)
    qPCR_input$sample_date <- lubridate::ymd(qPCR_input$sample_date)
    qPCR_input <- qPCR_input %>% filter(sample_date >= first_date & sample_date <= last_date)
    #Filter results by species list.
    if (selected_species_list != "None") {
      qPCR_input <- qPCR_input %>% filter(target_organism %in% SpeciesList_df$name)
    }
    if(nrow(qPCR_input) == 0 || ncol(qPCR_input) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    num_filteredSamples <- nrow(qPCR_input)
    
    # Read in Taxonomy output and filter it.
    TaxonomyInput <- tbl(con, "Taxonomy")
    # Get taxonomic rank
    taxonomic_rank <- tolower(unique(na.omit(qPCR_input$target_taxonomic_rank_of_organism)))
    taxonomic_num <- which(taxonomic_ranks == taxonomic_rank)
    if (nrow(qPCR_input) > 0) {
      TaxaList <- na.omit(unique(qPCR_input$target_organism))
    }
    if (nrow(qPCR_input) == 0) {
      TaxaList <- c()
    }
    if (taxonomic_rank != "kingdom") {
      TaxonomyInput <- TaxonomyInput %>%
        filter(rank==taxonomic_rank) %>%
        filter(!!sym(taxonomic_rank) %in% TaxaList) %>%
        select(taxonomic_ranks[2:taxonomic_num],taxonomic_key_ranks,Common_Name,Image_URL)
      TaxonomyDB <- as.data.frame(TaxonomyInput)
      #Figure out which taxonomy version is more complete.
      TaxonomyDB$rankCount <- rowSums(!is.na(TaxonomyDB[,colnames(TaxonomyDB) %in% taxonomic_ranks]))
      TaxonomyDB$rankKeyCount <- rowSums(!is.na(TaxonomyDB[,colnames(TaxonomyDB) %in% taxonomic_key_ranks]))
      TaxonomyDB <- TaxonomyDB %>%
        group_by(!!sym(taxonomic_rank)) %>%
        slice_max(order_by = rankCount, n = 1) %>%
        ungroup()
      TaxonomyDB <- TaxonomyDB %>%
        group_by(!!sym(taxonomic_rank)) %>%
        slice_max(order_by = rankKeyCount, n = 1) %>%
        ungroup()
      #Figure out which common_name is most common per taxon.
      TaxonomyDB <- TaxonomyDB %>%
        group_by(!!sym(taxonomic_rank)) %>%
        mutate(Most_Common_Name = ifelse(all(is.na(Common_Name)), NA, names(which.max(table(Common_Name[!is.na(Common_Name)]))))) %>%
        ungroup()
      TaxonomyDB$Common_Name <- TaxonomyDB$Most_Common_Name
      TaxonomyDB$Most_Common_Name <- NULL
      TaxonomyDB <- as.data.frame(TaxonomyDB)
      TaxonomyDB <- subset(TaxonomyDB, select = -grep("Key", colnames(TaxonomyDB)))
      TaxonomyDB$rankCount <- NULL
      TaxonomyDB <- TaxonomyDB[!duplicated(TaxonomyDB),]
    }
    if (taxonomic_rank == "kingdom") {
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
    
    if(nrow(qPCR_input)>0){
      qPCR_input <- qPCR_input %>%
        dplyr::group_by(target_organism) %>%
        dplyr::summarise(per = n_distinct(sample_id)/n_distinct(qPCR_input$sample_id))
      qPCR_input <- as.data.frame(qPCR_input)
      #Merge in taxonomy data.
      colnames(qPCR_input)[which(names(qPCR_input) == "target_organism")] <- taxonomic_rank
      qPCR_input <- dplyr::left_join(qPCR_input, TaxonomyDB,na_matches="never")
      qPCR_input$Image_URL <- ifelse(is.na(qPCR_input$Image_URL), 'https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png', qPCR_input$Image_URL)
      if (taxonomic_rank != "kingdom") {
        colnames(qPCR_input)[which(names(qPCR_input) == taxonomic_rank)] <- "Latin_Name"
      }
      if (taxonomic_rank == "kingdom") {
        qPCR_input$Latin_Name <- qPCR_input$kingdom
      }
      qPCR_input <- qPCR_input %>% mutate_all(~ifelse(is.na(.), " ", .))
      #De-duplicated dataframe for tronko table.
      qPCR_input <- qPCR_input[!duplicated(qPCR_input),]
      
      # Insert the number of samples and number of samples post-filtering as a return object.
      SampleDB <- data.frame(matrix(ncol = 2, nrow = 1))
      colnames(SampleDB) <- c("totalSamples", "filteredSamples")
      SampleDB$totalSamples <- total_Samples
      SampleDB$filteredSamples <- num_filteredSamples
      datasets <- list(datasets = list(results=plotly_json(p, FALSE),metadata=toJSON(SampleDB)))
      
      #Save plot.
      filename <- paste("Prevalence_qPCR_FirstDate",first_date,"LastDate",last_date,"species_list",selected_species_list,".json",sep="_")
      filename <- gsub("_.json",".json",filename)
      filename <- tolower(filename)
      write(toJSON(datasets), filename)
      system(paste("aws s3 cp ", filename, " s3://",bucket,"/projects/", project_id, "/plots/", filename, " --endpoint-url ",ENDPOINT_URL, sep = ""), intern = TRUE)
      system(paste("rm ", filename, sep = ""))
    }
  }
)