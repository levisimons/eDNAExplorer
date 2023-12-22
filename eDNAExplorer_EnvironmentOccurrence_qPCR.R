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
bucket <- Sys.getenv("S3_BUCKET")
Database_Driver <- dbDriver("PostgreSQL")
ENDPOINT_URL <- Sys.getenv("ENDPOINT_URL")

if (length(args) != 6) {
  stop("Need the following inputs: ProjectID, First_Date, Last_Date, SpeciesList, EnvironmentalParameter, Taxon_name.", call. = FALSE)
} else if (length(args) == 6) {
  ProjectID <- args[1]
  First_Date <- args[2]
  Last_Date <- args[3]
  SpeciesList <- args[4]
  EnvironmentalParameter <- args[5]
  Taxon_name <- args[6]
  
  First_Date <- lubridate::ymd(First_Date)
  Last_Date <- lubridate::ymd(Last_Date)
  SelectedSpeciesList <- as.character(SpeciesList)
  EnvironmentalVariable <- as.character(EnvironmentalParameter)
  ProjectID <- as.character(ProjectID)
  Taxon_name <- as.character(Taxon_name)
  
  CategoricalVariables <- c("site","grtgroup","biome_type","iucn_cat","eco_name","hybas_id")
  ContinuousVariables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
  FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
  #Save plot name.
  filename <- paste("EnvironmentOccurrence_qPCR_FirstDate",First_Date,"LastDate",Last_Date,"Variable",EnvironmentalVariable,"SpeciesList",SelectedSpeciesList,"Taxon",Taxon_name,".json",sep="_")
  filename <- gsub("_.json",".json",filename)
  filename <- tolower(filename)
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
    paste("s3://",bucket,"/errors/environment_occurrence/", new_filename, " --endpoint-url ",ENDPOINT_URL,sep = "")
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
# SpeciesList:string Name of csv file containing selected species list.
# EnvironmentalParameter:string Environmental variable to analyze against alpha diversity
# Taxon_name:string Scientific taxon name
# Rscript --vanilla ednaexplorer_EnvironmentOccurrence_qPCR.R "ProjectID" "First_Date" "Last_Date""SpeciesList" "EnvironmentalParameter" "Taxon"

# Generate the output filename for cached plots.
tryCatch(
  {
    # Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
    data_to_write <- list(generating = TRUE, lastRanAt = Sys.time())
    write(toJSON(data_to_write), filename)
    dest_filename <- sub("\\.json$", ".build", filename) # Write to a temporary file first as .build
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/plots/",dest_filename," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
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
    
    #Read in qPCR data and filter it.
    qPCR_input <- tbl(con,"QPCRSample_test")
    Keep_Vars <- c(CategoricalVariables, ContinuousVariables, FieldVars)[c(CategoricalVariables, ContinuousVariables, FieldVars) %in% dbListFields(con, "QPCRSample_test")]
    # Get the number of samples in a project before filtering.
    qPCR_file <- paste("s3://",bucket,"/projects/",ProjectID,"/QPCR.csv",sep="")
    total_Samples <- as.numeric(system(paste("aws s3 cp ",qPCR_file," --endpoint-url ",ENDPOINT_URL," - | wc -l",sep=""),intern=T))-1
    #Filter samples
    qPCR_input <- qPCR_input %>%
      filter(projectid == ProjectID) %>%
      filter(!is.na(latitude) & !is.na(longitude))
    qPCR_input <- as.data.frame(qPCR_input)
    qPCR_input$sample_date <- lubridate::ymd(qPCR_input$sample_date)
    qPCR_input <- qPCR_input %>% filter(sample_date >= First_Date & sample_date <= Last_Date) %>% 
      filter(!is.na(!!sym(EnvironmentalParameter))) %>% filter(target_organism == Taxon_name)
    #Filter results by species list.
    if (SelectedSpeciesList != "None") {
      qPCR_input <- qPCR_input %>% filter(target_organism %in% SpeciesList_df$name)
    }
    if(nrow(qPCR_input) == 0 || ncol(qPCR_input) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    filteredSamples <- nrow(qPCR_input)
    
    #Read in information to map categorical labels for certain variables.
    category_file <- paste("Categories_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://",bucket,"/analysis/Categories.csv ",category_file," --endpoint-url ",ENDPOINT_URL,sep=""))
    categories <- as.data.frame(fread(file=category_file,header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A")))
    system(paste("rm ",category_file,sep=""))
    #Read in information for legends and labels
    legends_file <- paste("LabelsAndLegends_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://",bucket,"/analysis/LabelsAndLegends.csv ",legends_file," --endpoint-url ",ENDPOINT_URL,sep=""))
    legends_and_labels <- as.data.frame(fread(file=legends_file,header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A")))
    system(paste("rm ",legends_file,sep=""))
    #Set up new legends and x-axis labels.
    new_legend <- legends_and_labels[legends_and_labels$Environmental_Variable==EnvironmentalParameter,"Legend"]
    
    #Plot qPCR presence/absence for a taxon against an environmental variable
    qPCR_plot <- qPCR_input[,c(EnvironmentalParameter,"target_was_organism_detected")]
    qPCR_plot$Detected <- as.factor(ifelse(qPCR_plot$target_was_organism_detected == 0, 'No', 'Yes'))
    qPCR_plot$x <- qPCR_plot$Detected

    if(EnvironmentalVariable %in% CategoricalVariables){
      if(EnvironmentalParameter %in% unique(categories$Environmental_Variable)){
        qPCR_plot$value <- as.character(qPCR_plot[,c(EnvironmentalParameter)])
        qPCR_plot <- dplyr::left_join(qPCR_plot,categories[categories$Environmental_Variable==EnvironmentalParameter,])
        qPCR_plot["description"][is.na(qPCR_plot["description"])] <- "no data available"
        qPCR_plot$y <- as.factor(qPCR_plot$description)
      }
      if(length(unique(na.omit(qPCR_plot$x)))>1 & length(unique(na.omit(qPCR_plot$y)))>1){
        test <- suppressWarnings(kruskal.test(qPCR_plot$y ~ qPCR_plot$x, data = qPCR_plot))
        Stats_Message <- paste("chi-squared = ",round(test$statistic,digits=3)," df = ",test$parameter," p-value = ",round(test$p.value,digits=3))
      } else {
        Stats_Message <- "Not enough variation in environmental variable to run Kruskal-Wallis test."
      }
      #Generate a scatterplot of detection versus an environmental variable.
      p <- ggplot(qPCR_plot, aes(x=Detected, y=EnvironmentalParameter))+
        labs(title=paste("Organism detection versus ",gsub("_"," ",EnvironmentalParameter),".\nSamples collected between: ",First_Date," and ",Last_Date,"\n",Stats_Message,sep=""),x="Detected", y = new_legend)+
        geom_violin()+theme_bw()+geom_point(position = position_jitter(seed = 1, width = 0.2))+guides(fill=guide_legend(title=new_legend))
    }
    
    if(EnvironmentalVariable %in% ContinuousVariables){
      qPCR_plot$y <- as.numeric(qPCR_plot[,EnvironmentalParameter])
      #Run a Kruskal-Wallis test between alpha diversity and selected environmental variable.
      if(length(unique(na.omit(qPCR_plot$x)))>1 & length(unique(na.omit(qPCR_plot$y)))>1){
        test <- suppressWarnings(kruskal.test(qPCR_plot$y ~ qPCR_plot$x, data = qPCR_plot))
        Stats_Message <- paste("chi-squared = ",round(test$statistic,digits=3)," df = ",test$parameter," p-value = ",round(test$p.value,digits=3))
      } else {
        Stats_Message <- "Not enough variation in environmental variable to run Kruskal-Wallis test."
      }
      #Generate a scatterplot of detection versus an environmental variable.
      p <- ggplot(qPCR_plot, aes(x=Detected, y=EnvironmentalParameter))+
        labs(title=paste("Organism detection versus ",gsub("_"," ",EnvironmentalParameter),".\nSamples collected between: ",First_Date," and ",Last_Date,"\n",Stats_Message,sep=""),x="Detected", y = new_legend)+
        geom_violin()+theme_bw()+geom_point(position = position_jitter(seed = 1, width = 0.2))+guides(fill=guide_legend(title=new_legend))
    }
    
    # Insert the number of samples and number of samples post-filtering as a return object.
    SampleDB <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(SampleDB) <- c("totalSamples", "filteredSamples")
    SampleDB$totalSamples <- total_Samples
    SampleDB$filteredSamples <- num_filteredSamples
    datasets <- list(datasets = list(results=plotly_json(p, FALSE),metadata=toJSON(SampleDB)))
    
    #Save plot.
    filename <- paste("EnvironmentOccurrence_qPCR_FirstDate",First_Date,"LastDate",Last_Date,"Variable",EnvironmentalVariable,"SpeciesList",SelectedSpeciesList,"Taxon",Taxon_name,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    write(toJSON(datasets), filename)
    system(paste("aws s3 cp ", filename, " s3://",bucket,"/projects/", Project_ID, "/plots/", filename, " --endpoint-url ",ENDPOINT_URL, sep = ""), intern = TRUE)
    system(paste("rm ", filename, sep = ""))
    sapply(dbListConnections(Database_Driver), dbDisconnect)
  },
  error = function(e) {
    process_error(e)
  }
)