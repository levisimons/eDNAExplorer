#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
#require(aws.s3)
require(tidyr)
require(dplyr)
require(vegan)
require(ggplot2)
require(lubridate)
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
}
require(phyloseq)
require(htmlwidgets)
require(plotly)
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
    paste("s3://ednaexplorer/errors/alpha/", new_filename, sep = "")
  } else {
    paste("s3://ednaexplorer/projects/", ProjectID, "/plots/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
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
# EnvironmentalParameter:string Environmental variable to analyze against alpha diversity
# AlphaDiversity:string Alpha diversity metric
# Rscript --vanilla eDNAExplorer_Alpha_Metabarcoding.R "ProjectID" "First_Date" "Last_Date" "Marker" "Num_Mismatch" "TaxonomicRank" "CountThreshold" "FilterThreshold" "SpeciesList" "EnvironmentalParameter" "AlphaDiversity"

tryCatch(
  {
    if (length(args) != 11) {
      stop("Need the following inputs: ProjectID, First_Date, Last_Date, Marker, Num_Mismatch, TaxonomicRank, CountThreshold, FilterThreshold, SpeciesList, EnvironmentalParameter, AlphaDiversity.", call. = FALSE)
    } else if (length(args) == 11) {
      ProjectID <- args[1]
      First_Date <- args[2]
      Last_Date <- args[3]
      Marker <- args[4]
      Num_Mismatch <- args[5]
      TaxonomicRank <- args[6]
      CountThreshold <- args[7]
      FilterThreshold <- args[8]
      SpeciesList <- args[9]
      EnvironmentalParameter <- args[10]
      AlphaDiversity <- args[11]
      #Define filters in Phyloseq as global parameters.
      sample_ProjectID <<- as.character(ProjectID)
      sample_First_Date <<- lubridate::ymd(First_Date)
      sample_Last_Date <<- lubridate::ymd(Last_Date)
      sample_Primer <<- as.character(Marker)
      sample_TaxonomicRank <<- as.character(TaxonomicRank)
      sample_Num_Mismatch <<- as.numeric(Num_Mismatch)
      sample_CountThreshold <<- as.numeric(CountThreshold)
      sample_FilterThreshold <<- as.numeric(FilterThreshold)
      EnvironmentalVariable <<- as.character(EnvironmentalParameter)
      AlphaDiversityMetric <<- as.character(AlphaDiversity)
      SelectedSpeciesList <<- as.character(SpeciesList)
      
      CategoricalVariables <- c("site","grtgroup","biome_type","iucn_cat","eco_name","hybas_id")
      ContinuousVariables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
      FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
      TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
      
      #Save plot name.
      filename <- paste("Alpha_Metabarcoding_FirstDate",sample_First_Date,"LastDate",sample_Last_Date,"Marker",sample_Primer,"Rank",sample_TaxonomicRank,"Mismatch",sample_Num_Mismatch,"CountThreshold",sample_CountThreshold,"AbundanceThreshold",format(sample_FilterThreshold,scientific=F),"Variable",EnvironmentalVariable,"Diversity",AlphaDiversityMetric,"SpeciesList",SelectedSpeciesList,".json",sep="_")
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
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",sample_ProjectID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
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
    tmp <- Metadata %>% filter(projectid == sample_ProjectID)
    tmp <- as.data.frame(tmp)
    total_Samples <- nrow(tmp)
    Metadata <- Metadata %>%
      filter(projectid == sample_ProjectID) %>%
      filter(!is.na(latitude) & !is.na(longitude))
    Metadata <- as.data.frame(Metadata)
    Metadata$sample_date <- lubridate::ymd(Metadata$sample_date)
    Metadata <- Metadata %>% filter(sample_date >= sample_First_Date & sample_date <= sample_Last_Date)
    Metadata$fastqid <- gsub("_", "-", Metadata$fastqid)
    
    #Create sample metadata matrix
    if(nrow(Metadata) == 0 || ncol(Metadata) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    Sample <- Metadata[!is.na(Metadata$fastqid),]
    rownames(Sample) <- Sample$fastqid
    Sample$fastqid <- NULL
    Sample <- sample_data(Sample)
    remaining_Samples <- rownames(Sample)
    
    # Read in Tronko output and filter it.
    TronkoFile <- paste(Marker, ".csv", sep = "")
    TronkoFile_tmp <- paste(Marker,"_prevalence_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://ednaexplorer/tronko_output/", sample_ProjectID, "/", TronkoFile, " ", TronkoFile_tmp, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""))
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- paste("subset_prevalence_",UUIDgenerate(),".csv",sep="")
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
    if (SelectedSpeciesList != "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)) & TronkoDB$species %in% SpeciesList_df$name, ]
    }
    if (SelectedSpeciesList == "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)), ]
    }
    system(paste("rm ",TronkoFile_tmp,sep=""))
    system(paste("rm ",SubsetFile,sep=""))
    
    if(nrow(TronkoDB) > 1){
      #Create OTU matrix
      otumat <- as.data.frame(pivot_wider(as.data.frame(table(TronkoDB[,c("SampleID",sample_TaxonomicRank)])), names_from = SampleID, values_from = Freq))
      rownames(otumat) <- otumat[,sample_TaxonomicRank]
      otumat <- otumat[,colnames(otumat) %in% unique(TronkoDB$SampleID)]
      otumat[sapply(otumat, is.character)] <- lapply(otumat[sapply(otumat, is.character)], as.numeric)
      OTU <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)
      
      #Create merged Phyloseq object.
      physeq <- phyloseq(OTU,Sample)
      CountFilter <- physeq
      
      #Filter on read abundance per sample.
      tmpFilter  = transform_sample_counts(CountFilter, function(x) x / sum(x) )
      tmpFiltered = filter_taxa(tmpFilter, function(x) sum(x) > sample_FilterThreshold, TRUE)
      keeptaxa <- taxa_names(tmpFiltered)
      AbundanceFiltered <- prune_taxa(keeptaxa,CountFilter)
      
      if(EnvironmentalVariable %in% CategoricalVariables){
        #Store diversity versus variable data.
        tmp <- ggplot_build(plot_richness(AbundanceFiltered,x=EnvironmentalVariable,measures=AlphaDiversityMetric))
        tmp <- tmp$data[[1]]
        tmp$x <- as.factor(tmp$x)
        #Run a Kruskal-Wallis test between alpha diversity and selected environmental variable.
        if(length(unique(tmp$x))>1){
          test <- suppressWarnings(kruskal.test(tmp$y ~ tmp$x, data = tmp))
          Stats_Message <- paste("chi-squared = ",round(test$statistic,digits=3)," df = ",test$parameter," p-value = ",round(test$p.value,digits=3))
        } else {
          Stats_Message <- "Not enough variation in environmental variable to run Kruskal-Wallis test."
        }
        #General a violin plot of alpha diversity versus an environmental variable.
        p <- ggplot(tmp, aes(x=x, y=y))+
          labs(title=paste(AlphaDiversityMetric," versus ",gsub("_"," ",EnvironmentalVariable),".\nSamples collected between: ",sample_First_Date," and ",sample_Last_Date,"\nRelative abundance minimum of ",100*sample_FilterThreshold,"%.\nReads per sample minimum: ",sample_CountThreshold,"\n",Stats_Message,sep=""),x=gsub("_"," ",EnvironmentalVariable), y = AlphaDiversityMetric)+
          geom_violin()+theme_bw()+geom_point(position = position_jitter(seed = 1, width = 0.2))
      }
      if(EnvironmentalVariable %in% ContinuousVariables){
        #Store diversity versus variable data.
        tmp <- ggplot_build(plot_richness(AbundanceFiltered,x=EnvironmentalVariable,measures=AlphaDiversityMetric))
        tmp <- tmp$data[[1]]
        tmp$x <- as.numeric(tmp$x)
        #Calculate summary statistics using a Kendall correlation test between alpha diversity and environmental variable.
        if(length(unique(tmp$x))>1){
          test <- suppressWarnings(cor.test(tmp$x,tmp$y,alternative="two.sided",method="kendall"))
          Stats_Message <- paste("z = ",round(test$statistic,digits=3)," tau = ",round(test$estimate,digits=3)," p-value = ",round(test$p.value,digits=3))
        } else {
          Stats_Message <- "Not enough variation in environmental variable to run Kruskal-Wallis test."
        }
        #Generate a scatterplot of alpha diversity versus an environmental variable.
        p <- ggplot(tmp, aes(x=x, y=y))+
          labs(title=paste(AlphaDiversityMetric," versus ",gsub("_"," ",EnvironmentalVariable),".\nSamples collected between: ",sample_First_Date," and ",sample_Last_Date,"\nRelative abundance minimum of ",100*sample_FilterThreshold,"%.\nReads per sample minimum: ",sample_CountThreshold,"\n",Stats_Message,sep=""),x=gsub("_"," ",EnvironmentalVariable), y = AlphaDiversityMetric)+
          theme_bw()+geom_point()+geom_smooth()
      }
    } else {
      Stat_test <- "Not enough data to perform a Kruskal-Wallis test on alpha diversity."
      p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
    }
    
    #Insert the number of samples and number of samples post-filtering as a return object.
    SampleDB <- data.frame(matrix(ncol=2,nrow=1))
    colnames(SampleDB) <- c("totalSamples","filteredSamples")
    SampleDB$totalSamples <- nrow(Metadata)
    SampleDB$filteredSamples <- nsamples(AbundanceFiltered)
    datasets <- list(datasets = list(results=plotly_json(p, FALSE),metadata=toJSON(SampleDB)))
    
    #Save plot as json object
    filename <- paste("Alpha_Metabarcoding_FirstDate",sample_First_Date,"LastDate",sample_Last_Date,"Marker",sample_Primer,"Rank",sample_TaxonomicRank,"Mismatch",sample_Num_Mismatch,"CountThreshold",sample_CountThreshold,"AbundanceThreshold",format(sample_FilterThreshold,scientific=F),"Variable",EnvironmentalVariable,"Diversity",AlphaDiversityMetric,"SpeciesList",SelectedSpeciesList,".json",sep="_")
    filename <- gsub("_.json",".json",filename)
    filename <- tolower(filename)
    write(toJSON(datasets),filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",sample_ProjectID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
  },
  error = function(e) {
    process_error(e, filename)
  }
)
