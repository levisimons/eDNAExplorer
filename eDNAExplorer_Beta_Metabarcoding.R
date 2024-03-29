#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
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
  stack_trace <- paste(capture.output(traceback()), collapse = "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, error = error_message, stack_trace = stack_trace))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  
  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://ednaexplorer/errors/beta/", new_filename, sep = "")
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
# BetaDiversity:string Alpha diversity metric
# Rscript --vanilla eDNAExplorer_Beta_Metabarcoding.R "ProjectID" "First_Date" "Last_Date" "Marker" "Num_Mismatch" "TaxonomicRank" "CountThreshold" "FilterThreshold" "SpeciesList" "EnvironmentalParameter" "BetaDiversity"

tryCatch(
  {
    if (length(args) != 11) {
      stop("Need the following inputs: ProjectID, First_Date, Last_Date, Marker, Num_Mismatch, TaxonomicRank, CountThreshold, FilterThreshold, SpeciesList, EnvironmentalParameter, BetaDiversity", call. = FALSE)
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
      BetaDiversity <- args[11]
      
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
      BetaDiversityMetric <<- as.character(BetaDiversity)
      SelectedSpeciesList <<- as.character(SpeciesList)
      
      CategoricalVariables <- c("site","grtgroup","biome_type","iucn_cat","eco_name","hybas_id")
      ContinuousVariables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
      FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
      TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
      TaxonomicNum <<- as.numeric(which(TaxonomicRanks == sample_TaxonomicRank))
      
      #Save plot name.
      filename <- paste("Beta_Metabarcoding_FirstDate",sample_First_Date,"LastDate",sample_Last_Date,"Marker",sample_Primer,"Rank",sample_TaxonomicRank,"Mismatch",sample_Num_Mismatch,"CountThreshold",sample_CountThreshold,"AbundanceThreshold",format(sample_FilterThreshold,scientific=F),"Variable",EnvironmentalVariable,"Diversity",BetaDiversityMetric,"SpeciesList",SelectedSpeciesList,",json",sep="_")
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
    Metadata_Unfiltered <- Metadata %>% filter(projectid == sample_ProjectID)
    Metadata_Unfiltered <- as.data.frame(Metadata_Unfiltered)
    total_Samples <- nrow(Metadata_Unfiltered)
    Metadata <- Metadata %>%
      filter(projectid == sample_ProjectID) %>%
      filter(!is.na(latitude) & !is.na(longitude)) %>%
      filter(!is.na(!!sym(EnvironmentalVariable)))
    Metadata <- as.data.frame(Metadata)
    Metadata$sample_date <- lubridate::ymd(Metadata$sample_date)
    Metadata <- Metadata %>% filter(sample_date >= sample_First_Date & sample_date <= sample_Last_Date)
    if(nrow(Metadata) == 0 || ncol(Metadata) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
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
    TronkoFile <- paste(sample_Primer, ".csv", sep = "")
    TronkoFile_tmp <- paste(sample_Primer,"_beta_",UUIDgenerate(),".csv",sep="")
    system(paste("aws s3 cp s3://ednaexplorer/tronko_output/", sample_ProjectID, "/", TronkoFile, " ", TronkoFile_tmp, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = ""))
    #Check if file exists.
    if(file.info(TronkoFile_tmp)$size== 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    SubsetFile <- paste("subset_beta_",UUIDgenerate(),".csv",sep="")
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
      filter(n() > sample_CountThreshold) %>% filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% select(-Mismatch)
    #Merege relative abudance results into data filtered by reads per sample and mismatches.
    #Then filtered on relative abundances.
    TronkoInput <- TronkoInput %>%
      left_join(Tronko_Unfiltered,na_matches="never") %>%
      filter(freq >= FilterThreshold) %>% select(-freq)
    TronkoDB <- as.data.frame(TronkoInput)
    #Remove taxa which are unknown at a given rank.
    TronkoDB <- TronkoDB[,c("SampleID",TaxonomicRanks[1:TaxonomicNum])]
    TronkoDB <- TronkoDB[!is.na(TronkoDB[, sample_TaxonomicRank]), ]
    #Filter results by species list.
    if (SelectedSpeciesList != "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)) & TronkoDB$species %in% SpeciesList_df$name, ]
    }
    if (SelectedSpeciesList == "None") {
      TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)), ]
    }
    system(paste("rm ",TronkoFile_tmp,sep=""))
    system(paste("rm ",SubsetFile,sep=""))
    
    if(nrow(TronkoDB)>1){
      #Create OTU matrix
      otumat <- as.data.frame(pivot_wider(as.data.frame(table(TronkoDB[,c("SampleID",sample_TaxonomicRank)])), names_from = SampleID, values_from = Freq))
      rownames(otumat) <- otumat[,sample_TaxonomicRank]
      otumat <- otumat[,colnames(otumat) %in% unique(TronkoDB$SampleID)]
      otumat[sapply(otumat, is.character)] <- lapply(otumat[sapply(otumat, is.character)], as.numeric)
      OTU <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)
            
      #Create merged Phyloseq object.
      physeq <- phyloseq(OTU,Sample)
      AbundanceFiltered <- physeq
            
      #Plot and analyze beta diversity versus an environmental variables.
      if(nsamples(AbundanceFiltered)>1 & ntaxa(AbundanceFiltered)>1){
        if(BetaDiversityMetric!="jaccard"){BetaDist = phyloseq::distance(AbundanceFiltered, method=BetaDiversityMetric, weighted=F)}
        if(BetaDiversityMetric=="jaccard"){BetaDist = phyloseq::distance(AbundanceFiltered, method=BetaDiversityMetric, weighted=F,binary=T)}
        if(sum(!is.nan(BetaDist))>1){
          ordination = ordinate(AbundanceFiltered, method="PCoA", distance=BetaDist)
          AbundanceFiltered_df <- data.frame(sample_data(AbundanceFiltered))
          if(length(unique(AbundanceFiltered_df[,EnvironmentalVariable]))>1){
            BetaExpression = paste("adonis2(BetaDist ~ sample_data(AbundanceFiltered)$",EnvironmentalVariable,")",sep="")
            test <- eval(parse(text=BetaExpression))
            Stat_test <- paste("PCA plot.  Results of PERMANOVA, using 999 permutations.\n",BetaDiversityMetric," beta diversity and ",EnvironmentalVariable,"\nDegrees of freedom: ",round(test$Df[1],3),". Sum of squares: ",round(test$SumOfSqs[1],3),". R-squared: ",round(test$R2[1],3),". F-statistic: ",round(test$F[1],3),". p: ",round(test$`Pr(>F)`[1],3),sep="")
          } else {
            Stat_test <- paste("PCA plot.  Not enough variation in ",EnvironmentalVariable," to perform a PERMANOVA on beta diversity.",sep="")
          }
          p <- plot_ordination(AbundanceFiltered, ordination, color=EnvironmentalVariable) + theme(aspect.ratio=1) + labs(title = Stat_test, color = EnvironmentalVariable)
          p <- p+theme_bw()
        } else{
          Stat_test <- "PCA plot.  Not enough remaining data after filters to perform a PERMANOVA on beta diversity."
          p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
        }
      } else {
        Stat_test <- "PCA plot.  Not enough data to perform a PERMANOVA on beta diversity."
        p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
      }
    } else {
      Stat_test <- "PCA plot.  Not enough data to perform a PERMANOVA on beta diversity."
      p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
    }

    #Insert the number of samples and number of samples post-filtering as a return object.
    SampleDB <- data.frame(matrix(ncol=2,nrow=1))
    colnames(SampleDB) <- c("totalSamples","filteredSamples")
    SampleDB$totalSamples <- total_Samples
    SampleDB$filteredSamples <- nsamples(AbundanceFiltered)
    datasets <- list(datasets = list(results=plotly_json(p, FALSE),metadata=toJSON(SampleDB)))
    
    #Save plot as json object
    filename <- paste("Beta_Metabarcoding_FirstDate",sample_First_Date,"LastDate",sample_Last_Date,"Marker",sample_Primer,"Rank",sample_TaxonomicRank,"Mismatch",sample_Num_Mismatch,"CountThreshold",sample_CountThreshold,"AbundanceThreshold",format(sample_FilterThreshold,scientific=F),"Variable",EnvironmentalVariable,"Diversity",BetaDiversityMetric,"SpeciesList",SelectedSpeciesList,",json",sep="_")
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
