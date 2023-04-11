# plumber.R
library(plumber)
require(aws.s3)
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
require(RSQLite)
require(digest)

Sys.setenv("AWS_ACCESS_KEY_ID" = "",
           "AWS_SECRET_ACCESS_KEY" = "")

#* Echo the parameter that was sent in
#* @param ProjectID:string
#* @param First_Date:string YYYY-MM-DD
#* @param Last_Date:string YYYY-MM-DD
#* @param Marker:string Target marker name
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @param EnvironmentalParameter:string Environmental variable to analyze against alpha diversity
#* @param AlphaDiversity:string Alpha diversity metric
#* @get /alpha
alpha <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity){
  
  #Define filters in Phyloseq as global parameters.
  sample_ProjectID <<- as.character(ProjectID)
  sample_First_Date <<- as.numeric(as.POSIXct(First_Date))
  sample_Last_Date <<- as.numeric(as.POSIXct(Last_Date))
  sample_Primer <<- as.character(Marker)
  sample_TaxonomicRank <<- as.character(TaxonomicRank)
  sample_Num_Mismatch <<- as.numeric(Num_Mismatch)
  sample_CountThreshold <<- as.numeric(CountThreshold)
  sample_FilterThreshold <<- as.numeric(FilterThreshold)
  EnvironmentalVariable <<- as.character(EnvironmentalParameter)
  AlphaDiversityMetric <<- as.character(AlphaDiversity)
  
  CategoricalVariables <- c("grtgroup","biome_type","IUCN_CAT","ECO_NAME","HYBAS_ID")
  ContinuousVariables <- c("bio01","bio12","gHM","elevation","NDVI","Average_Radiance")
  FieldVars <- c("FastqID","Sample Date","Latitude","Longitude","Spatial Uncertainty")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
  #Establish sql connection
  Database_Driver <- dbDriver("SQLite")
  db_host <- ""
  db_port <- 
  db_name <- ""
  db_user <- ""
  db_pass <- ""
  
  #Read in metadata and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  Metadata <- tbl(con,"TronkoMetadata")
  Keep_Vars <- c(CategoricalVariables,ContinuousVariables,FieldVars)[c(CategoricalVariables,ContinuousVariables,FieldVars) %in% dbListFields(con,"TronkoMetadata")]
  Metadata <- Metadata %>% filter(`Sample Date` >= sample_First_Date & `Sample Date` <= sample_Last_Date) %>%
    filter(`ProjectID` == sample_ProjectID) %>% filter(!is.na(Latitude) & !is.na(Longitude)) %>% select(Keep_Vars)
  Metadata <- as.data.frame(Metadata)
  dbDisconnect(con)
  
  #Create sample metadata matrix
  Sample <- Metadata[!is.na(Metadata$FastqID),]
  rownames(Sample) <- Sample$FastqID
  Sample$FastqID <- NULL
  Sample <- sample_data(Sample)
  remaining_Samples <- rownames(Sample)
  
  #Read in Tronko output and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  TronkoInput <- tbl(con,"TronkoOutput")
  TronkoInput <- TronkoInput %>% filter(`ProjectID` == sample_ProjectID) %>% filter(Primer == sample_Primer) %>% 
    filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(sample_TaxonomicRank))) %>%
    group_by(SampleID) %>% filter(n() > sample_CountThreshold) %>% 
    select(SampleID,sample_TaxonomicRank)
  TronkoDB <- as.data.frame(TronkoInput)
  TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample),]
  dbDisconnect(con)
  
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
  #Save plot as json object
  jfig <- plotly_json(p, FALSE)
  return(jfig)
}
