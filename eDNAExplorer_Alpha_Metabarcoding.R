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
require(RPostgreSQL)
require(digest)

readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")

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
#* @param SpeciesList:string Name of csv file containing selected species list.
#* @get /alpha
alpha <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity,SpeciesList){
  
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
  
  CategoricalVariables <- c("grtgroup","biome_type","iucn_cat","eco_name","hybas_id")
  ContinuousVariables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
  FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
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
  
  #Save default blank plot if not enough data is available.
  Stat_test <- "Not enough data to perform a Kruskal-Wallis test on alpha diversity."
  p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
  jfig <- plotly_json(p, FALSE)
  filename <- paste("Alpha_Metabarcoding_FirstDate",sample_First_Date,"LastDate",sample_Last_Date,"Marker",sample_Primer,"Rank",sample_TaxonomicRank,"Mismatch",sample_Num_Mismatch,"CountThreshold",sample_CountThreshold,"AbundanceThreshold",format(sample_FilterThreshold,scientific=F),"Variable",EnvironmentalVariable,"Diversity",AlphaDiversityMetric,"SpeciesList",SelectedSpeciesList,".json",sep="_")
  filename <- gsub("_.json",".json",filename)
  filename <- tolower(filename)
  write(jfig,filename)
  system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",sample_ProjectID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  system(paste("rm ",filename,sep=""))
  
  #Read in metadata and filter it.
  Metadata <- tbl(con,"TronkoMetadata")
  Keep_Vars <- c(CategoricalVariables,ContinuousVariables,FieldVars)[c(CategoricalVariables,ContinuousVariables,FieldVars) %in% dbListFields(con,"TronkoMetadata")]
  Metadata <- Metadata %>% filter(sample_date >= sample_First_Date & sample_date <= sample_Last_Date) %>%
    filter(projectid == sample_ProjectID) %>% filter(!is.na(latitude) & !is.na(longitude)) %>% select(Keep_Vars)
  Metadata <- as.data.frame(Metadata)
  Metadata$fastqid <- gsub("_","-",Metadata$fastqid)
  
  #Create sample metadata matrix
  Sample <- Metadata[!is.na(Metadata$fastqid),]
  rownames(Sample) <- Sample$fastqid
  Sample$fastqid <- NULL
  Sample <- sample_data(Sample)
  remaining_Samples <- rownames(Sample)
  
  #Read in Tronko output and filter it.
  TronkoFile <- paste(Marker,".csv",sep="")
  system(paste("aws s3 cp s3://ednaexplorer/tronko_output/",sample_ProjectID,"/",TronkoFile," ",TronkoFile," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
  system(paste("cut -d ',' -f 6,7,8,9,10,11,12,13,14,16 ",TronkoFile," > subset.csv",sep=""))
  TronkoInput <- fread(file="subset.csv",header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
  if(sample_TaxonomicRank != "species"){
    TronkoInput <- TronkoInput %>% filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(sample_TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > sample_CountThreshold) %>% 
      select(SampleID,species,sample_TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    if(SelectedSpeciesList != "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample) & TronkoDB$species %in% SpeciesList_df$name,]}
    if(SelectedSpeciesList == "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample),]}
    TronkoDB$species <- NULL
  } else{
    TronkoInput <- TronkoInput %>% filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(sample_TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > sample_CountThreshold) %>% 
      select(SampleID,sample_TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    if(SelectedSpeciesList != "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample) & TronkoDB$species %in% SpeciesList_df$name,]}
    if(SelectedSpeciesList == "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample),]}
  }
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  system(paste("rm",TronkoFile,sep=" "))
  system("rm subset.csv")
  
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
  #Save plot as json object
  jfig <- plotly_json(p, FALSE)
  filename <- paste("Alpha_Metabarcoding_FirstDate",sample_First_Date,"LastDate",sample_Last_Date,"Marker",sample_Primer,"Rank",sample_TaxonomicRank,"Mismatch",sample_Num_Mismatch,"CountThreshold",sample_CountThreshold,"AbundanceThreshold",format(sample_FilterThreshold,scientific=F),"Variable",EnvironmentalVariable,"Diversity",AlphaDiversityMetric,"SpeciesList",SelectedSpeciesList,".json",sep="_")
  filename <- gsub("_.json",".json",filename)
  filename <- tolower(filename)
  write(jfig,filename)
  system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",sample_ProjectID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  system(paste("rm ",filename,sep=""))  
}
