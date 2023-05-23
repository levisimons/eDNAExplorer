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
#* @param BetaDiversity:string Beta diversity metric. Options are chao, bray, or jaccard.
#* @get /beta
beta <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,BetaDiversity){
  
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
  
  CategoricalVariables <- c("grtgroup","biome_type","iucn_Cat","eco_name","hybas_id")
  ContinuousVariables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
  FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
  #Establish sql connection
  Database_Driver <- dbDriver("PostgreSQL")
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Read in metadata and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  Metadata <- tbl(con,"TronkoMetadata")
  Keep_Vars <- c(CategoricalVariables,ContinuousVariables,FieldVars)[c(CategoricalVariables,ContinuousVariables,FieldVars) %in% dbListFields(con,"TronkoMetadata")]
  Metadata <- Metadata %>% filter(sample_date >= sample_First_Date & sample_date <= sample_Last_Date) %>%
    filter(ProjectID == sample_ProjectID) %>% filter(!is.na(latitude) & !is.na(longitude)) %>% select(Keep_Vars)
  Metadata <- as.data.frame(Metadata)
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
  #Create sample metadata matrix
  Sample <- Metadata[!is.na(Metadata$fastqid),]
  rownames(Sample) <- Sample$fastqid
  Sample$fastqid <- NULL
  Sample <- sample_data(Sample)
  remaining_Samples <- rownames(Sample)
  
  #Read in Tronko output and filter it.
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  TronkoInput <- tbl(con,"TronkoOutput")
  TronkoInput <- TronkoInput %>% filter(ProjectID == sample_ProjectID) %>% filter(Primer == sample_Primer) %>% 
    filter(Mismatch <= sample_Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(sample_TaxonomicRank))) %>%
    group_by(SampleID) %>% filter(n() > sample_CountThreshold) %>% 
    select(SampleID,sample_TaxonomicRank)
  TronkoDB <- as.data.frame(TronkoInput)
  TronkoDB <- TronkoDB[TronkoDB$SampleID %in% rownames(Sample),]
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  
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
  
  #Plot and analyze beta diversity versus an environmental variables.
  if(nsamples(AbundanceFiltered)>1){
    if(BetaDiversityMetric!="jaccard"){BetaDist = phyloseq::distance(AbundanceFiltered, method=BetaDiversityMetric, weighted=F)}
    if(BetaDiversityMetric=="jaccard"){BetaDist = phyloseq::distance(AbundanceFiltered, method=BetaDiversityMetric, weighted=F,binary=T)}
    if(sum(!is.nan(BetaDist))>1){
      ordination = ordinate(AbundanceFiltered, method="PCoA", distance=BetaDist)
      if(length(unique(sample_data(AbundanceFiltered)[,EnvironmentalVariable]))>1){
        BetaExpression = paste("adonis2(BetaDist ~ sample_data(AbundanceFiltered)$",EnvironmentalVariable,")",sep="")
        test <- eval(parse(text=BetaExpression))
        Stat_test <- paste("PERMANOVA results, using 999 permutations.\n",DiversityMetric," beta diversity and ",EnvironmentalVariable,"\nDegrees of freedom: ",round(test$Df[1],3),". Sum of squares: ",round(test$SumOfSqs[1],3),". R-squared: ",round(test$R2[1],3),". F-statistic: ",round(test$F[1],3),". p: ",round(test$`Pr(>F)`[1],3),sep="")
      } else {
        Stat_test <- "Not enough variation in environmental data to analyze against beta diversity."
      }
      p <- plot_ordination(AbundanceFiltered, ordination, color=EnvironmentalVariable) + theme(aspect.ratio=1) + labs(title = Stat_test, color = EnvironmentalVariable)
      p <- p+theme_bw()
    } else{
      Stat_test <- "Not enough remaining data after filters to analyze beta diversity."
      p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
    }
  } else {
    Stat_test <- "Not enough data to analyze beta diversity."
    p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
  }
  #Save plot as json object
  jfig <- plotly_json(p, FALSE)
  return(jfig)
}
