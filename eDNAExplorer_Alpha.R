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

Sys.setenv("AWS_ACCESS_KEY_ID" = "e9190baae65b40a38bf43ade883b04a6",
           "AWS_SECRET_ACCESS_KEY" = "7aa129a7d84744efa76183cc9cf4b0a5")

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
  sample_First_Date <<- ymd(First_Date)
  sample_Last_Date <<- ymd(Last_Date)
  sample_Primer <<- as.character(Marker)
  sample_TaxonomicRank <<- as.character(TaxonomicRank)
  sample_Num_Mismatch <<- as.numeric(Num_Mismatch)
  sample_CountThreshold <<- as.numeric(CountThreshold)
  sample_FilterThreshold <<- as.numeric(FilterThreshold)
  EnvironmentalVariable <<- as.character(EnvironmentalParameter)
  AlphaDiversityMetric <<- as.character(AlphaDiversity)
  
  CategoricalVariables <- c("grtgroup","biome_type","IUCN_CAT","ECO_NAME","HYBAS_ID")
  ContinuousVariables <- c("bio01","bio12","gHM","elevation","NDVI","Average_Radiance")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
  #Read in parsed Tronko database.
  TronkoInput <- system(paste("aws s3 cp s3://ednaexplorer/projects/",sample_ProjectID,"/",sample_Primer,"/Taxa_Parsed.tsv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  TronkoInput <- read.table(text = TronkoInput,header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  
  #Coerce data type.
  TronkoDB <- as.data.frame(TronkoInput)
  
  #Filter by maximum number of read mismatches.
  TronkoDB <- TronkoDB[TronkoDB$Mismatch <= sample_Num_Mismatch & !is.na(TronkoDB$Mismatch),]
  
  #Read in Metadata
  Metadata <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/Metadata.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  Metadata <- read.table(text = Metadata,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  #Set date type
  Metadata$sample_date <- as.Date(lubridate::ymd(Metadata$sample_date))
    
  #Create Phyloseq object
  #Create a taxonomic id matrix.
  taxmat <- TronkoDB[,c("Taxonomic_Path",TaxonomicRanks)]
  taxmat <- taxmat[!duplicated(taxmat),]
  rownames(taxmat) <- taxmat$Taxonomic_Path
  taxmat <- as.matrix(taxmat[,TaxonomicRanks])
  TAX = suppressWarnings(tax_table(taxmat))
  #Create OTU matrix
  otumat <- as.data.frame(pivot_wider(as.data.frame(table(TronkoDB[,c("SampleID","Taxonomic_Path")])), names_from = SampleID, values_from = Freq))
  rownames(otumat) <- otumat$Taxonomic_Path
  otumat <- otumat[,colnames(otumat) %in% unique(TronkoDB$SampleID)]
  otumat[sapply(otumat, is.character)] <- lapply(otumat[sapply(otumat, is.character)], as.numeric)
  OTU <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)
  #Create sample metadata matrix
  Sample <- Metadata[,colnames(Metadata) %in% c("sample_id","sample_date",CategoricalVariables,ContinuousVariables)]
  Sample <- Sample[!duplicated(Sample),]
  rownames(Sample) <- Sample$sample_id
  Sample$SampleID <- NULL
  Sample <- sample_data(Sample)
  #Create merged Phyloseq object.
  physeq <- phyloseq(OTU,TAX,Sample)
  
  #Filter eDNA data by date range and primers
  physeq <- subset_samples(physeq,sample_date >= sample_First_Date & sample_date <= sample_Last_Date)
  
  #Aggregate reads to a particular taxonomic level.
  physeq <- tax_glom(physeq,taxrank=sample_TaxonomicRank)
  
  #Filter on read counts.
  if(sum(sample_sums(physeq)>sample_CountThreshold)==0){
    sample_CountThreshold <<- min(sample_sums(physeq))
  }
  CountFilter <- prune_samples(sample_sums(physeq)>sample_CountThreshold, physeq)
  
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
