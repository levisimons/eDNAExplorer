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

#* Echo the parameter that was sent in
#* @param First_Date:string YYYY-MM-DD
#* @param Last_Date:string YYYY-MM-DD
#* @param Marker:string Target marker name
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @param EnvironmentalParameter:string Environmental variable to analyze against alpha diversity
#* @param BetaDiversity:string Beta diversity metric. Options are chao, bray, or jaccard.
#* @get /Tronko_Input
beta <- function(First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,BetaDiversity){
  
  #Define filters in Phyloseq as global parameters.
  sample_First_Date <<- ymd(First_Date)
  sample_Last_Date <<- ymd(Last_Date)
  sample_Primer <<- as.character(Marker)
  sample_TaxonomicRank <<- as.character(TaxonomicRank)
  sample_Num_Mismatch <<- as.numeric(Num_Mismatch)
  sample_CountThreshold <<- as.numeric(CountThreshold)
  sample_FilterThreshold <<- as.numeric(FilterThreshold)
  EnvironmentalVariable <<- as.character(EnvironmentalParameter)
  BetaDiversityMetric <<- as.character(BetaDiversity)
  
  CategoricalVariables <- c("grtgroup","biome_type","IUCN_CAT","ECO_NAME","HYBAS_ID")
  ContinuousVariables <- c("bio01","bio12","gHM","elevation","NDVI","Average_Radiance")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  
  #Read in parsed Tronko database.
  Tronko_Input <- get_object("test-taxa.tsv",bucket="ednaexplorer",region="",as="text",base_url="js2.jetstream-cloud.org:8001") %>% data.table::fread(encoding="UTF-8")
  
  #Coerce date type.
  TronkoDB <- as.data.frame(Tronko_Input)
  TronkoDB$sample_date <- as.Date(TronkoDB$sample_date)
  
  #Filter by maximum number of read mismatches.
  TronkoDB <- TronkoDB[TronkoDB$Mismatch <= sample_Num_Mismatch & !is.na(TronkoDB$Mismatch),]
  
  #Convert all character variables to being categorical. Default is that all numeric columns are continuous.
  TronkoDB[sapply(TronkoDB, is.character)] <- lapply(TronkoDB[sapply(TronkoDB, is.character)], as.factor)
  TronkoDB$SampleID <- as.character(TronkoDB$SampleID)
  
  #Create Phyloseq object
  #Create a taxonomic id matrix.
  taxmat <- TronkoDB[,c("usageKey",TaxonomicRanks)]
  taxmat <- taxmat[!duplicated(taxmat),]
  rownames(taxmat) <- taxmat$usageKey
  taxmat <- as.matrix(taxmat[,TaxonomicRanks])
  TAX = suppressWarnings(tax_table(taxmat))
  #Create OTU matrix
  otumat <- as.data.frame(pivot_wider(as.data.frame(table(TronkoDB[,c("SampleID","usageKey")])), names_from = SampleID, values_from = Freq))
  rownames(otumat) <- otumat$usageKey
  otumat <- otumat[,colnames(otumat) %in% unique(TronkoDB$SampleID)]
  otumat[sapply(otumat, is.character)] <- lapply(otumat[sapply(otumat, is.character)], as.numeric)
  OTU <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)
  #Create sample metadata matrix
  Sample <- TronkoDB[,colnames(TronkoDB) %in% c("SampleID","sample_date","Primer",CategoricalVariables,ContinuousVariables)]
  Sample <- Sample[!duplicated(Sample),]
  rownames(Sample) <- Sample$SampleID
  Sample$SampleID <- NULL
  Sample <- sample_data(Sample)
  #Create merged Phyloseq object.
  physeq <- phyloseq(OTU,TAX,Sample)
  
  #Filter eDNA data by date range and primers
  physeq <- subset_samples(physeq,sample_date >= sample_First_Date & sample_date <= sample_Last_Date & Primer==sample_Primer)
  
  #Aggregate reads to a particular taxonomic level.
  physeq <- tax_glom(physeq,taxrank=sample_TaxonomicRank)
  
  #Filter on read counts.
  if(sum(sample_sums(physeq)>sample_CountThreshold)==0){
    sample_CountThreshold <<- min(sample_sums(physeq))
  }
  CountFilter <- prune_samples(sample_sums(physeq)>sample_CountThreshold, physeq)
  
  #Filter on read abundance per sample.
  AbundanceFilter  = transform_sample_counts(CountFilter, function(x) x / sum(x) )
  AbundanceFiltered = filter_taxa(AbundanceFilter, function(x) sum(x) > sample_FilterThreshold, TRUE)
  
  #Plot and analyze beta diversity versus an environmental variables.
  if(nsamples(AbundanceFiltered)>1){
    BetaDist = phyloseq::distance(AbundanceFiltered, method=BetaDiversityMetric, weighted=F)
    ordination = ordinate(AbundanceFiltered, method="PCoA", distance=BetaDist)
    if(length(unique(sample_data(AbundanceFiltered)[,EnvironmentalVariable]))>1){
      BetaExpression = paste("adonis2(BetaDist ~ sample_data(AbundanceFiltered)$",EnvironmentalVariable,")",sep="")
      test <- eval(parse(text=BetaExpression))
      Stat_test <- paste("PERMANOVA results, using 999 permutations.\n",DiversityMetric," beta diversity and ",EnvironmentalVariable,"\nDegrees of freedom: ",round(test$Df[1],3),". Sum of squares: ",round(test$SumOfSqs[1],3),". R-squared: ",round(test$R2[1],3),". F-statistic: ",round(test$F[1],3),". p: ",round(test$`Pr(>F)`[1],3),sep="")
    } else {
      Stat_test <- "Not enough variation in environmental data to analyze against beta diversity."
    }
    p <- plot_ordination(AbundanceFiltered, ordination, color=EnvironmentalVariable,title=Stat_test) + theme(aspect.ratio=1)
    p <- p+theme_bw()
    #saveWidget(ggplotly(p, tooltip = c(EnvironmentalVariable,"Axis.1","Axis.2")),file="Beta_Diversity.html")
  } else {
    Stat_test <- "Not enough data to analyze beta diversity."
    p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=Stat_test)
    #saveWidget(ggplotly(p, tooltip = c(EnvironmentalVariable,"Axis.1","Axis.2")),file="Beta_Diversity.html")
  }
  #Save plot as json object
  jfig <- plotly_json(p, FALSE)
  return(jfig)
}
