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
  TronkoInput <- system(paste("aws s3 cp s3://ednaexplorer/projects/",sample_ProjectID,"/",sample_Primer,"/Taxa_Parsed.tsv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  TronkoInput <- read.table(text = TronkoInput,header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  
  #Coerce data type.
  TronkoDB <- as.data.frame(TronkoInput)
  
  #Filter by maximum number of read mismatches.
  TronkoDB <- TronkoDB[TronkoDB$Mismatch <= sample_Num_Mismatch & !is.na(TronkoDB$Mismatch),]
  
  #Read in Metadata
  Metadata <- system(paste("aws s3 cp s3://ednaexplorer/projects/",sample_ProjectID,"/Metadata.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
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
  
  #Plot and analyze beta diversity versus an environmental variables.
  if(nsamples(AbundanceFiltered)>1){
    BetaDist = phyloseq::distance(AbundanceFiltered, method=BetaDiversityMetric, weighted=F)
    if(sum(!is.nan(BetaDist))>1){
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
