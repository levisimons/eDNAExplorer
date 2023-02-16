rm(list=ls())
require(tidyr)
require(dplyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(htmlwidgets)
require(plotly)

wd <- "~/Desktop/eDNAExplorer/Tronko/"

setwd(wd)

#Choose a threshold for filtering ASVs prior to analysis. The options are:
#1. No filtering.
#2. Filtering on ASVs with a read count of at least 0.003% the project total found in a sample: https://doi.org/10.3389/fenvs.2017.00011
#3. Filtering on ASVs with a read count of at least 0.01% the project total found in a sample: https://doi.org/10.1111/2041-210X.12849
#User input.
FilterThreshold <- 0.00003 
#Use total reads per sample histogram to decide on a read count threshold for retaining samples.
#User input
CountThreshold <- 1000

#Read in parsed Tronko-assign output files.  Generated from eDNAExplorer_Initializer.R.
TronkoDB <- suppressWarnings(read.table("Taxa_Parsed.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))

#Set taxonomic rank column names.  Select taxonomic level to aggregate data on.
TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
TaxonomicRank <- "class" #User input to parse Venn diagrams on.

#Read in complete metadata.
Metadata <- suppressWarnings(read.table("Metadata.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))
Metadata$sample_date <- as.Date(Metadata$sample_date)

#Convert all character variables to being categorical. Default is that all numeric columns are continuous.
Metadata[sapply(Metadata, is.character)] <- lapply(Metadata[sapply(Metadata, is.character)], as.factor)
Metadata$sample_id <- as.character(Metadata$sample_id)

#Select date range to filter data on.
First_Date <- as.Date("2021-06-12") #User input
Last_Date <- as.Date("2021-09-12") #User input

#Select primer(s) to filter data on.
Primers <- c("MiFish_12S_U") #User input
#Filter eDNA data by primer
TronkoDB <- TronkoDB[TronkoDB$Primer %in% Primers,]

#Designate a subset of categorical variables for users to select for plotting.
#grtgroup = Soil taxonomy great groups, a classification number from 0 to 433
#biome_type = Potential Natural Vegetation biomes global predictions of classes (based on predictions using the BIOMES 6000 dataset's 'current biomes' category.)
#IUCN_CAT = IUCN management category, one of: Ia (strict nature reserve), Ib (wilderness area), II (national park), III (natural monument or feature), IV (habitat/species management area), V (protected landscape/seascape), VI (PA with sustainable use of natural resources), not applicable, not assigned, or not reported.
#ECO_NAME = one of the 846 terrestrial ecoregions that represent our living planet
#HYBAS_ID = watershed ID number.
CategoricalVariables <- c("grtgroup","biome_type","IUCN_CAT","ECO_NAME","HYBAS_ID")
#Designate a subset of continuous variables for users to select for plotting.
#bio01 = Annual Mean Temperature, bio12 = Annual Precipitation
#gHM = global human modification index
#NDVI = normalized difference vegetation index, a proxy for photosynthetic activity.
#Average Radiance = Monthly average radiance composite images using nighttime data from the Visible Infrared Imaging Radiometer Suite (VIIRS) Day/Night Band
ContinuousVariables <- c("bio01","bio12","gHM","elevation","NDVI","Average_Radiance")

#Create a taxonomic id matrix.
taxmat <- as.matrix(TronkoDB[,TaxonomicRanks])

#Create OTU matrix
otumat <- as.matrix(TronkoDB[,Metadata$sample_id])

#Create phyloseq object for taxa and samples.
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)

#Create sample metadata matrix and merge into phyloseq object.
Sample <- Metadata[,colnames(Metadata) %in% c(CategoricalVariables,ContinuousVariables)]
rownames(Sample) <- Metadata$sample_id
Sample <- sample_data(Sample)
physeq <- merge_phyloseq(physeq,Sample)

#Filter eDNA data by date range.
physeq <- subset_samples(physeq,Metadata$sample_date >= First_Date & Metadata$sample_date <= Last_Date)

#Filter phyloseq object by relative abundance.
AbundanceFilter <- genefilter_sample(physeq, function(x) x/sum(x) > FilterThreshold)
AbundanceFiltered <- prune_taxa(AbundanceFilter,physeq)
#Filter phyloseq object by the number of remaining reads per sample.
#If threshold is too low, force it to lowest sequence per sample level.
if(sum(sample_sums(AbundanceFiltered)>CountThreshold)==0){
  CountThreshold <- min(sample_sums(AbundanceFiltered))
}
CountFilter <- prune_samples(sample_sums(AbundanceFiltered)>CountThreshold, AbundanceFiltered)

#Plot and analyze alpha diversity versus environmental variation.
#Enter an environmental variable to compare alpha diversity against
#User input
EnvironmentalVariable <- "Average_Radiance"
#Enter an alpha diversity metric.
#User input
DiversityMetrics <-  c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
DiversityMetric <- "Shannon"
if(EnvironmentalVariable %in% CategoricalVariables){
  #Store diversity versus variable data.
  tmp <- ggplot_build(plot_richness(CountFilter,x=EnvironmentalVariable,measures=DiversityMetric))
  tmp <- tmp$data[[1]]
  tmp$x <- as.factor(tmp$x)
  #General a violin plot of alpha diversity versus an environmental variable.
  p <- ggplot(tmp, aes(x=x, y=y))+
    labs(title=paste(DiversityMetric," versus ",gsub("_"," ",EnvironmentalVariable),".  Samples collected between: ",First_Date," and ",Last_Date,"\nRelative abundance minimum of ",100*FilterThreshold,"%. Reads per sample minimum: ",CountThreshold,sep=""),x=gsub("_"," ",EnvironmentalVariable), y = DiversityMetric)+
    geom_violin()+theme_bw()+geom_point(position = position_jitter(seed = 1, width = 0.2))
  #Run a Kruskal-Wallis test between alpha diversity and selected environmental variable.
  test <- kruskal.test(tmp$y ~ tmp$x, data = tmp)
  print(paste("chi-squared = ",round(test$statistic,digits=3)," df = ",test$parameter," p-value = ",round(test$p.value,digits=3)))
}
if(EnvironmentalVariable %in% ContinuousVariables){
  #Store diversity versus variable data.
  tmp <- ggplot_build(plot_richness(CountFilter,x=EnvironmentalVariable,measures=DiversityMetric))
  tmp <- tmp$data[[1]]
  tmp$x <- as.numeric(tmp$x)
  #Generate a scatterplot of alpha diversity versus an environmental variable.
  p <- ggplot(tmp, aes(x=x, y=y))+
    labs(title=paste(DiversityMetric," versus ",gsub("_"," ",EnvironmentalVariable),".  Samples collected between: ",First_Date," and ",Last_Date,"\nRelative abundance minimum of ",100*FilterThreshold,"%. Reads per sample minimum: ",CountThreshold,sep=""),x=gsub("_"," ",EnvironmentalVariable), y = DiversityMetric)+
    theme_bw()+geom_point()+geom_smooth()
  #Calculate summary statistics using a Kendall correlation test between alpha diversity and environmental variable.
  test <- cor.test(tmp$x,tmp$y,alternative="two.sided",method="kendall")
  print(paste("z = ",round(test$statistic,digits=3)," tau = ",round(test$estimate,digits=3)," p-value = ",round(test$p.value,digits=3)))
}
saveWidget(ggplotly(p, tooltip = "all"),file="Alpha_Diversity.html")
