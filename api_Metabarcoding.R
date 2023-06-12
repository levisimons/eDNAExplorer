source("eDNAExplorer_Alpha_Metabarcoding.R")
source("eDNAExplorer_Beta_Metabarcoding.R")
source("eDNAExplorer_Prevalence_Metabarcoding.R")
source("eDNAExplorer_Venn_Metabarcoding.R")
source("eDNAExplorer_Timeline_Metabarcoding.R")
source("eDNAExplorer_Map_Metabarcoding.R")
#' @param ProjectID
#' @param First_Date
#' @param Last_Date
#' @param Marker
#' @param Num_Mismatch
#' @param TaxonomicRank
#' @param CountThreshold
#' @param FilterThreshold
#' @param EnvironmentalParameter
#' @param AlphaDiversity
#' @param BetaDiversity
#' @param Geographic_Scale
#' @param Taxon_name
#' @get /alpha
#' @get /beta
#' @get /prevalence
#' @get /venn
#' @get /timeline
#' @get /map
function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity,BetaDiversity,Geographic_Scale,SpeciesList){
  alpha(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity,SpeciesList)
  beta(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,BetaDiversity,SpeciesList)
  prevalence(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,SpeciesList)
  venn(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,Geographic_Scale,SpeciesList)
  timeline(Marker,Taxon_name,TaxonomicRank,Num_Mismatch,CountThreshold,FilterThreshold)
  map(Marker,Taxon_name,TaxonomicRank,Num_Mismatch,CountThreshold,FilterThreshold)
}
