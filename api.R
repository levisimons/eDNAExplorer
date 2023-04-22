source("eDNAExplorer_Alpha.R")
source("eDNAExplorer_Beta.R")
source("eDNAExplorer_Prevalence.R")
source("eDNAExplorer_Venn.R")
source("eDNAExplorer_Timeline.R")
source("eDNAExplorer_Map.R")
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
function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity,BetaDiversity,Geographic_Scale){
  alpha(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity)
  beta(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,BetaDiversity)
  prevalence(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold)
  venn(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,Geographic_Scale)
  timeline(Marker,Taxon_name,TaxonomicRank,Num_Mismatch,CountThreshold,FilterThreshold)
  map(Marker,Taxon_name,TaxonomicRank,Num_Mismatch,CountThreshold,FilterThreshold)
}
