source("eDNAExplorer_Alpha.R")
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
#' @get /alpha
function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity){
  alpha(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity)
}
