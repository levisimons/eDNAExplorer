source("eDNAExplorer_Alpha.R")
#' @param Object_Name
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
function(Object_Name,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity){
  alpha(Object_Name,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity)
}
