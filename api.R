source("eDNAExplorer_Alpha.R")
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
function(First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity){
  alpha("test-taxa.tsv",First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,EnvironmentalParameter,AlphaDiversity)
}
