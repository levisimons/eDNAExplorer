# plumber.R
library(plumber)
require(aws.s3)
require(tidyr)
require(dplyr)
require(ggplot2)
require(rgbif)
require(gbifdb)
require(plotly)
require(jsonlite)

#* Echo the parameter that was sent in
#* @param Taxon_name:string Scientific taxon name
#* @get /timeline

timeline <- function(Taxon_name){
  
  #Select taxon to map.
  #User input
  Taxon <- Taxon_name
  #Get GBIF taxonomy key for taxon.
  Taxon_GBIF <- name_backbone(Taxon)$usageKey
  
  #Read in GBIF occurrences.
  gbif <- gbif_local()
  
  #Filter GBIF occurrences to a particular taxon.
  TaxonDB <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                occurrencestatus=="PRESENT",taxonkey==Taxon_GBIF) %>% select(year)
  TaxonDB <- as.data.frame(TaxonDB)
  
  #Get GBIF taxon occurrences over time.
  Taxa_Time <- as.data.frame(table(TaxonDB))
  colnames(Taxa_Time) <- c("year","Occurrences")
  Taxa_Time <- toJSON(Taxa_Time)
  return(Taxa_Time)
}
