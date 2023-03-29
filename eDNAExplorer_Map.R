# plumber.R
library(plumber)
require(tidyr)
require(dplyr)
require(rgbif)
require(jsonlite)
require(gbifdb)

#* Echo the parameter that was sent in
#* @param Taxon_name:string Scientific taxon name
#* @get /map

map <- function(Taxon_name){
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
                             occurrencestatus=="PRESENT",taxonkey==Taxon_GBIF) %>% select(decimallongitude,decimallatitude)
  TaxonDB <- as.data.frame(TaxonDB)
  #Get unique taxon locations
  TaxonDB <- TaxonDB[!duplicated(TaxonDB),]
  TaxonDB <- TaxonDB[complete.cases(TaxonDB),]
  colnames(TaxonDB) <- c("longitude","latitude")
  TaxonDB$Source <- "GBIF"
  
  #Return results
  Taxon_Map_Data <- jsonlite::toJSON(TaxonDB)
  return(Taxon_Map_Data)
}
