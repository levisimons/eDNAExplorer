# plumber.R
library(plumber)
require(aws.s3)
require(tidyr)
require(dplyr)
require(rgbif)
require(jsonlite)
require(gbifdb)

#* Echo the parameter that was sent in
#* @param Taxon_name:string Scientific taxon name
#* @get /Tronko_Input

GBIF_Map <- function(Taxon_name){
  #Select taxon to map.
  #User input
  Taxon <- Taxon_name
  #Get GBIF taxonomy key for taxon.
  Taxon_GBIF <- name_backbone(Taxon)$usageKey
  
  #Read in parsed Tronko database.
  Tronko_Input <- get_object("test-taxa.tsv",bucket="ednaexplorer",region="",as="text",base_url="js2.jetstream-cloud.org:8001") %>% data.table::fread(encoding="UTF-8")
  
  #Coerce date type.
  TronkoDB <- as.data.frame(Tronko_Input)
  
  #Taxonomic rank column names
  TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
  
  #Find where taxon occurs in Tronko output.
  TaxonDB <- dplyr::filter_all(TronkoDB[,c(TaxonomicRanks,"longitude","latitude")], any_vars(. ==Taxon))
  
  #Get taxon locations by eDNA.
  Taxon_Locations <- TaxonDB[,c("longitude","latitude")]
  Taxon_Locations <- Taxon_Locations[!duplicated(Taxon_Locations),]
  Taxon_Locations$Sample_Source <- "eDNA"
  
  #Read in GBIF occurrences.
  gbif <- gbif_local()
  
  #Filter GBIF occurrence locations by taxon.
  Taxon_Map <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                               coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                               occurrencestatus=="PRESENT",
                               taxonkey==Taxon_GBIF) %>% select(decimallongitude,decimallatitude)
  
  Taxon_Map <- as.data.frame(Taxon_Map)
  colnames(Taxon_Map) <- c("longitude","latitude")
  Taxon_Map$Sample_Source <- "GBIF"
  
  #Merge eDNA and GBIF taxon locations
  Taxon_Map <- rbind(Taxon_Locations,Taxon_Map)
  print(Taxon_Map)
  #Return results
  Taxon_Map_Data <- jsonlite::toJSON(Taxon_Map)
  return(Taxon_Map_Data)
}
