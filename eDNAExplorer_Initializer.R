require(tidyr)
require(dplyr)
require(sf)
require(sp)
require(lubridate)
require(httr)
require(gbifdb)
require(jsonlite)
require(rgbif)
require(duckdb)

## This script should run automatically once a project's metadata and sequence data are uploaded.
#1. Parse Tronko output into a taxa by sample dataframe.
#Calculate the traditional observation score between eDNA and GBIF data per taxon. Merge in environmental metadata.
#2. Update taxon to icon database for web graphics.
##

#Read in Tronko-assign output files.  Standardize sample IDs within them.
TronkoFiles <- list.files(pattern=".out.txt")
TronkoInputs <- data.frame()
for(TronkoFile in TronkoFiles){
  TronkoInput <- read.table(TronkoFile, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  TronkoInput$SampleID <- sub("_[^_]+$", "", TronkoFile)
  TronkoInput$Category <- "Tronko-assign"
  TronkoInputs <- rbind(TronkoInputs,TronkoInput)
}
TronkoDB <- TronkoInputs

#Get project primers.  Hard coded for now.
Primer <- "MiFish_12S_U"
TronkoDB$Primer <- Primer

#Standardize taxonomy naming schema.
names(TronkoDB)[names(TronkoDB) == 'Taxonomic_Path'] <- 'sum.taxonomy'

#Set taxonomic rank column names.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")

#Split taxonomy names into their component ranks.
TronkoDB <- suppressWarnings(tidyr::separate(TronkoDB,'sum.taxonomy',TaxonomicRanks,sep=";", extra="drop"))

#Assign total number of mismatches.
TronkoDB$Mismatch <- TronkoDB$Forward_Mismatch+TronkoDB$Reverse_Mismatch

#Read in initial metadata.
Metadata_Initial <- suppressWarnings(read.table("InputMetadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))
#Get field variables from initial metadata.
Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% c("sample_id","longitude","latitude","sample_date","spatial_uncertainty"))]
#Read in extracted metadata.
Metadata_Extracted <- suppressWarnings(read.table("MetadataOutput.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))
names(Metadata_Extracted)[names(Metadata_Extracted) == 'name'] <- 'sample_id'
#Merge metadata
Metadata <- dplyr::left_join(Metadata_Initial[,c("sample_id",Field_Variables)],Metadata_Extracted)
#Clean up date format.
Metadata$sample_date <- as.Date(lubridate::ymd_hms(Metadata$sample_date))

#Merge sample metadata with parsed Tronko output.
TronkoDB <- dplyr::left_join(TronkoDB,Metadata,by=c("SampleID" = "sample_id"))

#Get kingdom data for phyla from GBIF, if available.
Phylum_to_Kingdom <- data.frame()
for(phylum in na.omit(unique(TronkoDB$phylum))){
  Taxon_GBIF <- name_backbone(phylum,verbose=T,strict=F)
  Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
  tmp <- data.frame(matrix(nrow=1,ncol=2))
  colnames(tmp) <- c("phylum","kingdom")
  tmp$phylum <- phylum
  if(nrow(Taxon_GBIF)>0){
    tmp$kingdom <- Taxon_GBIF$kingdom
  } else {
    tmp$kingdom <- NA
  }
  Phylum_to_Kingdom <- rbind(Phylum_to_Kingdom,tmp)
}
TronkoDB <- dplyr::left_join(TronkoDB,Phylum_to_Kingdom)

#Get taxa which have already been added to the Phylopic database.
if(file.exists("PhylopicDB.txt")==TRUE){
  PhylopicDB_Initial <- read.table("PhylopicDB.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
} else {
  PhylopicDB_Initial <- data.frame()
}
#Get unique taxa list.
Unique_Taxa_df <- as.data.frame(unlist(TronkoDB[,TaxonomicRanks]))
colnames(Unique_Taxa_df) <- c("Taxon")
Unique_Taxa_df$TaxonomicRankImage <- rownames(Unique_Taxa_df)
Unique_Taxa_df$TaxonomicRankImage <- gsub('[[:digit:]]+', '', Unique_Taxa_df$TaxonomicRankImage)
Unique_Taxa_df <- Unique_Taxa_df[!duplicated(Unique_Taxa_df),]
Unique_Taxa_df <- Unique_Taxa_df[complete.cases(Unique_Taxa_df),]
Unique_Taxa_df <- Unique_Taxa_df[Unique_Taxa_df$Taxon!="unassigned",]
Unique_Taxa <- unique(Unique_Taxa_df$Taxon)

#Get unique Phylopic icons for each taxon.  If one does not exist, keep moving up a taxonomic level until an icon does exist.
PhylopicDB <- PhylopicDB_Initial
for(Unique_Taxon in Unique_Taxa){
  Taxon_GBIF <- name_backbone(Unique_Taxon,verbose=T,strict=F)
  Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
  if(nrow(Taxon_GBIF)>0){
    Taxon_Backbone <- na.omit(as.numeric(rev(unlist(Taxon_GBIF[1,unique(grep(paste(c("kingdom",TaxonomicRanks),collapse="|"),colnames(Taxon_GBIF[,grepl("Key",names(Taxon_GBIF))]), value=TRUE))]))))
    Taxon_Keys <- as.data.frame(Taxon_GBIF[1,!(colnames(Taxon_GBIF) %in% c("kingdom",TaxonomicRanks))])
    res <- httr::POST(url="https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true",body=jsonlite::toJSON(as.character(Taxon_Backbone), auto_unbox=TRUE),content_type("application/json"),encode="json")
    test <- fromJSON(rawToChar(res$content))
    Taxon_Image <- test[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
  } else{
    Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
    Taxon_Keys <- data.frame(matrix(nrow=1,ncol=1))
    colnames(Taxon_Keys) <- c("usageKey")
    Taxon_Keys$usageKey <- 0
  }
  tmp <- data.frame(matrix(nrow=1,ncol=2))
  colnames(tmp) <- c("Taxon","Image_URL")
  tmp$Taxon <- Unique_Taxon
  tmp$Image_URL <- Taxon_Image
  tmp <- cbind(tmp,Taxon_Keys)
  PhylopicDB <- dplyr::bind_rows(PhylopicDB,tmp)
}
#Save Phylopic database to table.
write.table(PhylopicDB,"PhylopicDB.txt",quote=FALSE,sep="\t",row.names = FALSE)
PhylopicDB <- dplyr::left_join(PhylopicDB,Unique_Taxa_df)
#Get common names and merge them into Phylopic database.
Common_Names <- data.frame()
for(key in unique(PhylopicDB$usageKey)){
  Common_Name <- as.data.frame(name_usage(key=key, data="vernacularNames")$data)
  if(nrow(Common_Name)>0){
    Common_Name <- Common_Name[Common_Name$language=="eng",]
    Common_Name <- Common_Name[1,"vernacularName"]
  } else {
    Common_Name <- NA
  }
  tmp <- data.frame(matrix(nrow=1,ncol=2))
  colnames(tmp) <- c("usageKey","Common_Name")
  tmp$usageKey <- key
  tmp$Common_Name <- Common_Name
  Common_Names <- rbind(Common_Names,tmp)
}
PhylopicDB <- dplyr::left_join(PhylopicDB,Common_Names)

#Read in GBIF occurrences.
gbif <- gbif_local()

#Read in state/province boundaries.
#Boundaries are from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
sf_use_s2(FALSE)
GADM_1_Boundaries <- sf::st_read("ne_10m_admin_1_states_provinces.shp")
#Determine the unique list of national and state/proving boundaries sample locations cover.
country_list <- st_join(st_as_sf(Metadata, coords = c("longitude", "latitude"), crs = 4326), GADM_1_Boundaries['iso_a2'],join = st_intersects)
country_list <- na.omit(unique(country_list$iso_a2))
state_province_list <- st_join(st_as_sf(Metadata, coords = c("longitude", "latitude"), crs = 4326), GADM_1_Boundaries['woe_name'],join = st_intersects)
state_province_list <- na.omit(unique(state_province_list$woe_name))

#Get local bounds for sample locations, add 0.5 degree buffer.
Local_East <- max(Metadata$longitude)+0.5
Local_West <- min(Metadata$longitude)-0.5
Local_South <- min(Metadata$latitude)-0.5
Local_North <- max(Metadata$latitude)+0.5

#Clip GBIF occurrence locations by local boundaries.
Taxa_Local <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                              coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                              occurrencestatus=="PRESENT",
                              decimallongitude >= Local_West & decimallongitude <= Local_East & decimallatitude >= Local_South & decimallatitude <= Local_North) %>% 
  select(decimallongitude,decimallatitude,taxonkey,taxonrank,species,genus,family,order,class,phylum,kingdom)
Taxa_Local <- as.data.frame(Taxa_Local)
#Assign weight value for species, genera, and families.
Taxa_Local <- Taxa_Local %>% dplyr:: mutate(Local_GBIFWeight = dplyr::case_when(taxonrank=="SPECIES" ~ 4, taxonrank=="SUBSPECIES" ~ 4, taxonrank=="GENUS" ~ 2, taxonrank=="FAMILY" ~ 1, !(taxonrank %in% c("SPECIES","GENUS","FAMILY"))~0))

#Clip GBIF occurrence locations by state/province boundaries.
Taxa_State <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                              coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                              occurrencestatus=="PRESENT", stateprovince %in% state_province_list) %>% select(decimallongitude,decimallatitude,taxonkey,taxonrank,species,genus,family,order,class,phylum,kingdom)
Taxa_State <- as.data.frame(Taxa_State)
#Assign weight value for species, genera, and families.
Taxa_State <- Taxa_State %>% dplyr:: mutate(State_GBIFWeight = dplyr::case_when(taxonrank=="SPECIES" ~ 4, taxonrank=="SUBSPECIES" ~ 4, taxonrank=="GENUS" ~ 2, taxonrank=="FAMILY" ~ 1, !(taxonrank %in% c("SPECIES","GENUS","FAMILY"))~0))

#Clip GBIF occurrence locations by national boundaries.
Taxa_Nation <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                              coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                              occurrencestatus=="PRESENT", countrycode %in% country_list) %>% select(decimallongitude,decimallatitude,taxonkey,taxonrank,species,genus,family,order,class,phylum,kingdom)
Taxa_Nation <- as.data.frame(Taxa_Nation)
#Assign weight value for species, genera, and families.
Taxa_Nation <- Taxa_Nation %>% dplyr:: mutate(Ecoregion_GBIFWeight = dplyr::case_when(taxonrank=="SPECIES" ~ 4, taxonrank=="SUBSPECIES" ~ 4, taxonrank=="GENUS" ~ 2, taxonrank=="FAMILY" ~ 1, !(taxonrank %in% c("SPECIES","GENUS","FAMILY"))~0))

#Count eDNA taxonomic resolution and weigh them.
#Species = 4, Genus = 2, Family = 1.  Everything else = 0.
TronkoDB$eDNAWeight <- 1*as.numeric(!is.na(TronkoDB$family))+2*as.numeric(!is.na(TronkoDB$genus))+4*as.numeric(!is.na(TronkoDB$species))

#Check if any of the eDNA reads show up in the local set of GBIF family observations.
TronkoDB$LocalFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Local$family))))
#Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
TronkoDB$EcoregionFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_State$family))))
#Check if any of the eDNA reads show up in the realm set of GBIF family observations.
TronkoDB$RealmFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Nation$family))))
#Check if any of the eDNA reads show up in the local set of GBIF family observations.
TronkoDB$LocalGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Local$genus))))
#Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
TronkoDB$EcoregionGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_State$genus))))
#Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
TronkoDB$RealmGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Nation$genus))))
#Check if any of the eDNA reads show up in the local set of GBIF family observations.
TronkoDB$LocalSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Local$species))))
#Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
TronkoDB$EcoregionSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_State$species))))
#Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
TronkoDB$RealmSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Nation$species))))

#Assign TOS scores for GBIF results.
TronkoDB$TOS_Local <- (1*TronkoDB$LocalFamilyPresentGBIF+2*TronkoDB$LocalGenusPresentGBIF+4*TronkoDB$LocalSpeciesPresentGBIF)/TronkoDB$eDNAWeight
TronkoDB$TOS_Ecoregion <- (1*TronkoDB$EcoregionFamilyPresentGBIF+2*TronkoDB$EcoregionGenusPresentGBIF+4*TronkoDB$EcoregionSpeciesPresentGBIF)/TronkoDB$eDNAWeight
TronkoDB$TOS_Realm <- (1*TronkoDB$RealmFamilyPresentGBIF+2*TronkoDB$RealmGenusPresentGBIF+4*TronkoDB$RealmSpeciesPresentGBIF)/TronkoDB$eDNAWeight
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
TronkoDB[is.nan(TronkoDB)] <- 0

#Merge phylopic image links into ecological dataframe.
#Update taxonomic rank column names.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
tmp <- TronkoDB[,c("kingdom",TaxonomicRanks)]

#Get most finely resolved taxon per row.
TronkoDB$LeafTaxa <- tmp[cbind(1:nrow(tmp), max.col(!is.na(tmp), ties.method = 'last'))]

#Merge in Phylopic urls corresponding to most finely resolved taxon in each row.
TronkoDB <- dplyr::left_join(TronkoDB,PhylopicDB,by=c("LeafTaxa"="Taxon"))

#Export dataframe to file for downstream use.
write.table(TronkoDB,"Taxa_Parsed.txt",quote=FALSE,sep="\t",row.names = FALSE)
