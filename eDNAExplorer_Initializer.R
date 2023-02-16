rm(list=ls())
require(tidyr)
require(dplyr)
require(sf)
requie(sp)
require(lubridate)

wd <- "~/Desktop/eDNAExplorer/Tronko/"

setwd(wd)
## This script should run automatically once a project's metadata and sequence data are uploaded.
#1. Parse Tronko output into a taxa by sample dataframe.  Calculate the traditional observation score between eDNA and GBIF data per taxon.
#2. Update taxon to icon database for web graphics.
##

#1. Parse Tronko output into a taxa by sample dataframe.  Calculate the traditional observation score between eDNA and GBIF data per taxon.
#Read in Tronko-assign output files.  Standardize sample IDs within them.
TronkoFiles <- list.files(pattern=".out.txt")
TronkoInputs <- data.frame()
for(TronkoFile in TronkoFiles){
  TronkoInput <- read.table(TronkoFile, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  TronkoInput$SampleID <- sub("_[^_]+$", "", TronkoFile)
  TronkoInput$Category <- "Tronko-assign"
  TronkoInputs <- rbind(TronkoInputs,TronkoInput)
}

#Get project primers.  Hard coded for now.
Primer <- "MiFish_12S_U"

#Retain sample ID and taxonomic ID.
#Create a read count matrix with rows as taxa and columns as samples.
TronkoDB <- as.data.frame(table(TronkoInputs[,c("SampleID","Taxonomic_Path")]))
TronkoDB <- pivot_wider(TronkoDB,names_from=SampleID,values_from=Freq)
names(TronkoDB)[names(TronkoDB) == 'Taxonomic_Path'] <- 'sum.taxonomy'

#Set taxonomic rank column names.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
#Select taxonomic rank to proceed with analysis.
TaxonomicRank <- "species" #User input

#Split taxonomy names into their component ranks.
TronkoDB <- suppressWarnings(tidyr::separate(TronkoDB,'sum.taxonomy',TaxonomicRanks,sep=";", extra="drop"))

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
#Save merged metadata.
write.table(Metadata,"Metadata.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Get local bounds for sample locations, add 0.5 degree buffer.
Local_East <- min(Metadata$longitude)-0.5
Local_West <- max(Metadata$longitude)+0.5
Local_South <- min(Metadata$latitude)-0.5
Local_North <- max(Metadata$latitude)+0.5
ylims <- c(Local_South,Local_North)
xlims <- c(Local_West,Local_East)
box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(4326))
#get the bounding box of the two x & y coordintates, make sfc
Local_Bounds <- st_bbox(box_coords) %>% st_as_sfc()

#Read in GBIF occurrences.
GBIF_Input <- suppressWarnings(read.table("GBIF.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))
GBIF_Input <- GBIF_Input[!is.na(GBIF_Input$decimalLatitude) & !is.na(GBIF_Input$decimalLongitude),]

GBIF_Input <- sf::st_as_sf(GBIF_Input,coords=c("decimalLongitude","decimalLatitude"),crs=4326)
#Read in the boundaries for ecoregions and realms
EcoBoundaries <- sf::st_read("Ecoregions2017.shp")
#Get realms present in sample data.
Realms <- EcoBoundaries[EcoBoundaries$REALM %in% unique(Metadata$REALM),]
#Get ecoregions present in sample data.
Eco_Regions <- EcoBoundaries[EcoBoundaries$ECO_NAME %in% unique(Metadata$ECO_NAME),]

#Clip GBIF occurrence locations by local boundaries.
Taxa_Local <- suppressWarnings(sf::st_intersection(GBIF_Input,Local_Bounds))
#Assign weight value for species, genera, and families.
Taxa_Local <- Taxa_Local %>% dplyr:: mutate(Local_GBIFWeight = dplyr::case_when(taxonRank=="SPECIES" ~ 4, taxonRank=="SUBSPECIES" ~ 4, taxonRank=="GENUS" ~ 2, taxonRank=="FAMILY" ~ 1, !(taxonRank %in% c("SPECIES","GENUS","FAMILY"))~0))

#Clip GBIF occurrence locations by realms.
Taxa_Realm <- suppressWarnings(sf::st_intersection(GBIF_Input,Realms[7:10,]))
#Assign weight value for species, genera, and families.
Taxa_Realm <- Taxa_Realm %>% dplyr:: mutate(Realm_GBIFWeight = dplyr::case_when(taxonRank=="SPECIES" ~ 4, taxonRank=="SUBSPECIES" ~ 4, taxonRank=="GENUS" ~ 2, taxonRank=="FAMILY" ~ 1, !(taxonRank %in% c("SPECIES","GENUS","FAMILY"))~0))

#Clip GBIF occurrence locations by ecoregions.
Taxa_Ecoregion <- suppressWarnings(sf::st_intersection(GBIF_Input,Eco_Regions))
#Assign weight value for species, genera, and families.
Taxa_Ecoregion <- Taxa_Ecoregion %>% dplyr:: mutate(Ecoregion_GBIFWeight = dplyr::case_when(taxonRank=="SPECIES" ~ 4, taxonRank=="SUBSPECIES" ~ 4, taxonRank=="GENUS" ~ 2, taxonRank=="FAMILY" ~ 1, !(taxonRank %in% c("SPECIES","GENUS","FAMILY"))~0))

#Count eDNA taxonomic resolution and weigh them.
#Species = 4, Genus = 2, Family = 1.  Everything else = 0.
TronkoDB$eDNAWeight <- 1*as.numeric(!is.na(TronkoDB$family))+2*as.numeric(!is.na(TronkoDB$genus))+4*as.numeric(!is.na(TronkoDB$species))

#Check if any of the eDNA reads show up in the local set of GBIF family observations.
TronkoDB$LocalFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Local$family))))
#Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
TronkoDB$EcoregionFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Ecoregion$family))))
#Check if any of the eDNA reads show up in the realm set of GBIF family observations.
TronkoDB$RealmFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Realm$family))))
#Check if any of the eDNA reads show up in the local set of GBIF family observations.
TronkoDB$LocalGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Local$genus))))
#Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
TronkoDB$EcoregionGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Ecoregion$genus))))
#Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
TronkoDB$RealmGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Realm$genus))))
#Check if any of the eDNA reads show up in the local set of GBIF family observations.
TronkoDB$LocalSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Local$species))))
#Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
TronkoDB$EcoregionSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Ecoregion$species))))
#Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
TronkoDB$RealmSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Realm$species))))

#Assign TOS scores for GBIF results.
TronkoDB$TOS_Local <- (1*TronkoDB$LocalFamilyPresentGBIF+2*TronkoDB$LocalGenusPresentGBIF+4*TronkoDB$LocalSpeciesPresentGBIF)/TronkoDB$eDNAWeight
TronkoDB$TOS_Ecoregion <- (1*TronkoDB$EcoregionFamilyPresentGBIF+2*TronkoDB$EcoregionGenusPresentGBIF+4*TronkoDB$EcoregionSpeciesPresentGBIF)/TronkoDB$eDNAWeight
TronkoDB$TOS_Realm <- (1*TronkoDB$RealmFamilyPresentGBIF+2*TronkoDB$RealmGenusPresentGBIF+4*TronkoDB$RealmSpeciesPresentGBIF)/TronkoDB$eDNAWeight
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
TronkoDB[is.nan(TronkoDB)] <- 0
#Add primer information.
TronkoDB$Primer <- Primer

#Remove entries where superkingdom is unassigned or unknown.
TronkoDB <- TronkoDB[TronkoDB$superkingdom!="unassigned" & !is.na(TronkoDB$superkingdom),]

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
#Export dataframe to file for downstream use.
write.table(TronkoDB,"Taxa_Parsed.txt",quote=FALSE,sep="\t",row.names = FALSE)

#2. Update taxon to icon database for web graphics.
#Get unique taxa list.
Unique_Taxa_df <- as.data.frame(unlist(TronkoDB[,TaxonomicRanks]))
colnames(Unique_Taxa_df) <- c("Taxon")
Unique_Taxa_df$TaxonomicRank <- rownames(Unique_Taxa_df)
Unique_Taxa_df$TaxonomicRank <- gsub('[[:digit:]]+', '', Unique_Taxa_df$TaxonomicRank)
Unique_Taxa_df <- Unique_Taxa_df[!duplicated(Unique_Taxa_df),]
Unique_Taxa_df <- Unique_Taxa_df[complete.cases(Unique_Taxa_df),]
Unique_Taxa <- unique(Unique_Taxa_df$Taxon)

#Get taxa which have already been added to the Phylopic database.
PhylopicDB_Initial <- read.table("PhylopicDB.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")

#Get taxa which are not yet in the Phylopic database.
Unique_Taxa <- Unique_Taxa[!(Unique_Taxa %in% PhylopicDB_Initial$Taxon)]

#Get unique Phylopic icons for each taxon.  If one does not exist, keep moving up a taxonomic level until an icon does exist.
PhylopicDB <- PhylopicDB_Initial
for(Unique_Taxon in Unique_Taxa){
  Taxon_GBIF <- name_backbone(Unique_Taxon,verbose=T,strict=F)
  Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
  if(nrow(Taxon_GBIF)>0){
    Taxon_Backbone <- na.omit(as.numeric(rev(unlist(Taxon_GBIF[1,unique(grep(paste(c("kingdom",TaxonomicRanks),collapse="|"),colnames(Taxon_GBIF[,grepl("Key",names(Taxon_GBIF))]), value=TRUE))]))))
    res <- httr::POST(url="https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true",body=jsonlite::toJSON(as.character(Taxon_Backbone), auto_unbox=TRUE),content_type("application/json"),encode="json")
    test <- fromJSON(rawToChar(res$content))
    Taxon_Image <- test[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
  } else{
    Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
  }
  tmp <- data.frame(matrix(nrow=1,ncol=2))
  colnames(tmp) <- c("Taxon","Image_URL")
  tmp$Taxon <- Unique_Taxon
  tmp$Image_URL <- Taxon_Image
  print(paste(Unique_Taxon,Taxon_Image))
  PhylopicDB <- rbind(PhylopicDB,tmp)
}
PhylopicDB <- dplyr::left_join(PhylopicDB,Unique_Taxa_df)
#Save Phylopic database to table.
write.table(PhylopicDB,"PhylopicDB.txt",quote=FALSE,sep="\t",row.names = FALSE)
