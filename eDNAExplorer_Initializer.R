rm(list=ls())
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
require(aws.s3)
require(data.table)

## This script should run automatically once a project's metadata and sequence data are uploaded.
#1. Parse Tronko output into a taxa by sample dataframe.
#Calculate the traditional observation score between eDNA and GBIF data per taxon. Merge in environmental metadata.
#2. Update taxon to icon database for web graphics.
##

ProjectID <- "LARiverRound1" #This is hard-coded for now.

#Read in initial metadata.
Metadata_Initial <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/InputMetadata.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Initial <- read.table(text = Metadata_Initial,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
Metadata_Initial$sample_date <- lubridate::ymd(Metadata_Initial$sample_date)
#Get field variables from initial metadata.
Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% c("sample_id","longitude","latitude","sample_date","spatial_uncertainty"))]
#Read in extracted metadata.
Metadata_Extracted <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/MetadataOutput.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Extracted <- read.table(text = Metadata_Extracted,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
Metadata_Extracted$sample_date <- lubridate::ymd_hms(Metadata_Extracted$sample_date)
names(Metadata_Extracted)[names(Metadata_Extracted) == 'name'] <- 'sample_id'
#Merge metadata
Metadata <- dplyr::left_join(Metadata_Initial[,c("sample_id","sample_date","latitude","longitude","spatial_uncertainty",Field_Variables)],Metadata_Extracted)
#Clean up date format.
Metadata$sample_date <- as.Date(lubridate::ymd(Metadata$sample_date))

#Read in state/province boundaries.
#Boundaries are from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
sf_use_s2(FALSE)
SpatialBucket <- system("aws s3 ls s3://ednaexplorer/spatial --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",intern=TRUE)
SpatialBucket <- read.table(text = paste(SpatialBucket,sep = ""),header = FALSE)
colnames(SpatialBucket) <- c("Date", "Time", "Size","Filename")
SpatialFiles <- unique(SpatialBucket$Filename)
for(SpatialFile in SpatialFiles){
  system(paste("aws s3 cp s3://ednaexplorer/",SpatialFile," . --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
}
GADM_1_Boundaries <- sf::st_read("ne_10m_admin_1_states_provinces.shp")
#Determine the unique list of national and state/proving boundaries sample locations cover.
GADM_Boundaries <- st_join(st_as_sf(Metadata[!is.na(Metadata$latitude) & !is.na(Metadata$latitude),], coords = c("longitude", "latitude"), crs = 4326), GADM_1_Boundaries[c('iso_a2','woe_name')],join = st_intersects)
GADM_Boundaries <- GADM_Boundaries %>% st_drop_geometry()
GADM_Boundaries <- as.data.frame(GADM_Boundaries[,c("sample_id","sample_date","iso_a2","woe_name")])
names(GADM_Boundaries)[names(GADM_Boundaries) == "woe_name"] <- "State"
names(GADM_Boundaries)[names(GADM_Boundaries) == "iso_a2"] <- "Nation"
names(GADM_Boundaries)[names(GADM_Boundaries) == "sample_id"] <- "SampleID"
Metadata <- dplyr::left_join(Metadata,GADM_Boundaries,by=c("sample_id"="SampleID"))
country_list <- na.omit(unique(GADM_Boundaries$Nation))
state_province_list <- na.omit(unique(GADM_Boundaries$State))

#Save metadata.
write.table(Metadata,"Metadata.csv",quote=FALSE,sep=",",row.names = FALSE)
system(paste("aws s3 cp Metadata.csv s3://ednaexplorer/projects/",ProjectID,"/Metadata.csv --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))

#Read in GBIF occurrences.
gbif <- gbif_local()

#Get local bounds for sample locations, add 0.5 degree buffer.
Local_East <- max(na.omit(Metadata$longitude))+0.5
Local_West <- min(na.omit(Metadata$longitude))-0.5
Local_South <- min(na.omit(Metadata$latitude))-0.5
Local_North <- max(na.omit(Metadata$latitude))+0.5

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

#Get primers
Markers <- setdiff(colnames(Metadata[,grepl("Marker_",colnames(Metadata))]),colnames(Metadata[,grepl("Marker_(.*?)_",colnames(Metadata))]))
Primers <- unique(unlist(Metadata[,Markers]))
j=1
TronkoProject <- list()
for(Primer in Primers){
  #Read in Tronko-assign output files.  Standardize sample IDs within them.
  TronkoBucket <- system(paste("aws s3 ls s3://ednaexplorer/projects/",ProjectID," --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  TronkoBucket <- read.table(text = paste(TronkoBucket,sep = ""),header = FALSE)
  colnames(TronkoBucket) <- c("Date", "Time", "Size","Filename")
  TronkoFiles <- unique(TronkoBucket$Filename)
  TronkoFiles <- TronkoFiles[grepl(paste("projects",ProjectID,Primer,sep="/"),TronkoFiles)]
  TronkoFiles <- TronkoFiles[grepl("*.txt$",TronkoFiles)]
  TronkoInputs <- list()
  i=1
  TronkoHeaders <- c("Readname","Taxonomic_Path","Score","Forward_Mismatch","Reverse_Mismatch","Tree_Number","Node_Number")
  for(TronkoFile in TronkoFiles){
    TronkoInput <- system(paste("aws s3 cp s3://ednaexplorer/",TronkoFile," - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    if(length(TronkoInput)>0){
      TronkoInput <- read.table(text = paste(TronkoInput,sep = "\t"),header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
      if(nrow(TronkoInput)>0){
        print(paste(Primer,i,length(TronkoFiles)))
        if(identical(colnames(TronkoInput),TronkoHeaders)==FALSE){
          header_row <- as.data.frame(t(colnames(TronkoInput)))
          colnames(header_row) <- TronkoHeaders
          colnames(TronkoInput) <- TronkoHeaders
        }
        TronkoInput$SampleID <- gsub("\\..*","",basename(TronkoFile))
        TronkoInputs[[i]] <- TronkoInput
        i=i+1
      }
    }
  }
  TronkoInputs <- rbindlist(TronkoInputs, use.names=TRUE, fill=TRUE)
  
  TronkoDB <- as.data.frame(TronkoInputs)
  
  #Standardize taxonomy naming schema.
  TronkoDB$sum.taxonomy <- TronkoDB$Taxonomic_Path
  
  TronkoDB$ProjectID <- ProjectID
  
  #Set taxonomic rank column names.
  TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
  
  #Split taxonomy names into their component ranks.
  TronkoDB <- suppressWarnings(tidyr::separate(TronkoDB,'sum.taxonomy',TaxonomicRanks,sep=";", extra="drop"))
  
  TronkoDB[TronkoDB=="NA"] <- NA
  TronkoDB[TronkoDB==""] <- NA
  
  #Assign total number of mismatches.
  TronkoDB$Forward_Mismatch <- as.numeric(TronkoDB$Forward_Mismatch)
  TronkoDB$Reverse_Mismatch <- as.numeric(TronkoDB$Reverse_Mismatch)
  TronkoDB <- TronkoDB %>% mutate(Mismatch = rowSums(select(., Forward_Mismatch, Reverse_Mismatch), na.rm = TRUE))
  print(paste(Primer,"Mismatches calculated"))
  #Get kingdom data for phyla from GBIF, if available.
  Phylum_to_Kingdom <- list()
  i=1
  for(phylum in na.omit(unique(TronkoDB$phylum))){
    Taxon_GBIF <- name_backbone(phylum,verbose=T,strict=F)
    Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
    if(!("kingdom" %in% colnames(Taxon_GBIF))){
      test_class <- names(sort(table(TronkoDB[TronkoDB$phylum==phylum,"class"]),decreasing=TRUE)[1])
      if(!is.null(test_class)){
        Taxon_GBIF <- name_backbone(test_class,verbose=T,strict=F)
        Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
      }
    }
    tmp <- data.frame(matrix(nrow=1,ncol=2))
    colnames(tmp) <- c("phylum","kingdom")
    tmp$phylum <- phylum
    if(nrow(Taxon_GBIF)>0){
      tmp$kingdom <- names(sort(table(Taxon_GBIF$kingdom),decreasing=TRUE)[1])
    } else {
      tmp$kingdom <- NA
    }
    Phylum_to_Kingdom[[i]] <- tmp
    print(paste(Primer,i,phylum))
    i=i+1
  }
  Phylum_to_Kingdom <- rbindlist(Phylum_to_Kingdom, use.names=TRUE, fill=TRUE)
  Phylum_to_Kingdom <- as.data.frame(Phylum_to_Kingdom)
  TronkoDB <- dplyr::left_join(TronkoDB,Phylum_to_Kingdom)
  print(paste(Primer,"Kingdoms added"))
  
  #Count eDNA taxonomic resolution and weigh them.
  #Species = 4, Genus = 2, Family = 1.  Everything else = 0.
  TronkoDB$eDNAWeight <- 1*as.numeric(!is.na(TronkoDB$family))+2*as.numeric(!is.na(TronkoDB$genus))+4*as.numeric(!is.na(TronkoDB$species))
  
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  TronkoDB$LocalFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Local$family))))
  #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
  TronkoDB$StateFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_State$family))))
  #Check if any of the eDNA reads show up in the realm set of GBIF family observations.
  TronkoDB$NationFamilyPresentGBIF <- as.numeric(lapply(TronkoDB$family,is.element,unique(na.omit(Taxa_Nation$family))))
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  TronkoDB$LocalGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Local$genus))))
  #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
  TronkoDB$StateGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_State$genus))))
  #Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
  TronkoDB$NationGenusPresentGBIF <- as.numeric(lapply(TronkoDB$genus,is.element,unique(na.omit(Taxa_Nation$genus))))
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  TronkoDB$LocalSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Local$species))))
  #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
  TronkoDB$StateSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_State$species))))
  #Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
  TronkoDB$NationSpeciesPresentGBIF <- as.numeric(lapply(TronkoDB$species,is.element,unique(na.omit(Taxa_Nation$species))))
  
  #Assign TOS scores for GBIF results.
  TronkoDB$TOS_Local <- (1*TronkoDB$LocalFamilyPresentGBIF+2*TronkoDB$LocalGenusPresentGBIF+4*TronkoDB$LocalSpeciesPresentGBIF)/TronkoDB$eDNAWeight
  TronkoDB$TOS_State <- (1*TronkoDB$StateFamilyPresentGBIF+2*TronkoDB$StateGenusPresentGBIF+4*TronkoDB$StateSpeciesPresentGBIF)/TronkoDB$eDNAWeight
  TronkoDB$TOS_Nation <- (1*TronkoDB$NationFamilyPresentGBIF+2*TronkoDB$NationGenusPresentGBIF+4*TronkoDB$NationSpeciesPresentGBIF)/TronkoDB$eDNAWeight
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  TronkoDB[is.nan(TronkoDB)] <- 0
  print(paste(Primer,"TOS scores added"))
  TronkoDB$Primer <- Primer
  #Save Tronko-assign data set.
  Retained_Columns <- c("Taxonomic_Path","Score","SampleID","superkingdom","kingdom",
                        "phylum","class","order","family","genus","species","Primer",
                        "Mismatch","TOS_Local","TOS_State","TOS_Nation")
  write.table(TronkoDB[Retained_Columns],"Taxa_Parsed.tsv",quote=FALSE,sep="\t",row.names = FALSE)
  system(paste("aws s3 cp Taxa_Parsed.tsv s3://ednaexplorer/projects/",ProjectID,"/",Primer,"/Taxa_Parsed.tsv --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
  #Save all of the Tronko-assign data from a project.
  TronkoProject[[j]] <- TronkoDB[,c("superkingdom","kingdom","phylum","class","order","family","genus","species")]
  j=j+1
}
#Get all Tronko-assign taxa information for a project.
TronkoProject <- data.table::rbindlist(TronkoProject,use.names=TRUE, fill=TRUE)
TronkoProject <- as.data.frame(TronkoProject)

#Check existing Phylopic database
check_PhylopicDB <- suppressWarnings(system("aws s3 ls s3://ednaexplorer/PhylopicDB.tsv --endpoint-url https://js2.jetstream-cloud.org:8001/",intern=T))
#Get taxa which have already been added to the Phylopic database.
if(length(check_PhylopicDB>0)){
  PhylopicDB_Initial <- system("aws s3 cp s3://ednaexplorer/PhylopicDB.tsv - --endpoint-url https://js2.jetstream-cloud.org:8001/",intern=TRUE)
  PhylopicDB_Initial <- read.table(text = paste(PhylopicDB_Initial,sep = "\t"),header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  Taxa_Completed <- unique(na.omit(PhylopicDB_Initial$Taxon))
  usageKey_Completed <- unique(na.omit(PhylopicDB_Initial$usageKey))
  Common_Names_Initial <- PhylopicDB_Initial[,c("usageKey","Common_Name")]
  Common_Names_Initial <- as.data.frame(Common_Names_Initial[!duplicated(Common_Names_Initial),])
} else {
  PhylopicDB_Initial <- data.frame()
  Taxa_Completed <- c()
  usageKey_Completed <- c()
  Common_Names_Initial <- data.frame()
}

#Get unique taxa list.
Unique_Taxa_df <- as.data.frame(unlist(TronkoProject))
colnames(Unique_Taxa_df) <- c("Taxon")
Unique_Taxa_df$TaxonomicRankImage <- rownames(Unique_Taxa_df)
Unique_Taxa_df$TaxonomicRankImage <- gsub('[[:digit:]]+', '', Unique_Taxa_df$TaxonomicRankImage)
Unique_Taxa_df$rank <- toupper(Unique_Taxa_df$TaxonomicRankImage)
Unique_Taxa_df$rank <- gsub('SUPERKINGDOM', 'KINGDOM', Unique_Taxa_df$rank)
Unique_Taxa_df <- Unique_Taxa_df[complete.cases(Unique_Taxa_df),]
Unique_Taxa_df <- Unique_Taxa_df[Unique_Taxa_df$Taxon!="unassigned",]
Unique_Taxa_df <- Unique_Taxa_df[!duplicated(Unique_Taxa_df),]
Unique_Taxa <- unique(na.omit(Unique_Taxa_df$Taxon))
Unique_Taxa <- Unique_Taxa[!(Unique_Taxa %in% Taxa_Completed)]

#Get unique Phylopic icons for each taxon.  If one does not exist, keep moving up a taxonomic level until an icon does exist.
TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
PhylopicDB <- list()
PhylopicDB[[1]] <- PhylopicDB_Initial
i=1
for(Unique_Taxon in Unique_Taxa){
  Taxon_GBIF <- name_backbone(Unique_Taxon,verbose=T,strict=F)
  Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
  if(nrow(Taxon_GBIF)>0){
    Taxon_Backbone <- as.numeric(na.omit(rev(unlist(Taxon_GBIF[1,unique(grep(paste(TaxonomicRanks,collapse="|"),colnames(Taxon_GBIF[,grepl("Key",names(Taxon_GBIF))]), value=TRUE))]))))
    Taxon_Keys <- as.data.frame(Taxon_GBIF[1,!(colnames(Taxon_GBIF) %in% TaxonomicRanks)])
    if(length(Taxon_Backbone)>1){
      res <- httr::POST(url="https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true",body=jsonlite::toJSON(as.character(Taxon_Backbone), auto_unbox=TRUE),content_type("application/json"),encode="json",http_version=2)
    } else{
      res <- httr::POST(url="https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true",body=jsonlite::toJSON(as.character(Taxon_Backbone), auto_unbox=FALSE),content_type("application/json"),encode="json",http_version=2)
    }
    test <- fromJSON(rawToChar(res$content))
    Taxon_Image <- test[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
    if(is.null(Taxon_Image)){
      Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
    }
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
  print(paste(Primer,i,length(Unique_Taxa)))
  i=i+1
  PhylopicDB[[i]] <- tmp
}
PhylopicDB <- data.table::rbindlist(PhylopicDB,use.names=TRUE, fill=TRUE)
PhylopicDB <- as.data.frame(PhylopicDB)
PhylopicDB <- dplyr::left_join(PhylopicDB,Unique_Taxa_df)
print(paste(Primer,"Phylopics added"))
#Get common names and merge them into Phylopic database.
usageKey_All <- unique(na.omit(PhylopicDB$usageKey))
usageKey_Remaining <- usageKey_All[!(usageKey_All %in% usageKey_Completed)]
Common_Names <- list()
Common_Names[[1]] <- Common_Names_Initial
i=1
for(key in usageKey_Remaining){
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
  Common_Names[[i]] <- tmp
  print(paste(Primer,i,length(unique(PhylopicDB$usageKey))))
  i=i+1
}
Common_Names <- data.table::rbindlist(Common_Names,use.names=TRUE, fill=TRUE)
Common_Names <- as.data.frame(Common_Names)
PhylopicDB <- dplyr::left_join(PhylopicDB,Common_Names)
print(paste(Primer,"Common names added"))
#Save updated Phylopic database.
write.table(PhylopicDB,"PhylopicDB.tsv",quote=FALSE,sep="\t",row.names = FALSE)
system("aws s3 cp PhylopicDB.tsv s3://ednaexplorer/PhylopicDB.tsv --endpoint-url https://js2.jetstream-cloud.org:8001/")

system("rm ne_10m_admin_1_states_provinces.*")
system("rm Metadata.csv")
system("rm PhylopicDB.tsv")
system("rm Taxa_Parsed.tsv")
