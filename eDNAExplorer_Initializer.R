rm(list=ls())
require(tidyr)
require(sf)
require(sp)
require(lubridate)
require(httr)
require(curl)
httr::set_config(httr::config(http_version = 2))
curl::handle_setopt(new_handle(),http_version=2)
require(gbifdb)
require(jsonlite)
require(rgbif)
require(data.table)
require(dplyr)
require(DBI)
require(RPostgreSQL)
require(digest)

Sys.setenv("AWS_ACCESS_KEY_ID" = "e9190baae65b40a38bf43ade883b04a6","AWS_SECRET_ACCESS_KEY" = "7aa129a7d84744efa76183cc9cf4b0a5")

ProjectID <- "LARiverRound1" #This is hard-coded for now.

#Read in initial metadata.
Metadata_Initial <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/InputMetadata.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Initial <- read.table(text = Metadata_Initial,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
Metadata_Initial$`Sample Date` <- lubridate::mdy(Metadata_Initial$`Sample Date`)
#Get field variables from initial metadata.
Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% c("Sample ID","Longitude","Latitude","Sample Date","Spatial Uncertainty"))]
#Read in extracted metadata.
Metadata_Extracted <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/MetadataOutput.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Extracted <- read.table(text = Metadata_Extracted,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
Metadata_Extracted$Sample_Date <- lubridate::ymd_hms(Metadata_Extracted$Sample_Date)

#Merge metadata
Metadata <- dplyr::left_join(Metadata_Initial[,c("Sample ID","Sample Date","Latitude","Longitude","Spatial Uncertainty",Field_Variables)],Metadata_Extracted,by=c("Sample ID"="name","Sample Date"="Sample_Date","Latitude","Longitude","Spatial Uncertainty"="Spatial_Uncertainty"))

#Add project ID
Metadata$ProjectID <- ProjectID

#Add Fastq ID to make sure metadata and Tronko-assign output lines up.
Metadata$FastqID <- gsub("_R.*","",Metadata$`Fastq Forward Reads Filename`)

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
GADM_Boundaries <- st_join(st_as_sf(Metadata[!is.na(Metadata$Latitude) & !is.na(Metadata$Latitude),], coords = c("Longitude", "Latitude"), crs = 4326), GADM_1_Boundaries[c('iso_a2','woe_name')],join = st_intersects)
GADM_Boundaries <- GADM_Boundaries %>% st_drop_geometry()
GADM_Boundaries <- as.data.frame(GADM_Boundaries[,c("Sample ID","Sample Date","iso_a2","woe_name")])
names(GADM_Boundaries)[names(GADM_Boundaries) == "woe_name"] <- "State"
names(GADM_Boundaries)[names(GADM_Boundaries) == "iso_a2"] <- "Nation"
Metadata <- dplyr::left_join(Metadata,GADM_Boundaries,by=c("Sample ID","Sample Date"))
country_list <- na.omit(unique(GADM_Boundaries$Nation))
state_province_list <- na.omit(unique(GADM_Boundaries$State))

#Remove rows without an associated sequence file.
Metadata <- Metadata[!is.na(Metadata$`Fastq Forward Reads Filename`) & !is.na(Metadata$`Fastq Reverse Reads Filename`),]

#Generate unique code for each Tronko-assign output
Metadata$UniqueID <- sapply(paste(Metadata$ProjectID,Metadata$FastqID,Metadata$`Sample Date`,Metadata$Latitude,Metadata$Longitude,Metadata$`Spatial Uncertainty`),digest,algo="md5")

#Match metadata column names to format in SQL database.
colnames(Metadata) <- gsub(" ","_",tolower(colnames(Metadata)))

#Establish database credentials.
db_host <- "db.jfsudbghjnulaznuohbj.supabase.co"
db_port <- 5432
db_name <- "postgres"
db_user <- "postgres"
db_pass <- "Bazjic-xanbog-hokma2"
Database_Driver <- dbDriver("PostgreSQL")
#Force close any possible postgreSQL connections.
sapply(dbListConnections(Database_Driver), dbDisconnect)

#Create Metadata database.
con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)

#Check for redundant data.
#Add new metadata.
if(dbExistsTable(con,"TronkoMetadata")){
  Metadata_Check <-  tbl(con,"TronkoMetadata")
  Metadata_IDs <- Metadata$uniqueid
  Metadata_Check <- Metadata_Check %>% filter(uniqueid %in% Metadata_IDs)
  Metadata_Check <- as.data.frame(Metadata_Check)
  Metadata_Check_IDs <- Metadata_Check$uniqueid
  Metadata_Append <- Metadata[!(Metadata_IDs %in% Metadata_Check_IDs),]
  dbWriteTable(con,"TronkoMetadata",Metadata_Append,row.names=FALSE,append=TRUE)
} else{
  dbWriteTable(con,"TronkoMetadata",Metadata,row.names=FALSE,append=TRUE)
}

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
Markers <- setdiff(colnames(Metadata[,grepl("^marker_[[:digit:]]$",colnames(Metadata))]),colnames(Metadata[,grepl("marker_(.*?)_",colnames(Metadata))]))
Primers <- unique(unlist(Metadata[,Markers]))
#Loop over primers to add Tronko-assign data to database, along with associate Phylopic metadata.
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
        TronkoInput$SampleID <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(TronkoFile))
        TronkoInputs[[i]] <- TronkoInput
        i=i+1
      }
    }
  }
  TronkoInputs <- rbindlist(TronkoInputs, use.names=TRUE, fill=TRUE)
  
  TronkoDB <- as.data.frame(TronkoInputs)
  #Remove reads with non-assigned scores.
  TronkoDB$Score <- as.numeric(TronkoDB$Score)
  TronkoDB <- TronkoDB[!is.na(TronkoDB$Score),]
  
  #Standardize taxonomy naming schema.
  TronkoDB$sum.taxonomy <- TronkoDB$Taxonomic_Path
  
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
    Taxon_GBIF <- name_backbone(phylum,verbose=T,strict=F,curlopts=list(http_version=2))
    print(paste(Primer,i,phylum))
    Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
    if(!("kingdom" %in% colnames(Taxon_GBIF))){
      test_class <- names(sort(table(TronkoDB[TronkoDB$phylum==phylum,"class"]),decreasing=TRUE)[1])
      if(!is.null(test_class)){
        Taxon_GBIF <- name_backbone(test_class,verbose=T,strict=F,curlopts=list(http_version=2))
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
    i=i+1
  }
  Phylum_to_Kingdom <- rbindlist(Phylum_to_Kingdom, use.names=TRUE, fill=TRUE)
  Phylum_to_Kingdom <- as.data.frame(Phylum_to_Kingdom)
  TronkoDB <- dplyr::left_join(TronkoDB,Phylum_to_Kingdom)
  print(paste(Primer,"Kingdoms added"))
  
  #Count eDNA taxonomic resolution and weigh them.
  #Species = 4, Genus = 2, Family = 1.  Everything else = 0.
  tmp <- TronkoDB[,c("family","genus","species")]
  tmp <- tmp[!duplicated(tmp),]
  tmp$eDNAWeight <- 1*as.numeric(!is.na(tmp$family))+2*as.numeric(!is.na(tmp$genus))+4*as.numeric(!is.na(tmp$species))
  
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  tmp$LocalFamilyPresentGBIF <- as.numeric(lapply(tmp$family,is.element,unique(na.omit(Taxa_Local$family))))
  #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
  tmp$StateFamilyPresentGBIF <- as.numeric(lapply(tmp$family,is.element,unique(na.omit(Taxa_State$family))))
  #Check if any of the eDNA reads show up in the realm set of GBIF family observations.
  tmp$NationFamilyPresentGBIF <- as.numeric(lapply(tmp$family,is.element,unique(na.omit(Taxa_Nation$family))))
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  tmp$LocalGenusPresentGBIF <- as.numeric(lapply(tmp$genus,is.element,unique(na.omit(Taxa_Local$genus))))
  #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
  tmp$StateGenusPresentGBIF <- as.numeric(lapply(tmp$genus,is.element,unique(na.omit(Taxa_State$genus))))
  #Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
  tmp$NationGenusPresentGBIF <- as.numeric(lapply(tmp$genus,is.element,unique(na.omit(Taxa_Nation$genus))))
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  tmp$LocalSpeciesPresentGBIF <- as.numeric(lapply(tmp$species,is.element,unique(na.omit(Taxa_Local$species))))
  #Check if any of the eDNA reads show up in the ecoregional set of GBIF family observations.
  tmp$StateSpeciesPresentGBIF <- as.numeric(lapply(tmp$species,is.element,unique(na.omit(Taxa_State$species))))
  #Check if any of the eDNA reads show up in the realm set of GBIF genus observations.
  tmp$NationSpeciesPresentGBIF <- as.numeric(lapply(tmp$species,is.element,unique(na.omit(Taxa_Nation$species))))
  
  #Assign TOS scores for GBIF results.
  tmp$TOS_Local <- (1*tmp$LocalFamilyPresentGBIF+2*tmp$LocalGenusPresentGBIF+4*tmp$LocalSpeciesPresentGBIF)/tmp$eDNAWeight
  tmp$TOS_State <- (1*tmp$StateFamilyPresentGBIF+2*tmp$StateGenusPresentGBIF+4*tmp$StateSpeciesPresentGBIF)/tmp$eDNAWeight
  tmp$TOS_Nation <- (1*tmp$NationFamilyPresentGBIF+2*tmp$NationGenusPresentGBIF+4*tmp$NationSpeciesPresentGBIF)/tmp$eDNAWeight
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  tmp[is.nan(tmp)] <- 0
  print(paste(Primer,"TOS scores added"))
  
  #Merge in TOS results.
  TronkoDB <- dplyr::left_join(TronkoDB,tmp,by=c("family","genus","species"))
  
  TronkoDB$Primer <- Primer
  TronkoDB$ProjectID <- ProjectID
  #Save Tronko-assign data set on a per project/primer basis.
  Retained_Columns <- c("ProjectID","Primer","Taxonomic_Path","Score","SampleID","superkingdom","kingdom",
                        "phylum","class","order","family","genus","species","Readname",
                        "Mismatch","TOS_Local","TOS_State","TOS_Nation")
  
  TronkoProject <- as.data.frame(TronkoDB[,Retained_Columns])
  
  #Generate unique code for each Tronko-assign output
  TronkoProject$UniqueID <- sapply(paste(TronkoProject$ProjectID,TronkoProject$Primer,TronkoProject$Taxonomic_Path,TronkoProject$SampleID,TronkoProject$Score,TronkoProject$Mismatch,TronkoProject$Readname),digest,algo="md5")
  
  #Create Tronko output database.
  #Check for redundant data.
  #Add only new data.
  chunk <- 10000
  chunks <- split(1:nrow(TronkoProject), ceiling(seq_along(1:nrow(TronkoProject))/chunk))
  TronkoInput <-  tbl(con,"TronkoOutput")
  if(dbExistsTable(con,"TronkoOutput")==TRUE){
    for(i in 1:ceiling(nrow(TronkoProject)/chunk)){
      TronkoProject_Subset <- TronkoProject[min(chunks[[i]]):max(chunks[[i]]),]
      Tronko_IDs <- TronkoProject_Subset$UniqueID
      Tronko_Check <- TronkoInput %>% filter(UniqueID %in% Tronko_IDs)
      Tronko_Check <- as.data.frame(Tronko_Check)
      Tronko_Check_IDs <- Tronko_Check$UniqueID
      TronkoProject_Subset <- TronkoProject_Subset[!(Tronko_IDs %in% Tronko_Check_IDs),]
      dbWriteTable(con,"TronkoOutput",TronkoProject_Subset,row.names=FALSE,append=TRUE)
    }
  } 
  if(dbExistsTable(con,"TronkoOutput")==FALSE){
    for(i in 1:ceiling(nrow(TronkoProject)/chunk)){
      dbWriteTable(con,"TronkoOutput",TronkoProject[min(chunks[[i]]):max(chunks[[i]]),],row.names=FALSE,append=TRUE)
    }
  }
  
  #Create database of Phylopic images and common names for taxa.
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  TaxonomicKeyRanks <- c("speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")
  TaxaDB <- TronkoDB[,TaxonomicRanks]
  TaxaDB <- TaxaDB[!duplicated(TaxaDB),]
  TaxaDB$Taxon <- TaxaDB[cbind(1:nrow(TaxaDB), max.col(!is.na(TaxaDB), ties.method = 'last'))]
  TaxaDB$rank <- TaxonomicRanks[max.col(!is.na(TaxaDB[TaxonomicRanks]), ties.method="last")]
  #Create unique ID for the Phylopic database.
  TaxaDB$UniqueID <- sapply(paste(TaxaDB$Taxon,TaxaDB$rank),digest,algo="md5")
  
  #Check for pre-existing Phylopic database entries.  Only leave new and unique entries to append.
  if(dbExistsTable(con,"Taxonomy")==TRUE){
    Phylopic_Check <-  tbl(con,"Taxonomy")
    Phylopic_IDs <- TaxaDB$UniqueID
    Phylopic_Check <- Phylopic_Check %>% filter(UniqueID %in% Phylopic_IDs)
    Phylopic_Check <- as.data.frame(Phylopic_Check)
    Phylopic_Check_IDs <- Phylopic_Check$UniqueID
    TaxaDB <- TaxaDB[!(Phylopic_IDs %in% Phylopic_Check_IDs),]
  } 
  if(dbExistsTable(con,"Taxonomy")==FALSE){
    TaxaDB <- TaxaDB
  }
  
  #Get Phylopic urls and common names for each unique taxon to append to database.
  if(nrow(TaxaDB)>0){
    GBIF_Keys <- c()
    for(i in 1:nrow(TaxaDB)){
      GBIF_Key <- TaxaDB[i,]
      if(is.na(GBIF_Key$kingdom)){GBIF_Key$kingdom <- GBIF_Key$superkingdom}
      tmp <- name_backbone(name=GBIF_Key$Taxon,rank=GBIF_Key$rank,genus=GBIF_Key$genus,
                           family=GBIF_Key$family,order=GBIF_Key$order,
                           class=GBIF_Key$class,phylum=GBIF_Key$phylum,kingdom=GBIF_Key$kingdom,curlopts=list(http_version=2))
      GBIF_Key <- cbind(GBIF_Key,as.data.frame(tmp[,colnames(tmp) %in% TaxonomicKeyRanks]))
      #Get GBIF backbone for querying Phylopic images
      Taxon_Backbone <- as.numeric(na.omit(rev(unlist(GBIF_Key[,colnames(GBIF_Key) %in% TaxonomicKeyRanks]))))
      #Get Phylopic images for each taxon
      res <- httr::GET(url=paste("https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true&objectIDs=",paste(Taxon_Backbone,collapse=","),sep=""),config = httr::config(connecttimeout = 100))
      Sys.sleep(0.5)
      test <- fromJSON(rawToChar(res$content))
      Taxon_Image <- test[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
      if(is.null(Taxon_Image)){
        Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
      }
      print(paste(i,Taxon_Image))
      #Get common names if available
      if(length(na.omit(Taxon_Backbone))==0){Common_Name <- NA}
      if(length(na.omit(Taxon_Backbone))>0){
        Common_Name <- as.data.frame(name_usage(key=max(na.omit(Taxon_Backbone)),rank=GBIF_Key[1,"rank"], data="vernacularNames",curlopts=list(http_version=2))$data)
        if(nrow(Common_Name)>0){
          Common_Name <- Common_Name[Common_Name$language=="eng",]
          Common_Name <- Common_Name[1,"vernacularName"]
        } else {
          Common_Name <- NA
        }
      }
      GBIF_Key$Common_Name <- Common_Name
      GBIF_Key$Image_URL <- Taxon_Image
      GBIF_Keys[[i]] <- GBIF_Key
    }
    
    GBIF_Keys <- rbindlist(GBIF_Keys, use.names=TRUE, fill=TRUE)
    GBIF_Keys <- as.data.frame(GBIF_Keys)
    #Create unique ID for the Phylopic database.
    GBIF_Keys$UniqueID <- sapply(paste(GBIF_Keys$Taxon,GBIF_Keys$rank),digest,algo="md5")
    
    #Check for redundant data.
    #Add new Phylopic data.
    if(dbExistsTable(con,"Taxonomy")==TRUE){
      Phylopic_Check <-  tbl(con,"Taxonomy")
      Phylopic_IDs <- GBIF_Keys$UniqueID
      Phylopic_Check <- Phylopic_Check %>% filter(UniqueID %in% Phylopic_IDs)
      Phylopic_Check <- as.data.frame(Phylopic_Check)
      Phylopic_Check_IDs <- Phylopic_Check$UniqueID
      Phylopic_Append <- GBIF_Keys[!(Phylopic_IDs %in% Phylopic_Check_IDs),]
      dbWriteTable(con,"Taxonomy",Phylopic_Append,row.names=FALSE,append=TRUE)
    } 
    if(dbExistsTable(con,"Taxonomy")==FALSE){
      dbWriteTable(con,"Taxonomy",GBIF_Keys,row.names=FALSE,append=TRUE)
    } 
  }
}
RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
system("rm ne_10m_admin_1_states_provinces.*")
