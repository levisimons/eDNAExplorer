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
require(anytime)

#Establish database credentials.
readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
Database_Driver <- dbDriver("PostgreSQL")
#Force close any possible postgreSQL connections.
sapply(dbListConnections(Database_Driver), dbDisconnect)

#Get project ID.
#Rscript --vanilla eDNAExplorer_Metabarcoding_qPCR_Initializer.R "project ID string"
if (length(args)<1) {
  stop("Need a project ID", call.=FALSE)
} else if (length(args)==1) {
  ProjectID <- args[1]
}

#Find qPCR project data file and read it into a dataframe.
Project_Scan <- system(paste("aws s3 ls s3://ednaexplorer/projects/",ProjectID," --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=T)
Project_Scan <- read.table(text = paste(Project_Scan,sep = ""),header = FALSE)
colnames(Project_Scan) <- c("Date", "Time", "Size","Filename")
Project_Scan <- Project_Scan[grep(".csv$",Project_Scan$Filename),]
for(csv_file in unique(Project_Scan$Filename)){
  Project_Data <- system(paste("aws s3 cp s3://ednaexplorer/",csv_file," - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  #Check if file is qPCR input metadata.
  if(length(grep("Target 1",Project_Data))==1){
    #Read in qPCR project data.
    Project_Data <- system(paste("aws s3 cp s3://ednaexplorer/",csv_file," - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    Project_Data <- gsub("[\r\n]", "", Project_Data)
    Project_Data <- read.table(text = Project_Data,header=FALSE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
    colnames(Project_Data) <- Project_Data[5,]
    Project_Data <- Project_Data[6:nrow(Project_Data),]
    addFormats(c("%m/%d/%y","%m-%d-%y","%d/%m/%y","%y/%m/%d"))
    Project_Data$`Sample Date` <- anytime::anydate(Project_Data$`Sample Date`)
    Project_Data$`Data type` <- NULL
    Project_Data$`Additional environmental metadata....` <- NULL
    gsub('^Target [[:digit:]] qPCR Probe Fluorophore (dye)$','^Target [[:digit:]] qPCR Probe Fluorophore$',colnames(Project_Data))
    gsub('^Target [[:digit:]] Cycle Threshold (ct)$','^Target [[:digit:]] Cycle Threshold (ct)$',colnames(Project_Data))
    Project_Data <- Project_Data %>% dplyr::mutate_at(c("Latitude","Longitude","Spatial Uncertainty"),as.numeric)
  }
}

Field_Variables <- colnames(Project_Data)[!(colnames(Project_Data) %in% c("Sample ID","Longitude","Latitude","Sample Date","Spatial Uncertainty"))]
#Read in extracted metadata.
Metadata_Extracted <- system(paste("aws s3 cp s3://ednaexplorer/projects/",ProjectID,"/MetadataOutput_qPCR.csv - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
Metadata_Extracted <- read.table(text = Metadata_Extracted,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
Metadata_Extracted$Sample_Date <- as.Date(as.POSIXct(Metadata_Extracted$Sample_Date))
#Shorten some variable names for downstream database storage.
colnames(Metadata_Extracted) <- gsub("Sea water salinity in practical salinity units at a depth of ","salinity units at depth of ",colnames(Metadata_Extracted))

#Remove duplicate rows in extracted metadata
Metadata_Extracted <- Metadata_Extracted[!duplicated(Metadata_Extracted),]

#Merge metadata and project data
MergedData <- dplyr::left_join(Project_Data[,c("Sample ID","Sample Date","Latitude","Longitude","Spatial Uncertainty",Field_Variables)],Metadata_Extracted,by=c("Sample ID"="name","Sample Date"="Sample_Date","Latitude","Longitude","Spatial Uncertainty"="Spatial_Uncertainty"),na_matches = "never")

#Add project ID
MergedData$ProjectID <- ProjectID

#Get target organism variables from initial metadata.
Target_Variables <- colnames(MergedData)[grep("^Target [[:digit:]]",colnames(MergedData))]

#Get non-target variables from initial metadata.
Non_Target_Variables <- colnames(MergedData)[!(colnames(MergedData) %in% Target_Variables)]

#Get number of target organisms.
Target_Numbers <- unique(as.numeric(gsub("\\D", "", Target_Variables)))

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
GADM_Boundaries <- st_join(st_as_sf(MergedData[!is.na(MergedData$Latitude) & !is.na(MergedData$Latitude),], coords = c("Longitude", "Latitude"), crs = 4326), GADM_1_Boundaries[c('iso_a2','woe_name')],join = st_intersects)
GADM_Boundaries <- GADM_Boundaries %>% st_drop_geometry()
GADM_Boundaries <- as.data.frame(GADM_Boundaries[,c("Sample ID","Sample Date","iso_a2","woe_name")])
names(GADM_Boundaries)[names(GADM_Boundaries) == "woe_name"] <- "State"
names(GADM_Boundaries)[names(GADM_Boundaries) == "iso_a2"] <- "Nation"
GADM_Boundaries <- GADM_Boundaries[!duplicated(GADM_Boundaries),]
MergedData <- dplyr::left_join(MergedData,GADM_Boundaries,by=c("Sample ID","Sample Date"),na_matches = "never")
country_list <- na.omit(unique(GADM_Boundaries$Nation))
state_province_list <- na.omit(unique(GADM_Boundaries$State))

#Generate unique code for each sample
MergedData <- MergedData[!duplicated(MergedData),]
Target_Variables <- colnames(MergedData)[grep("^Target [[:digit:]]",colnames(MergedData))]
MergedData$UniqueID <- sapply(apply(MergedData[,c("Sample ID","Sample Date","Latitude","Longitude","Spatial Uncertainty",Target_Variables)],1,paste,collapse = "" ),digest,algo="md5")

#Open qPCR sample database.
con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)

#Read in GBIF occurrences.
gbif <- gbif_local()

#Get local bounds for sample locations, add 0.5 degree buffer.
Local_East <- max(na.omit(MergedData$Longitude))+0.5
Local_West <- min(na.omit(MergedData$Longitude))-0.5
Local_South <- min(na.omit(MergedData$Latitude))-0.5
Local_North <- max(na.omit(MergedData$Latitude))+0.5

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

#Get unique taxa with full taxonomy
TaxaList <- unique(unlist(MergedData[,grep("^Target [[:digit:]] Organism$",colnames(MergedData))]))
Taxa <- list()
i=1
TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
TaxonomicKeyRanks <- c("speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")
for(Taxon in TaxaList){
  Taxon_GBIF <- name_backbone(Taxon,verbose=T,strict=F,curlopts=list(http_version=2))
  Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",colnames(Taxon_GBIF) %in% c(TaxonomicRanks,TaxonomicKeyRanks)]
  if(nrow(Taxon_GBIF)>0){
    Taxon_GBIF <- Taxon_GBIF[1,]
    Taxa[[i]] <- Taxon_GBIF
    i=i+1
  }
}
TaxaDB <- rbindlist(Taxa, use.names=TRUE, fill=TRUE)
TaxaDB <- as.data.frame(TaxaDB[!duplicated(TaxaDB),])
TaxaDB <- TaxaDB %>% dplyr::mutate(superkingdom = case_when(
  kingdom=="Bacteria" ~ "Bacteria",
  kingdom=="Archaea" ~ "Archaea",
  !is.na(kingdom) & kingdom!="Bacteria" & kingdom!="Archaea" ~ "Eukaryota"
))
tmp <- TaxaDB[,colnames(TaxaDB) %in% TaxonomicRanks]
tmp <- setcolorder(tmp, intersect(TaxonomicRanks, names(tmp)))
tmp$Taxon <- tmp[cbind(1:nrow(tmp), max.col(!is.na(tmp), ties.method = 'last'))]
tmp$rank <- TaxonomicRanks[max.col(!is.na(tmp), ties.method="last")]
TaxaDB <- dplyr::left_join(TaxaDB,tmp)

if(nrow(TaxaDB)>0){
  GBIF_Keys <- c()
  for(i in 1:nrow(TaxaDB)){
    GBIF_Key <- TaxaDB[i,]
    if(is.na(GBIF_Key$kingdom)){GBIF_Key$kingdom <- GBIF_Key$superkingdom}
    tmp <- name_backbone(name=GBIF_Key$Taxon,rank=GBIF_Key$rank,genus=GBIF_Key$genus,
                         family=GBIF_Key$family,order=GBIF_Key$order,
                         class=GBIF_Key$class,phylum=GBIF_Key$phylum,kingdom=GBIF_Key$kingdom,curlopts=list(http_version=2))
    GBIF_Key <- dplyr::left_join(GBIF_Key,as.data.frame(tmp[,colnames(tmp) %in% TaxonomicKeyRanks]))
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
}
GBIF_Keys <- rbindlist(GBIF_Keys, use.names=TRUE, fill=TRUE)
GBIF_Keys <- as.data.frame(GBIF_Keys)
#Create unique ID for the Phylopic database.
GBIF_Keys$UniqueID <- sapply(paste(GBIF_Keys$Taxon,GBIF_Keys$rank),digest,algo="md5")

#Count eDNA taxonomic resolution and weigh them.
#Species = 4, Genus = 2, Family = 1.  Everything else = 0.
tmp <- GBIF_Keys[,c("family","genus","species")]
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

#Merge in TOS results.
tmp <- dplyr::left_join(GBIF_Keys[c("Taxon","family","genus","species")],tmp,by=c("family","genus","species"))
for(TargetVar in unique(grep("^Target [[:digit:]] Organism$",colnames(MergedData),value=T))){
  MergedData <- MergedData %>% dplyr::left_join(tmp[,c("Taxon","TOS_Local","TOS_State","TOS_Nation")],by=setNames("Taxon",TargetVar))
}
MergedData <- dplyr::left_join(MergedData,tmp[,!(colnames(tmp) %in% c("family","genus","species"))],by=c("Target Organism"="Taxon"))

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

#Match qPCR project column names to format in SQL database.
colnames(MergedData) <- gsub(" ","_",tolower(colnames(MergedData)))

#Remove non-standard columns
Required_Variables <- unique(c("site","sample_id","sample_type","longitude","latitude","sample_date","spatial_uncertainty","sample_replicate_number",grep("^target_[[:digit:]]",colnames(MergedData),value=T),gsub(" ","_",tolower(colnames(Metadata_Extracted)))))
Required_Variables <- Required_Variables[Required_Variables != "name"]
MergedData <- MergedData[,colnames(MergedData) %in% Required_Variables]

#Check for redundant data.
#Add new project data.
if(dbExistsTable(con,"QPCRSample")){
  MergedData_Check <-  tbl(con,"QPCRSample")
  MergedData_IDs <- MergedData$uniqueid
  MergedData_Check <- MergedData_Check %>% filter(uniqueid %in% MergedData_IDs)
  MergedData_Check <- as.data.frame(MergedData_Check)
  MergedData_Check_IDs <- MergedData_Check$uniqueid
  MergedData_Append <- MergedData[!(MergedData_IDs %in% MergedData_Check_IDs),]
  dbWriteTable(con,"QPCRSample",MergedData_Append,row.names=FALSE,append=TRUE)
} else{
  dbWriteTable(con,"QPCRSample",MergedData,row.names=FALSE,append=TRUE)
}

RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
system("rm ne_10m_admin_1_states_provinces.*")


###
Metadata <- read.table(file="~/Desktop/MetabarcodingJeffreyMiller.csv",header=FALSE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
colnames(Metadata) <- Metadata[5,]
Metadata <- Metadata[6:nrow(Metadata),]
colnames(Metadata) <- gsub(" ","_",tolower(colnames(Metadata)))
Markers <- grep("^marker_[[:digit:]]$",colnames(Metadata),value=T)
Metadata <- Metadata[,!(colnames(Metadata) %in% grep("^marker_[[:alpha:]]",colnames(Metadata),value=T))]
