#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
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

# Fetch project ID early so we can use it for error output when possible.
ProjectID <- args[1]

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, error = error_message))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  
  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://ednaexplorer/errors/taxonomy/", new_filename, sep = "")
  } else {
    paste("s3://ednaexplorer/tronko_output/", ProjectID, "/", filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  stop(error_message)
}

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

tryCatch(
  {
    #Get project ID.
    #Rscript --vanilla eDNAExplorer_Metabarcoding_Taxonomy_Initializer.R "project ID string"
    if (length(args)==0) {
      stop("Need a project ID", call.=FALSE)
    } else if (length(args)==1) {
      ProjectID <- args[1]
    }
  },
  error = function(e) {
    process_error(e)
  }
)

# Update taxonomy tables and cache processed Tronko-assign output.
tryCatch(
  {
    #Read in project metadata.
    con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)
    Metadata <- tbl(con,"TronkoMetadata")
    Metadata <- Metadata %>% filter(projectid == ProjectID)
    Metadata <- as.data.frame(Metadata)
    #Get state and nation lists.
    country_list <- na.omit(unique(Metadata$nation))
    state_province_list <- na.omit(unique(Metadata$state))
    
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
    Markers <- grep("^marker_[[:digit:]]$",colnames(Metadata),value=T)
    Primers <- na.omit(unique(unlist(Metadata[,Markers])))
    j=1
    TaxaCount <- list()
    #Loop over primers to add Tronko-assign data to database, along with associate Phylopic metadata.
    for(Primer in Primers){
      #Read in Tronko-assign output files.  Standardize sample IDs within them.
      TronkoBucket <- system(paste("aws s3 ls s3://ednaexplorer/projects/",ProjectID," --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
      TronkoBucket <- read.table(text = paste(TronkoBucket,sep = ""),header = FALSE)
      colnames(TronkoBucket) <- c("Date", "Time", "Size","Filename")
      TronkoFiles <- unique(TronkoBucket$Filename)
      TronkoFiles <- TronkoFiles[grepl(paste("projects",ProjectID,"assign",Primer,sep="/"),TronkoFiles)]
      TronkoFiles <- TronkoFiles[grepl("*.txt$",TronkoFiles)]
      if(length(TronkoFiles)>0){
        TronkoInputs <- list()
        i=1
        #TronkoHeaders <- c("Readname","Taxonomic_Path","Score","Forward_Mismatch","Reverse_Mismatch","Tree_Number","Node_Number")
        for(TronkoFile in TronkoFiles){
          TronkoInput <- system(paste("aws s3 cp s3://ednaexplorer/",TronkoFile," - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
          if(length(TronkoInput)>0){
            TronkoInput <- read.table(text = paste(TronkoInput,sep = "\t"),header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
            if(nrow(TronkoInput)>0){
              print(paste(Primer,i,length(TronkoFiles)))
              #if(identical(colnames(TronkoInput),TronkoHeaders)==FALSE){
              #  header_row <- as.data.frame(t(colnames(TronkoInput)))
              #  colnames(header_row) <- TronkoHeaders
              #  colnames(TronkoInput) <- TronkoHeaders
              #}
              #TronkoInput$SampleID <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(TronkoFile))
              TronkoInputs[[i]] <- TronkoInput
              i=i+1
            }
          }
        }
        TronkoInputs <- rbindlist(TronkoInputs, use.names=TRUE, fill=TRUE)
        #Get ASV to sampleID information
        TronkoASVs <- unique(TronkoBucket$Filename)
        TronkoASVs <- TronkoASVs[grepl(paste("projects",ProjectID,"assign",Primer,sep="/"),TronkoASVs)]
        TronkoASVs <- TronkoASVs[grepl("*.asv$",TronkoASVs)]
        TronkoASVs <- TronkoASVs[!grepl("*-paired_R.asv$",TronkoASVs)]
        if(length(TronkoASVs)>0){
          ASVInputs <- list()
          m=1
          for(TronkoASV in TronkoASVs){
            ASVInput <- system(paste("aws s3 cp s3://ednaexplorer/",TronkoASV," - --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
            if(length(ASVInput)>0){
              ASVInput <- read.table(text = paste(ASVInput,sep = "\t"),header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
              if(nrow(ASVInput)>0){
                print(paste(Primer,j,length(TronkoASVs)))
                #TronkoInput$SampleID <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(TronkoFile))
                ASVInputs[[m]] <- ASVInput
                m=m+1
              }
            }
          }
        }
        ASVInputs <- rbindlist(ASVInputs, use.names=TRUE, fill=TRUE)
        colnames(ASVInputs) <- gsub(paste(Primer,"_",sep=""),"",colnames(ASVInputs))
        ASVtoSample <- ASVInputs
        ASVtoSample$sequence <- NULL
        ASVtoSample <- ASVtoSample %>%
          pivot_longer(cols = -seq_number, names_to = "SampleID", values_to = "observations") %>%
          filter(observations > 0)
        ASVtoSample$observations <- NULL
        TronkoWithASVs <- dplyr::left_join(TronkoInputs,ASVtoSample,by=c("Readname"="seq_number"),multiple="all")
        
        TronkoDB <- as.data.frame(TronkoWithASVs)
        #Remove reads with non-assigned scores.
        TronkoDB$Score <- as.numeric(TronkoDB$Score)
        TronkoDB <- TronkoDB[!is.na(TronkoDB$Score),]
        
        #Get the count of unique organisms per sample per primer.
        TaxaCount_Primer <- TronkoDB[,c("SampleID","Taxonomic_Path")]
        TaxaCount_Primer <- TaxaCount_Primer[!duplicated(TaxaCount_Primer),]
        TaxaCount[[j]] <- TaxaCount_Primer
        j=j+1
        
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
        k=1
        for(phylum in na.omit(unique(TronkoDB$phylum))){
          Taxon_GBIF <- name_backbone(phylum,verbose=T,strict=F,curlopts=list(http_version=2))
          print(paste(Primer,k,phylum))
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
          Phylum_to_Kingdom[[k]] <- tmp
          k=k+1
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
        #Calculate a geographically weighted TOS.
        tmp$gw_TOS <- (4*tmp$TOS_Local+2*tmp$TOS_State+1*tmp$TOS_Nation)/7
        print(paste(Primer,"TOS scores added"))
        
        #Remove rows without known family, genus, species data.
        tmp <- tmp[!is.na(tmp$family) & !is.na(tmp$genus) & !is.na(tmp$species),]
        
        #Merge in TOS results.
        TronkoDB <- dplyr::left_join(TronkoDB,tmp,by=c("family","genus","species"),na_matches = "never")
        
        TronkoDB$Primer <- Primer
        TronkoDB$ProjectID <- ProjectID
        #Save Tronko-assign data set on a per project/primer basis.
        Retained_Columns <- c("ProjectID","Primer","Taxonomic_Path","Score","SampleID","superkingdom","kingdom",
                              "phylum","class","order","family","genus","species","Readname",
                              "Mismatch","TOS_Local","TOS_State","TOS_Nation","gw_TOS")
        
        TronkoProject <- as.data.frame(TronkoDB[,Retained_Columns])
        
        #Save Tronko output.
        TronkoOutput_Filename <- paste(Primer,".csv",sep="")
        write.table(x=TronkoProject,file=TronkoOutput_Filename,quote=FALSE,sep=",",row.names = FALSE)
        system(paste("aws s3 cp ",TronkoOutput_Filename," s3://ednaexplorer/tronko_output/",ProjectID,"/",TronkoOutput_Filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
        system(paste("rm ",TronkoOutput_Filename,sep=""))
        
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
          for(n in 1:nrow(TaxaDB)){
            GBIF_Key <- TaxaDB[n,]
            tmp <- name_backbone(name=GBIF_Key$Taxon,rank=GBIF_Key$rank,genus=GBIF_Key$genus,
                                 family=GBIF_Key$family,order=GBIF_Key$order,
                                 class=GBIF_Key$class,phylum=GBIF_Key$phylum,kingdom=GBIF_Key$kingdom,curlopts=list(http_version=2))
            #Resolve taxon name/rank mismatches.
            if("rank" %in% colnames(tmp)){
              if(tolower(tmp$rank)!=GBIF_Key$rank & GBIF_Key$rank!="superkingdom"){
                GBIF_Key$rank <- tolower(tmp$rank)
                GBIF_Key[,tolower(tmp$rank)] <- as.character(tmp[,tolower(tmp$rank)])
                GBIF_Key[,paste(tolower(tmp$rank),"Key",sep="")] <- as.numeric(tmp[,paste(tolower(tmp$rank),"Key",sep="")])
                
                #rank_check <- name_lookup(GBIF_Key$Taxon,GBIF_Key$rank)
                #rank_check <- as.data.frame(rank_check$data)
                #if(nrow(rank_check)!=0 | !is.null(nrow(rank_check))){
                #  #Resolve missing taxon key for relevant taxonomic rank.
                #  GBIF_Key[,paste(GBIF_Key$rank,"Key",sep="")] <- as.numeric(names(sort(-table(rank_check[!is.na(rank_check$acceptedKey),"acceptedKey"])))[1])
                #}
              }
            }
            if(nrow(as.data.frame(tmp[,colnames(tmp) %in% TaxonomicKeyRanks & !(colnames(tmp) %in% colnames(GBIF_Key))]))>0){
              GBIF_Key <- cbind(GBIF_Key,as.data.frame(tmp[,colnames(tmp) %in% TaxonomicKeyRanks & !(colnames(tmp) %in% colnames(GBIF_Key))]))
            }
            
            #Get GBIF backbone for querying Phylopic images
            match_ranks <- GBIF_Key[,colnames(GBIF_Key)[colnames(GBIF_Key) %in% TaxonomicKeyRanks]]
            if(!is.null(dim(match_ranks))){
              match_ranks <- match_ranks[,colSums(is.na(match_ranks))<nrow(match_ranks)]
              match_ranks <- colnames(match_ranks)
            }
            if(length(match_ranks)>0){
              match_ranks <- match_ranks[order(match(match_ranks,TaxonomicKeyRanks[length(TaxonomicKeyRanks):1]))]
              Taxon_Backbone <- as.numeric(GBIF_Key[,match_ranks[length(match_ranks):1]])
              #Get Phylopic images for each taxon
              res <- httr::GET(url=paste("https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true&objectIDs=",paste(Taxon_Backbone,collapse=","),sep=""),config = httr::config(connecttimeout = 100))
              Sys.sleep(0.5)
              test <- fromJSON(rawToChar(res$content))
              Taxon_Image <- test[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
              if(is.null(Taxon_Image)){
                Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
              }
              print(paste(n,GBIF_Key$rank,GBIF_Key$Taxon,Taxon_Image))
              #Get common names if available
              if(length(na.omit(Taxon_Backbone))==0){Common_Name <- NA}
              if(length(na.omit(Taxon_Backbone))>0){
                Taxon_Key <- GBIF_Key[1,paste(tolower(tmp[1,"rank"]),"Key",sep="")]
                Common_Name <- as.data.frame(name_usage(key=Taxon_Key,rank=tolower(tmp[1,"rank"]), data="vernacularNames",curlopts=list(http_version=2))$data)
                if(nrow(Common_Name)>0){
                  Common_Name <- Common_Name[Common_Name$language=="eng","vernacularName"]
                  if(length(Common_Name)>0){
                    Common_Name <- names(sort(-table(Common_Name)))[1]
                    Common_Name <- gsub("[[:punct:]]", " ", Common_Name)
                  } else {
                    Common_Name <- NA
                  }
                } else {
                  Common_Name <- NA
                }
              }
              GBIF_Key$Common_Name <- Common_Name
              GBIF_Key$Image_URL <- Taxon_Image
            }
            if(length(match_ranks)==0){
              GBIF_Key$Common_Name <- NA
              GBIF_Key$Image_URL <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
            }
            GBIF_Keys[[n]] <- GBIF_Key
          }
          
          GBIF_Keys <- rbindlist(GBIF_Keys, use.names=TRUE, fill=TRUE)
          GBIF_Keys <- as.data.frame(GBIF_Keys)
          #Create unique ID for the Phylopic database.
          GBIF_Keys$UniqueID <- sapply(apply(GBIF_Keys, 1, function(row) paste(row, collapse = "")),digest,algo="md5")
          
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
    }
    #Get the number of unique organisms per sample per project.
    TaxaCounts <- rbindlist(TaxaCount, use.names=TRUE, fill=TRUE)
    TaxaCount_Project <- as.data.frame(table(TaxaCounts[,c("SampleID")]))
    colnames(TaxaCount_Project) <- c("SampleID","Unique_Taxa")
    TaxaCount_Filename <- paste("taxa_counts_",ProjectID,".json",sep="")
    write(toJSON(TaxaCount_Project),TaxaCount_Filename)
    system(paste("aws s3 cp ",TaxaCount_Filename," s3://ednaexplorer/projects/",ProjectID,"/plots/",TaxaCount_Filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
    system(paste("rm ",TaxaCount_Filename,sep=""))
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  },
  error = function(e) {
    process_error(e, filename)
  }
)
