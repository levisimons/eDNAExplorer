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

#Establish database credentials.
readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
gbif_dir <- Sys.getenv("GBIF_HOME")
bucket <- Sys.getenv("S3_BUCKET")
taxonomy_home <- Sys.getenv("taxonomy_home")
Database_Driver <- dbDriver("PostgreSQL")

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, lastRanAt = Sys.time(), error = error_message))
  write(json_content, filename)
  
  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  dest_filename <- sub("\\.json$", ".build", filename)
  
  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://",bucket,"/errors/taxonomy/", new_filename, sep = "")
  } else {
    paste("s3://",bucket,"/projects/", ProjectID, "/plots/", dest_filename, " --endpoint-url https://js2.jetstream-cloud.org:8001/", sep = "")
  }
  
  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ",filename,sep=""))
  stop(error_message)
}

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
      TronkoBucket <- system(paste("aws s3 ls s3://",bucket,"/projects/",ProjectID," --recursive --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
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
          system(paste("aws s3 cp s3://",bucket,"/",TronkoFile," ",basename(TronkoFile)," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
          TronkoInput <- fread(file=basename(TronkoFile),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
          if(nrow(TronkoInput)>0){
            TronkoInput <- as.data.frame(TronkoInput)
            print(paste(Primer,i,length(TronkoFiles)))
            TronkoInputs[[i]] <- TronkoInput
            i=i+1
          }
          system(paste("rm ",basename(TronkoFile),sep=""))
        }
        TronkoInputs <- rbindlist(TronkoInputs, use.names=TRUE, fill=TRUE)
        TronkoInputs$Taxonomic_Path <- gsub(","," ",TronkoInputs$Taxonomic_Path)
        #Get ASV to sampleID information
        TronkoASVs <- unique(TronkoBucket$Filename)
        TronkoASVs <- TronkoASVs[grepl(paste("projects",ProjectID,"assign",Primer,sep="/"),TronkoASVs)]
        TronkoASVs <- TronkoASVs[grepl("*.asv$",TronkoASVs)]
        TronkoASVs <- TronkoASVs[!grepl("*-paired_R.asv$",TronkoASVs)]
        if(length(TronkoASVs)>0){
          ASVInputs <- list()
          m=1
          for(TronkoASV in TronkoASVs){
            system(paste("aws s3 cp s3://",bucket,"/",TronkoASV," ",basename(TronkoASV)," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
            ASVInput <- fread(file=basename(TronkoASV),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
            if(nrow(ASVInput)>0){
              print(paste(Primer,j,length(TronkoASVs)))
              ASVInput <- as.data.frame(ASVInput)
              ASVInputs[[m]] <- ASVInput
              m=m+1
            }
            system(paste("rm ",basename(TronkoASV),sep=""))
          }
        }
        ASVInputs <- rbindlist(ASVInputs, use.names=TRUE, fill=TRUE)
        TronkoInputs$Taxonomic_Path <- gsub(","," ",TronkoInputs$Taxonomic_Path)
        colnames(ASVInputs) <- gsub(paste(Primer,"_",sep=""),"",colnames(ASVInputs))
        names(ASVInputs)[grep('seq_number', names(ASVInputs))] <- 'seq_number'
        
        max_rows_per_split <- 1000
        asv_table_filename <- paste(ProjectID,Primer,"asv_table.txt",sep="_")
        # Calculate the number of splits needed
        num_splits <- ceiling(nrow(ASVInputs) / max_rows_per_split)
        # Create an empty list to store the splits
        split_list <- list()
        # Split the data frame and store in the list
        for (p in 1:num_splits) {
          start_row <- (p - 1) * max_rows_per_split + 1
          end_row <- min(p * max_rows_per_split, nrow(ASVInputs))
          asv_subset <- ASVInputs[start_row:end_row,!c("sequence")]
          asv_subset <- asv_subset %>%
            pivot_longer(cols = -seq_number, names_to = "SampleID", values_to = "observations") %>%
            filter(observations > 0) %>% select("seq_number","SampleID")
          if (p==1) {
            write.table(asv_subset, asv_table_filename, sep = "\t", col.names = TRUE, row.names = FALSE)
          } else {
            write.table(asv_subset, asv_table_filename, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
          }
          print(paste(p,num_splits))
        }
        
        ASVtoSample <- fread(file=asv_table_filename,header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
        system(paste("rm ",asv_table_filename,sep=""))
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
          Taxon_GBIF <- name_backbone(phylum,verbose=T,strict=F,curlopts=list(http_version=2,timeout_ms=5000))
          print(paste(Primer,k,phylum))
          Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$matchType!="NONE",]
          if(!("kingdom" %in% colnames(Taxon_GBIF))){
            test_class <- names(sort(table(TronkoDB[TronkoDB$phylum==phylum,"class"]),decreasing=TRUE)[1])
            if(!is.null(test_class)){
              Taxon_GBIF <- name_backbone(test_class,verbose=T,strict=F,curlopts=list(http_version=2,timeout_ms=5000))
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
        system(paste("aws s3 cp ",TronkoOutput_Filename," s3://",bucket,"/tronko_output/",ProjectID,"/",TronkoOutput_Filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
        system(paste("rm ",TronkoOutput_Filename,sep=""))
      }
    }
    
    #Get the number of unique organisms per sample per project.
    TaxaCounts <- rbindlist(TaxaCount, use.names=TRUE, fill=TRUE)
    TaxaCount_Project <- as.data.frame(table(TaxaCounts[,c("SampleID")]))
    colnames(TaxaCount_Project) <- c("SampleID","Unique_Taxa")
    TaxaCount_Filename <- paste("taxa_counts_",ProjectID,".json",sep="")
    write(toJSON(TaxaCount_Project),TaxaCount_Filename)
    system(paste("aws s3 cp ",TaxaCount_Filename," s3://",bucket,"/projects/",ProjectID,"/plots/",TaxaCount_Filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
    system(paste("rm ",TaxaCount_Filename,sep=""))
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
    
    #Save log file.
    filename <- paste(gsub(" ","_",date()),"eDNAExplorer_Metabarcoding_Taxonomy_Initializer.R.log",sep="_")
    system(paste("echo > ",filename,sep=""))
    system(paste("aws s3 cp ",filename," s3://",bucket,"/projects/",ProjectID,"/log/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
    system(paste("rm ",filename))
  },
  error = function(e) {
    process_error(e, filename)
  }
)
