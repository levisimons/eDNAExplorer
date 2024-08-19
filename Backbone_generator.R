rm(list=ls())
require(dplyr)
require(taxonbridge)#Make sure taxonkit is installed: conda install -c bioconda taxonkit
require(httr)
require(curl)
httr::set_config(httr::config(http_version = 2))
curl::handle_setopt(new_handle(),http_version=2)
require(jsonlite)
require(data.table)
require(digest)
require(DBI)
require(RPostgreSQL)
readRenviron(".env")
taxonomy_home <- Sys.getenv("taxonomy_home")
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
Database_Driver <- dbDriver("PostgreSQL")

setwd(taxonomy_home)

#Download GBIF taxonomic backbone
download.file(url="https://hosted-datasets.gbif.org/datasets/backbone/current/backbone.zip",destfile=paste(taxonomy_home,"backbone.zip",sep="/"))
system("unzip -o backbone.zip")
system("rm backbone.zip")
#Download NCBI taxonomic backbone
download.file(url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",destfile=paste(taxonomy_home,"taxdmp.zip",sep="/"))
system("unzip -o taxdmp.zip")
system("rm taxdmp.zip")
#Load NCBI taxonomic backbone into All.lineages.tsv.gz
system(paste("taxonkit list --data-dir=",taxonomy_home," --ids 1 | taxonkit lineage --show-lineage-taxids --show-lineage-ranks --show-rank --show-name --data-dir=",taxonomy_home," | taxonkit reformat --taxid-field 1 --data-dir=",taxonomy_home," -o All.lineages.tsv.gz",sep=""))
#Load combined NCBI and GBIF taxonomies.
custom_taxonomy <- load_taxonomies(paste(taxonomy_home,"Taxon.tsv",sep="/"), paste(taxonomy_home,"All.lineages.tsv.gz",sep="/"))
#Export initial GBIF-NCBI Taxonbridge
write.table(custom_taxonomy, paste(taxonomy_home,"ncbi_gbif_backbone_full.tsv",sep="/"), sep = "\t", col.names = TRUE, row.names = FALSE)

#NCBI taxonomic ranks.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
#GBIF taxonomic ranks.
GBIF_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

#Read in NCBI-GBIF Taxonbridge 
Taxon_Bridge <- fread(input="ncbi_gbif_backbone_full.tsv",sep="\t")
Taxon_Bridge <- Taxon_Bridge[Taxon_Bridge$taxonomicStatus!="doubtful"]
#Set genus
Taxon_Bridge$genus <- ifelse(Taxon_Bridge$taxonRank %in% c("species","genus"),Taxon_Bridge$genericName,NA_character_)
#Set species
Taxon_Bridge$species <- ifelse(Taxon_Bridge$taxonRank=="species" & !is.na(Taxon_Bridge$specificEpithet),paste(Taxon_Bridge$genus,Taxon_Bridge$specificEpithet),NA_character_) 
#Use the most accepted GBIF taxon ID
Taxon_Bridge$taxonID <- ifelse(Taxon_Bridge$taxonomicStatus!="accepted",Taxon_Bridge$acceptedNameUsageID,Taxon_Bridge$taxonID)
#Preferentially using GBIF entries with an accepted taxonomy for the GBIF-NCBI backbone
Taxon_Bridge_Accepted <- Taxon_Bridge %>%
  group_by(taxonID,ncbi_id) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()
Taxon_Bridge_Accepted <- Taxon_Bridge_Accepted[Taxon_Bridge_Accepted$taxonRank %in% GBIF_ranks,]
#Export
write.table(Taxon_Bridge_Accepted,"Taxon_Bridge_Accepted.tsv",quote=FALSE,sep="\t",row.names = FALSE)
rm(Taxon_Bridge)

Taxon_Bridge_Accepted <- fread(input="Taxon_Bridge_Accepted.tsv",sep="\t")

#Read in Open Tree of Life (OTL) backbone https://tree.opentreeoflife.org/about/taxonomy-version/ott3.5
OTL_taxonomy <- fread(input="otl_taxonomy.tsv",sep="\t")
OTL_taxonomy <- OTL_taxonomy[,c("uid","parent_uid","name","rank","sourceinfo")]

#Get NCBI and GBIF taxonomic ranks and names from OTL
OTL_db <- OTL_taxonomy
OTL_db <- OTL_db %>% mutate(rank = ifelse(rank == "domain", "superkingdom", rank))
#Retain entries with at least a NCBI ID
OTL_db <- OTL_db %>% filter(str_detect(sourceinfo, "ncbi:\\d+"))
#Extract the first returns for GBIF and NCBI IDs, along with how many of them are associated with a particular OTL ID.
OTL_standardized <- OTL_db[OTL_db$rank %in% TaxonomicRanks,]
OTL_standardized <- OTL_standardized %>% mutate(ncbi_count = str_count(sourceinfo, "ncbi:\\d+"))
OTL_standardized <- OTL_standardized %>% mutate(ncbi_id = str_extract(sourceinfo, "ncbi:\\d+") %>% str_extract("\\d+"))
OTL_standardized <- OTL_standardized %>% mutate(gbif_count = str_count(sourceinfo, "gbif:\\d+"))
OTL_standardized <- OTL_standardized %>% mutate(gbif_id = str_extract(sourceinfo, "gbif:\\d+") %>% str_extract("\\d+"))
OTL_standardized$gbif_id <- as.integer(OTL_standardized$gbif_id)
OTL_standardized$ncbi_id <- as.integer(OTL_standardized$ncbi_id)
#Add in GBIF IDs and their taxonomic status from Taxon Bridge to check for missing entries downstream.
OTL_standardized <- dplyr::left_join(OTL_standardized,Taxon_Bridge[,c("taxonID","taxonRank","taxonomicStatus")],by=c("gbif_id"="taxonID"))
OTL_standardized <- OTL_standardized[!duplicated(OTL_standardized),]
OTL_standardized$gbif_id <- as.integer(OTL_standardized$gbif_id)
OTL_standardized$ncbi_id <- as.integer(OTL_standardized$ncbi_id)

write.table(OTL_standardized,"OTL_standardized.tsv",quote=FALSE,sep="\t",row.names = FALSE)
OTL_standardized <- fread(input="OTL_standardized.tsv",sep="\t")

#Function to find the most taxonomically resolved unique link between NCBI and GBIF
#taxonomic IDs within the Open Tree of Life database
select_rows <- function(row_num) {
  #result <- c()
  i <- row_num
  df <- OTL_standardized
  current_row <- df[i, ]
  original_name <- current_row$name
  original_rank <- current_row$rank
  original_ncbi_id <- current_row$ncbi_id
  final_rank <- current_row$rank
  # Check if ncbi_count and gbif_count both are equal to 1
  if (current_row$ncbi_count == 1 & current_row$gbif_count == 1) {
    current_row$name <- original_name
    current_row$rank <- original_rank
    current_row$ncbi_id <- original_ncbi_id
    #result[[i]] <- current_row
    result <- current_row
  } else {
    # Keep checking until ncbi_count and gbif_count both are equal to 1
    while (!(current_row$ncbi_count == 1 & current_row$gbif_count == 1)) {
      current_row <- df[df$uid == current_row$parent_uid, ]
      if (nrow(current_row) == 0) break
      current_row$name <- original_name
      current_row$rank <- original_rank
      current_row$ncbi_id <- original_ncbi_id
      current_row <- current_row %>%
        group_by(uid, parent_uid) %>%
        slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
        ungroup()
      current_row <- as.data.frame(current_row)
    }
    #result[[i]] <- current_row
    result <- current_row
    #print(paste(i,"of",nrow(df)))
  }
  return(result)
}

#Convert OTL database to only have IDs and names for entries with standard taxonomic ranks.
k_min <- 1
k_delta <- 10000
k <- 1
k_max <- k*k_delta
otl_taxa_reformatted <- c()
while(k_max < nrow(OTL_standardized)){
  #Get full taxonomic paths for each taxon in NCBI
  tmp <- pbmclapply(k_min:k_max, select_rows,mc.cores=detectCores())
  tmp <- rbind.fill(tmp)
  otl_taxa_reformatted[[k]] <- tmp
  print(paste(k,k_min,k_max))
  k <- k+1
  k_min <- k_max+1
  k_max <- k*k_delta
}
k_max <- nrow(OTL_standardized)
print(paste(k,k_min,k_max))
tmp <- pbmclapply(k_min:k_max, select_rows,mc.cores=detectCores())
tmp <- rbind.fill(tmp)
otl_taxa_reformatted[[k]] <- tmp
#Get full taxonomic paths for each taxon in OTL
OTL_GBIF_updated <- rbind.fill(otl_taxa_reformatted)
#Retain most accepted taxa from OTL
OTL_GBIF_updated <- OTL_GBIF_updated %>%
  group_by(ncbi_id,gbif_id) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()
#Clear doubtful and depreciated GBIF taxonomic entries.
OTL_GBIF_updated <- OTL_GBIF_updated[!is.na(OTL_GBIF_updated$taxonRank),]
#Rename columns
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "name"] <- "ncbi_name"
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "rank"] <- "ncbi_rank"
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "taxonRank"] <- "gbif_rank"
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "name"] <- "ncbi_name"
#Store results
write.table(OTL_GBIF_updated,"OTL_GBIF_updated.tsv",quote=FALSE,sep="\t",row.names = FALSE)
OTL_GBIF_updated <- fread(input="OTL_GBIF_updated.tsv",sep="\t")

#NCBI node and name descriptions: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
#Download https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip (2024-06-21 version)
#Read in ncbi nodes
ncbi_nodes <- fread("nodes.dmp",sep="\t")
ncbi_nodes <- ncbi_nodes[,c("V1","V3","V5")]
ncbi_nodes <- ncbi_nodes[!duplicated(ncbi_nodes),]
colnames(ncbi_nodes) <- c("ncbi_id","parent_id","ncbi_rank")
#Read in ncbi names
ncbi_names <- fread("names.dmp",sep="\t")
ncbi_names <- ncbi_names[!is.na(ncbi_names$V5) & ncbi_names$V7=="scientific name",]
ncbi_names <- ncbi_names[,c("V1","V3")]
ncbi_names <- ncbi_names[!duplicated(ncbi_names),]
colnames(ncbi_names) <- c("ncbi_id","ncbi_name")
#Merge names and nodes
ncbi_taxa_named_nodes <- dplyr::left_join(ncbi_nodes,ncbi_names)

#Take NCBI database and reformat it so that all of the taxonomic ranks are arranged in
#order columns
ncbi_taxa <- c()
ncbi_taxa_rows <- nrow(ncbi_taxa_named_nodes)
i=1
ncbi_taxa[[i]] <- ncbi_taxa_named_nodes
while(ncbi_taxa_rows>1){
  i=i+1
  ncbi_taxa[[i]] <- ncbi_taxa[[i-1]][(ncbi_id) %in% (parent_id), ]
  ncbi_taxa_rows <- nrow(ncbi_taxa[[i]])
  print(paste(i,ncbi_taxa_rows))
}
i_max <- i
for(i in 1:i_max){
  ncbi_taxa[[i]] <- setnames(ncbi_taxa[[i]], old = names(ncbi_taxa[[i]]), new = c(paste("ncbi_id_",(i),sep=""),paste("ncbi_id_",(i+1),sep=""),paste("ncbi_rank_",(i),sep=""),paste("ncbi_name_",(i),sep="")))
  if(i==1){ncbi_taxa_full <- ncbi_taxa[[i]]}
  if(i>1){
    ncbi_taxa_full <- dplyr::left_join(ncbi_taxa_full,ncbi_taxa[[i]])
  }
}

#Define a function to only retain standard taxonomic ranks.
rank_cols <- colnames(ncbi_taxa_full)[grep("ncbi_rank",colnames(ncbi_taxa_full))]
is_valid_column <- function(column) {
  any(column %in% c("species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"))
}

#Define a function to only retain columns related to the names and IDs of
#NCBI entries with standard taxonomic ranks.
ncbi_fix <- function(row_num) {
  j <- row_num
  valid_ranks <- rank_cols[sapply(ncbi_taxa_full[j, ..rank_cols], is_valid_column)]
  if(length(valid_ranks)>0){
    valid_indices <- as.integer(gsub("[^0-9]", "", valid_ranks))
    valid_names <- paste0("ncbi_name_",valid_indices,sep="")
    valid_ids <- paste0("ncbi_id_",valid_indices,sep="")
    tmp1 <- ncbi_taxa_full[j, ..valid_names]
    colnames(tmp1) <- paste0("ncbi_",unlist(ncbi_taxa_full[j, ..valid_ranks]),sep="")
    tmp2 <- ncbi_taxa_full[j, ..valid_ids]
    colnames(tmp2) <- paste0("ncbi_",unlist(ncbi_taxa_full[j, ..valid_ranks]),"_id",sep="")
    tmp <- cbind(tmp1,tmp2)
    ncbi_taxa_reformatted <- tmp
    return(ncbi_taxa_reformatted)
  }
}

#Convert NCBI database to only have IDs and names for entries with standard taxonomic ranks.
k_min <- 1
k_delta <- 10000
k <- 1
k_max <- k*k_delta
ncbi_taxa_reformatted <- c()
while(k_max < nrow(ncbi_taxa_full)){
  #Get full taxonomic paths for each taxon in NCBI
  tmp <- pbmclapply(k_min:k_max, ncbi_fix,mc.cores=detectCores())
  tmp <- rbind.fill(tmp)
  ncbi_taxa_reformatted[[k]] <- tmp
  print(paste(k,k_min,k_max))
  k <- k+1
  k_min <- k_max+1
  k_max <- k*k_delta
}
k_max <- nrow(ncbi_taxa_full)
print(paste(k,k_min,k_max))
tmp <- pbmclapply(k_min:k_max, ncbi_fix,mc.cores=detectCores())
tmp <- rbind.fill(tmp)
ncbi_taxa_reformatted[[k]] <- tmp
#Get full taxonomic paths for each taxon in NCBI
ncbi_taxa_reformatted <- rbind.fill(ncbi_taxa_reformatted)

#Generate a NCBI ID column
ncbi_taxa_reformatted <- ncbi_taxa_reformatted %>% mutate(ncbi_id = coalesce(ncbi_species_id, ncbi_genus_id, ncbi_family_id, ncbi_order_id, ncbi_class_id, ncbi_phylum_id, ncbi_kingdom_id, ncbi_superkingdom_id))
#Export results
write.table(ncbi_taxa_reformatted,"ncbi_taxa_full.tsv",quote=FALSE,sep="\t",row.names = FALSE)
ncbi_taxa_full <- fread(input="ncbi_taxa_full.tsv",sep="\t")

#Merge Open Tree of Life database entries to the full GBIF taxonomic paths.
GBIF_NCBI <- dplyr::left_join(OTL_GBIF_updated[,c("ncbi_name","ncbi_rank","ncbi_id","gbif_id","gbif_rank")],
                              Taxon_Bridge_Accepted[,c("taxonID","taxonomicStatus","canonicalName","taxonRank","species","genus","family","order","class","phylum","kingdom")],
                              by=c("gbif_id"="taxonID"))
rm(Taxon_Bridge_Accepted)

#Mege in NCBI backbone.
GBIF_NCBI <- dplyr::left_join(GBIF_NCBI,ncbi_taxa_full,by=c("ncbi_id"="ncbi_id"))
#Preferentially using GBIF entries with an accepted taxonomy for the GBIF-NCBI backbone
GBIF_NCBI <- GBIF_NCBI %>%
  group_by(gbif_id,ncbi_id) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()

#Remove extraneous columns
drop_cols <- c("canonicalName","taxonRank","taxonomicStatus")
GBIF_NCBI <- GBIF_NCBI %>% dplyr::select(-one_of(drop_cols))
GBIF_NCBI <- GBIF_NCBI[!duplicated(GBIF_NCBI),]
GBIF_NCBI <- setDT(GBIF_NCBI)
#Clear out entries with NCBI entries with no corresponding GBIF entries
GBIF_NCBI <- GBIF_NCBI[rowSums(is.na(GBIF_NCBI[, ..GBIF_ranks])) != length(GBIF_ranks)]

#Export
write.table(GBIF_NCBI,"GBIF_NCBI_export.tsv",quote=FALSE,sep="\t",row.names = FALSE)
GBIF_NCBI <- fread(input="GBIF_NCBI_export.tsv",sep="\t")

#Pull in full GBIF taxonomic data set and determine the most accepted taxon keys for each
#name and rank.
Taxon_GBIF <- fread(input="Taxon.tsv",sep="\t",select=c("taxonID","taxonRank","taxonomicStatus","parentNameUsageID","canonicalName","species","genus","family","order","class","phylum","kingdom"),na.strings=c("NA",""))
#Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$taxonomicStatus!="doubtful",]
#Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$taxonomicStatus!="doubtful" & Taxon_GBIF$taxonRank %in% c("species","genus","family","order","class","phylum","kingdom"),]
#Preferentially select accepted taxonomic assignments.
Taxon_GBIF <- Taxon_GBIF[, species := ifelse(taxonRank == "species", canonicalName, NA)]
Taxon_GBIF <- Taxon_GBIF %>% mutate(gbif_name = coalesce(species, genus, family, order, class, phylum,kingdom))

#Take GBIF database and reformat it so that all of the taxonomic ranks are arranged in
#order columns
gbif_taxa <- c()
gbif_taxa_rows <- nrow(Taxon_GBIF)
i=1
gbif_taxa[[i]] <- Taxon_GBIF[,.(taxonID,parentNameUsageID,taxonRank,gbif_name,taxonomicStatus)]
while(gbif_taxa_rows>1){
  i=i+1
  gbif_taxa[[i]] <- gbif_taxa[[i-1]][(taxonID) %in% (parentNameUsageID),]
  gbif_taxa_rows <- nrow(gbif_taxa[[i]])
  print(paste(i,gbif_taxa_rows))
}
i_max <- i
for(i in 1:i_max){
  gbif_taxa[[i]] <- setnames(gbif_taxa[[i]], old = names(gbif_taxa[[i]]), new = c(paste("gbif_id_",(i),sep=""),paste("gbif_id_",(i+1),sep=""),paste("gbif_rank_",(i),sep=""),paste("gbif_name_",(i),sep=""),paste("gbif_status_",(i),sep="")))
  if(i==1){gbif_taxa_full <- gbif_taxa[[i]]}
  if(i>1){
    gbif_taxa_full <- dplyr::left_join(gbif_taxa_full,gbif_taxa[[i]])
  }
}

#Define a function to only retain standard taxonomic ranks.
rank_cols <- colnames(gbif_taxa_full)[grep("gbif_rank",colnames(gbif_taxa_full))]
is_valid_column <- function(column) {
  any(column %in% c("species", "genus", "family", "order", "class", "phylum", "kingdom"))
}

#Define a function to only retain columns related to the names and IDs of
#GBIF entries with standard taxonomic ranks.
gbif_fix <- function(row_num) {
  j <- row_num
  valid_ranks <- rank_cols[sapply(gbif_taxa_full[j, ..rank_cols], is_valid_column)]
  if(length(valid_ranks)>0){
    valid_indices <- as.integer(gsub("[^0-9]", "", valid_ranks))
    valid_names <- paste0("gbif_name_",valid_indices,sep="")
    valid_ids <- paste0("gbif_id_",valid_indices,sep="")
    tmp1 <- gbif_taxa_full[j, ..valid_names]
    colnames(tmp1) <- paste0("gbif_",unlist(gbif_taxa_full[j, ..valid_ranks]),sep="")
    tmp2 <- gbif_taxa_full[j, ..valid_ids]
    colnames(tmp2) <- paste0("gbif_",unlist(gbif_taxa_full[j, ..valid_ranks]),"_id",sep="")
    tmp <- cbind(tmp1,tmp2)
    gbif_taxa_reformatted <- tmp
    return(gbif_taxa_reformatted)
  }
}

#Convert GBIF database to only have IDs and names for entries with standard taxonomic ranks.
k_min <- 1
k_delta <- 10000
k <- 1
k_max <- k*k_delta
gbif_taxa_reformatted <- c()
while(k_max < nrow(gbif_taxa_full)){
  #Get full taxonomic paths for each taxon in GBIF
  tmp <- pbmclapply(k_min:k_max, gbif_fix,mc.cores=detectCores())
  tmp <- rbind.fill(tmp)
  gbif_taxa_reformatted[[k]] <- tmp
  print(paste(k,k_min,k_max))
  k <- k+1
  k_min <- k_max+1
  k_max <- k*k_delta
}
k_max <- nrow(gbif_taxa_full)
print(paste(k,k_min,k_max))
tmp <- pbmclapply(k_min:k_max, gbif_fix,mc.cores=detectCores())
tmp <- rbind.fill(tmp)
gbif_taxa_reformatted[[k]] <- tmp
#Get full taxonomic paths for each taxon in GBIF
gbif_taxa_reformatted <- rbind.fill(gbif_taxa_reformatted)
#Generate a GBIF ID column
gbif_taxa_reformatted <- gbif_taxa_reformatted %>% mutate(taxonID = coalesce(gbif_species_id, gbif_genus_id, gbif_family_id, gbif_order_id, gbif_class_id, gbif_phylum_id, gbif_kingdom_id))
#Rename GBIF columns
gbif_taxa_reformatted <- gbif_taxa_reformatted %>% dplyr::rename(species=gbif_species, genus=gbif_genus, family=gbif_family, order=gbif_order, class=gbif_class, phylum=gbif_phylum, kingdom=gbif_kingdom, speciesKey=gbif_species_id, genusKey=gbif_genus_id, familyKey=gbif_family_id, orderKey=gbif_order_id, classKey=gbif_class_id, phylumKey=gbif_phylum_id, kingdomKey=gbif_kingdom_id)
gbif_taxa_reformatted <- gbif_taxa_reformatted[!duplicated(gbif_taxa_reformatted),]
#Export results
write.table(gbif_taxa_reformatted,"gbif_taxa_reformatted.tsv",quote=FALSE,sep="\t",row.names = FALSE)
rm(Taxon_GBIF)

#Merge in reformmated GBIF backbone
gbif_taxa_full <- dplyr::left_join(Taxon_GBIF,gbif_taxa_reformatted)
#Retain standard taxonomic ranks
gbif_taxa_full <- gbif_taxa_full[gbif_taxa_full$taxonRank %in% GBIF_ranks,]

#Preferentially retain entries with an accepted taxonomic status
gbif_taxa_full <- gbif_taxa_full %>%
  group_by(gbif_name, taxonRank) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()
gbif_taxa_full <- setDT(gbif_taxa_full)

#Remove entries where both the taxon name and taxon key are not either both present or both missing.
gbif_taxa_full <- gbif_taxa_full %>%
  filter(!(is.na(species) & !is.na(speciesKey)) & 
           !(is.na(speciesKey) & !is.na(species)) &
           !(is.na(genus) & !is.na(genusKey)) &
           !(is.na(genusKey) & !is.na(genus)) &
           !(is.na(family) & !is.na(familyKey)) &
           !(is.na(familyKey) & !is.na(family)) &
           !(is.na(order) & !is.na(orderKey)) &
           !(is.na(orderKey) & !is.na(order)) &
           !(is.na(class) & !is.na(classKey)) &
           !(is.na(classKey) & !is.na(class)) &
           !(is.na(phylum) & !is.na(phylumKey)) &
           !(is.na(phylumKey) & !is.na(phylum)) &
           !(is.na(kingdom) & !is.na(kingdomKey)) &
           !(is.na(kingdomKey) & !is.na(kingdom)))
#Export GBIF taxonomy table with keys.
write.table(gbif_taxa_full,"gbif_taxa_full.tsv",quote=FALSE,sep="\t",row.names = FALSE)
gbif_taxa_full <- fread(input="gbif_taxa_full.tsv",sep="\t")

#Merge in GBIF keys with GBIF-NCBI backbone.
GBIF_NCBI_WithKeys <- dplyr::left_join(GBIF_NCBI,gbif_taxa_full)
rm(gbif_taxa_full)

#Get the most common English common names per GBIF taxon ID
Vernacular_Input <- fread(file="VernacularName.tsv",sep="\t",na.strings=c("NA",""))
common_names <- Vernacular_Input[Vernacular_Input$language=="en" & !is.na(Vernacular_Input$language),c("taxonID","vernacularName")]
# Group by 'taxonID' and 'vernacularName', count occurrences
grouped <- common_names %>% dplyr::group_by(taxonID, vernacularName) %>% dplyr::summarize(count = n())
# Find the most common 'vernacularName' for each 'taxonID'
most_common <- grouped %>% group_by(taxonID) %>% slice_max(order_by = count, n = 1) %>% ungroup()
#Keep only the first value of 'vernacularName' for each 'taxonID'
first_common <- most_common %>% group_by(taxonID) %>% slice(1)

#Get available common names for all taxa and ranks. Merge in common names.
GBIF_NCBI_WithKeys_WithCommonNames <- GBIF_NCBI_WithKeys
for(taxonomicRank in GBIF_ranks){
  tmp <- first_common[,c("taxonID","vernacularName")]
  names(tmp)[names(tmp) == "vernacularName"] <- paste("common",taxonomicRank,sep="_")
  names(tmp)[names(tmp) == "taxonID"] <- paste(taxonomicRank,"Key",sep="")
  GBIF_NCBI_WithKeys_WithCommonNames <- dplyr::left_join(GBIF_NCBI_WithKeys_WithCommonNames,tmp)
}

#Get the most resolved common name.
GBIF_NCBI_WithKeys_WithCommonNames <- GBIF_NCBI_WithKeys_WithCommonNames %>%
  mutate(Common_Name = coalesce(
    common_species,
    common_genus,
    common_family,
    common_order,
    common_class,
    common_phylum,
    common_kingdom
  ))

#Get Phylopic images for all of the taxa in the GBIF-NCBI backbone.
phylopic_export <- GBIF_NCBI_WithKeys_WithCommonNames[,c("gbif_name","Common_Name","speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")]
phylopic_export <- phylopic_export[!duplicated(phylopic_export),]
phylopic_assigner <- function(row_num){
  i <- row_num
  Taxon_Backbone <- as.numeric(phylopic_export[i,c("speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")])
  #res <- httr::GET(url=paste("https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true&objectIDs=",paste(Taxon_Backbone,collapse=","),sep=""),config = httr::config(connecttimeout = 100))
  res <- httr::GET(url=paste("https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true&?filter_license_nc=false&objectIDs=",paste(Taxon_Backbone,collapse=","),sep=""),config = httr::config(connecttimeout = 100))
  phylopic_query <- fromJSON(rawToChar(res$content))
  Taxon_Image <- phylopic_query[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
  if(is.null(Taxon_Image)){
    Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
  }
  phylopic_export[i,"Image_URL"] <- Taxon_Image
  phylopic_assigned <- phylopic_export[i,]
  return(phylopic_assigned)
}
phylopic_assigned <- c()
tmp <- pbmclapply(1:nrow(phylopic_export), phylopic_assigner,mc.cores=detectCores())
phylopic_export <- rbindlist(tmp)
#Export Phylopic image table
write.table(phylopic_export, "phylopic_export.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

phylopic_export <- fread(input="phylopic_export.tsv", sep = "\t")
#Merge Phylopic image URLs in with other taxonomic data.
GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic <-  dplyr::left_join(GBIF_NCBI_WithKeys_WithCommonNames,phylopic_export)

#Read in GBIF file containing protected status information.
#https://www.gbif.org/dataset/19491596-35ae-4a91-9a98-85cf505f1bd3
#This file is a download of presence points with the IucnRedListCategory values of EX, NE, DD, LC, NT, VU, EN, CR, EW
#Get IUCN status information
IUCN_distribution <- fread(input="distribution.txt",sep="\t",na.strings=c("NA",""))
IUCN_distribution <- IUCN_distribution[,.(V1,V4)]
colnames(IUCN_distribution) <- c("iucnKey","iucnStatus")
IUCN_distribution$iucnStatus <- str_to_title(IUCN_distribution$iucnStatus)

#Create map for IUCN status values.
IUCN_categories <- data.frame(iucnRedListCategory  = c("EX","NE","DD","LC","NT","VU","EN","CR","EW"),
                              iucnStatus = c("Extinct","Not Evaluated","Data Deficient","Least Concern","Near Threatened","Vulnerable","Endangered","Critically Endangered","Extinct in the Wild"))

#Merge IUCN categories into status values.
IUCN_distribution <- dplyr::left_join(IUCN_distribution,IUCN_categories)

#Get taxonomic data to merge into IUCN status data.
IUCN_taxa <- fread(input="taxon.txt",sep="\t",na.strings=c("NA",""))
IUCN_taxa <- IUCN_taxa[,.(V14,V3,V4,V5,V6,V7,V8,V9,V13)]
colnames(IUCN_taxa) <- c("iucnKey","kingdom","phylum","class","order","family","genus","species_suffix","taxonomicStatus")
#Clean up taxonomic data
IUCN_taxa[,c("kingdom","phylum","class","order","family")] <- lapply(IUCN_taxa[,c("kingdom","phylum","class","order","family")], str_to_title)
IUCN_taxa$species <- ifelse(!is.na(IUCN_taxa$species_suffix),paste(IUCN_taxa$genus,IUCN_taxa$species_suffix),IUCN_taxa$species_suffix)
IUCN_taxa$species_suffix <- NULL
IUCN_taxa <- IUCN_taxa %>% mutate(gbif_name = coalesce(species,genus,family,order,class,phylum,kingdom))
#Preferentially using GBIF entries with an accepted taxonomy
IUCN_taxa <- IUCN_taxa %>%
  group_by(iucnKey) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()

#Merge taxonomies and IUCN status information.
IUCN <- dplyr::left_join(IUCN_taxa,IUCN_distribution)
IUCN <- IUCN[,c("species","genus","family","order","class","phylum","kingdom","gbif_name","iucnStatus","iucnRedListCategory")]
IUCN <- IUCN[!duplicated(IUCN),]

#Merge IUCN status into unified NCBI-GBIF backbone
GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN <- dplyr::left_join(GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic,IUCN)

#Replace NA values in iucn_status with Not Evaluated
GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN$iucnStatus[is.na(GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN$iucnStatus)] <- "Not Evaluated"
GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN$iucnRedListCategory[is.na(GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN$iucnRedListCategory)] <- "NE"
GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN <- GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN[!duplicated(GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN),]

#Export table
write.table(GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN, "full_backbone_export.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN <- fread(input="full_backbone_export.tsv",sep="\t")

#Export selected columns to postgresql database.
require(digest)
require(DBI)
require(RPostgreSQL)
readRenviron(".env")
taxonomy_home <- Sys.getenv("taxonomy_home")
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
Database_Driver <- dbDriver("PostgreSQL")

#Select and rename columns
GBIF_NCBI_eDNAExplorer_Export <- GBIF_NCBI_WithKeys_WithCommonNames_WithPhylopic_WithIUCN[,c("ncbi_superkingdom","ncbi_kingdom","ncbi_phylum","ncbi_class","ncbi_order","ncbi_family","ncbi_genus","ncbi_species","ncbi_name","ncbi_rank","kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey","Common_Name","iucnStatus","Image_URL")]
GBIF_NCBI_eDNAExplorer_Export <- GBIF_NCBI_eDNAExplorer_Export %>% dplyr::rename(superkingdom=ncbi_superkingdom, kingdom=ncbi_kingdom, phylum=ncbi_phylum, class=ncbi_class, order=ncbi_order, family=ncbi_family, genus=ncbi_genus, species=ncbi_species, Taxon=ncbi_name, rank=ncbi_rank)


#Columns for exporting to database.
#export_columns <- c("id","superkingdom","kingdom","phylum","class","order","family","genus","species","Taxon","rank","kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey","Common_Name","iucnStatus","Image_URL","UniqueID")

#Create unique ID for the Taxonomy database.
GBIF_NCBI_eDNAExplorer_Export$UniqueID <- sapply(paste(GBIF_NCBI_eDNAExplorer_Export$species,GBIF_NCBI_eDNAExplorer_Export$genus,GBIF_NCBI_eDNAExplorer_Export$family,GBIF_NCBI_eDNAExplorer_Export$order,GBIF_NCBI_eDNAExplorer_Export$class,GBIF_NCBI_eDNAExplorer_Export$phylum,GBIF_NCBI_eDNAExplorer_Export$kingdom,GBIF_NCBI_eDNAExplorer_Export$Taxon,GBIF_NCBI_eDNAExplorer_Export$rank,GBIF_NCBI_eDNAExplorer_Export$Image_URL,GBIF_NCBI_eDNAExplorer_Export$Common_Name,GBIF_NCBI_eDNAExplorer_Export$iucnStatus),digest,algo="md5")

#Export file
write.table(GBIF_NCBI_eDNAExplorer_Export, "GBIF_NCBI_eDNAExplorer_Export.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

#Export to posgresql database
con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)
#Clear old taxonomy entries.
dbExecute(con,'DELETE FROM "Taxonomy"')
#Write out new taxonomy table
dbWriteTable(con,"Taxonomy",GBIF_NCBI_eDNAExplorer_Export,row.names=FALSE,append=TRUE)
sapply(dbListConnections(Database_Driver), dbDisconnect)
