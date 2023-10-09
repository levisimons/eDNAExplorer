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
#taxonomy_home <- "~/Desktop/backbone/"

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

#Save combined NCBI and GBIF taxonomies.
write.table(custom_taxonomy, paste(taxonomy_home,"ncbi_gbif_backbone_full.tsv",sep="/"), sep = "\t", col.names = TRUE, row.names = FALSE)

#Filter taxonomies to only include entries with both accepted NCBI and GBIF entries.
custom_taxonomy <- read.table(file=paste(taxonomy_home,"ncbi_gbif_backbone_full.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"",as.is=TRUE, encoding = "UTF-8",na = c("", "NA", "N/A"))
#Remove doubtful taxonomies
filtered_taxonomy <- custom_taxonomy[custom_taxonomy$taxonomicStatus!="doubtful",]
#Keep only the first unique entry of taxonID.
filtered_taxonomy <- filtered_taxonomy %>% distinct(taxonID, .keep_all = TRUE)
#Standardize GBIF species naming to be genus+species.
taxonomicRanks <- c("species","genus","family","order","class","phylum","kingdom")
names(filtered_taxonomy)[names(filtered_taxonomy) == 'genericName'] <- 'genus'
names(filtered_taxonomy)[names(filtered_taxonomy) == 'specificEpithet'] <- 'species'
tmp <- filtered_taxonomy[,c("taxonID","genus","species")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp$species <- paste(tmp$genus,tmp$species)
filtered_taxonomy$species <- NULL
filtered_taxonomy <- dplyr::left_join(filtered_taxonomy,tmp[,c("taxonID","species")])

#Read in taxonomy NCBI - GBIF synonyms from OTL https://tree.opentreeoflife.org/about/taxonomy-version/ott3.5
OTL_taxonomy <- read.table(file=paste(taxonomy_home,"otl_taxonomy.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"",as.is=TRUE, encoding = "UTF-8",na = c("", "NA", "N/A"))
OTL_bridge <- as.data.frame(OTL_taxonomy$sourceinfo)
colnames(OTL_bridge) <- c("sourceinfo")
#Get rows which have both gbif and ncbi entries.
OTL_bridge <- OTL_bridge %>% 
  filter(grepl("gbif",sourceinfo)) %>%
  filter(grepl("ncbi",sourceinfo))
#Create a data frame linking GBIF and NCBI ids.
combinations <- c()
for(i in 1:nrow(OTL_bridge)){
  OTL_bridge_subset <- OTL_bridge[i,]
  OTL_bridge_subset_list <- as.list(unlist(strsplit(OTL_bridge_subset, ",")))
  gbif_list <- c()
  for(j in 1:length(OTL_bridge_subset_list)){
    if(grepl("gbif",OTL_bridge_subset_list[j])){
      gbif_list <- c(gbif_list,gsub("gbif:","",OTL_bridge_subset_list[j]))
    }
  }
  ncbi_list <- c()
  for(j in 1:length(OTL_bridge_subset_list)){
    if(grepl("ncbi",OTL_bridge_subset_list[j])){
      ncbi_list <- c(ncbi_list,gsub("ncbi:","",OTL_bridge_subset_list[j]))
    }
  }
  combinations[[i]] <- as.data.frame(expand.grid(taxonID = gbif_list, ncbi_id = ncbi_list))
  print(paste(i,nrow(OTL_bridge)))
}
OTL_GBIF_NCBI <- rbindlist(combinations, use.names=TRUE, fill=TRUE)
OTL_GBIF_NCBI <- as.data.frame(OTL_GBIF_NCBI)
OTL_GBIF_NCBI$taxonID <- as.integer(as.character(OTL_GBIF_NCBI$taxonID))
OTL_GBIF_NCBI$ncbi_id <- as.integer(as.character(OTL_GBIF_NCBI$ncbi_id))

#Find the gbif ids for where a ncbi one is missing
gbif_missing_ncbi <- na.omit(filtered_taxonomy[is.na(filtered_taxonomy$ncbi_id),"taxonID"])
#Find additional ncbi ids from the OTL backbone.
additional_ncbi <- as.numeric(as.character(na.omit(OTL_GBIF_NCBI[na.omit(OTL_GBIF_NCBI$taxonID) %in% gbif_missing_ncbi,"ncbi_id"])))
ncbi_to_add <- filtered_taxonomy[na.omit(filtered_taxonomy$ncbi_id),]
ncbi_to_add <- ncbi_to_add[ncbi_to_add$ncbi_id %in% additional_ncbi,colnames(ncbi_to_add)[sapply(colnames(ncbi_to_add), function(x) grepl("ncbi", x, ignore.case = TRUE))]]
#Add in additional gbif ids.
ncbi_to_add <- dplyr::left_join(ncbi_to_add,OTL_GBIF_NCBI,multiple="all",relationship = "many-to-many")
#Combine back with entries with corresponding gbif ids.
ncbi_to_add <- dplyr::left_join(ncbi_to_add,filtered_taxonomy[filtered_taxonomy$taxonID %in% ncbi_to_add$taxonID,colnames(filtered_taxonomy)[sapply(colnames(filtered_taxonomy), function(x) !grepl("ncbi", x, ignore.case = TRUE))]]) 
ncbi_to_add <- ncbi_to_add[,colnames(filtered_taxonomy)]
#Add additional entries into taxonomy table.
filtered_taxonomy_expanded <- rbind(filtered_taxonomy,ncbi_to_add)

#Create taxonomic keys for GBIF taxonomies.
filtered_taxonomy_withKeys <- filtered_taxonomy_expanded[filtered_taxonomy_expanded$taxonRank %in% taxonomicRanks,]
for(taxonomicRank in taxonomicRanks){
  taxon_tmp <- filtered_taxonomy_withKeys[filtered_taxonomy_withKeys$taxonRank==taxonomicRank,c("taxonID","taxonomicStatus",taxonomicRank)]
  taxon_tmp <- taxon_tmp[!duplicated(taxon_tmp),]
  taxon_tmp <- taxon_tmp[complete.cases(taxon_tmp),]
  #Preferentially retain taxa with accepted taxonomies.
  taxon_tmp <- taxon_tmp %>% group_by(!!sym(taxonomicRank)) %>% arrange(desc(taxonomicStatus == 'accepted')) %>% slice(1) %>% ungroup()
  # Group by 'taxonID' and taxon, count occurrences
  grouped <- taxon_tmp %>% group_by(!!sym(taxonomicRank),taxonID) %>% summarize(count = n())
  # Find the most common 'vernacularName' for each 'taxonID'
  most_common <- grouped %>% group_by(taxonID) %>% slice_max(order_by = count, n = 1) %>% ungroup()
  #Keep only the first value of 'vernacularName' for each 'taxonID'
  first_common <- most_common %>% group_by(!!sym(taxonomicRank)) %>% slice(1) %>% select(-count)
  first_common$taxonID <- as.numeric(first_common$taxonID)
  taxon_tmp <- first_common
  colnames(taxon_tmp)[which(names(taxon_tmp) == "taxonID")] <- paste(taxonomicRank,"Key",sep="")
  filtered_taxonomy_withKeys <- dplyr::left_join(filtered_taxonomy_withKeys,taxon_tmp,multiple="all")
}

#Get the most common English common names per GBIF taxon ID
Vernacular_Input <- read.table(file=paste(taxonomy_home,"VernacularName.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "", encoding = "UTF-8",na = c("", "NA", "N/A"),allowEscapes=TRUE)
common_names <- Vernacular_Input[Vernacular_Input$language=="en" & !is.na(Vernacular_Input$language),c("taxonID","vernacularName")]
remove <- common_names[grep("^[0-9@#$%^&*()_+{}\":;']+$", common_names$vernacularName), ]
common_names <- common_names[!(common_names$vernacularName %in% remove$vernacularName),]
# Group by 'taxonID' and 'vernacularName', count occurrences
grouped <- common_names %>% group_by(taxonID, vernacularName) %>% summarize(count = n())
# Find the most common 'vernacularName' for each 'taxonID'
most_common <- grouped %>% group_by(taxonID) %>% slice_max(order_by = count, n = 1) %>% ungroup()
#Keep only the first value of 'vernacularName' for each 'taxonID'
first_common <- most_common %>% group_by(taxonID) %>% slice(1)
first_common$taxonID <- as.numeric(first_common$taxonID)

#Get available common names for all taxa and ranks. Merge in common names.
filtered_taxonomy_withNames <- filtered_taxonomy_withKeys
for(taxonomicRank in taxonomicRanks){
  tmp <- first_common[,c("taxonID","vernacularName")]
  taxon_list <- unique(na.omit(filtered_taxonomy_withNames[filtered_taxonomy_withNames$taxonRank %in% taxonomicRank,c("taxonID")]))
  tmp <- tmp[tmp$taxonID %in% taxon_list,]
  names(tmp)[names(tmp) == "vernacularName"] <- paste("common",taxonomicRank,sep="_")
  names(tmp)[names(tmp) == "taxonID"] <- paste(taxonomicRank,"Key",sep="")
  #filtered_taxonomy <- dplyr::left_join(filtered_taxonomy,tmp)
  filtered_taxonomy_withNames <- dplyr::left_join(filtered_taxonomy_withNames,tmp)
}

#Get the most resolved common name.
filtered_taxonomy_withNames <- filtered_taxonomy_withNames %>%
  mutate(Common_Name = coalesce(
    common_species,
    common_genus,
    common_family,
    common_order,
    common_class,
    common_phylum,
    common_kingdom
  ))

#Retain certain columns and then remove duplicate rows.
taxonomy_export <- filtered_taxonomy_withNames[,c("taxonID","canonicalName","taxonRank","kingdom","phylum","class","order","family","genus","species","ncbi_id","ncbi_rank","ncbi_kingdom","ncbi_phylum","ncbi_class","ncbi_order","ncbi_family","ncbi_genus","ncbi_species","speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey","Common_Name")]
taxonomy_export <- taxonomy_export[!is.na(taxonomy_export$ncbi_id),]
taxonomy_export <- taxonomy_export %>% distinct()
#Retain the most common value of the gbif id for a given ncbi id.
taxonomy_export <- taxonomy_export %>%
  group_by(ncbi_id, taxonID) %>%
  summarise(count = n()) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  filter(rank(desc(count), ties.method = "first") == 1) %>%
  ungroup() %>%
  left_join(taxonomy_export, by = c("ncbi_id", "taxonID")) %>%
  select(-count)

#Get phylopic images.
phylopic_export <- taxonomy_export[,c("canonicalName","Common_Name","speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")]
phylopic_export <- phylopic_export[!duplicated(phylopic_export),]
for(i in 1:nrow(phylopic_export)){
  Taxon_Backbone <- as.numeric(phylopic_export[i,c("speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")])
  res <- httr::GET(url=paste("https://api.phylopic.org/resolve/gbif.org/species?embed_primaryImage=true&objectIDs=",paste(Taxon_Backbone,collapse=","),sep=""),config = httr::config(connecttimeout = 100))
  phylopic_query <- fromJSON(rawToChar(res$content))
  Taxon_Image <- phylopic_query[["_embedded"]][["primaryImage"]][["_links"]][["rasterFiles"]][["href"]][[1]]
  if(is.null(Taxon_Image)){
    Taxon_Image <- "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png"
  }
  #print(paste(i,phylopic_export[i,"canonicalName"],phylopic_export[i,"Common_Name"],Taxon_Image))
  phylopic_export[i,"Image_URL"] <- Taxon_Image
  if (i==1) {
    write.table(phylopic_export[i,], "phylopic_export.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    write.table(phylopic_export[i,], "phylopic_export.tsv", sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}
#Read in table containing canonical names, common names, GBIF keys, and phylopic images.
phylopic_export <- read.table(file=paste(taxonomy_home,"phylopic_export.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"",as.is=TRUE, encoding = "UTF-8",na = c("", "NA", "N/A"))
#Merge back in with other taxonomic data.
filtered_taxonomy_export <- dplyr::left_join(phylopic_export,filtered_taxonomy_withNames)

#Get superkingdom data.
filtered_taxonomy_export <- filtered_taxonomy_export %>%
  mutate(superkingdom = case_when(
    kingdom %in% c("Bacteria", "Archaea", "Viruses") ~ kingdom,
    kingdom %in% c("Plantae", "Animalia", "Chromista", "Fungi", "Protozoa") ~ "Eukaryota",
    kingdom == "incertae sedis" ~ NA_character_,
    TRUE ~ NA_character_
  ))

#Prepare columns for export.
colnames(filtered_taxonomy_export)[colnames(filtered_taxonomy_export) == 'canonicalName'] <- 'Taxon'
colnames(filtered_taxonomy_export)[colnames(filtered_taxonomy_export) == 'taxonRank'] <- 'rank'
#Create unique ID for the Phylopic database.
filtered_taxonomy_export$UniqueID <- sapply(paste(filtered_taxonomy_export$species,filtered_taxonomy_export$genus,filtered_taxonomy_export$family,filtered_taxonomy_export$order,filtered_taxonomy_export$class,filtered_taxonomy_export$phylum,filtered_taxonomy_export$kingdom,filtered_taxonomy_export$Taxon,filtered_taxonomy_export$rank,filtered_taxonomy_export$Image_URL,filtered_taxonomy_export$Common_Name),digest,algo="md5")
#Columns to export.
retained_columns <- c("superkingdom","kingdom","phylum","class","order","family","genus","species","Taxon","rank","Common_Name","Image_URL","kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey","UniqueID")
filtered_taxonomy_export <- filtered_taxonomy_export[,retained_columns]
filtered_taxonomy_export <- filtered_taxonomy_export[!duplicated(filtered_taxonomy_export),]

#Save combined NCBI and GBIF taxonomies with common names merged in.
write.table(filtered_taxonomy_export, paste(taxonomy_home,"filtered_taxonomy_export.csv",sep="/"), sep = ",", col.names = TRUE, row.names = FALSE)
con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name, user = db_user, password = db_pass)
#Clear old taxonomy entries.
dbExecute(con,'DELETE FROM "Taxonomy"')
#Write out new taxonomy table
dbWriteTable(con,"Taxonomy",filtered_taxonomy_export,row.names=FALSE,append=TRUE)
sapply(dbListConnections(Database_Driver), dbDisconnect)
