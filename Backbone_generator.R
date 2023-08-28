rm(list=ls())
require(dplyr)
readRenviron(".env")
taxonomy_home <- Sys.getenv("taxonomy_home")
#taxonomy_home <- "~/Desktop/backbone/"

#Get taxon backbone files from the most recent backbone.zip file listed here: https://hosted-datasets.gbif.org/datasets/backbone/current/
Backbone_Input <- read.table(file=paste(taxonomy_home,"Taxon.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "", encoding = "UTF-8",na = c("", "NA", "N/A"))

taxonomy <- Backbone_Input[Backbone_Input$taxonRank %in% c("kingdom","phylum","class","order","family","genus","species") & Backbone_Input$taxonomicStatus=="accepted" & !is.na(Backbone_Input$taxonomicStatus),c("taxonID","taxonRank","kingdom","phylum","class","order","family","genus","canonicalName")]
#taxonomy <- Backbone_Input[Backbone_Input$taxonRank %in% c("kingdom","phylum","class","order","family","genus","species"),c("taxonID","taxonRank","kingdom","phylum","class","order","family","genus","canonicalName")]

Vernacular_Input <- read.table(file=paste(taxonomy_home,"VernacularName.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "", encoding = "UTF-8",na = c("", "NA", "N/A"))

common_names <- Vernacular_Input[Vernacular_Input$language=="en" & !is.na(Vernacular_Input$language),c("taxonID","vernacularName")]

Merged_Taxonomy <- dplyr::left_join(taxonomy,common_names,by=c("taxonID"),multiple="all")

# Group by 'taxonID' and 'vernacularName', count occurrences
grouped <- Merged_Taxonomy %>%
  group_by(taxonID, vernacularName) %>%
  summarize(count = n())

# Find the most common 'vernacularName' for each 'taxonID'
most_common <- grouped %>%
  group_by(taxonID) %>%
  slice_max(order_by = count, n = 1) %>%
  ungroup()

#Keep only the first value of 'vernacularName' for each 'taxonID'
first_common <- most_common %>%
  group_by(taxonID) %>%
  slice(1)

# Extract relevant columns and create the final data frame
result <- first_common %>%
  select(taxonID, vernacularName)

result <- as.data.frame(result)
result$taxonID <- as.numeric(result$taxonID)
result$vernacularName <- gsub("[[:punct:]]", " ", result$vernacularName)

Merged_Taxonomy$vernacularName <- NULL
Merged_Taxonomy <- Merged_Taxonomy[!duplicated(Merged_Taxonomy),]
Merged_Taxonomy <- dplyr::left_join(Merged_Taxonomy,result,by=c("taxonID"))
colnames(Merged_Taxonomy)[colnames(Merged_Taxonomy) == "vernacularName"] ="Common_Name"
colnames(Merged_Taxonomy)[colnames(Merged_Taxonomy) == "taxonRank"] ="rank"
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(kingdomKey = ifelse(rank == "kingdom", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(phylumKey = ifelse(rank == "phylum", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(classKey = ifelse(rank == "class", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(orderKey = ifelse(rank == "order", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(familyKey = ifelse(rank == "family", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(genusKey = ifelse(rank == "genus", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(speciesKey = ifelse(rank == "species", taxonID, NA))
Merged_Taxonomy <- Merged_Taxonomy %>% mutate(species = ifelse(rank == "species", canonicalName, NA))

write.table(Merged_Taxonomy, paste(taxonomy_home,"gbif_backbone.tsv",sep="/"), sep = "\t", col.names = TRUE, row.names = FALSE)

#test <- read.table(file=paste(taxonomy_home,"gbif_backbone.tsv",sep="/"),header=TRUE, sep="\t",skip=0,fill=TRUE,check.names=FALSE,quote = "\"",as.is=TRUE, encoding = "UTF-8",na = c("", "NA", "N/A"))
