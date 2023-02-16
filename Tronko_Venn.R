rm(list=ls())
require(tidyr)
require(dplyr)
require(sf)
require(sp)
require(VennDiagram)
require(RColorBrewer)

wd <- "~/Desktop/eDNAExplorer/Tronko/"

setwd(wd)

#Read in parsed Tronko-assign output files.  Generated from eDNAExplorer_Initializer.R.
TronkoDB <- suppressWarnings(read.table("Taxa_Parsed.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))

#Set taxonomic rank column names.  Select taxonomic level to aggregate data on.
TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
TaxonomicRank <- "class" #User input to parse Venn diagrams on.

#Select geographic scale to compare eDNA with GBIF.
Geographic_Scales <- c("Local","Ecoregion","Realm")
Geographic_Scale <- "Local" #User input

#Select primer(s) to filter data on.
Primers <- c("MiFish_12S_U") #User input

#Read in complete metadata.
Metadata <- suppressWarnings(read.table("Metadata.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))
Metadata$sample_date <- as.Date(Metadata$sample_date)

#Select date range to filter data on.
First_Date <- as.Date("2021-06-12") #User input
Last_Date <- as.Date("2021-09-12") #User input

#Filter eDNA data by primer and date range.
Samples_Remove <- Metadata[!(Metadata$sample_date >= First_Date & Metadata$sample_date <= Last_Date),"sample_id"]
TronkoDB <- TronkoDB[TronkoDB$Primer %in% Primers,!(colnames(TronkoDB) %in% Samples_Remove)]

#Read in GBIF occurrences.
GBIF_Input <- suppressWarnings(read.table("GBIF.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))
GBIF_Input <- GBIF_Input[!is.na(GBIF_Input$decimalLatitude) & !is.na(GBIF_Input$decimalLongitude),]

#Convert GBIF data to spatial object.
GBIF_Input <- sf::st_as_sf(GBIF_Input,coords=c("decimalLongitude","decimalLatitude"),crs=4326)
#Read in the boundaries for ecoregions and realms
EcoBoundaries <- sf::st_read("Ecoregions2017.shp")
#Get realms present in sample data.
Realms <- EcoBoundaries[EcoBoundaries$REALM %in% unique(Metadata$REALM),]
#Get ecoregions present in sample data.
Eco_Regions <- EcoBoundaries[EcoBoundaries$ECO_NAME %in% unique(Metadata$ECO_NAME),]

if(Geographic_Scale=="Local"){
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
  #Clip GBIF occurrence locations by local boundaries.
  Taxa_Local <- as.data.frame(suppressWarnings(sf::st_intersection(GBIF_Input,Local_Bounds)))
  #Make Venn diagram of GBIF vs. eDNA data.
  p <- ggVennDiagram(list(GBIF=unique(na.omit(Taxa_Local[,TaxonomicRank])),eDNA=unique(na.omit(TronkoDB[,TaxonomicRank]))))
  p <- p + scale_fill_distiller(palette="YlGn",direction=1) + theme(legend.position = "none")
  ggsave(file="Tronko_GBIF_Venn.svg", plot=p, width=8, height=8)
  #Text output on species found in either data set at the local scale.
  GBIF_Only <- unique(na.omit(Taxa_Local[,TaxonomicRank]))
  eDNA_Only <- unique(na.omit(TronkoDB[,TaxonomicRank]))
  Shared_Taxa <- intersect(GBIF_Only,eDNA_Only)
}

if(Geographic_Scale=="Ecoregion"){
  #Clip GBIF occurrence locations by ecoregions.
  Taxa_Ecoregion <- as.data.frame(suppressWarnings(sf::st_intersection(GBIF_Input,Eco_Regions)))
  #Make Venn diagram of GBIF vs. eDNA data.
  p <- ggVennDiagram(list(GBIF=unique(na.omit(Taxa_Ecoregion[,TaxonomicRank])),eDNA=unique(na.omit(TronkoDB[,TaxonomicRank]))))
  p <- p + scale_fill_distiller(palette="YlGn",direction=1) + theme(legend.position = "none")
  ggsave(file="Tronko_GBIF_Venn.svg", plot=p, width=8, height=8)
  #Text output on species found in either data set at the ecoregion scale.
  GBIF_Only <- unique(na.omit(Taxa_Ecoregion[,TaxonomicRank]))
  eDNA_Only <- unique(na.omit(TronkoDBp,TaxonomicRank))
  Shared_Taxa <- intersect(GBIF_Only,eDNA_Only)
}

if(Geographic_Scale=="Realm"){
  #Clip GBIF occurrence locations by realms.
  Taxa_Realm <- as.data.frame(suppressWarnings(sf::st_intersection(GBIF_Input,Realms[7:10,])))
  #Make Venn diagram of GBIF vs. eDNA data.
  p <- ggVennDiagram(list(GBIF=unique(na.omit(Taxa_Realm[,TaxonomicRank])),eDNA=unique(na.omit(TronkoDB[,TaxonomicRank]))))
  p <- p + scale_fill_distiller(palette="YlGn",direction=1) + theme(legend.position = "none")
  ggsave(file="Tronko_GBIF_Venn.svg", plot=p, width=8, height=8)
  #Text output on species found in either data set at the ecoregion scale.
  GBIF_Only <- unique(na.omit(Taxa_Realm[,TaxonomicRank]))
  eDNA_Only <- unique(na.omit(TronkoDB[,TaxonomicRank]))
  Shared_Taxa <- intersect(GBIF_Only,eDNA_Only)
}
