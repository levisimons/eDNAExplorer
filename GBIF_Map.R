rm(list=ls())
require(tidyr)
require(dplyr)
require(ggplot2)
require(plotly)
require(rgbif)
require(leaflet)
require(htmlwidgets)

wd <- "~/Desktop/eDNAExplorer/Tronko/"

setwd(wd)

#Read in parsed Tronko-assign output files.  Generated from eDNAExplorer_Initializer.R.
TronkoDB <- suppressWarnings(read.table("Taxa_Parsed.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))

#Set taxonomic rank column names.
TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")

#Select taxon to map.
#User input
Taxon <- "Bos taurus"
#Get GBIF taxonomy key for taxon.
Taxon_GBIF <- name_backbone(Taxon)$usageKey

#Read in full metadata.
Metadata <- suppressWarnings(read.table("Metadata.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8"))

#Find where taxon occurs in Tronko output.
TaxonDB <- dplyr::filter_all(TronkoDB, any_vars(. ==Taxon))
TaxonDB <- TaxonDB %>% dplyr::mutate_at(Metadata$sample_id, as.numeric)

#Get taxon locations
Sample_Occurrences <- TaxonDB[,Metadata$sample_id]
Sample_Occurrences <- Sample_Occurrences %>% select_if(~any(. > 0))
Taxon_Locations <- Metadata[Metadata$sample_id %in% colnames(Sample_Occurrences),c("longitude","latitude")]

# create style raster layer 
projection = '3857' # projection code
style = 'style=osm-bright' # map style
tileRaster = paste0('https://tile.gbif.org/',projection,'/omt/{z}/{x}/{y}@4x.png?',style)
# create our polygons layer 
prefix = 'https://api.gbif.org/v2/map/occurrence/density/{z}/{x}/{y}@4x.png?'
polygons = 'style=classic.poly&bin=hex&hexPerTile=70' # ploygon styles 
taxonKey = paste("taxonKey=",Taxon_GBIF,sep="")
tilePolygons = paste0(prefix,polygons,'&',taxonKey)
# plot the styled map
map <- leaflet() %>%
  setView(lng = 5.4265362, lat = 43.4200248, zoom = 01) %>%
  addTiles(urlTemplate=tileRaster) %>%
  addTiles(urlTemplate=tilePolygons) %>%
  addMarkers(lng=Taxon_Locations$longitude,lat=Taxon_Locations$latitude)
saveWidget(map,file="Taxa_Map.html")

#Get GBIF occurrences over time.
Taxa_Time <- data.frame()
for(GBIF_Year in 1950:as.integer(format(Sys.Date(), "%Y"))){
  tmp <- data.frame(matrix(nrow=1,ncol=2))
  colnames(tmp) <- c("Year","Occurrences")
  tmp$Year <- GBIF_Year
  tmp$Occurrences <- occ_count(taxonKey=Taxon_GBIF,georeferenced=TRUE,year=GBIF_Year)
  Taxa_Time <- rbind(Taxa_Time,tmp)
}
#Plot GBIF occurrences over time.
p <- ggplot(data=Taxa_Time,aes(x=Year,y=Occurrences))+geom_bar(stat="identity",fill="gold")
#Save timeline output.
write.table(Taxa_Time,"Taxa_Timeline.txt",quote=FALSE,sep="\t",row.names = FALSE)
saveWidget(ggplotly(p, tooltip = "all"),file="Taxa_Timeline.html")
