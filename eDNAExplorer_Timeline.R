# plumber.R
library(plumber)
require(aws.s3)
require(tidyr)
require(dplyr)
require(ggplot2)
require(rgbif)
require(plotly)
require(jsonlite)

#* Echo the parameter that was sent in
#* @param Taxon_name:string Scientific taxon name
#* @get /Tronko_Input

timeline <- function(Taxon_name){
  
  #Select taxon to map.
  #User input
  Taxon <- Taxon_name
  #Get GBIF taxonomy key for taxon.
  Taxon_GBIF <- name_backbone(Taxon)$usageKey
  
  #Read in parsed Tronko database.
  Tronko_Input <- get_object("test-taxa.tsv",bucket="ednaexplorer",region="",as="text",base_url="js2.jetstream-cloud.org:8001") %>% data.table::fread(encoding="UTF-8")
    
  #Coerce date type.
  TronkoDB <- as.data.frame(Tronko_Input)
  TronkoDB$sample_date <- as.Date(TronkoDB$sample_date)
  
  #Convert all character variables to being categorical. Default is that all numeric columns are continuous.
  TronkoDB[sapply(TronkoDB, is.character)] <- lapply(TronkoDB[sapply(TronkoDB, is.character)], as.factor)
  TronkoDB$SampleID <- as.character(TronkoDB$SampleID)
  
  #Find where taxon occurs in Tronko output.
  TaxonDB <- TronkoDB[TronkoDB$usageKey==Taxon_GBIF & !is.na(TronkoDB$usageKey),]
  
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
  #Save plot as json object
  jfig <- plotly_json(p, FALSE)
  return(jfig)
}
