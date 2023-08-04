# plumber.R
library(plumber)
require(tidyr)
require(dplyr)
require(lubridate)
require(jsonlite)
require(data.table)
require(DBI)
require(RPostgreSQL)
require(digest)

#Establish database credentials.
readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")

#* Echo the parameter that was sent in
#* @param ProjectID:string
#* @param First_Date:string YYYY-MM-DD
#* @param Last_Date:string YYYY-MM-DD
#* @param Marker:string Target marker name
#* @param Num_Mismatch:numeric Maximum number of sequence mismatches allowed with Tronko-assign output
#* @param TaxonomicRank:string Taxonomic level to aggregate results to
#* @param CountThreshold:numeric Read count threshold for retaining samples
#* @param FilterThreshold:numeric Choose a threshold for filtering ASVs prior to analysis
#* @param SpeciesList:string Name of csv file containing selected species list.
#* @get /prevalence

prevalence <- function(ProjectID,First_Date,Last_Date,Marker,Num_Mismatch,TaxonomicRank,CountThreshold,FilterThreshold,SpeciesList){
  CategoricalVariables <- c("grtgroup","biome_type","iucn_Cat","eco_name","hybas_id")
  ContinuousVariables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
  FieldVars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
  TaxonomicRanks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  First_Date <- lubridate::ymd(First_Date)
  Last_Date <- lubridate::ymd(Last_Date)
  Num_Mismatch <- as.numeric(Num_Mismatch)
  CountThreshold <- as.numeric(CountThreshold)
  FilterThreshold <- as.numeric(FilterThreshold)
  SelectedSpeciesList <- as.character(SpeciesList)
  Project_ID <- as.character(ProjectID)
  
  #Generate the output filename for cached plots.
  filename <- paste("Prevalence_Metabarcoding_FirstDate",First_Date,"LastDate",Last_Date,"Marker",Marker,"Rank",TaxonomicRank,"Mismatch",Num_Mismatch,"CountThreshold",CountThreshold,"AbundanceThreshold",format(FilterThreshold,scientific=F),"SpeciesList",SelectedSpeciesList,sep="_")
  filename <- paste(filename,".json",sep="")
  filename <- tolower(filename)
  #Output a blank json output for plots as a default.  This gets overwritten is actual plot material exists.
  write(toJSON(data.frame(error=c("No results found"))),filename)
  system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
  system(paste("rm ",filename,sep=""))
  
  #Establish sql connection
  Database_Driver <- dbDriver("PostgreSQL")
  sapply(dbListConnections(Database_Driver), dbDisconnect)
  con <- dbConnect(Database_Driver,host = db_host,port = db_port,dbname = db_name,user = db_user,password = db_pass)
  
  #Read in species list
  if(SelectedSpeciesList != "None"){
    SpeciesList_df <- tbl(con,"SpeciesListItem")
    SpeciesList_df <- SpeciesList_df %>% filter(species_list == SelectedSpeciesList)
    SpeciesList_df <- as.data.frame(SpeciesList_df)
  }
  
  #Read in metadata and filter it.
  Metadata <- tbl(con,"TronkoMetadata")
  Keep_Vars <- c(CategoricalVariables,ContinuousVariables,FieldVars)[c(CategoricalVariables,ContinuousVariables,FieldVars) %in% dbListFields(con,"TronkoMetadata")]
  Metadata <- Metadata %>% filter(sample_date >= First_Date & sample_date <= Last_Date) %>%
    filter(projectid == Project_ID) %>% filter(!is.na(latitude) & !is.na(longitude))
  Metadata <- as.data.frame(Metadata)
  Metadata$fastqid <- gsub("_","-",Metadata$fastqid)
  
  #Read in Tronko output and filter it.
  TronkoFile <- paste(Marker,".csv",sep="")
  system(paste("aws s3 cp s3://ednaexplorer/tronko_output/",Project_ID,"/",TronkoFile," ",TronkoFile," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""))
  #Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
  system(paste("cut -d ',' -f 6,7,8,9,10,11,12,13,14,16 ",TronkoFile," > subset1.csv",sep=""))
  #Filter on the number of mismatches.  Remove entries with NA for mismatches and for the selected taxonomic rank.
  TaxonomicNum <- which(TaxonomicRanks==TaxonomicRank)+1
  system(paste("awk -F',' 'NR == 1 || ( $10 != \"\" && $",TaxonomicNum," != \"\")' subset1.csv > subset2.csv",sep=""))
  TronkoInput <- fread(file="subset2.csv",header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A"))
  if(TaxonomicRank != "species"){
    TronkoInput <- TronkoInput %>% filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > CountThreshold) %>% 
      select(SampleID,species,kingdom,TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    if(SelectedSpeciesList != "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)) & TronkoDB$species %in% SpeciesList_df$name,]}
    if(SelectedSpeciesList == "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)),]}
    TronkoDB$species <- NULL
  } else{
    TronkoInput <- TronkoInput %>% filter(Mismatch <= Num_Mismatch & !is.na(Mismatch)) %>% filter(!is.na(!!sym(TaxonomicRank))) %>%
      group_by(SampleID) %>% filter(n() > CountThreshold) %>% 
      select(SampleID,kingdom,TaxonomicRank)
    TronkoDB <- as.data.frame(TronkoInput)
    if(SelectedSpeciesList != "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)) & TronkoDB$species %in% SpeciesList_df$name,]}
    if(SelectedSpeciesList == "None"){TronkoDB <- TronkoDB[TronkoDB$SampleID %in% unique(na.omit(Metadata$fastqid)),]}
  }
  system(paste("rm",TronkoFile,sep=" "))
  system("rm subset*.csv")
  
  #Read in Taxonomy output and filter it.
  TaxonomyInput <- tbl(con,"Taxonomy")
  if(nrow(TronkoDB)>0){TaxaList <- na.omit(unique(TronkoDB[,TaxonomicRank]))}
  if(nrow(TronkoDB)==0){TaxaList <- c()}
  if(TaxonomicRank!="kingdom"){TaxonomyInput <- TaxonomyInput %>% filter(rank == TaxonomicRank) %>% filter(Taxon %in% TaxaList) %>% select(Taxon,Common_Name,Image_URL)}
  TaxonomyDB <- as.data.frame(TaxonomyInput)
  colnames(TaxonomyDB) <- c(TaxonomicRank,"Common_Name","Image_URL")
  if(TaxonomicRank=="kingdom"){
    TaxonomyDB <- data.frame(kingdom=c("Fungi","Plantae","Animalia","Bacteria","Archaea","Protista","Monera","Chromista"),
                             Image_URL=c("https://images.phylopic.org/images/7ebbf05d-2084-4204-ad4c-2c0d6cbcdde1/raster/958x1536.png",
                                         "https://images.phylopic.org/images/573bc422-3b14-4ac7-9df0-27d7814c099d/raster/1052x1536.png",
                                         "https://images.phylopic.org/images/0313dc90-c1e2-467e-aacf-0f7508c92940/raster/681x1536.png",
                                         "https://images.phylopic.org/images/d8c9f603-8930-4973-9a37-e9d0bc913a6b/raster/1536x1128.png",
                                         "https://images.phylopic.org/images/7ccfe198-154b-4a2f-a7bf-60390cfe6135/raster/1177x1536.png",
                                         "https://images.phylopic.org/images/4641171f-e9a6-4696-bdda-e29bc4508538/raster/336x1536.png",
                                         "https://images.phylopic.org/images/018ee72f-fde6-4bc3-9b2e-087d060ee62d/raster/872x872.png",
                                         "https://images.phylopic.org/images/1fd55f6f-553c-4838-94b4-259c16f90c31/raster/1054x1536.png"))
  }
  #Filter by relative abundance per taxon per sample.
  if(nrow(TronkoDB)>0){
    TronkoDB <- TronkoDB[!is.na(TronkoDB[,TaxonomicRank]),]
    KingdomMatch <- TronkoDB[,c("kingdom",TaxonomicRank)]
    KingdomMatch <- KingdomMatch[!duplicated(KingdomMatch),]
    KingdomMatch <- KingdomMatch[!is.na(KingdomMatch[,"kingdom"]),]
    KingdomMatch <- as.data.frame(KingdomMatch)
    TronkoDB <- TronkoDB %>% dplyr::group_by(SampleID,!!sym(TaxonomicRank)) %>% 
      dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n/sum(n)) %>% 
      dplyr::ungroup() %>% dplyr::filter(freq > FilterThreshold) %>% select(-n,-freq)
    num_filteredSamples <- length(unique(TronkoDB$SampleID))
    TronkoDB <- TronkoDB %>% dplyr::group_by(!!sym(TaxonomicRank)) %>% dplyr::summarise(per=n()/length(unique(TronkoDB$SampleID)))
    TronkoDB <- as.data.frame(TronkoDB)
    if(TaxonomicRank!="kingdom"){TronkoDB <- dplyr::left_join(TronkoDB,KingdomMatch)}
    TronkoDB <- dplyr::left_join(TronkoDB,TaxonomyDB)
    if(TaxonomicRank!="kingdom"){
      colnames(TronkoDB)[which(names(TronkoDB)==TaxonomicRank)] <- "Latin_Name"
    }
    if(TaxonomicRank=="kingdom"){
      TronkoDB$Latin_Name <- TronkoDB$kingdom
    }
    #Insert the number of samples and number of samples post-filtering as a return object.
    SampleDB <- data.frame(matrix(ncol=2,nrow=1))
    colnames(SampleDB) <- c("totalSamples","filteredSamples")
    SampleDB$totalSamples <- nrow(Metadata)
    SampleDB$filteredSamples <- num_filteredSamples
    datasets <- list(datasets = list(results=TronkoDB,metadata=SampleDB))
    write(toJSON(datasets),filename)
    system(paste("aws s3 cp ",filename," s3://ednaexplorer/projects/",Project_ID,"/plots/",filename," --endpoint-url https://js2.jetstream-cloud.org:8001/",sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
  }
  sapply(dbListConnections(Database_Driver), dbDisconnect)
}
