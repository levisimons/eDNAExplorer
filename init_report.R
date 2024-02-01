require(DBI)
require(RPostgreSQL)
require(lubridate)
require(uuid)
require(tidyr)
require(dplyr)

readRenviron(".env")
Sys.setenv("AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
           "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY"))
db_host <- Sys.getenv("db_host")
db_port <- Sys.getenv("db_port")
db_name <- Sys.getenv("db_name")
db_user <- Sys.getenv("db_user")
db_pass <- Sys.getenv("db_pass")
bucket <- Sys.getenv("S3_BUCKET")
home_dir <- Sys.getenv("home_dir")
ENDPOINT_URL <- Sys.getenv("ENDPOINT_URL")

if (length(args) < 1) {
  stop("Need the following inputs: report_id", call. = FALSE)
} else {
  database_driver <- dbDriver("PostgreSQL")
  sapply(dbListConnections(database_driver), dbDisconnect)
  con <- dbConnect(database_driver, host = db_host, port = db_port, dbname = db_name, user = db_user, password = db_pass)

  # Extract report_id from args.
  report_id <- args[1]
  # Update report with QUEUED state.
  updateReport(report_id, "QUEUED", con, reset = TRUE)

  # Fetch report data
  sql_query <- sprintf("SELECT * FROM \"Report\" WHERE id = '%s'", report_id)
  report_data <- dbGetQuery(con, sql_query)

  project_id <- report_data$projectId
  first_date <- report_data$firstDate
  last_date <- report_data$lastDate
  marker <- report_data$marker
  num_mismatch <- report_data$numMismatch
  taxonomic_rank <- tolower(report_data$taxonomicRank)
  count_threshold <- report_data$countThreshold
  filter_threshold <- report_data$filterThreshold
  species_list <- report_data$speciesList
  environmental_parameter <- report_data$environmentalParameter
  beta_diversity <- report_data$betaDiversity
  alpha_diversity <- report_data$alphaDiversity
  sites <- report_data$sites
  geographic_scale <- report_data$geographicScale
  
  # Define filters in Phyloseq as global parameters.
  sample_project_id <<- as.character(project_id)
  sample_first_date <<- lubridate::ymd(first_date)
  sample_last_date <<- lubridate::ymd(last_date)
  sample_primer <<- as.character(marker)
  sample_taxonomic_rank <<- as.character(taxonomic_rank)
  sample_num_mismatch <<- as.numeric(num_mismatch)
  sample_count_threshold <<- as.numeric(count_threshold)
  sample_filter_threshold <<- as.numeric(filter_threshold)
  environmental_variable <<- as.character(environmental_parameter)
  beta_diversity_metric <<- as.character(beta_diversity)
  alpha_diversity_metric <<- as.character(alpha_diversity)
  selected_species_list <<- as.character(species_list)
  
  has_sites = !is.null(sites) && sites != "{}" && length(sites) > 0 
  if (has_sites) {
    # Get alphabetized site list for a given project.
    selected_site_list <- sites
    project_sites <- tbl(con, "ProjectSite")
    project_sites <- project_sites %>% filter(projectId == sample_project_id) %>% select(id,name)
    project_sites <- as.data.frame(project_sites)
    filter_sites <- project_sites[project_sites$id %in% selected_site_list,]
    
    # Sort sites alphabetically
    filter_sites <- filter_sites[order(filter_sites$id),]

    # Get site IDs
    # Get site names corresponding to selected site IDs.
    filter_site_names <<- filter_sites$name
  } else {
    print("No sites were selected for filtering.")
  }
  
  categorical_variables <- c("site","grtgroup","biome_type","iucn_cat","eco_name","hybas_id")
  continuous_variables <- c("bio01","bio12","ghm","elevation","ndvi","average_radiance")
  field_vars <- c("fastqid","sample_date","latitude","longitude","spatial_uncertainty")
  taxonomic_ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  taxonomic_key_ranks <- c("kingdomKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")
  taxonomic_num <<- as.numeric(which(taxonomic_ranks == sample_taxonomic_rank))

  # Save plot name.
  filename <- paste(report_id,"json",sep=".")
}

print("Script initialized.")