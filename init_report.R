require(DBI)
require(RPostgreSQL)
require(lubridate)
require(uuid)
require(tidyr)
require(dplyr)
library(sentryR)
library(jsonlite)


tryCatch(
  {
    readRenviron(".env")

    configure_sentry(
      dsn = Sys.getenv("SENTRY_DSN"),
      app_name = "r-report-service", app_version = "1.1.0",
      environment = Sys.getenv("APP_ENV"),
      runtime = NULL
    )
    cat("Configured Sentry: ", Sys.getenv("SENTRY_DSN"), "\n")

    Sys.setenv(
      "AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
      "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY")
    )

    # Define the command and arguments to retrieve the secret
    command <- "aws"
    aws_args <- c(
      "secretsmanager",
      "get-secret-value",
      "--secret-id",
      "prod/ednaExplorer/postgres",
      "--query",
      "SecretString",
      "--output",
      "text"
    )

    # Execute the command and capture the output
    output <- system2(command, aws_args, stdout = TRUE)

    # Parse the JSON output
    parsed_output <- fromJSON(output)

    # Extract the parameter value
    conn_str <- parsed_output$DATABASE_URL

    # Extract the main part of the connection string without the
    # protocol and parameters
    main_conn_str <- sub("postgresql://(.*?)\\?.*", "\\1", conn_str)

    # Split the main part into user:password and the rest
    credentials_rest <- strsplit(main_conn_str, "@")[[1]]

    # Extract user and password
    user_pass <- strsplit(credentials_rest[1], ":")[[1]]
    db_user <- user_pass[1]
    db_pass <- user_pass[2]

    # Extract host and port
    host_port <- strsplit(credentials_rest[2], ":")[[1]]
    db_host <- host_port[1]
    db_port <- as.numeric(strsplit(host_port[2], "/")[[1]][1])

    # Check if db_host is NA
    db_host <- ifelse(is.na(db_port), strsplit(credentials_rest[2], "/")[[1]], db_host)

    # Check if db_port is NA, and default to 5432 if it is
    db_port <- ifelse(is.na(db_port), 5432, db_port)

    # Extract database name
    db_name <- sub(".*/([^/?]+).*", "\\1", conn_str)

    bucket <- Sys.getenv("S3_BUCKET")
    home_dir <- Sys.getenv("home_dir")
    ENDPOINT_URL <- Sys.getenv("ENDPOINT_URL")

    if (length(args) < 1) {
      capture("Need the following inputs: report_id")
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

      has_sites <- !is.null(sites) && sites != "{}" && length(sites) > 0
      if (has_sites) {
        # Get alphabetized site list for a given project.
        selected_site_list <- sites
        project_sites <- tbl(con, "ProjectSite")
        project_sites <- project_sites %>%
          filter(projectId == sample_project_id) %>%
          select(id, name)
        project_sites <- as.data.frame(project_sites)
        filter_sites <- project_sites[project_sites$id %in% selected_site_list, ]

        # Sort sites alphabetically
        filter_sites <- filter_sites[order(filter_sites$id), ]

        # Get site IDs
        # Get site names corresponding to selected site IDs.
        filter_site_names <<- filter_sites$name
      } else {
        print("No sites were selected for filtering.")
      }

      categorical_variables <- c("site", "grtgroup", "biome_type", "iucn_cat", "eco_name", "hybas_id")
      continuous_variables <- c("bio01", "bio12", "ghm", "elevation", "ndvi", "average_radiance")
      field_vars <- c("fastqid", "sample_date", "latitude", "longitude", "spatial_uncertainty")
      taxonomic_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
      taxonomic_key_ranks <- c("kingdomKey", "phylumKey", "classKey", "orderKey", "familyKey", "genusKey", "speciesKey")
      taxonomic_num <<- as.numeric(which(taxonomic_ranks == sample_taxonomic_rank))

      # Save plot name.
      filename <- paste(report_id, "json", sep = ".")
    }

    print("Script initialized.")
  },
  error = capture_exception
)
