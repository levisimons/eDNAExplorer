updateReport <- function(report_id, newBuildState, con, reset = FALSE) {
  library(DBI)
  # Get the current Unix timestamp
  currentTimestamp <- as.integer(Sys.time())

  # Initialize sql_query variable
  sql_query <- ""

  if (reset) {
    # If reset is TRUE, reset executionTimes with current timestamp
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"buildState\" = '%s', \"executionTimes\" = ARRAY[%d],  \"error\" = NULL, \"errorType\" = NULL, \"updatedAt\" = NOW() WHERE \"id\" = '%s'",
      newBuildState, currentTimestamp, report_id
    )
  } else {
    # If reset is FALSE, append current timestamp to executionTimes
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"buildState\" = '%s', \"executionTimes\" = array_append(\"executionTimes\", %d),  \"updatedAt\" = NOW() WHERE \"id\" = '%s'",
      newBuildState, currentTimestamp, report_id
    )
  }
  print(paste("Updating report", report_id, "with build state", newBuildState, "and timestamp", currentTimestamp, sep = " "))

  # Execute the query
  dbExecute(con, sql_query)
}

# Write error output to our json file.
process_error <- function(e, report_id, project_id, con, errorType = "error") {
  library(DBI)
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  updateReport(report_id, "FAILED", con)

  # Assuming `conn` is your database connection
  sql_query <- "UPDATE \"Report\" SET \"error\" = $1 WHERE \"id\" = $2"

  # Execute the query with parameters
  dbExecute(con, sql_query, list(e$message, report_id))

  stop(error_message)
}

download_primer <- function(primer, project_id, bucket) {
  library(uuid)

  # Constructing filenames
  original_file_name <- paste(project_id, "_", primer, ".csv", sep = "")
  temp_file_name <- paste(primer, "_", UUIDgenerate(), ".tmp", sep = "")

  # Check if the file already exists
  if (!file.exists(original_file_name)) {
    # File doesn't exist, proceed with download
    s3_download_cmd <- paste("aws s3 cp s3://", bucket, "/tronko_output/", project_id, "/", primer, ".csv ", temp_file_name, " --endpoint-url ", ENDPOINT_URL, sep = "")
    system(s3_download_cmd)

    # Rename the downloaded file
    file.rename(temp_file_name, original_file_name)
  } else {
    # File already exists
    print(paste("File", original_file_name, "already exists. Skipping download."))
  }

  return(original_file_name)
}

process_metadata <- function(con, project_id, has_sites, filter_site_names, sample_first_date, sample_last_date, environmental_variable) {
  # Read in metadata and filter it
  metadata <- tbl(con, "TronkoMetadata")

  # Get the number of samples in a project before filtering
  sql_query <- sprintf("SELECT COUNT(*) AS \"total_samples\" FROM \"TronkoMetadata\" WHERE \"projectid\" = '%s'", project_id)
  result <- dbGetQuery(con, sql_query)
  total_samples <- result$total_samples

  if (is.null(environmental_variable) ||
    (is.character(environmental_variable) && (length(environmental_variable) == 0 || any(environmental_variable == "")))) {
    environmental_variable <- "grtgroup"
  }

  # Metadata filtering based on sites
  # Common operations: filtering by project_id and selecting specific columns
  metadata <- metadata %>%
    filter(projectid == project_id) %>%
    select(projectid, all_of(environmental_variable), sample_date, sample_id, sample_type, fastqid) %>%
    as.data.frame()

  # Further processing based on additional conditions (e.g., has_sites)
  if (has_sites) {
    metadata <- metadata %>%
      filter(site %in% filter_site_names)
  }

  # Transform sample_date to a date object
  metadata$sample_date <- lubridate::ymd(metadata$sample_date)

  # Additional filtering based on sample dates, if required
  if (!is.null(sample_first_date) && !is.null(sample_last_date)) {
    metadata <- metadata %>%
      filter(sample_date >= sample_first_date & sample_date <= sample_last_date)
  }

  if (nrow(metadata) == 0 || ncol(metadata) == 0) {
    stop("Error: Sample data frame is empty. Cannot proceed.")
  }
  metadata$fastqid <- gsub("_", "-", metadata$fastqid)

  # Return processed metadata and total samples
  list(metadata = metadata, total_samples = total_samples)
}

profile_statement <- (function() {
  last_time <<- NULL

  function(statement, restart_interval = FALSE) {
    current_time <- Sys.time()

    if (!is.null(last_time) && !restart_interval) {
      elapsed <- difftime(current_time, last_time, units = "secs")
      cat(sprintf("%s, %.2fs\n", statement, elapsed))
    } else {
      cat(sprintf("%s\n", statement))
    }

    last_time <<- current_time
  }
})()
