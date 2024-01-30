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

process_metadata <- function(con, project_id, has_sites, filter_site_names, sample_first_date, sample_last_date) {
  # Read in metadata and filter it
  metadata <- tbl(con, "TronkoMetadata")

  # Get the number of samples in a project before filtering
  metadata_unfiltered <- metadata %>% filter(projectid == project_id)
  metadata_unfiltered <- as.data.frame(metadata_unfiltered)
  total_samples <- nrow(metadata_unfiltered)

  # Metadata filtering based on sites
  if(has_sites) {
    metadata <- metadata %>%
      filter(projectid == project_id) %>%
      as.data.frame() %>%
      filter(site %in% filter_site_names)
  } else {
    metadata <- metadata %>%
      filter(projectid == project_id) %>%
      as.data.frame()
  }

  # Further processing
  metadata$sample_date <- lubridate::ymd(metadata$sample_date)
  metadata <- metadata %>% filter(sample_date >= sample_first_date & sample_date <= sample_last_date)
  if(nrow(metadata) == 0 || ncol(metadata) == 0) {
    stop("Error: Sample data frame is empty. Cannot proceed.")
  }
  metadata$fastqid <- gsub("_", "-", metadata$fastqid)

  # Return processed metadata and total samples
  list(metadata = metadata, total_samples = total_samples)
}