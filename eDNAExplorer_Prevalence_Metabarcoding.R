#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(jsonlite)
require(data.table)
require(digest)

source("helpers.R")
source("init_report.R")

# Rscript --vanilla eDNAExplorer_Prevalence_Metabarcoding.R "report_id"

tryCatch(
  {
    print(paste("Prevalence script started."))

    # Read in species list
    if (selected_species_list != "None") {
      print(paste("Species list selected:", selected_species_list, sep = " "))
      species_list_df <- tbl(con, "SpeciesListItem")
      species_list_df <- species_list_df %>% filter(species_list == selected_species_list)
      species_list_df <- as.data.frame(species_list_df)
    }

    print(paste("Connecting to database", db_host, db_port, db_name, db_user, sep = ", "))
    result <- process_metadata(
      con = con,
      project_id = project_id,
      has_sites = has_sites,
      filter_site_names = filter_site_names,
      sample_first_date = sample_first_date,
      sample_last_date = sample_last_date,
      environmental_variable = environmental_parameter
    )
    print(paste("Metadata processed", result$total_samples, sep = ", "))
    metadata <- result$metadata
    total_samples <- result$total_samples

    # Read in Tronko output and filter it.
    updateReport(report_id, "LOADING", con)
    tronko_file_tmp <- download_primer(sample_primer, sample_project_id, bucket)

    # Throw if the tronko file is empty
    if (file.info(tronko_file_tmp)$size == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }

    # Primer has been downloaded loaded, update the report record to the BUILDING state.
    updateReport(report_id, "BUILDING", con)

    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    subset_file <- paste("subset_prevalence_", UUIDgenerate(), ".csv", sep = "")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s", tronko_file_tmp, subset_file)
    system(awk_command, intern = TRUE)
    tronko_input <- fread(file = subset_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    tronko_input$Mismatch <- as.numeric(as.character(tronko_input$Mismatch))

    # Remove samples with missing coordinates, and which are outside of the date filters.
    tronko_input <- tronko_input[tronko_input$SampleID %in% unique(na.omit(metadata$fastqid)), ]

    # Store the unfiltered reads.
    tronko_unfiltered <- tronko_input

    # Calculate relative abundance of taxa with a given rank in the unfiltered reads.
    tronko_unfiltered <- tronko_unfiltered %>%
      dplyr::group_by(SampleID, !!sym(taxonomic_rank)) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>%
      dplyr::ungroup() %>%
      select(-n)

    # Filter on the number of reads per sample, then a mismatch threshold.
    tronko_input <- tronko_input %>%
      group_by(SampleID) %>%
      filter(n() > count_threshold) %>%
      filter(Mismatch <= num_mismatch & !is.na(Mismatch)) %>%
      select(-Mismatch)

    # Merge relative abudance results into data filtered by reads per sample and mismatches.
    # and then filtered on relative abundances.
    tronko_input <- tronko_input %>%
      left_join(tronko_unfiltered, na_matches = "never") %>%
      filter(freq >= filter_threshold) %>%
      select(-freq)
    tronko_db <- as.data.frame(tronko_input)

    # Remove taxa which are unknown at a given rank.
    tronko_db <- tronko_db[, c("SampleID", taxonomic_ranks[1:taxonomic_num])]
    tronko_db <- tronko_db[!is.na(tronko_db[, taxonomic_rank]), ]

    # Filter results by species list.
    if (selected_species_list != "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)) & tronko_db$species %in% species_list_df$name, ]
    }
    if (selected_species_list == "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    }
    system(paste("rm ", subset_file, sep = ""))

    # Read in Taxonomy output and filter it.
    taxonomy_input <- tbl(con, "Taxonomy")
    if (nrow(tronko_db) > 0) {
      taxa_list <- na.omit(unique(tronko_db[, taxonomic_rank]))
    }
    if (nrow(tronko_db) == 0) {
      taxa_list <- c()
    }
    if (taxonomic_rank != "kingdom") {
      taxonomy_input <- taxonomy_input %>%
        filter(rank == taxonomic_rank) %>%
        filter(!!sym(taxonomic_rank) %in% taxa_list) %>%
        select(taxonomic_ranks[2:taxonomic_num], taxonomic_key_ranks, Common_Name, Image_URL)
      taxonomy_db <- as.data.frame(taxonomy_input)

      # Figure out which taxonomy version is more complete.
      taxonomy_db$rankCount <- rowSums(!is.na(taxonomy_db[, colnames(taxonomy_db) %in% taxonomic_ranks]))
      taxonomy_db$rankKeyCount <- rowSums(!is.na(taxonomy_db[, colnames(taxonomy_db) %in% taxonomic_key_ranks]))
      taxonomy_db <- taxonomy_db %>%
        group_by(!!sym(taxonomic_rank)) %>%
        slice_max(order_by = rankCount, n = 1) %>%
        ungroup()
      taxonomy_db <- taxonomy_db %>%
        group_by(!!sym(taxonomic_rank)) %>%
        slice_max(order_by = rankKeyCount, n = 1) %>%
        ungroup()

      # Figure out which common_name is most common per taxon.
      taxonomy_db <- taxonomy_db %>%
        group_by(!!sym(taxonomic_rank)) %>%
        mutate(Most_Common_Name = ifelse(all(is.na(Common_Name)), NA, names(which.max(table(Common_Name[!is.na(Common_Name)]))))) %>%
        ungroup()
      taxonomy_db$Common_Name <- taxonomy_db$Most_Common_Name
      taxonomy_db$Most_Common_Name <- NULL
      taxonomy_db <- as.data.frame(taxonomy_db)
      taxonomy_db <- subset(taxonomy_db, select = -grep("Key", colnames(taxonomy_db)))
      taxonomy_db$rankCount <- NULL
      taxonomy_db <- taxonomy_db[!duplicated(taxonomy_db), ]
    }
    if (taxonomic_rank == "kingdom") {
      taxonomy_db <- data.frame(
        kingdom = c("Fungi", "Plantae", "Animalia", "Bacteria", "Archaea", "Protista", "Monera", "Chromista"),
        Image_URL = c(
          "https://images.phylopic.org/images/7ebbf05d-2084-4204-ad4c-2c0d6cbcdde1/raster/958x1536.png",
          "https://images.phylopic.org/images/573bc422-3b14-4ac7-9df0-27d7814c099d/raster/1052x1536.png",
          "https://images.phylopic.org/images/0313dc90-c1e2-467e-aacf-0f7508c92940/raster/681x1536.png",
          "https://images.phylopic.org/images/d8c9f603-8930-4973-9a37-e9d0bc913a6b/raster/1536x1128.png",
          "https://images.phylopic.org/images/7ccfe198-154b-4a2f-a7bf-60390cfe6135/raster/1177x1536.png",
          "https://images.phylopic.org/images/4641171f-e9a6-4696-bdda-e29bc4508538/raster/336x1536.png",
          "https://images.phylopic.org/images/018ee72f-fde6-4bc3-9b2e-087d060ee62d/raster/872x872.png",
          "https://images.phylopic.org/images/1fd55f6f-553c-4838-94b4-259c16f90c31/raster/1054x1536.png"
        )
      )
    }

    # Calculate prevalence of taxa per sample and output results for plotting.
    if (nrow(tronko_db) > 0) {
      num_filtered_samples <- length(unique(tronko_db$SampleID))
      tronko_db <- tronko_db %>%
        dplyr::group_by(!!sym(taxonomic_rank)) %>%
        dplyr::summarise(per = n_distinct(SampleID) / n_distinct(tronko_db$SampleID))
      tronko_db <- as.data.frame(tronko_db)

      # Merge in taxonomy data.
      tronko_db <- dplyr::left_join(tronko_db, taxonomy_db, na_matches = "never")
      tronko_db$Image_URL <- ifelse(is.na(tronko_db$Image_URL), "https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png", tronko_db$Image_URL)
      if (taxonomic_rank != "kingdom") {
        colnames(tronko_db)[which(names(tronko_db) == taxonomic_rank)] <- "Latin_Name"
      }
      if (taxonomic_rank == "kingdom") {
        tronko_db$Latin_Name <- tronko_db$kingdom
      }
      tronko_db <- tronko_db %>% mutate_all(~ ifelse(is.na(.), " ", .))

      # De-duplicated dataframe for tronko table.
      tronko_db <- tronko_db[!duplicated(tronko_db), ]

      # Insert the number of samples and number of samples post-filtering as a return object.
      sample_db <- data.frame(matrix(ncol = 2, nrow = 1))
      colnames(sample_db) <- c("totalSamples", "filteredSamples")
      sample_db$totalSamples <- total_samples
      sample_db$filteredSamples <- num_filtered_samples

      datasets <- list(datasets = list(results = tronko_db, metadata = sample_db))
      write(toJSON(datasets), filename)
      file_key <- paste("projects/", project_id, "/plots/", filename, sep = "")
      system(paste("aws s3 cp ", filename, " s3://", bucket, "/", file_key, " --endpoint-url ", ENDPOINT_URL, sep = ""), intern = TRUE)
      system(paste("rm ", filename, sep = ""))
    }
    if (nrow(tronko_db) == 0) {
      stop("Error: Filters are too stringent. Cannot proceed.")
    }

    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"fileKey\" = '%s' WHERE \"id\" = '%s'",
      file_key, report_id
    )
    print(paste("Updating report", report_id, "with file key", file_key, sep = " "))
    dbExecute(con, sql_query)
    updateReport(report_id, "COMPLETED", con)
    RPostgreSQL::dbDisconnect(con, shutdown = TRUE)
  },
  error = function(e) {
    process_error(e, report_id, project_id, con)
  }
)
