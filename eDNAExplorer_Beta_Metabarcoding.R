#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(vegan)
require(ggplot2)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
}
require(phyloseq)
require(htmlwidgets)
require(plotly)
require(jsonlite)
require(data.table)
require(digest)

source("helpers.R")
source("init_report.R")

# Rscript --vanilla eDNAExplorer_Beta_Metabarcoding.R "report_id"
tryCatch(
  {
    # Read in species list
    if (selected_species_list != "None") {
      species_list_df <- tbl(con, "species_listItem")
      species_list_df <- species_list_df %>% filter(species_list == selected_species_list)
      species_list_df <- as.data.frame(species_list_df)
    }

    # Read in information to map categorical labels for certain variables.
    category_file <- paste("Categories_", UUIDgenerate(), ".csv", sep = "")
    system(paste("aws s3 cp s3://", bucket, "/analysis/Categories.csv ", category_file, " --endpoint-url ", ENDPOINT_URL, sep = ""))
    categories <- as.data.frame(fread(file = category_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A")))
    system(paste("rm ", category_file, sep = ""))

    # Read in information for legends and labels
    legends_file <- paste("LabelsAndLegends_", UUIDgenerate(), ".csv", sep = "")
    system(paste("aws s3 cp s3://", bucket, "/analysis/LabelsAndLegends.csv ", legends_file, " --endpoint-url ", ENDPOINT_URL, sep = ""))
    legends_and_labels <- as.data.frame(fread(file = legends_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A")))
    system(paste("rm ", legends_file, sep = ""))

    # Set up new legends and x-axis labels.
    new_legend <- legends_and_labels[legends_and_labels$Environmental_Variable == environmental_variable, "Legend"]
    new_axis_label <- legends_and_labels[legends_and_labels$Environmental_Variable == environmental_variable, "x_axis"]

    result <- process_metadata(
      con = con,
      project_id = project_id,
      has_sites = has_sites,
      filter_site_names = filter_site_names,
      sample_first_date = sample_first_date,
      sample_last_date = sample_last_date
    )
    metadata <- result$metadata
    total_samples <- result$total_samples

    # Change values of selected variable
    if (environmental_variable %in% unique(categories$Environmental_Variable)) {
      colnames(metadata)[colnames(metadata) == environmental_variable] <- "value"
      metadata <- dplyr::left_join(metadata, categories[categories$Environmental_Variable == environmental_variable, ])
      colnames(metadata)[colnames(metadata) == "description"] <- environmental_variable
    }

    # Create sample metadata matrix
    if (nrow(metadata) == 0 || ncol(metadata) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    sample <- metadata[!is.na(metadata$fastqid), ]
    rownames(sample) <- sample$fastqid
    sample$fastqid <- NULL
    sample <- sample_data(sample)

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
    subset_file <- paste("subset_beta_", UUIDgenerate(), ".csv", sep = "")
    print("Filtering Tronko output.")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s", tronko_file_tmp, subset_file)
    system(awk_command, intern = TRUE)
    print("Finished filtering Tronko output.")
    tronko_input <- fread(file = subset_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    print("Finished reading Tronko output.")
    tronko_input$Mismatch <- as.numeric(as.character(tronko_input$Mismatch))
    print("Finished converting Tronko output.")

    # Remove samples with missing coordinates, and which are outside of the date filters.
    tronko_input <- tronko_input[tronko_input$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    print("Finished removing samples with missing coordinates.")

    # Store the unfiltered reads.
    tronko_unfiltered <- tronko_input
    print("Finished storing unfiltered reads.")
    print(names(tronko_unfiltered))
    print(paste("Calculating relative abundance of taxa with a given rank in the unfiltered reads.", taxonomic_rank, sep = " "))

    # Calculate relative abundance of taxa with a given rank in the unfiltered reads.
    tronko_unfiltered <- tronko_unfiltered %>%
      dplyr::group_by(SampleID, !!sym(taxonomic_rank)) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>%
      dplyr::ungroup() %>%
      select(-n)
    print("Finished calculating relative abundance of taxa with a given rank in the unfiltered reads.")

    # Filter on the number of reads per sample, then a mismatch threshold.
    tronko_input <- tronko_input %>%
      group_by(SampleID) %>%
      filter(n() > sample_count_threshold) %>%
      filter(Mismatch <= sample_num_mismatch & !is.na(Mismatch)) %>%
      select(-Mismatch)
    print("Finished filtering on the number of reads per sample, then a mismatch threshold.")

    # Merge relative abudance results into data filtered by reads per sample and mismatches
    # and then filtered on relative abundances.
    tronko_db <- tronko_input %>%
      left_join(tronko_unfiltered, na_matches = "never") %>%
      filter(freq >= filter_threshold) %>%
      select(-freq)
    print("Finished merging relative abudance results into data filtered by reads per sample and mismatches.")
    tronko_db <- as.data.frame(tronko_db)
    print("Finished filtering on relative abundances.")

    # Remove duplicated
    tronko_db <- tronko_db[!duplicated(tronko_db), ]
    print("Finished removing duplicated.")

    # Remove taxa which are unknown at a given rank.
    tronko_db <- tronko_db[, c("SampleID", taxonomic_ranks[1:taxonomic_num])]
    tronko_db <- tronko_db[!is.na(tronko_db[, sample_taxonomic_rank]), ]
    print("Finished removing taxa which are unknown at a given rank.")

    # Filter results by species list.
    if (selected_species_list != "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)) & tronko_db$species %in% species_list_df$name, ]
    }
    if (selected_species_list == "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    }
    system(paste("rm ", subset_file, sep = ""))

    if (nrow(tronko_db) > 1) {
      print("Generating OTU matrix.")
      # Create OTU matrix
      otumat <- as.data.frame(pivot_wider(as.data.frame(table(tronko_db[, c("SampleID", sample_taxonomic_rank)])), names_from = SampleID, values_from = Freq))
      rownames(otumat) <- otumat[, sample_taxonomic_rank]
      otumat <- otumat[, colnames(otumat) %in% unique(tronko_db$SampleID)]
      otumat[sapply(otumat, is.character)] <- lapply(otumat[sapply(otumat, is.character)], as.numeric)
      otu <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)
      print("Finished generating OTU matrix.")

      # Create merged Phyloseq object.
      print("Creating merged Phyloseq object.")
      physeq <- phyloseq(otu, sample)
      abundance_filtered <- physeq
      valid_rows <- !is.na(sample_data(abundance_filtered)[[environmental_variable]])
      filtered_phyloseq <- prune_samples(valid_rows, abundance_filtered)
      print("Finished creating merged Phyloseq object.")

      # Plot and analyze beta diversity versus an environmental variables.
      if (nsamples(abundance_filtered) > 1 & ntaxa(abundance_filtered) > 1) {
        print("Calculating beta diversity.")
        if (beta_diversity_metric != "jaccard") {
          print("Calculating beta diversity. (non-jaccard)")
          BetaDist <- phyloseq::distance(filtered_phyloseq, method = beta_diversity_metric, weighted = F)
        }
        if (beta_diversity_metric == "jaccard") {
          print("Calculating beta diversity. (jaccard)")
          BetaDist <- phyloseq::distance(filtered_phyloseq, method = beta_diversity_metric, weighted = F, binary = T)
        }
        if (sum(BetaDist) <= 0) {
          stop(paste("There is not enough data in this project to calculate", beta_diversity_metric, "beta diversity for the current filter settings."))
        }
        if (sum(!is.nan(BetaDist)) > 1) {
          print("Calculating beta diversity. (non-nan)")
          ordination <- ordinate(abundance_filtered, method = "PCoA", distance = BetaDist)
          abundance_filtered_df <- data.frame(sample_data(abundance_filtered))
          print("Calculating beta diversity. (using abundance_filtered_df)")
          if (length(unique(abundance_filtered@sam_data[[environmental_variable]])) > 1) {
            print("Calculating beta diversity. (unique abundance_filtered@sam_data[[environmental_variable]] > 1)")
            print(paste("Calculating beta diversity. (adonis2)", environmental_variable, sep = ", "))
            filtered_data <- sample_data(abundance_filtered)[valid_rows, ]

            # Constructing the formula dynamically
            formula <- paste("BetaDist ~ ", environmental_variable, sep = "")
            browser()
            test <- adonis2(as.formula(formula), data = data.frame(filtered_data))
            stat_test <- paste("PCA plot.  Results of PERMANOVA, using 999 permutations.\n", beta_diversity_metric, " beta diversity and ", new_legend, "\nDegrees of freedom: ", round(test$Df[1], 3), ". Sum of squares: ", round(test$SumOfSqs[1], 3), ". R-squared: ", round(test$R2[1], 3), ". F-statistic: ", round(test$F[1], 3), ". p: ", round(test$`Pr(>F)`[1], 3), sep = "")
          } else {
            print("Calculating beta diversity. (unique abundance_filtered@sam_data[[environmental_variable]] <= 1)")
            stat_test <- paste("PCA plot.  Not enough variation in ", new_legend, " to perform a PERMANOVA on beta diversity.", sep = "")
          }
          print("Calculating beta diversity. (plot_ordination)")
          p <- plot_ordination(abundance_filtered, ordination, color = environmental_variable) + theme(aspect.ratio = 1) + labs(title = stat_test, color = environmental_variable)
          p <- p + theme_bw()
        } else {
          print("There is not enough remaining data after filters to perform a PERMANOVA on beta diversity.")
          stat_test <- "PCA plot.  Not enough remaining data after filters to perform a PERMANOVA on beta diversity."
          p <- ggplot(data.frame()) +
            geom_point() +
            xlim(0, 1) +
            ylim(0, 1) +
            labs(title = stat_test)
        }
      } else {
        print("There is not enough data in this project to calculate beta diversity for the current filter settings.")
        stat_test <- "PCA plot.  Not enough data to perform a PERMANOVA on beta diversity."
        p <- ggplot(data.frame()) +
          geom_point() +
          xlim(0, 1) +
          ylim(0, 1) +
          labs(title = stat_test)
      }
    } else {
      print("There is not enough data in this project to calculate beta diversity for the current filter settings.")
      stat_test <- "PCA plot.  Not enough data to perform a PERMANOVA on beta diversity."
      p <- ggplot(data.frame()) +
        geom_point() +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(title = stat_test)
    }

    print("Inserting the number of samples and number of samples post-filtering as a return object.")
    # Insert the number of samples and number of samples post-filtering as a return object.
    sample_db <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(sample_db) <- c("totalSamples", "filteredSamples")
    sample_db$totalSamples <- total_samples
    sample_db$filteredSamples <- nsamples(abundance_filtered)
    print("Finished inserting the number of samples and number of samples post-filtering as a return object.")
    print("Changing legend label.")
    # Change legend label
    datasets <- list(datasets = list(results = gsub(environmental_variable, new_legend, plotly_json(p, FALSE)), metadata = toJSON(sample_db)))
    print("Finished changing legend label.")

    # Save plot as json object.
    write(toJSON(datasets), filename)
    file_key <- paste("projects/", sample_project_id, "/plots/", filename, sep = "")
    system(paste("aws s3 cp ", filename, " s3://", bucket, "/", file_key, " --endpoint-url ", ENDPOINT_URL, sep = ""), intern = TRUE)
    system(paste("rm ", filename, sep = ""))

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
    print("Error in eDNAExplorer_Beta_Metabarcoding.R:")
    print(e$message)
    process_error(e, report_id, project_id, con)
  }
)
