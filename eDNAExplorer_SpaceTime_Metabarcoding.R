#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(jsonlite)
require(data.table)
require(digest)
require(zoo)

source("helpers.R")
source("init_report.R")

# Rscript --vanilla eDNAExplorer_Spacetime_Metabarcoding.R "report_id"

tryCatch(
  { 
    # Read in species list
    if(selected_species_list != "None"){
      SpeciesList_df <- tbl(con,"SpeciesListItem")
      SpeciesList_df <- SpeciesList_df %>% filter(species_list == selected_species_list)
      SpeciesList_df <- as.data.frame(SpeciesList_df)
    }
    
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
    
    # Read in Tronko output and filter it.
    updateReport(report_id, "LOADING", con)
    tronko_file_tmp = download_primer(sample_primer, sample_project_id, bucket)
 
    # Throw if the tronko file is empty
    if(file.info(tronko_file_tmp)$size== 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }

    # Primer has been downloaded loaded, update the report record to the BUILDING state.
    updateReport(report_id, "BUILDING", con)
    
    # Select relevant columns in bash (SampleID, taxonomic ranks, Mismatch)
    subset_file <- paste("subset_spacetime_",UUIDgenerate(),".csv",sep="")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",tronko_file_tmp, subset_file)
    system(awk_command, intern = TRUE)
    tronko_input <- fread(file=subset_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    tronko_input$Mismatch <- as.numeric(as.character(tronko_input$Mismatch))
    summarise(tronko_input)
    
    # Remove samples with missing coordinates, and which are outside of the date filters.
    tronko_input <- tronko_input <- tronko_input[tronko_input$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    summarise(tronko_input)

    # Store the unfiltered reads.
    print("Storing the unfiltered reads.")
    tronko_unfiltered <- tronko_input
    summarise(tronko_unfiltered)

    # Calculate relative abundance of taxa with a given rank in the unfiltered reads.
    print("Calculating relative abundance of taxa with a given rank in the unfiltered reads.")
    tronko_unfiltered <- tronko_unfiltered %>%
      dplyr::group_by(SampleID, !!sym(taxonomic_rank)) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>%
      dplyr::ungroup() %>% select(-n)
    
    # Filter on the number of reads per sample, then a mismatch threshold.
    print("Filtering on number of reads per sample and mismatch threshold")
    tronko_input <- tronko_input %>% group_by(SampleID) %>%
      filter(n() > count_threshold) %>% filter(Mismatch <= num_mismatch & !is.na(Mismatch)) %>% select(-Mismatch)
    
    # Merege relative abudance results into data filtered by reads per sample and mismatches
    # and then filtered on relative abundances.
    tronko_input <- tronko_input %>%
      left_join(tronko_unfiltered,na_matches="never") %>%
      filter(freq >= filter_threshold) %>% select(-freq)
    tronko_db <- as.data.frame(tronko_input)
    
    # Remove taxa which are unknown at a given rank.
    tronko_db <- tronko_db[,c("SampleID",taxonomic_ranks[1:taxonomic_num])]
    tronko_db <- tronko_db[!is.na(tronko_db[, taxonomic_rank]), ]
    
    # Filter results by species list.
    if (selected_species_list != "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)) & tronko_db$species %in% SpeciesList_df$name, ]
    }
    if (selected_species_list == "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    }
    system(paste("rm ",subset_file,sep=""))
    
    # Read in Taxonomy output and filter it.
    TaxonomyInput <- tbl(con, "Taxonomy")
    if (nrow(tronko_db) > 0) {
      TaxaList <- na.omit(unique(tronko_db[, taxonomic_rank]))
    }
    if (nrow(tronko_db) == 0) {
      TaxaList <- c()
    }
    if (taxonomic_rank != "kingdom") {
      TaxonomyInput <- TaxonomyInput %>%
        filter(rank==taxonomic_rank) %>%
        filter(!!sym(taxonomic_rank) %in% TaxaList) %>%
        select(taxonomic_ranks[2:taxonomic_num],taxonomic_key_ranks,Common_Name,Image_URL)
      taxonomy_db <- as.data.frame(TaxonomyInput)
      
      # Figure out which taxonomy version is more complete.
      taxonomy_db$rankCount <- rowSums(!is.na(taxonomy_db[,colnames(taxonomy_db) %in% taxonomic_ranks]))
      taxonomy_db$rankKeyCount <- rowSums(!is.na(taxonomy_db[,colnames(taxonomy_db) %in% taxonomic_key_ranks]))
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
      taxonomy_db <- taxonomy_db[!duplicated(taxonomy_db),]
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
    
    # Determine the number of samples passing the filter.
    # Merge in taxonomic data.
    if(nrow(tronko_db)>0){
      num_filteredSamples <- length(unique(tronko_db$SampleID))
      
      # Merge in taxonomy data.
      tronko_db <- dplyr::left_join(tronko_db, taxonomy_db[,c(taxonomic_rank,"Common_Name","Image_URL")],na_matches="never")
      tronko_db$Image_URL <- ifelse(is.na(tronko_db$Image_URL), 'https://images.phylopic.org/images/5d646d5a-b2dd-49cd-b450-4132827ef25e/raster/487x1024.png', tronko_db$Image_URL)
      if (taxonomic_rank != "kingdom") {
        colnames(tronko_db)[which(names(tronko_db) == taxonomic_rank)] <- "Latin_Name"
      }
      if (taxonomic_rank == "kingdom") {
        tronko_db$Latin_Name <- tronko_db$kingdom
      }
      tronko_db <- tronko_db %>% mutate_all(~ifelse(is.na(.), " ", .))
      
      # De-duplicated dataframe for tronko table.
      tronko_db <- tronko_db[!duplicated(tronko_db),]
    }
    if(nrow(tronko_db)==0){
      stop("Error: Filters are too stringent. Cannot proceed.")
    }
    
    # Merge Tronko output with sample metadata
    project_db <- dplyr::left_join(tronko_db,metadata[,c(categorical_variables,continuous_variables,field_vars,"sample_id")],by=c("SampleID"="fastqid"))
    print(paste("project_db",nrow(project_db),ncol(project_db)))
    
    # Generate the number of samples and number of samples post-filtering as a return object,
    # to append to output json objects.
    sample_db <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(sample_db) <- c("totalSamples", "filteredSamples")
    sample_db$totalSamples <- total_samples
    sample_db$filteredSamples <- num_filteredSamples
    
    # Aggregate merged data by the appropriate time interval to find taxa presence by time.
    project_db_by_time <- project_db %>% dplyr::distinct(site,SampleID,sample_date,Latin_Name,.keep_all=T) %>% 
      dplyr::group_by(sample_date,Latin_Name) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    
    # Export taxa presence by time.
    project_db_by_time <- as.data.frame(project_db_by_time)
    
    # Store date range information
    unique_dates <- unique(as.character(project_db_by_time$sample_date))
    
    # Merge back in taxonomic information.
    tmp <- project_db[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    project_db_by_time <- dplyr::left_join(project_db_by_time,tmp)
    
    # Convert to a wide data frame.
    project_db_by_time <- tidyr::pivot_wider(project_db_by_time[c("Latin_Name","Common_Name","Image_URL","sample_date","freq")], names_from = sample_date, values_from = freq, values_fill = 0)
    project_db_by_time <- as.data.frame(project_db_by_time)
    
    # Convert data set to JSON object.
    json_list <- list()
    
    # Iterate through rows of the data frame
    df <- project_db_by_time
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
    
      # Create a list for the date ranges and Presence values
      sample_data <- list()
      for(unique_date in unique_dates){
        freq <- df[i, unique_date]
        sample_data[[unique_date]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name = Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        Sample_Date = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }

    # Initialize an empty list to accumulate file keys
    labeled_file_keys <- list()

    output_filename <- paste("presence_by_time",filename,sep="_")
    datasets <- list(datasets = list(results = json_list, metadata = sample_db))
    write(toJSON(datasets), output_filename)
    
    # Add the file key to the named list with a label
    file_key = paste("projects/",project_id,"/plots/",output_filename,sep="")
    labeled_file_keys[["presence_by_time"]] <- file_key
    system(paste("aws s3 cp ",output_filename," s3://",bucket,"/",file_key," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",output_filename,sep=""))

    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"labeledFileKeys\" = '%s' WHERE \"id\" = '%s'",
      toJSON(labeled_file_keys, auto_unbox = TRUE), report_id
    )
    print(paste("Updating report", report_id, "with file labeledFileKeys", file_key, sep = " "))
    dbExecute(con, sql_query)
    
    # Aggregate merged data by site to find taxa presence by site.
    project_db_by_site <- project_db %>% dplyr::distinct(site,SampleID,Latin_Name,.keep_all=T) %>% 
      dplyr::group_by(site,Latin_Name) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    
    # Export taxa presence by site.
    project_db_by_site <- as.data.frame(project_db_by_site)
    
    # Get unique sites.
    unique_sites <- unique(project_db_by_site$site)
    
    # Merge back in taxonomic information.
    tmp <- project_db[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    project_db_by_site <- dplyr::left_join(project_db_by_site,tmp)
    
    # Convert to a wide data frame.
    project_db_by_site <- tidyr::pivot_wider(project_db_by_site, names_from = site, values_from = freq, values_fill = 0)
    project_db_by_site <- as.data.frame(project_db_by_site)
    
    # Convert data set to JSON object.
    json_list <- list()
    
    # Iterate through rows of the data frame
    df <- project_db_by_site
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      
      # Create a list for the Site and Presence values
      sample_data <- list()
      for(unique_site in unique_sites){
        freq <- df[i, unique_site]
        sample_data[[unique_site]] <- freq
      }
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name= Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        Site = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    output_filename <- paste("presence_by_site",filename,sep="_")
    datasets <- list(datasets = list(results = json_list, metadata = sample_db))
    write(toJSON(datasets), output_filename)

    # Add the file key to the named list with a label
    file_key = paste("projects/",project_id,"/plots/",output_filename,sep="")
    labeled_file_keys[["presence_by_site"]] <- file_key
    system(paste("aws s3 cp ",output_filename," s3://",bucket,"/",file_key," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",output_filename,sep=""))

    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"labeledFileKeys\" = '%s' WHERE \"id\" = '%s'",
      toJSON(labeled_file_keys, auto_unbox = TRUE), report_id
    )
    print(paste("Updating report", report_id, "with file labeledFileKeys", file_key, sep = " "))
    dbExecute(con, sql_query)
    
    # Find taxa presence/absence by sample.
    project_db_by_sample <- project_db[,c("Latin_Name","sample_id")]
    project_db_by_sample <- project_db_by_sample[!duplicated(project_db_by_sample),]
    project_db_by_sample <- project_db_by_sample[!is.na(project_db_by_sample[,"Latin_Name"]),]
    colnames(project_db_by_sample) <- c("taxa","SampleID")
    project_db_by_sample$Presence <- 1
    
    # Get sample names.
    unique_samples <- unique(project_db$sample_id)
    
    # Create a reference data frame with all possible combinations of taxa and samples.
    all_combinations <- expand.grid(
      taxa = na.omit(unique(project_db_by_sample$taxa)),
      SampleID = unique(project_db_by_sample$SampleID)
    )
    
    # Merge the original data frame with the reference data frame
    project_db_by_sample <- merge(all_combinations, project_db_by_sample, by = c("taxa", "SampleID"), all.x = TRUE)
    
    # Convert NA values to 0 in Presence.
    project_db_by_sample["Presence"][is.na(project_db_by_sample["Presence"])] <- 0
    
    # Merge back in taxonomic information.
    tmp <- project_db[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    project_db_by_sample <- dplyr::left_join(project_db_by_sample,tmp,by=c("taxa"="Latin_Name"))
    
    # Convert to a presence/absence data frame.
    project_db_by_sample <- tidyr::pivot_wider(project_db_by_sample, names_from = SampleID, values_from = Presence, values_fill = 0)
    project_db_by_sample <- as.data.frame(project_db_by_sample)
    
    # Re-insert taxonomic rank.
    colnames(project_db_by_sample)[which(names(project_db_by_sample) == "taxa")] <- "Latin_Name"
    
    # Convert data set to JSON object.
    json_list <- list()
    
    # Iterate through rows of the data frame
    df <- project_db_by_sample
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      # Create a list for the SampleID and Presence values
      sample_data <- list()
      for(unique_sample in unique_samples){
        presence <- df[i, unique_sample]
        sample_data[[unique_sample]] <- presence
      }
      
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name = Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        SampleID = sample_data
      )
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }

    output_filename <- paste("presence_by_sample",filename,sep="_")
    datasets <- list(datasets = list(results = json_list, metadata = sample_db))

    # Add the file key to the named list with a label
    file_key = paste("projects/",project_id,"/plots/",output_filename,sep="")
    labeled_file_keys[["presence_by_sample"]] <- file_key
    write(toJSON(datasets), output_filename)
    system(paste("aws s3 cp ",output_filename," s3://",bucket,"/",file_key," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",output_filename,sep=""))

    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"labeledFileKeys\" = '%s' WHERE \"id\" = '%s'",
      toJSON(labeled_file_keys, auto_unbox = TRUE), report_id
    )
    print(paste("Updating report", report_id, "with file labeledFileKeys", file_key, sep = " "))
    dbExecute(con, sql_query)
    
    # Aggregate merged data by the appropriate time interval to find taxa presence by time.
    project_db_by_site_time <- project_db %>% dplyr::distinct(site,SampleID,sample_date,Latin_Name,.keep_all=T) %>% 
      dplyr::group_by(site,sample_date,Latin_Name) %>% 
      dplyr::summarise(n=n_distinct(SampleID)) %>% dplyr::mutate(freq=n/max(n)) %>% select(-n)
    
    # Export taxa presence by time.
    project_db_by_site_time <- as.data.frame(project_db_by_site_time)
    
    # Get unique sites.
    unique_sites <- unique(project_db_by_site_time$site)
    
    # Get unique date ranges.
    unique_dates <- unique(as.character(project_db_by_site_time$sample_date))
    
    # Merge back in taxonomic information.
    tmp <- project_db[,c("Latin_Name","Common_Name","Image_URL")]
    tmp <- tmp[!duplicated(tmp),]
    project_db_by_site_time <- dplyr::left_join(project_db_by_site_time,tmp)
    
    # Convert to a wide data frame.
    project_db_by_site_time <- tidyr::pivot_wider(project_db_by_site_time[,c("site","Latin_Name","freq","sample_date","Common_Name","Image_URL")], names_from = site, values_from = freq, values_fill = 0)
    project_db_by_site_time <- as.data.frame(project_db_by_site_time)
    
    # Convert data set to JSON object.
    json_list <- list()
    
    # Iterate through rows of the data frame
    df <- project_db_by_site_time
    for (i in 1:nrow(df)){
      Latin_Name <- df[i,"Latin_Name"]
      Common_Name <- df[i,"Common_Name"]
      Image_URL <- df[i,"Image_URL"]
      Sample_Date <- df[i,"sample_date"]
      
      # Create a list for the site/date/frequency values
      sample_data <- list()
      for(unique_site in unique_sites){
        freq <- df[i, unique_site]
        sample_data[[unique_site]] <- freq
      }
      
      # Create the JSON object for the current row
      json_obj <- list(
        latin_name = Latin_Name,
        common_name = Common_Name,
        image_url = Image_URL,
        sample_date = Sample_Date,
        Site = sample_data
      )
      
      # Append the JSON object to the list
      json_list[[i]] <- json_obj
    }
    
    output_filename <- paste("presence_by_site_and_time",filename,sep="_")
    datasets <- list(datasets = list(results = json_list, metadata = sample_db))
    write(toJSON(datasets),output_filename)
    
    # Add the file key to the named list with a label
    file_key = paste("projects/",project_id,"/plots/",output_filename,sep="")
    labeled_file_keys[["presence_by_site_and_time"]] <- file_key
    system(paste("aws s3 cp ",output_filename," s3://",bucket,"/",file_key," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",output_filename,sep=""))
    
    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"labeledFileKeys\" = '%s' WHERE \"id\" = '%s'",
      toJSON(labeled_file_keys, auto_unbox = TRUE), report_id
    )
    print(paste("Updating report", report_id, "with file labeledFileKeys", file_key, sep = " "))
    dbExecute(con, sql_query)

    # Export filtered taxonomy table.
    tronko_table <- project_db[,c("SampleID","Latin_Name")]
    tronko_table <- as.data.frame(table(tronko_table[,c("SampleID","Latin_Name")]))
    tronko_table <- as.data.frame(pivot_wider(tronko_table, names_from = SampleID, values_from = Freq))
    
    output_filename <- paste("filtered_taxonomy",filename,sep="_")
    file_key = paste("projects/",project_id,"/tables/",output_filename,sep="")
    labeled_file_keys[["filtered_taxonomy"]] <- file_key
    write.table(tronko_table,output_filename,quote=FALSE,sep=",",row.names = FALSE)
    system(paste("aws s3 cp ",output_filename," s3://",bucket,"/",file_key," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",output_filename,sep=""))

    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"labeledFileKeys\" = '%s' WHERE \"id\" = '%s'",
      toJSON(labeled_file_keys, auto_unbox = TRUE), report_id
    )
    print(paste("Updating report", report_id, "with file labeledFileKeys", file_key, sep = " "))
    dbExecute(con, sql_query)

    updateReport(report_id, "COMPLETED", con)
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  },
  error = function(e) {
    print("Error in eDNAExplorer_Alpha_Metabarcoding.R:")
    print(e$message)
    process_error(e, report_id, project_id, con)
  }
)
