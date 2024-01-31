#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(sf)
require(sp)
require(ggVennDiagram)
require(RColorBrewer)
require(ggplot2)
require(gbifdb)
require(jsonlite)
require(data.table)
require(digest)

source("helpers.R")
source("init_report.R")

# Rscript --vanilla ednaexplorer_staging_Venn_Metabarcoding.R "report_id"
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
    subset_file <- paste("subset_venn_",UUIDgenerate(),".csv",sep="")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",tronko_file_tmp, subset_file)
    system(awk_command, intern = TRUE)
    tronko_input <- fread(file=subset_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
    tronko_input$Mismatch <- as.numeric(as.character(tronko_input$Mismatch))
    
    # Remove samples with missing coordinates, and which are outside of the date filters.
    tronko_input <- tronko_input <- tronko_input[tronko_input$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    
    # Store the unfiltered reads.
    tronko_unfiltered <- tronko_input

    # Calculate relative abundance of taxa with a given rank in the unfiltered reads.
    tronko_unfiltered <- tronko_unfiltered %>%
      dplyr::group_by(SampleID, !!sym(taxonomic_rank)) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>%
      dplyr::ungroup() %>% select(-n)
    
    # Filter on the number of reads per sample, then a mismatch threshold.
    print(tronko_input)
    str(tronko_input)
    summary(tronko_input)

    tronko_input <- tronko_input %>% group_by(SampleID) %>%
      filter(n() > count_threshold) %>% filter(Mismatch <= num_mismatch & !is.na(Mismatch)) %>% select(-Mismatch)

    # Merge relative abudance results into data filtered by reads per sample and mismatches
    # and then filtered on relative abundances.
    tronko_input <- tronko_input %>%
      left_join(tronko_unfiltered,na_matches="never") %>%
      filter(freq >= filter_threshold) %>% select(-freq)
    
    # Remove taxa which are unknown at a given rank.
    tronko_db <- as.data.frame(tronko_input)
    tronko_db <- tronko_db[,c("SampleID",taxonomic_ranks[1:taxonomic_num])]
    tronko_db <- tronko_db[!is.na(tronko_db[, taxonomic_rank]), ]
    
    # Filter results by species list.
    if (selected_species_list != "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)) & tronko_db$species %in% SpeciesList_df$name, ]
    }
    if (selected_species_list == "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    }
    system(paste("rm ",tronko_file_tmp,sep=""))
    system(paste("rm ",subset_file,sep=""))
    
    # Get eDNA results summarized.
    if(nrow(tronko_db) >= 1){
      num_filtered_samples <- length(unique(tronko_db$SampleID))
      
      # Get unique taxa list from Tronko-assign
      tronko_taxa <- na.omit(unique(tronko_db[,taxonomic_rank]))
      tronko_taxa <- as.data.frame(tronko_taxa)
      colnames(tronko_taxa) <- c("eDNA")
    } else {
      tronko_taxa <- data.frame(matrix(ncol=1,nrow=1))
      colnames(tronko_taxa) <- c("eDNA")
      tronko_taxa$eDNA <- NA
      num_filtered_samples <- 0
    }
    
    # Read in GBIF occurrences.
    gbif <- gbif_local()
    
    # Get unique states and nations in project.
    country_list <- na.omit(unique(metadata$nation))
    state_province_list <- na.omit(unique(metadata$state))
    
    if(geographic_scale=="Local"){
      # Get local bounds for sample locations, add 0.5 degree buffer.
      local_east <- max(na.omit(metadata$longitude))+0.5
      local_west <- min(na.omit(metadata$longitude))-0.5
      local_south <- min(na.omit(metadata$latitude))-0.5
      local_north <- max(na.omit(metadata$latitude))+0.5
      
      # Clip GBIF occurrence locations by local boundaries.
      Taxa_Local <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                    coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                    occurrencestatus=="PRESENT",
                                    decimallongitude >= local_west & decimallongitude <= local_east & decimallatitude >= local_south & decimallatitude <= local_north) %>% 
        select(!!sym(taxonomic_rank))
      taxa_gbif <- as.data.frame(Taxa_Local)
      taxa_gbif <- na.omit(unique(taxa_gbif[,taxonomic_rank]))
      taxa_gbif <- as.data.frame(taxa_gbif)
      colnames(taxa_gbif) <- c("Taxa_Local")
      
      # Define results shared between eDNA and GBIF.
      Both <- as.data.frame(intersect(tronko_taxa$eDNA,taxa_gbif$Taxa_Local))
      colnames(Both) <- c("Both")
      
      # Define eDNA results as those unique from GBIF
      eDNA_only <- as.data.frame(tronko_taxa$eDNA[!(tronko_taxa$eDNA %in% taxa_gbif$Taxa_Local)])
      colnames(eDNA_only) <- c("eDNA")
      
      # Define GBIF results those unique from eDNA
      gbif_only <- as.data.frame(taxa_gbif$Taxa_Local[!(taxa_gbif$Taxa_Local %in% tronko_taxa$eDNA)])
      colnames(gbif_only) <- c("Taxa_Local")
      
      rm(taxa_gbif)
      taxa_gbif <- gbif_only
      rm(tronko_taxa)
      tronko_taxa <- eDNA_only
    }
    if(geographic_scale=="State"){
      # Clip GBIF occurrence locations by state/province boundaries.
      Taxa_State <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                    coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                    occurrencestatus=="PRESENT", stateprovince %in% state_province_list) %>% select(!!sym(taxonomic_rank))
      taxa_gbif <- as.data.frame(Taxa_State)
      taxa_gbif <- na.omit(unique(taxa_gbif[,taxonomic_rank]))
      taxa_gbif <- as.data.frame(taxa_gbif)
      colnames(taxa_gbif) <- c("Taxa_State")
      
      # Define results shared between eDNA and GBIF.
      Both <- as.data.frame(intersect(tronko_taxa$eDNA,taxa_gbif$Taxa_State))
      colnames(Both) <- c("Both")
      
      # Define eDNA results as those unique from GBIF
      eDNA_only <- as.data.frame(tronko_taxa$eDNA[!(tronko_taxa$eDNA %in% taxa_gbif$Taxa_State)])
      colnames(eDNA_only) <- c("eDNA")
      
      # Define GBIF results those unique from eDNA
      gbif_only <- as.data.frame(taxa_gbif$Taxa_State[!(taxa_gbif$Taxa_State %in% tronko_taxa$eDNA)])
      colnames(gbif_only) <- c("Taxa_State")
      
      rm(taxa_gbif)
      taxa_gbif <- gbif_only
      rm(tronko_taxa)
      tronko_taxa <- eDNA_only
    }
    if(geographic_scale=="Nation"){
      # Clip GBIF occurrence locations by national boundaries.
      Taxa_Nation <- gbif %>% filter(basisofrecord %in% c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION"),
                                     coordinateuncertaintyinmeters <= 100 & !is.na(coordinateuncertaintyinmeters),
                                     occurrencestatus=="PRESENT", countrycode %in% country_list) %>% select(!!sym(taxonomic_rank))
      taxa_gbif <- as.data.frame(Taxa_Nation)
      taxa_gbif <- na.omit(unique(taxa_gbif[,taxonomic_rank]))
      taxa_gbif <- as.data.frame(taxa_gbif)
      colnames(taxa_gbif) <- c("Taxa_Nation")
      
      # Define results shared between eDNA and GBIF.
      Both <- as.data.frame(intersect(tronko_taxa$eDNA,taxa_gbif$Taxa_Nation))
      colnames(Both) <- c("Both")
      
      # Define eDNA results as those unique from GBIF
      eDNA_only <- as.data.frame(tronko_taxa$eDNA[!(tronko_taxa$eDNA %in% taxa_gbif$Taxa_Nation)])
      colnames(eDNA_only) <- c("eDNA")
      
      # Define GBIF results those unique from eDNA
      gbif_only <- as.data.frame(taxa_gbif$Taxa_Nation[!(taxa_gbif$Taxa_Nation %in% tronko_taxa$eDNA)])
      colnames(gbif_only) <- c("Taxa_Nation")
      
      rm(taxa_gbif)
      taxa_gbif <- gbif_only
      rm(tronko_taxa)
      tronko_taxa <- eDNA_only
    }
    
    # Insert the number of samples and number of samples post-filtering as a return object.
    SampleDB <- data.frame(matrix(ncol=2,nrow=1))
    colnames(SampleDB) <- c("totalSamples","filteredSamples")
    SampleDB$totalSamples <- total_samples
    SampleDB$filteredSamples <- num_filtered_samples
    datasets <- list(datasets = list(eDNA=tronko_taxa[,1],Both=Both[,1],GBIF=taxa_gbif[,1],metadata=SampleDB))
    
    # Save plot as json object.
    write(toJSON(datasets),filename)
    file_key <- paste("projects/",sample_project_id,"/plots/",filename,sep="")
    system(paste("aws s3 cp ",filename," s3://",bucket,"/",file_key," --endpoint-url ",ENDPOINT_URL,sep=""),intern=TRUE)
    system(paste("rm ",filename,sep=""))
    
    # Update report with file key.
    sql_query <- sprintf(
      "UPDATE \"Report\" SET \"fileKey\" = '%s' WHERE \"id\" = '%s'",
      file_key, report_id
    )
    print(paste("Updating report", report_id, "with file key", file_key, sep = " "))
    dbExecute(con, sql_query)
    updateReport(report_id, "COMPLETED", con)
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  },
  error = function(e) {
    print("Error in eDNAExplorer_Venn_Metabarcoding.R:")
    print(e$message)
    process_error(e, report_id, project_id, con)
  }
)