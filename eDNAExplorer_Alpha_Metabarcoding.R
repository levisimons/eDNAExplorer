#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(vegan)
require(ggplot2)
if (!require("BiocManager", quietly = TRUE)){
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


# Rscript --vanilla eDNAExplorer_Alpha_Metabarcoding.R "report_id"
tryCatch(
  {     
    # Read in species list
    if(selected_species_list != "None"){
      species_list_df <- tbl(con,"SpeciesListItem")
      species_list_df <- species_list_df %>% filter(species_list == selected_species_list)
      species_list_df <- as.data.frame(species_list_df)
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
    
    # Create sample metadata matrix
    if(nrow(metadata) == 0 || ncol(metadata) == 0) {
      stop("Error: Sample data frame is empty. Cannot proceed.")
    }
    sample <- metadata[!is.na(metadata$fastqid),]
    rownames(sample) <- sample$fastqid
    sample$fastqid <- NULL
    sample <- sample_data(sample)
    
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
    subset_file <- paste("subset_alpha_",UUIDgenerate(),".csv",sep="")
    awk_command <- sprintf("awk -F, 'BEGIN {OFS=\",\"} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} {print $col[\"SampleID\"], $col[\"superkingdom\"], $col[\"kingdom\"], $col[\"phylum\"], $col[\"class\"], $col[\"order\"], $col[\"family\"], $col[\"genus\"], $col[\"species\"], $col[\"Mismatch\"]}' %s > %s",tronko_file_tmp, subset_file)
    system(awk_command, intern = TRUE)
    tronko_input <- fread(file=subset_file, header = TRUE, sep = ",", skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A"))
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
     dplyr::ungroup() %>% select(-n)

    # Filter on the number of reads per sample, then a mismatch threshold.
    tronko_input <- tronko_input %>% group_by(SampleID) %>%
      filter(n() > sample_count_threshold) %>% filter(Mismatch <= sample_num_mismatch & !is.na(Mismatch)) %>% select(-Mismatch)
    
    # Merge relative abudance results into data filtered by reads per sample and mismatches
    # and then filtered on relative abundances.
    tronko_input <- tronko_input %>%
      left_join(tronko_unfiltered,na_matches="never") %>%
      filter(freq >= filter_threshold) %>% select(-freq)
    tronko_db <- as.data.frame(tronko_input)
    
    # Remove taxa which are unknown at a given rank.
    tronko_db <- tronko_db[,c("SampleID",taxonomic_ranks[1:taxonomic_num])]
    tronko_db <- tronko_db[!is.na(tronko_db[, sample_taxonomic_rank]), ]
    
    # Filter results by species list.
    if (selected_species_list != "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)) & tronko_db$species %in% species_list_df$name, ]
    }
    if (selected_species_list == "None") {
      tronko_db <- tronko_db[tronko_db$SampleID %in% unique(na.omit(metadata$fastqid)), ]
    }
    system(paste("rm ",subset_file,sep=""))
    
    if(nrow(tronko_db) > 1){
      # Read in information to map categorical labels for certain variables.
      category_file <- paste("Categories_",UUIDgenerate(),".csv",sep="")
      system(paste("aws s3 cp s3://",bucket,"/analysis/Categories.csv ",category_file," --endpoint-url ",ENDPOINT_URL,sep=""))
      categories <- as.data.frame(fread(file=category_file,header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A")))
      system(paste("rm ",category_file,sep=""))

      # Read in information for legends and labels
      legends_file <- paste("LabelsAndLegends_",UUIDgenerate(),".csv",sep="")
      system(paste("aws s3 cp s3://",bucket,"/analysis/LabelsAndLegends.csv ",legends_file," --endpoint-url ",ENDPOINT_URL,sep=""))
      legends_and_labels <- as.data.frame(fread(file=legends_file,header=TRUE, sep=",",skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na = c("", "NA", "N/A")))
      system(paste("rm ",legends_file,sep=""))

      # Set up new legends and x-axis labels.
      new_legend <- legends_and_labels[legends_and_labels$Environmental_Variable==environmental_variable,"Legend"]
      
      # Create OTU matrix
      otumat <- as.data.frame(pivot_wider(as.data.frame(table(tronko_db[,c("SampleID",sample_taxonomic_rank)])), names_from = SampleID, values_from = Freq))
      rownames(otumat) <- otumat[,sample_taxonomic_rank]
      otumat <- otumat[,colnames(otumat) %in% unique(tronko_db$SampleID)]
      otumat[sapply(otumat, is.character)] <- lapply(otumat[sapply(otumat, is.character)], as.numeric)
      otu <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)
      
      # Create merged Phyloseq object.
      physeq <- phyloseq(otu,sample)
      abundance_filtered <- physeq
      
      if(environmental_variable %in% categorical_variables){
        # Store diversity versus variable data.
        tmp <- ggplot_build(plot_richness(abundance_filtered,x=environmental_variable,measures=alpha_diversity_metric))
        tmp <- tmp$data[[1]]

        # Store variable values for categorical data.
        if(environmental_variable %in% unique(categories$Environmental_Variable)){
          tmp$x <- as.character(sample_data(abundance_filtered)[[environmental_variable]])
          tmp <- dplyr::left_join(tmp,categories[categories$Environmental_Variable==environmental_variable,],by=c("x"="value"))
          tmp["description"][is.na(tmp["description"])] <- "no data available"
          tmp$x <- as.factor(tmp$description)
        } else {
          
          # Rename missing data
          x_label <- sample_data(abundance_filtered)[[environmental_variable]]
          x_label[is.na(x_label)] <- "no data available"
          tmp$x <- as.factor(x_label)
        }
        
        # Determine variable frequency.
        tmp <- dplyr::left_join(tmp,as.data.frame(table(tmp$x)),by=c("x"="Var1"))
        
        # Run a Kruskal-Wallis test between alpha diversity and selected environmental variable.
        if(length(unique(na.omit(tmp$x)))>1){
          test <- suppressWarnings(kruskal.test(tmp$y ~ tmp$x, data = tmp))
          Stats_Message <- paste("chi-squared = ",round(test$statistic,digits=3)," df = ",test$parameter," p-value = ",round(test$p.value,digits=3))
        } else {
          Stats_Message <- "Not enough variation in environmental variable to run Kruskal-Wallis test."
        }

        # General a violin plot of alpha diversity versus an environmental variable.
        if(max(tmp$Freq) > 1 & min(tmp$Freq)==1){
          p <- ggplot()+geom_violin(data = tmp %>% filter(Freq > 1),aes(x = x, y = y))+ 
            geom_point(data = tmp %>% filter(Freq==1),aes(x = x, y = y),position = position_jitter(seed = 1, width = 0.2))+
            labs(title=paste(alpha_diversity_metric," versus ",gsub("_"," ",environmental_variable),".\nSamples collected between: ",sample_first_date," and ",sample_last_date,"\nRelative abundance minimum of ",100*sample_filter_threshold,"%.\nReads per sample minimum: ",sample_count_threshold,"\n",Stats_Message,sep=""),x=new_legend, y = alpha_diversity_metric)+
            theme_bw()+guides(fill=guide_legend(title=new_legend))
        }

        if(max(tmp$Freq) > 1){
          p <- ggplot()+geom_violin(data = tmp,aes(x = x, y = y))+ 
            geom_point(data = tmp,aes(x = x, y = y),position = position_jitter(seed = 1, width = 0.2))+
            labs(title=paste(alpha_diversity_metric," versus ",gsub("_"," ",environmental_variable),".\nSamples collected between: ",sample_first_date," and ",sample_last_date,"\nRelative abundance minimum of ",100*sample_filter_threshold,"%.\nReads per sample minimum: ",sample_count_threshold,"\n",Stats_Message,sep=""),x=new_legend, y = alpha_diversity_metric)+
            theme_bw()+guides(fill=guide_legend(title=new_legend))
        }
        
        if(max(tmp$Freq) <= 1){
          p <- ggplot()+ 
            geom_point(data = tmp %>% filter(Freq==1),aes(x = x, y = y),position = position_jitter(seed = 1, width = 0.2))+
            labs(title=paste(alpha_diversity_metric," versus ",gsub("_"," ",environmental_variable),".\nSamples collected between: ",sample_first_date," and ",sample_last_date,"\nRelative abundance minimum of ",100*sample_filter_threshold,"%.\nReads per sample minimum: ",sample_count_threshold,"\n",Stats_Message,sep=""),x=new_legend, y = alpha_diversity_metric)+
            theme_bw()+guides(fill=guide_legend(title=new_legend))
        }
      }
      
      # Store diversity versus variable data.
      if(environmental_variable %in% continuous_variables){
        tmp <- ggplot_build(plot_richness(abundance_filtered,x=environmental_variable,measures=alpha_diversity_metric))
        tmp <- tmp$data[[1]]
        tmp$x <- as.numeric(tmp$x)
        
        # Calculate summary statistics using a Kendall correlation test between alpha diversity and environmental variable.
        if(length(unique(tmp$x))>1){
          test <- suppressWarnings(cor.test(tmp$x,tmp$y,alternative="two.sided",method="kendall"))
          Stats_Message <- paste("z = ",round(test$statistic,digits=3)," tau = ",round(test$estimate,digits=3)," p-value = ",round(test$p.value,digits=3))
        } else {
          Stats_Message <- "Not enough variation in environmental variable to run Kruskal-Wallis test."
        }
        
        # Generate a scatterplot of alpha diversity versus an environmental variable.
        p <- ggplot(tmp, aes(x=x, y=y))+
          labs(title=paste(alpha_diversity_metric," versus ",gsub("_"," ",environmental_variable),".\nSamples collected between: ",sample_first_date," and ",sample_last_date,"\nRelative abundance minimum of ",100*sample_filter_threshold,"%.\nReads per sample minimum: ",sample_count_threshold,"\n",Stats_Message,sep=""),x=new_legend, y = alpha_diversity_metric)+
          theme_bw()+geom_point()+geom_smooth()+guides(fill=guide_legend(title=new_legend))
      }
    }
    if(nrow(tronko_db) <= 1){
      stat_test <- "Not enough data to perform a Kruskal-Wallis test on alpha diversity."
      p <- ggplot(data.frame())+geom_point()+xlim(0, 1)+ylim(0, 1)+labs(title=stat_test)
    }
    
    # Insert the number of samples and number of samples post-filtering as a return object.
    sample_db <- data.frame(matrix(ncol=2,nrow=1))
    colnames(sample_db) <- c("totalSamples","filteredSamples")
    sample_db$totalSamples <- total_samples
    sample_db$filteredSamples <- nsamples(abundance_filtered)
    datasets <- list(datasets = list(results=plotly_json(p, FALSE),metadata=toJSON(sample_db)))
        
    updateReport(report_id, "COMPLETED", con)
    
    # Save plot as json object
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
    
    
    RPostgreSQL::dbDisconnect(con, shutdown=TRUE)
  },
  error = function(e) {
    print("Error in eDNAExplorer_Alpha_Metabarcoding.R:")
    print(e$message)
    process_error(e, report_id, project_id, con)
  }
)