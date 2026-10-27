###Metanalysis R Data

26/04/2025

##### 0 Packages #####

# detach loaded packages and data
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()
rm(list = ls())

# --- Core data science packages ---
library(tidyverse)      # Loads dplyr, ggplot2, tidyr, readr, etc.

# --- Table formatting and plotting extensions ---
library(ggpubr)         # advanced ggplot2 plots
library(patchwork)      # plot arrangement
library(viridis)        # color scales for ggplot2
library(knitr)          # for tables in reports

# --- Microbiome & phylogenetics ---
library(phyloseq)       # core microbiome object manipulation
library(microbiome)     # extra microbiome utilities
library(vegan)          # diversity, adonis
library(ape)            # phylogenetic trees

# --- Specialized matrix/data tools ---       # wide/long data reshaping
library(splitstackshape) # core microbiome matrix, cSplit

# --- Differential abundance/statistics ---
library(DESeq2)
library(MMUPHin)        # batch correction
library(metamicrobiomeR) # meta-analysis for microbiome
library(metafor)        # meta-analysis (general)
library(meta)           # forest plots etc
library(gamlss)         # GAMLSS models

library(pROC)

#### Constants and Functions #####

### Constants
# FDR cutoff for differential abundances
fdr.threshold <- 0.05
set.seed(123)
# working dir
dir.create(here::here("results/meta_analysis"), recursive = TRUE, showWarnings = FALSE)
setwd(here::here("results/meta_analysis"))

#### Load phyloseq #####
# Define the mapping between filenames and your preferred study names
study_mapping <- c(
  "ps_PRJEB13256_clean.rds" = "Study2",
  "ps_PRJNA717340_clean.rds" = "Study5", 
  "ps_PRJNA436876_clean.rds" = "Study7"
)

# Load each phyloseq object with your custom naming
ps_list <- list()
for(file_name in names(study_mapping)) {
  file_path <- here::here("data/processed", file_name)
  if(file.exists(file_path)) {
    ps <- readRDS(file_path)
    study_name <- study_mapping[file_name]
    
    # Add StudyID to the sample data if not already present
    if(!"StudyID" %in% colnames(sample_data(ps))) {
      sample_data(ps)$StudyID <- study_name
    } else {
      # Override existing StudyID with your preferred naming
      sample_data(ps)$StudyID <- study_name
    }
    
    ps_list[[study_name]] <- ps
    cat(paste0("Loaded ", file_name, " as ", study_name, "\n"))
  } else {
    warning(paste0("File not found: ", file_path))
  }
}

# Name the list elements with the study IDs
names(ps_list) <- sapply(ps_list, function(ps) {
  unique(sample_data(ps)$StudyID)[1]
})

# Optional: Examine the studies to ensure they're loaded correctly
lapply(ps_list, function(ps) {
  cat("Samples:", nsamples(ps), "\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Ranks available:", rank_names(ps), "\n\n")
})

#### Check metadata structure ####

# Function to examine metadata columns and unique values
check_metadata <- function(ps, study_name) {
  cat(paste0("\n=== Metadata for ", study_name, " ===\n"))
  
  meta <- data.frame(sample_data(ps))
  cat("Number of samples:", nrow(meta), "\n")
  cat("Metadata columns:", paste(colnames(meta), collapse=", "), "\n\n")
  
  # Print unique values for each column (up to 10 values)
  for(col in colnames(meta)) {
    unique_vals <- unique(meta[[col]])
    n_unique <- length(unique_vals)
    
    # Truncate if too many unique values
    if(n_unique > 10) {
      unique_vals <- unique_vals[1:10]
      truncated <- TRUE
    } else {
      truncated <- FALSE
    }
    
    cat(paste0(col, " (", n_unique, " unique values): ", 
               paste(unique_vals, collapse=", "),
               ifelse(truncated, "...", ""), "\n"))
  }
  
  # Identify potential treatment columns
  potential_treatment_cols <- c()
  for(col in colnames(meta)) {
    # Look for columns that might indicate treatment
    if(grepl("treat|rapa|dose|drug|condition|group", col, ignore.case=TRUE)) {
      cat("\nPotential treatment column:", col, "\n")
      cat("Values:", paste(unique(meta[[col]]), collapse=", "), "\n")
      potential_treatment_cols <- c(potential_treatment_cols, col)
    }
  }
  
  return(potential_treatment_cols)
}

# Check metadata for each study
treatment_columns <- list()
for(study_name in names(ps_list)) {
  treatment_columns[[study_name]] <- check_metadata(ps_list[[study_name]], study_name)
}

#### 1. Create filtered list first and then modify it ####

filter_study7 <- function(ps) {
  meta <- data.frame(sample_data(ps))
  
  # Filter using the correct column name "Group" (capital G)
  keep_samples <- meta$Group != "Probiotics"
  
  # Prune to keep only desired samples
  ps_filtered <- prune_samples(keep_samples, ps)
  
  cat("Study7 filtering:\n")
  cat("Original samples:", nsamples(ps), "\n")
  cat("Filtered samples:", nsamples(ps_filtered), "\n")
  
  return(ps_filtered)
}
# 1. Create filtered_ps_list first
filtered_ps_list <- ps_list
filtered_ps_list$Study7 <- filter_study7(ps_list$Study7)

# Alternatively, if the above method still fails:
# Get the metadata from Study 2
s2_meta <- as(sample_data(filtered_ps_list$Study2), "data.frame")
# Update the value in the data frame
s2_meta[rownames(s2_meta) == "ERR1428190", "application_method"] <- "Control_food"
# Convert back to sample_data and update
sample_data(filtered_ps_list$Study2) <- sample_data(s2_meta)
# Verify
cat("\nVerification after alternative method:\n")
print(sample_data(filtered_ps_list$Study2)$application_method[rownames(sample_data(filtered_ps_list$Study2)) == "ERR1428190"])

#### 2. Examine sample distribution after filtering ##################

# Simplified plot function with minimal conditionals/loops
plot_simplified <- function() {
  # Study 2 - Get data from application_method column
  study2_meta <- data.frame(sample_data(filtered_ps_list$Study2))
  study2_counts <- as.data.frame(table(study2_meta$application_method))
  names(study2_counts) <- c("Treatment", "Count")
  study2_counts$Study <- "Study2"
  
  # Study 5 - Try using most likely column
  study5_meta <- data.frame(sample_data(filtered_ps_list$Study5))
  study5_counts <- as.data.frame(table(study5_meta$Group))  # Assuming Group is the column
  names(study5_counts) <- c("Treatment", "Count")
  study5_counts$Study <- "Study5"
  
  # Study 7 - Using Group column
  study7_meta <- data.frame(sample_data(filtered_ps_list$Study7))
  study7_counts <- as.data.frame(table(study7_meta$Group))
  names(study7_counts) <- c("Treatment", "Count")
  study7_counts$Study <- "Study7"
  
  # Combine all data
  all_counts <- rbind(study2_counts, study5_counts, study7_counts)
  
  # Create and display plot
  p <- ggplot(all_counts, aes(x = Treatment, y = Count, fill = Study)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Sample Distribution in Filtered Data",
         x = "Treatment Group", y = "Number of Samples") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  return(p)
}

# Use the function
p <- plot_simplified()

plot.sample_distribution_before <- p

# Direct sample count verification
verify_study2_counts <- function() {
  # Get Study 2 metadata
  study2_meta <- data.frame(sample_data(filtered_ps_list$Study2))
  
  # Print column names to verify
  cat("Study 2 columns:", paste(colnames(study2_meta), collapse=", "), "\n\n")
  
  # Create a more direct count
  cat("Direct sample counts by application_method:\n")
  rapamycin_ip_count <- sum(study2_meta$application_method == "Rapamycin_IP")
  vehicle_ip_count <- sum(study2_meta$application_method == "Vehicle_IP")
  rapamycin_food_count <- sum(study2_meta$application_method == "Rapamycin_food")
  control_food_count <- sum(study2_meta$application_method == "Control_food")
  
  cat("Rapamycin_IP samples:", rapamycin_ip_count, "\n")
  cat("Vehicle_IP samples:", vehicle_ip_count, "\n")
  cat("Rapamycin_food samples:", rapamycin_food_count, "\n")
  cat("Control_food samples:", control_food_count, "\n\n")
  
  # Check for case sensitivity or spacing issues
  cat("Unique application_method values:\n")
  print(unique(study2_meta$application_method))
  
  # Create a data frame with the correct counts
  correct_counts <- data.frame(
    Treatment = c("Rapamycin_IP", "Vehicle_IP", "Rapamycin_food", "Control_food"),
    Count = c(rapamycin_ip_count, vehicle_ip_count, rapamycin_food_count, control_food_count),
    Study = rep("Study2", 4)
  )
  
  # Create a corrected plot
  tempplot2 <- ggplot(correct_counts, aes(x = Treatment, y = Count, fill = Study)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Corrected Study 2 Sample Counts",
         x = "Treatment Group", y = "Number of Samples") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(tempplot2)
  plot.sample_distribution_study2 <- tempplot2
  
  # Also check the other studies while we're at it
  cat("\nStudy 5 Group values:\n")
  print(table(sample_data(filtered_ps_list$Study5)$Group))
  
  cat("\nStudy 7 Group values:\n")
  print(table(sample_data(filtered_ps_list$Study7)$Group))
}

# Run the verification
verify_study2_counts()

#### 3. Standardize treatment names across studies ####
standardize_treatments_with_routes <- function(input_ps_list) {
  # Fixed: Use input_ps_list parameter name consistently
  # Create a copy of the list to modify
  std_ps_list <- input_ps_list
  
  # Process each study individually
  for(study_name in names(std_ps_list)) {
    ps <- std_ps_list[[study_name]]
    meta <- data.frame(sample_data(ps))
    
    # Create new standardized columns
    meta$StandardTreatment <- NA     # Basic Control/Rapamycin
    meta$AdminRoute <- NA           # Administration route (IP, Food, etc.)
    meta$DetailedTreatment <- NA    # Combined treatment+route
    
    # Study-specific standardization
    if(study_name == "Study2") {
      # For Study2, use the application_method column
      if("application_method" %in% colnames(meta)) {
        # Extract the basic treatment (Control/Rapamycin)
        meta$StandardTreatment <- ifelse(
          grepl("Control|Vehicle", meta$application_method, ignore.case=TRUE),
          "Control", "Rapamycin")
        
        # Extract the administration route (IP/Food)
        meta$AdminRoute <- ifelse(
          grepl("_IP$|Vehicle_IP", meta$application_method),
          "IP", "Food")
        
        # Create detailed treatment
        meta$DetailedTreatment <- paste(meta$StandardTreatment, meta$AdminRoute, sep="_")
      } else {
        warning("application_method column not found in Study2")
      }
    }
    
    if(study_name == "Study5") {
      # Check for relevant columns
      treatment_col <- NULL
      if("Group" %in% colnames(meta)) {
        treatment_col <- "Group"
      } else if("treatment" %in% colnames(meta)) {
        treatment_col <- "treatment"
      } else if("group" %in% colnames(meta)) {
        treatment_col <- "group"
      }
      
      if(!is.null(treatment_col)) {
        # Map to standard treatment
        meta$StandardTreatment <- ifelse(
          grepl("control|vehicle|WT$|wt$|wild", tolower(meta[[treatment_col]]), ignore.case=TRUE), 
          "Control", "Rapamycin")
        
        # Set all Study5 administration route to IP (assuming from previous visualizations)
        meta$AdminRoute <- "IP"
        
        # Create detailed treatment
        meta$DetailedTreatment <- paste(meta$StandardTreatment, meta$AdminRoute, sep = "_")
      }
    }
    
    if(study_name == "Study7") {
      # For Study7, we already know the column is "Group"
      if("Group" %in% colnames(meta)) {
        meta$StandardTreatment <- ifelse(meta$Group == "Control", "Control", "Rapamycin")
        
        # Set all Study7 administration route to IP (assuming from previous discussion)
        meta$AdminRoute <- "IP"
        
        # Create detailed treatment
        meta$DetailedTreatment <- paste(meta$StandardTreatment, meta$AdminRoute, sep="_")
      }
    }
    
    # Update the sample data
    sample_data(ps) <- sample_data(meta)
    std_ps_list[[study_name]] <- ps
  }
  
  # Report standardization summary
  cat("=== Standardized Treatment Summary ===\n")
  for(study_name in names(std_ps_list)) {
    ps <- std_ps_list[[study_name]]
    meta <- data.frame(sample_data(ps))
    cat(paste0(study_name, ":\n"))
    
    if("StandardTreatment" %in% colnames(meta)) {
      cat("Basic Treatment Groups:\n")
      print(table(meta$StandardTreatment, useNA = "ifany"))
    }
    
    if("AdminRoute" %in% colnames(meta)) {
      cat("\nAdministration Routes:\n")
      print(table(meta$AdminRoute, useNA = "ifany"))
    }
    
    if("DetailedTreatment" %in% colnames(meta)) {
      cat("\nDetailed Treatment Groups:\n")
      print(table(meta$DetailedTreatment, useNA = "ifany"))
    }
    
    cat("\n")
  }
  
  return(std_ps_list)
}

# Apply standardization to filtered_ps_list
standardized_ps_list <- standardize_treatments_with_routes(filtered_ps_list)

head(sample_data(standardized_ps_list$Study2))
head(sample_data(standardized_ps_list$Study5))
head(sample_data(standardized_ps_list$Study7))

#### Taxa Harmonisation ####

# Helper to strip SILVA‐style prefixes and trim whitespace
clean_genus <- function(x) {
  x <- gsub("^[kpcofgs]__", "", x, ignore.case=TRUE)
  str_trim(x)
}

# Apply to each study’s phyloseq
ps_genus <- lapply(standardized_ps_list, function(ps) {
  # 1) Collapse to Genus level
  ps_g <- tax_glom(ps, taxrank="Genus", NArm=FALSE)
  
  # 2) Extract and clean the Genus column
  tt <- as.data.frame(tax_table(ps_g))
  if ("Genus" %in% colnames(tt)) {
    tt$Genus <- clean_genus(tt$Genus)
    tax_table(ps_g) <- tax_table(as.matrix(tt))
    
    # 3) **Overwrite the taxa_names** so that downstream code uses the cleaned genus
    taxa_names(ps_g) <- make.unique(as.character(tt$Genus))
  }
  ps_g
})

#### Add Dosage information ####

# 2) Define a dosing function
add_dosing <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")
  # initialize
  meta$DailyDose_mgkg <- 0
  meta$TotalDose_mgkg <- 0
  meta$DurationDays   <- 0
  
  # assign per-study and route
  sid <- unique(meta$StudyID)
  if (sid == "Study2") {
    meta$DurationDays[meta$AdminRoute=="IP" ] <- 90
    meta$DailyDose_mgkg[meta$AdminRoute=="IP" ] <- 8
    meta$TotalDose_mgkg[meta$AdminRoute=="IP" ] <- 8 * 90
    meta$DailyDose_mgkg[meta$AdminRoute=="Food"] <- 11
    meta$TotalDose_mgkg[meta$AdminRoute=="Food"] <- 11 * 90
  }
  if (sid == "Study5") {
    meta$DurationDays      <- 14
    meta$DailyDose_mgkg    <- 2.5
    meta$TotalDose_mgkg    <- 2.5 * 14
  }
  if (sid == "Study7") {
    meta$DurationDays      <- 84
    meta$DailyDose_mgkg    <- 0.5
    meta$TotalDose_mgkg    <- 0.5 * 84
  }
  
  # zero-out controls
  ctrl <- meta$StandardTreatment=="Control"
  meta$DailyDose_mgkg[ctrl] <- 0
  meta$TotalDose_mgkg[ctrl] <- 0
  
  # log‐transform
  meta$LogTotalDose <- log2(meta$TotalDose_mgkg + 1)
  
  # write back
  sample_data(ps) <- sample_data(meta)
  ps
}

# 3) Apply to each genus‐level phyloseq
ps_genus_dose <- lapply(ps_genus, add_dosing)


head(sample_data(ps_genus_dose$Study2))
head(sample_data(ps_genus_dose$Study5))
head(sample_data(ps_genus_dose$Study7))

# Now ps_genus_dose$Study2, ps_genus_dose$Study5, ps_genus_dose$Study7
# each have cleaned genus taxonomy + dosage columns in sample_data().

#### Normalize and run Taxa.mat ####
# assume ps_genus_dose is your named list of phyloseqs
ps_genus_dose <- lapply(ps_genus_dose, function(ps) {
  # pull out the sample‐data as a data.frame
  md <- as(sample_data(ps), "data.frame")
  
  # coerce “StandardTreatment” into a two‐level factor
  md$StandardTreatment <- factor(md$StandardTreatment,
                                 levels = c("Control","Rapamycin"))
  
  # write it back
  sample_data(ps) <- sample_data(md)
  ps
})

# sanity check: show levels for each study
lapply(ps_genus_dose, function(ps) {
  levels(sample_data(ps)$StandardTreatment)
})

#### TAXA SANITIZATION  ##################

sanitize_taxa_names <- function(ps_list) {
  cat("Sanitizing taxa names in phyloseq objects...\n")
  
  # Define a function to clean a single taxon name
  clean_name <- function(name) {
    # FIRST: Collapse unknown/unclassified taxa before other cleaning
    name <- name %>%
      gsub("Tuzzerella", "Tyzzerella", .) %>%
      gsub("Erysipelatoclostridiaceae", "Erysipelatoclostridium", .)
    
    # THEN: Replace any non-alphanumeric character with underscore
    clean <- gsub("[^a-zA-Z0-9_]", "_", name)
    
    # Ensure it doesn't start with a number
    if(grepl("^[0-9]", clean)) {
      clean <- paste0("X", clean)
    }
    
    # Collapse multiple consecutive underscores and remove leading/trailing ones
    clean <- gsub("_+", "_", clean)
    clean <- gsub("^_|_$", "", clean)
    
    return(clean)
  }
  
  # Process each phyloseq object
  sanitized_list <- list()
  mapping_list <- list()  # To store all mappings
  
  for(study_name in names(ps_list)) {
    ps <- ps_list[[study_name]]
    
    
    # Remove them
    ps <- subset_taxa(ps, !Genus %in% c("Chloroplast", "Mitochondria"))
    ps <- subset_taxa(ps, !grepl("Chloroplast", Genus, ignore.case = TRUE))
    ps <- subset_taxa(ps, !grepl("Mitochondria", Genus, ignore.case = TRUE))
    
    cat("  Taxa after removal:", ntaxa(ps), "\n")
    
    
    # Get current taxa names
    old_names <- taxa_names(ps)
    
    # Create new sanitized names
    new_names <- sapply(old_names, clean_name)
    
    # Check for duplicates in new names (which would cause issues)
    if(length(unique(new_names)) < length(new_names)) {
      # Make them unique by adding a numeric suffix
      new_names <- make.unique(new_names, sep = "_")
    }
    
    # Store the mapping for reference
    name_mapping <- data.frame(
      original = old_names,
      sanitized = new_names,
      stringsAsFactors = FALSE
    )
    
    # Store mapping in list
    mapping_list[[study_name]] <- name_mapping
    
    # Replace taxa names in the phyloseq object
    taxa_names(ps) <- new_names
    
    # Store in the sanitized list
    sanitized_list[[study_name]] <- ps
    
    cat("Sanitized", length(old_names), "taxa names in", study_name, "\n")
  }
  
  # Return both the sanitized objects and the mappings
  return(list(
    sanitized_objects = sanitized_list,
    name_mappings = mapping_list
  ))
}

# Run sanitization
cat("Sanitizing taxa names to fix the $ operator error...\n")
sanitization_result <- sanitize_taxa_names(ps_genus_dose)

# Extract the sanitized objects
ps_genus_dose_sanitized_before <- sanitization_result$sanitized_objects

# Store the mappings in a variable for later reference
taxa_name_mappings <- sanitization_result$name_mappings

# Print summary of mappings
cat("\nSummary of taxa name mappings:\n")
for(study_name in names(taxa_name_mappings)) {
  mapping <- taxa_name_mappings[[study_name]]
  special_chars <- sum(mapping$original != mapping$sanitized)
  cat(study_name, ":", special_chars, "of", nrow(mapping), 
      "taxa names were changed (", 
      round(special_chars/nrow(mapping)*100, 1), "%)\n")
}

# Print the taxa names (genus names) from merged_genus
cat("=== GENERA IN MERGED_GENUS ===\n")
print(taxa_name_mappings)


#### POST-MERGE UNKNOWN TAXA COLLAPSING ####

# Function to collapse unknown taxa after merging
collapse_unknown_after_merge <- function(ps_merged) {
  cat("Collapsing unknown taxa after merge...\n")
  
  # Get OTU table and tax table
  otu_mat <- as.matrix(otu_table(ps_merged))
  tax_df <- as.data.frame(tax_table(ps_merged))
  
  # Apply comprehensive unknown taxa collapsing
  tax_df$Genus_clean <- tax_df$Genus %>%
    # Remove .1, .2, .3 suffixes first
    gsub("(\\.\\d+)$", "", .) %>%
    gsub("(_\\d+)$", "", .) %>%
    # Apply genus corrections - ADD THESE LINES
    gsub("^Tuzzerella$", "Tyzzerella", .) %>%
    gsub("^Erysipelatoclostridiaceae$", "Erysipelatoclostridium", .) %>%
    # Then collapse unknown variants
    gsub("^[Uu]ncultured.*$", "Uncultured", .) %>%
    gsub("^[Uu]nclassified.*$", "Unclassified", .) %>%
    gsub("^[Uu]nknown.*$", "Unknown", .) %>%
    gsub("^[Uu]nidentified.*$", "Unidentified", .)
  
  # Convert to data frame with genus info
  otu_df <- as.data.frame(otu_mat)
  otu_df$Genus_clean <- tax_df$Genus_clean
  
  # Aggregate by genus (this will combine all the variants)
  otu_collapsed <- otu_df %>%
    group_by(Genus_clean) %>%
    summarise(across(everything(), sum), .groups = 'drop')
  
  # Convert back to matrix format
  otu_collapsed_mat <- as.matrix(otu_collapsed[, -1])  # Remove genus column
  rownames(otu_collapsed_mat) <- otu_collapsed$Genus_clean
  
  # Create new tax_table with unique genera
  unique_genera <- unique(tax_df$Genus_clean)
  tax_collapsed <- data.frame(
    Kingdom = "Bacteria",
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = unique_genera,
    row.names = unique_genera
  )
  
  # Create new phyloseq object
  ps_collapsed <- phyloseq(
    otu_table(otu_collapsed_mat, taxa_are_rows = TRUE),
    tax_table(as.matrix(tax_collapsed)),
    sample_data(ps_merged)
  )
  
  # Report what was collapsed
  cat("Before collapsing:", ntaxa(ps_merged), "taxa\n")
  cat("After collapsing:", ntaxa(ps_collapsed), "taxa\n")
  cat("Taxa collapsed:", ntaxa(ps_merged) - ntaxa(ps_collapsed), "\n")
  
  # Show examples of what got collapsed
  original_genera <- sort(unique(tax_df$Genus))
  collapsed_genera <- sort(unique(tax_df$Genus_clean))
  
  uncultured_variants <- original_genera[grepl("[Uu]ncultured", original_genera)]
  unclassified_variants <- original_genera[grepl("[Uu]nclassified", original_genera)]
  
  if(length(uncultured_variants) > 0) {
    cat("\nUncultured variants found:", length(uncultured_variants), "\n")
    cat("Examples:", paste(head(uncultured_variants, 5), collapse = ", "), "\n")
  }
  
  if(length(unclassified_variants) > 0) {
    cat("Unclassified variants found:", length(unclassified_variants), "\n")
    cat("Examples:", paste(head(unclassified_variants, 5), collapse = ", "), "\n")
  }
  
  return(ps_collapsed)
}
# Remove trees from each study (but keep them separate)
ps_genus_dose_sanitized_before <- lapply(ps_genus_dose_sanitized_before, function(ps) {
  phyloseq(otu_table(ps), tax_table(ps), sample_data(ps))
})

# Apply unknown taxa collapsing to each study individually
ps_genus_final <- lapply(ps_genus_dose_sanitized_before, function(ps) {
  collapse_unknown_after_merge(ps)  # This function works on individual studies too
})

####FINAL##
ps_genus_final
ps_genus_dose_sanitized <- ps_genus_final


create_corrected_collapse_summary <- function(comparison_list) {
  summary_df <- data.frame(
    Study = character(),
    Total_Before = numeric(),
    Total_After = numeric(),
    Total_Collapsed = numeric(),
    Uncultured_Type = numeric(),
    Unclassified_Type = numeric(),
    Other_Collapsed = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(study_name in names(comparison_list)) {
    comp_data <- comparison_list[[study_name]]
    collapsed_taxa <- comp_data$Taxa_Name[comp_data$Status == "Collapsed"]
    kept_taxa <- comp_data$Taxa_Name[comp_data$Status == "Kept"]
    
    total_before <- length(collapsed_taxa) + length(kept_taxa)
    total_after <- length(kept_taxa)
    total_collapsed <- length(collapsed_taxa)
    
    # Categorize collapsed taxa (no double counting)
    uncultured_count <- sum(grepl("^[Uu]ncultured", collapsed_taxa))
    unclassified_count <- sum(grepl("^[Uu]nclassified", collapsed_taxa))
    
    # Everything else that was collapsed but isn't uncultured/unclassified
    other_collapsed <- total_collapsed - uncultured_count - unclassified_count
    
    summary_df <- rbind(summary_df, data.frame(
      Study = study_name,
      Total_Before = total_before,
      Total_After = total_after,
      Total_Collapsed = total_collapsed,
      Uncultured_Type = uncultured_count,
      Unclassified_Type = unclassified_count,
      Other_Collapsed = other_collapsed,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_df)
}


# Function to create detailed side-by-side comparison
create_detailed_comparison <- function(ps_before_list, ps_after_list) {
  
  all_comparisons <- list()
  
  for(study_name in names(ps_before_list)) {
    cat(paste0("\n=== ", study_name, " DETAILED COMPARISON ===\n"))
    
    # Get sorted taxa names
    taxa_before <- sort(taxa_names(ps_before_list[[study_name]]))
    taxa_after <- sort(taxa_names(ps_after_list[[study_name]]))
    
    # Create a comprehensive comparison dataframe
    all_taxa <- unique(c(taxa_before, taxa_after))
    
    comparison_df <- data.frame(
      Taxon = all_taxa,
      Before_Collapsing = ifelse(all_taxa %in% taxa_before, "Present", "Absent"),
      After_Collapsing = ifelse(all_taxa %in% taxa_after, "Present", "Absent"),
      Status = NA,
      stringsAsFactors = FALSE
    )
    
    # Assign status
    comparison_df$Status[comparison_df$Before_Collapsing == "Present" & comparison_df$After_Collapsing == "Present"] <- "Kept"
    comparison_df$Status[comparison_df$Before_Collapsing == "Present" & comparison_df$After_Collapsing == "Absent"] <- "Collapsed"
    comparison_df$Status[comparison_df$Before_Collapsing == "Absent" & comparison_df$After_Collapsing == "Present"] <- "New"
    
    # Sort by status and then by name
    comparison_df <- comparison_df[order(comparison_df$Status, comparison_df$Taxon), ]
    
    # Print summary
    cat("Total before:", length(taxa_before), "\n")
    cat("Total after:", length(taxa_after), "\n")
    cat("Kept:", sum(comparison_df$Status == "Kept"), "\n")
    cat("Collapsed:", sum(comparison_df$Status == "Collapsed"), "\n")
    if(any(comparison_df$Status == "New")) cat("New:", sum(comparison_df$Status == "New"), "\n")
    
    # Show collapsed taxa
    collapsed_taxa <- comparison_df$Taxon[comparison_df$Status == "Collapsed"]
    if(length(collapsed_taxa) > 0) {
      cat("\nCOLLAPSED TAXA:\n")
      for(i in 1:length(collapsed_taxa)) {
        cat(sprintf("%3d. %s\n", i, collapsed_taxa[i]))
      }
    }
    
    # Store comparison
    all_comparisons[[study_name]] <- comparison_df
    
    cat("\n", rep("-", 60), "\n")
  }
  
  return(all_comparisons)
}

# Function to get all taxa from ps_genus_dose_sanitized
get_all_taxa_lists <- function(ps_list) {
  cat("=== ALL TAXA IN PS_GENUS_DOSE_SANITIZED ===\n\n")
  
  all_taxa_lists <- list()
  
  for(study_name in names(ps_list)) {
    taxa <- sort(taxa_names(ps_list[[study_name]]))
    cat(paste0(study_name, " (", length(taxa), " taxa):\n"))
    
    # Print in columns for readability
    for(i in 1:length(taxa)) {
      cat(sprintf("%3d. %s\n", i, taxa[i]))
    }
    
    all_taxa_lists[[study_name]] <- taxa
    cat("\n", rep("=", 50), "\n\n")
  }
  
  return(all_taxa_lists)
}

# Run the detailed comparison
detailed_comparison <- create_detailed_comparison(ps_genus_dose_sanitized, ps_genus_final)

# Get all taxa lists
cat("\n\n")
all_taxa_lists <- get_all_taxa_lists(ps_genus_dose_sanitized)

# Optional: Create a master comparison table showing which taxa are in which studies
create_master_taxa_table <- function(taxa_lists) {
  all_unique_taxa <- unique(unlist(taxa_lists))
  
  master_df <- data.frame(
    Taxon = sort(all_unique_taxa),
    stringsAsFactors = FALSE
  )
  
  for(study_name in names(taxa_lists)) {
    master_df[[study_name]] <- ifelse(master_df$Taxon %in% taxa_lists[[study_name]], "Yes", "No")
  }
  
  # Add a count column
  master_df$Studies_Present <- rowSums(master_df[, names(taxa_lists)] == "Yes")
  
  return(master_df)
}

# Create master table
cat("\n=== MASTER TAXA TABLE (Which taxa are in which studies) ===\n")
master_taxa_table <- create_master_taxa_table(all_taxa_lists)
print(head(master_taxa_table, 20))  # Show first 20 rows

# Check the BEZI family directly
bezi_info <- BEZI()
cat("BEZI Family Information:\n")
cat("Link for mu (mean parameter):", bezi_info$mu.link, "\n")
cat("Link for sigma (precision):", bezi_info$sigma.link, "\n") 
cat("Link for nu (zero-inflation):", bezi_info$nu.link, "\n")

#### SPLIT STUDY2 IN THE LIST ####
cat("\n=== SPLITTING STUDY2 BY ADMINISTRATION ROUTE (LIST VERSION) ===\n")

# Check current list structure
cat("Studies before split:\n")
print(names(ps_genus_dose_sanitized))

# Check if Study2 exists
if("Study2" %in% names(ps_genus_dose_sanitized)) {
  
  # Get Study2 phyloseq object
  study2_ps <- ps_genus_dose_sanitized$Study2
  
  # Check AdminRoute distribution
  meta_study2 <- data.frame(sample_data(study2_ps))
  cat("\nStudy2 AdminRoute distribution:\n")
  print(table(meta_study2$AdminRoute))
  cat("Total samples in Study2:", nsamples(study2_ps), "\n")
  
  # Create Study2_Food: subset samples where AdminRoute == "Food"
  ps_genus_dose_sanitized$Study2_Food <- subset_samples(
    study2_ps, 
    AdminRoute == "Food"
  )
  
  # Create Study2_IP: subset samples where AdminRoute == "IP"  
  ps_genus_dose_sanitized$Study2_IP <- subset_samples(
    study2_ps, 
    AdminRoute == "IP"
  )
  
  # Remove original Study2 from the list
  ps_genus_dose_sanitized$Study2 <- NULL
  
  # Report results
  cat("\n✓ Split completed:\n")
  cat("  - Study2_Food:", nsamples(ps_genus_dose_sanitized$Study2_Food), "samples\n")
  cat("  - Study2_IP:", nsamples(ps_genus_dose_sanitized$Study2_IP), "samples\n")
  cat("  - Original Study2 removed from list\n")
  
} else {
  cat("⚠ Study2 not found in ps_genus_dose_sanitized\n")
}

# Show final list structure
cat("\n=== Studies after split ===\n")
print(names(ps_genus_dose_sanitized))
cat("\nSample counts per study:\n")
for(study_name in names(ps_genus_dose_sanitized)) {
  cat(study_name, ":", nsamples(ps_genus_dose_sanitized[[study_name]]), "samples\n")
}

###### GAMLSS ANALYZE ALL TAXA ACROSS STUDIES ######
 length(meta$StandardTreatment)  # Should equal ncol(rel_abund)
# GAMLSS function with only treatment effect (no AdminRoute)
run_gamlss_treatment_only <- function(ps_list) {
  # Create a list to store results
  all_results <- list()
  
  for(study_name in names(ps_list)) {
    cat("\n\n========== ANALYZING", study_name, "==========\n")
    
    # Extract phyloseq object
    ps <- ps_list[[study_name]]
    
    # Get OTU table
    otu <- as(otu_table(ps), "matrix")
    
    # Check if taxa are rows or columns
    if(!taxa_are_rows(ps)) {
      otu <- t(otu)  # This ensures taxa are in rows
    }
    
    # Get number of taxa for logging
    cat("Processing", nrow(otu), "taxa in", study_name, "\n")
    
    # Calculate relative abundance correctly (sample-wise)
    rel_abund_by_sample <- apply(otu, 2, function(x) x / sum(x) * 100)
    
    # Keep taxa as rows, samples as columns
    rel_abund <- rel_abund_by_sample
    
    # Get metadata
    meta <- data.frame(sample_data(ps))
    
    # Check dimensions match
    if(ncol(rel_abund) != nrow(meta)) {
      stop("Dimension mismatch in ", study_name, 
           ": rel_abund has ", ncol(rel_abund), 
           " columns but meta has ", nrow(meta), " rows")
    }
    
    # Get all taxa
    taxa_list <- rownames(rel_abund)
    total_taxa <- length(taxa_list)
    
    cat("Processing", total_taxa, "taxa in", study_name, "\n")
    
    # Create a list to store results for this study
    study_results <- list()
    
    # Processing counters
    success_count <- 0
    
    # Process each taxon
    for(i in 1:total_taxa) {
      taxon <- taxa_list[i]
      
      # Print progress
      if(i %% 10 == 0) {
        cat("Processing taxon", i, "of", total_taxa, "\n")
      }
      
      # Get abundance values
      abundance_values <- as.numeric(rel_abund[taxon, ])
      
      # Skip taxa with no variation
      if(length(unique(abundance_values)) <= 1) {
        next
      }
      
      # Store taxon_data locally in this environment
      data_for_model <- data.frame(
        abundance = abundance_values,
        StandardTreatment = factor(meta$StandardTreatment, levels = c("Control", "Rapamycin")),
        row.names = colnames(rel_abund)
      )
      
      # # Handle zero-inflation
      # epsilon <- 0.5 / 100  # 0.005 in percent units, or 0.00005 in proportion units
      # data_for_model$abundance[data_for_model$abundance == 0] <- epsilon
      # data_for_model$abundance[data_for_model$abundance == 100] <- 100 - epsilon
      # 
      # Scale to 0-1 (no epsilon replacement - let BEZI handle zeros)
      data_for_model$scaled <- data_for_model$abundance / 100
      
      # Run GAMLSS - ONLY include treatment effect, NOT AdminRoute
      model <- NULL
      tryCatch({
        # Use the same model for ALL studies - only treatment effect
        model <- gamlss(
          scaled ~ StandardTreatment,  # ONLY treatment effect
          family = BEZI(),  # Zero-inflated beta for microbiome data
          data = data_for_model,
          trace = FALSE
        )
      }, error = function(e) {
        # Skip silently
      })
      
      # Check if model was successfully fit
      if(!is.null(model)) {
        # Try to extract from summary output
        summ_output <- tryCatch({
          capture.output(summary(model))
        }, error = function(e) {
          return(NULL)
        })
        
        if(!is.null(summ_output)) {
          # Extract treatment effect line
          treatment_line <- grep("StandardTreatmentRapamycin", summ_output, value = TRUE)
          
          if(length(treatment_line) > 0) {
            # Find the line with actual numeric values
            for(line in treatment_line) {
              # Split the line
              values <- strsplit(trimws(line), "\\s+")[[1]]
              
              # Check if it has numeric values
              if(length(values) >= 4 && !is.na(suppressWarnings(as.numeric(values[2])))) {
                # First try the standard format
                effect <- as.numeric(values[2])
                std_error <- as.numeric(values[3])
                
                # Check if t-value and p-value are available
                if(length(values) >= 5) {
                  t_value <- as.numeric(values[4])
                  p_value <- as.numeric(values[5])
                } else {
                  # If not, calculate p-value from t-value
                  t_value <- effect / std_error
                  df <- nrow(data_for_model) - length(coef(model))
                  p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
                }
                
                # Store results
                result <- data.frame(
                  taxon = taxon,
                  effect = effect,
                  std_error = std_error,
                  t_value = t_value,
                  p_value = p_value,
                  stringsAsFactors = FALSE
                )
                
                # Add to results list
                study_results[[length(study_results) + 1]] <- result
                success_count <- success_count + 1
                
                # Break the loop once we found the values
                break
              }
            }
          }
        }
      }
    }
    
    cat("Successfully analyzed", success_count, "taxa out of", total_taxa, "\n")
    
    # Combine results for this study
    if(length(study_results) > 0) {
      study_df <- do.call(rbind, study_results)
      
      # Add FDR correction
      study_df$p_adjust <- p.adjust(study_df$p_value, method = "BH")
      
      # Store in all results
      all_results[[study_name]] <- study_df
      
      # Show summary
      cat("\nResults summary for", study_name, ":\n")
      cat("- Successfully analyzed", nrow(study_df), "taxa\n")
      
      # Report significant taxa
      sig_taxa <- sum(study_df$p_value < 0.05, na.rm = TRUE)
      cat("- Significant taxa (p < 0.05):", sig_taxa, "\n")
      
      # Show top significant taxa
      if(sig_taxa > 0) {
        top_taxa <- head(study_df[order(study_df$p_value), ], 5)
        cat("- Top significant taxa:\n")
        print(top_taxa[, c("taxon", "effect", "t_value", "p_value")])
      }
    } else {
      cat("No successful analyses for", study_name, "\n")
    }
  }
  
  return(all_results)
}

# Run the analysis without administration route
treatment_only_results <- run_gamlss_treatment_only(ps_genus_dose_sanitized)

treatment_only_results

ps_genus_dose_sanitized

####Check Results ####
# Get all taxa with p_value < 0.05 and show Study, taxon, and effect
all_sig_taxa <- list()

for(study_name in names(treatment_only_results)) {
  study_df <- treatment_only_results[[study_name]]
  
  # Filter for p_value < 0.05
  sig_taxa <- study_df[study_df$p_value < 0.05, ]
  
  # Add study name
  sig_taxa$Study <- study_name
  
  all_sig_taxa[[study_name]] <- sig_taxa
}


dash <- bind_rows(treatment_only_results, .id = "study") |>
  group_by(study) |>
  summarise(
    n_tested     = n(),
    n_raw_sig    = sum(p_value   < 0.05, na.rm = TRUE),
    n_fdr_sig    = sum(p_adjust  < 0.05, na.rm = TRUE),
    med_abs_beta = median(abs(effect), na.rm = TRUE),
    .groups = "drop"
  )

print(dash)

#Splitted # A tibble: 4 × 5
# study       n_tested n_raw_sig n_fdr_sig med_abs_beta
# <chr>          <int>     <int>     <int>        <dbl>
# 1 Study2_Food       97        12         5        0.488
# 2 Study2_IP         74         6         4        0.268
# 3 Study5            92         9         0        0.647
# 4 Study7            74        11        10        0.715
#

# Look at the structure of our treatment_only_results
str_treatment_results <- function() {
  cat("Structure of treatment_only_results:\n")
  for(study_name in names(treatment_only_results)) {
    cat("\nStudy:", study_name, "\n")
    print(str(treatment_only_results[[study_name]]))
    
    # Print first few rows
    cat("First few rows:\n")
    print(head(treatment_only_results[[study_name]]))
    break  # Just look at one study for brevity
  }
}

str_treatment_results()

# with BEZI no espilon 
# First few rows:
# taxon    effect    std_error  t_value  p_value    p_adjust
# A2 -0.323700 2.176000e-01   -1.487 1.47e-01 0.577111111
# ASF356  0.066110 2.069700e-01    0.319 7.51e-01 1.000000000
# Acetatifactor -0.348900 1.709000e-01   -2.041 4.96e-02 0.350506667
# Aerococcus  8.445000 2.425356e+04    0.000 1.00e+00 1.000000000
# Akkermansia  0.274200 4.181000e-01    0.656 5.17e-01 1.000000000
# Alistipes  0.204800 2.497000e-01    0.820 4.18e-01 1.000000000

# Count how many std_error values are missing in each data set
na_se_summary <- vapply(
  treatment_only_results,
  function(df) sum(is.na(df$std_error)),
  FUN.VALUE = integer(1)
)

cat("Missing standard-error counts:\n")
print(na_se_summary)

##### Filter out high SE taxa####
filtered_results <- lapply(treatment_only_results, function(study_data) {
  study_data[study_data$std_error <= 5, ]
})
# Extract list of high SE taxa from treatment_only_results
check_separation_patterns <- function(treatment_results, ps_list, se_threshold = 5) {
  # First identify taxa with high standard errors
  high_se_taxa_list <- lapply(treatment_results, function(x) {
    x$taxon[x$std_error > se_threshold]
  })
  
  # Print summary of high SE taxa counts
  cat("===== HIGH STANDARD ERROR TAXA (SE >", se_threshold, ") =====\n")
  for(study_name in names(high_se_taxa_list)) {
    cat(study_name, ":", length(high_se_taxa_list[[study_name]]), "taxa\n")
  }
  
  separation_results <- list()
  
  # Process each study
  for(study_name in names(high_se_taxa_list)) {
    # Get the list of high SE taxa for this study
    high_se_taxa <- high_se_taxa_list[[study_name]]
    
    # Skip if no high SE taxa in this study
    if(length(high_se_taxa) == 0) {
      cat(study_name, ": No high standard error taxa to check\n")
      next
    }
    
    # Get the phyloseq object
    ps <- ps_list[[study_name]]
    
    # Extract OTU table and metadata
    otu <- as(otu_table(ps), "matrix")
    if(!taxa_are_rows(ps)) otu <- t(otu)
    meta <- data.frame(sample_data(ps))
    
    # Create results data frame for this study
    study_results <- data.frame(
      taxon = character(),
      control_n = integer(),
      control_present = integer(),
      control_percent = numeric(),
      rapa_n = integer(),
      rapa_present = integer(),
      rapa_percent = numeric(),
      diff_percent = numeric(),
      separation_type = character(),
      fisher_p = numeric(),
      stringsAsFactors = FALSE
    )
    
    # For each high SE taxon, check presence patterns
    for(taxon in high_se_taxa) {
      # Check if taxon exists in OTU table
      if(taxon %in% rownames(otu)) {
        # Get counts for this taxon
        counts <- otu[taxon, ]
        
        # Calculate presence (binary) for each group
        control_idx <- meta$StandardTreatment == "Control"
        rapa_idx <- meta$StandardTreatment == "Rapamycin"
        
        control_n <- sum(control_idx)
        rapa_n <- sum(rapa_idx)
        
        control_present <- sum(counts[control_idx] > 0)
        rapa_present <- sum(counts[rapa_idx] > 0)
        
        control_percent <- control_present / control_n * 100
        rapa_percent <- rapa_present / rapa_n * 100
        
        # Calculate difference in presence percentages
        diff_percent <- rapa_percent - control_percent
        
        # Determine separation type
        if(control_present == 0 && rapa_present == 0) {
          separation_type <- "All Zero"
        } else if(control_present == 0 && rapa_present > 0) {
          if(rapa_present == rapa_n) {
            separation_type <- "Complete: Rapa Only"
          } else {
            separation_type <- "Quasi-complete: Rapa Sometimes"
          }
        } else if(control_present > 0 && rapa_present == 0) {
          if(control_present == control_n) {
            separation_type <- "Complete: Control Only"
          } else {
            separation_type <- "Quasi-complete: Control Sometimes"
          }
        } else if(control_present == control_n && rapa_present == rapa_n) {
          separation_type <- "All Present"
        } else if(control_present == control_n && rapa_present < rapa_n) {
          separation_type <- "Quasi-complete: Control Always, Rapa Sometimes"
        } else if(control_present < control_n && rapa_present == rapa_n) {
          separation_type <- "Quasi-complete: Control Sometimes, Rapa Always"
        } else {
          # Both groups have some present and some absent
          separation_type <- "No Separation"
        }
        
        # Run Fisher's exact test
        presence <- counts > 0
        cont_table <- table(
          Treatment = meta$StandardTreatment,
          Present = factor(presence, levels = c(FALSE, TRUE))
        )
        
        fisher_test <- tryCatch({
          fisher.test(cont_table)
        }, error = function(e) {
          list(p.value = NA)
        })
        
        # Add to results
        taxon_result <- data.frame(
          taxon = taxon,
          control_n = control_n,
          control_present = control_present,
          control_percent = control_percent,
          rapa_n = rapa_n,
          rapa_present = rapa_present,
          rapa_percent = rapa_percent,
          diff_percent = diff_percent,
          separation_type = separation_type,
          fisher_p = fisher_test$p.value,
          stringsAsFactors = FALSE
        )
        
        study_results <- rbind(study_results, taxon_result)
      }
    }
    
    # Add study name and count separation types
    if(nrow(study_results) > 0) {
      study_results$study <- study_name
      
      # Count separation types
      sep_counts <- table(study_results$separation_type)
      cat("\n===== SEPARATION PATTERNS IN", study_name, "=====\n")
      print(sep_counts)
      
      # Store results
      separation_results[[study_name]] <- study_results
    }
  }
  
  # Combined results across studies
  if(length(separation_results) > 0) {
    all_results <- do.call(rbind, separation_results)
    
    # Overall summary
    cat("\n===== OVERALL SEPARATION PATTERNS =====\n")
    print(table(all_results$separation_type))
    
    # Calculate correlation between separation patterns and standard errors
    cat("\n===== RELATIONSHIP BETWEEN SEPARATION AND STANDARD ERRORS =====\n")
    
    # For each study, calculate correlation between absolute difference in presence % and standard error
    for(study_name in names(treatment_results)) {
      if(study_name %in% names(separation_results)) {
        study_data <- separation_results[[study_name]]
        study_orig <- treatment_results[[study_name]]
        
        # Merge data
        combined <- merge(study_data, study_orig, by = "taxon")
        
        if(nrow(combined) > 0) {
          # Calculate correlation
          cor_test <- cor.test(abs(combined$diff_percent), combined$std_error, method = "spearman")
          
          cat(study_name, "- Spearman correlation between |∆Presence%| and SE:", 
              round(cor_test$estimate, 3), "(p =", format.pval(cor_test$p.value, digits = 3), ")\n")
          
          # Calculate mean SE by separation type
          se_by_type <- aggregate(std_error ~ separation_type, combined, mean)
          cat("Mean standard errors by separation type:\n")
          print(se_by_type)
        }
      }
    }
  }
  
  return(separation_results)
}

# Run the separation check for treatment_only_results
separation_analysis <- check_separation_patterns(treatment_only_results, ps_genus_dose_sanitized)

# Create visualizations
plot_separation_patterns <- function(separation_results) {
  # Check if any results exist
  if(length(separation_results) == 0) {
    cat("No separation results to plot\n")
    return(NULL)
  }
  
  # Combine results
  all_data <- do.call(rbind, separation_results)
  
  # Basic separation pattern plot
  cat("\nCreating plots of separation patterns...\n")
  
  # 1. Scatter plot of presence percentages
  plot(all_data$control_percent, all_data$rapa_percent, 
       xlab = "% Present in Control", ylab = "% Present in Rapamycin",
       main = "Presence Patterns for High SE Taxa",
       pch = 19, col = as.factor(all_data$separation_type),
       xlim = c(0, 100), ylim = c(0, 100))
  abline(0, 1, lty = 2, col = "gray")
  legend("topleft", legend = unique(all_data$separation_type), 
         col = 1:length(unique(all_data$separation_type)), pch = 19,
         cex = 0.8)
  
  # 2. Histograms of difference in presence
  hist(all_data$diff_percent, 
       main = "Distribution of Presence Differences", 
       xlab = "Difference in % Present (Rapa - Control)",
       breaks = 20)
  abline(v = 0, lty = 2, col = "red")
  
  # 3. Separation types by study
  study_counts <- table(all_data$study, all_data$separation_type)
  cat("\nSeparation types by study:\n")
  print(study_counts)
  
  # 4. Count taxa with various degrees of separation
  complete_separation <- sum(grepl("Complete:", all_data$separation_type))
  quasi_separation <- sum(grepl("Quasi-complete:", all_data$separation_type))
  no_separation <- sum(all_data$separation_type == "No Separation")
  
  cat("\nSummary of separation patterns:\n")
  cat("- Complete separation:", complete_separation, 
      "(", round(complete_separation/nrow(all_data)*100, 1), "%)\n")
  cat("- Quasi-complete separation:", quasi_separation, 
      "(", round(quasi_separation/nrow(all_data)*100, 1), "%)\n")
  cat("- No separation:", no_separation, 
      "(", round(no_separation/nrow(all_data)*100, 1), "%)\n")
  
  # Return the combined data
  return(all_data)
}

# Plot the separation patterns
all_separation_data <- plot_separation_patterns(separation_analysis)

# Additional analysis: Compare standard errors for taxa with different separation patterns
analyze_se_distribution <- function(treatment_results, separation_results) {
  cat("\n===== STANDARD ERROR ANALYSIS BY SEPARATION TYPE =====\n")
  
  # Combine all separation results
  all_separation <- do.call(rbind, separation_results)
  
  # Create a data frame to store all standard errors
  all_se_data <- data.frame()
  
  for(study_name in names(treatment_results)) {
    study_data <- treatment_results[[study_name]]
    
    # Add study name
    study_data$study <- study_name
    
    # Add to combined data
    all_se_data <- rbind(all_se_data, study_data[, c("taxon", "study", "std_error")])
  }
  
  # Merge with separation data
  merged_data <- merge(all_se_data, all_separation, by = c("taxon", "study"), all.x = TRUE)
  
  # Fill in separation_type for taxa that weren't in high_se_taxa_list
  merged_data$separation_type[is.na(merged_data$separation_type)] <- "Non-high SE"
  
  # Summary statistics by separation type
  cat("Number of taxa in each separation category:\n")
  print(table(merged_data$separation_type))
  
  cat("\nStandard error summary by separation type:\n")
  summary_stats <- aggregate(std_error ~ separation_type, merged_data, 
                             function(x) c(mean = mean(x), median = median(x), 
                                           min = min(x), max = max(x), n = length(x)))
  print(summary_stats)
  
  # Create box plot of standard errors by separation type
  boxplot(std_error ~ separation_type, data = merged_data,
          main = "Standard Errors by Separation Type",
          xlab = "Separation Type", ylab = "Standard Error",
          las = 2, cex.axis = 0.7)
  
  # Create histogram of standard errors
  hist(merged_data$std_error[merged_data$std_error < 100], 
       main = "Distribution of Standard Errors",
       xlab = "Standard Error", breaks = 50)
  abline(v = 5, col = "red", lty = 2)
  
  return(merged_data)
}

# Run standard error analysis
se_analysis <- analyze_se_distribution(treatment_only_results, separation_analysis)

#create table of excludetd taxa 
make_high_se_table <- function(treatment_results, separation_analysis, meta_list, se_threshold = 5) {
  all_rows <- list()
  
  for (study in names(treatment_results)) {
    res <- treatment_results[[study]]
    sep <- separation_analysis[[study]]
    # Only taxa with SE > threshold
    res_high_se <- res %>% filter(std_error > se_threshold)
    # For each high SE taxon, add separation type
    for (i in seq_len(nrow(res_high_se))) {
      taxon <- res_high_se$taxon[i]
      se_type <- sep$separation_type[match(taxon, sep$taxon)]
      # Get mean abundance per group if meta_list (list of phyloseq objects) is supplied
      if (!missing(meta_list)) {
        ps <- meta_list[[study]]
        otu <- as(otu_table(ps), "matrix")
        if (!taxa_are_rows(ps)) otu <- t(otu)
        meta <- data.frame(sample_data(ps))
        abund <- otu[taxon, ]
        means <- tapply(abund, meta$StandardTreatment, mean)
        mean_ctrl <- means["Control"]
        mean_rapa <- means["Rapamycin"]
      } else {
        mean_ctrl <- NA
        mean_rapa <- NA
      }
      row <- data.frame(
        Study = study,
        Taxon = taxon,
        SeparationType = se_type,
        MeanAbundanceControl = mean_ctrl,
        MeanAbundanceRapamycin = mean_rapa,
        StdError = res_high_se$std_error[i],
        Pvalue = res_high_se$p_value[i],
        stringsAsFactors = FALSE
      )
      all_rows[[length(all_rows) + 1]] <- row
    }
  }
  
  df <- do.call(rbind, all_rows)
  # Optional: sort by Study, then SeparationType, then Taxon
  df <- df %>% arrange(Study, SeparationType, Taxon)
  return(df)
}

# Example usage:
# If you have your phyloseq objects in a list called ps_genus_dose_sanitized:
high_se_summary <- make_high_se_table(
  treatment_only_results, separation_analysis, ps_genus_dose_sanitized, se_threshold = 5
)

# Show table

table.high_se_summary <- high_se_summary %>%
  mutate(
    MeanAbundanceControl = round(MeanAbundanceControl, 2),
    MeanAbundanceRapamycin = round(MeanAbundanceRapamycin, 2)
  )
kable(
  table.high_se_summary,
  caption = "Mean abundance (%) and separation type for taxa with quasi-complete separation"
)
# Tabulate presence by treatment and maybe by cage
table(Abundant = table.high_se_summary$MeanAbundanceControl > 0, Treatment = table.high_se_summary$MeanAbundanceRapamycin > 0)
# or
table(table.high_se_summary$Study, table.high_se_summary$SeparationType)

ggplot(table.high_se_summary, aes(x=SeparationType, y=MeanAbundanceControl)) + 
  geom_boxplot() + 
  labs(title="Abundance by Separation Type in Controls")


#### Filtered out HIGH SE CHECK ####

get_group_presence <- function(ps, group_var = "StandardTreatment") {
  # Get abundance table (as % if not already)
  abund <- as(otu_table(ps), "matrix")
  if(!taxa_are_rows(ps)) abund <- t(abund)
  # Convert to relative abundance if necessary
  abund <- apply(abund, 2, function(x) x / sum(x) * 100)
  
  meta <- data.frame(sample_data(ps))
  # Aggregate by group (e.g., Control vs. Rapamycin)
  by_group <- split(colnames(abund), meta[[group_var]])
  # For each genus, is mean abundance > 0 in each group?
  genus_names <- rownames(abund)
  result <- data.frame(
    Taxon = genus_names,
    MeanAbundanceControl = rowMeans(abund[, by_group$Control, drop=FALSE], na.rm=TRUE),
    MeanAbundanceRapamycin = rowMeans(abund[, by_group$Rapamycin, drop=FALSE], na.rm=TRUE)
  )
  result
}

# Loop over studies
all_overlap <- map2_dfr(
  ps_genus_dose_sanitized,
  names(ps_genus_dose_sanitized),
  ~ get_group_presence(.x) %>% mutate(Study = .y)
)
all_overlap <- all_overlap %>%
  mutate(
    InControl = MeanAbundanceControl > 0,
    InRapamycin = MeanAbundanceRapamycin > 0,
    OverlapType = case_when(
      InControl & InRapamycin ~ "Both",
      InControl & !InRapamycin ~ "Control only",
      !InControl & InRapamycin ~ "Rapamycin only",
      TRUE ~ "Neither"
    )
  )

table(all_overlap$OverlapType)
# Or by study:
table(all_overlap$Study, all_overlap$OverlapType)

# Which separated taxa appear in multiple studies?
separated_taxa <- all_overlap %>% 
  filter(OverlapType %in% c("Control only", "Rapamycin only"))

# Cross-study consistency
separated_cross_study <- separated_taxa %>%
  group_by(Taxon, OverlapType) %>%
  summarise(n_studies = n_distinct(Study), .groups = "drop") %>%
  filter(n_studies > 1) %>%
  arrange(desc(n_studies), OverlapType)

print(separated_cross_study)

# After the separated_cross_study analysis, add this:
separated_with_study_details <- separated_taxa %>%
  group_by(Taxon, OverlapType) %>%
  summarise(
    n_studies = n_distinct(Study),
    studies_list = paste(Study, collapse = ", "),
    mean_abundance_ctrl = mean(MeanAbundanceControl, na.rm = TRUE),
    mean_abundance_rapa = mean(MeanAbundanceRapamycin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_studies), OverlapType, Taxon)

print(separated_with_study_details)

# Example: Check Akkermansia across all studies
high_se_summary %>% 
  filter(Taxon == "Akkermansia") %>%
  select(Study, Taxon, SeparationType, MeanAbundanceControl, MeanAbundanceRapamycin, StdError)

# Or check all "Complete: Control Only" taxa
high_se_summary %>% 
  filter(SeparationType == "Complete: Control Only") %>%
  select(Study, Taxon, MeanAbundanceControl, MeanAbundanceRapamycin)

# Summary table
table(separated_cross_study$OverlapType, separated_cross_study$n_studies)

consistent_separated <- separated_taxa %>%
  group_by(Taxon, OverlapType) %>%
  summarise(n_studies = n_distinct(Study)) %>%
  filter(n_studies > 1)

# Which taxa are present in both groups in any study?
both_taxa <- all_overlap %>% filter(OverlapType == "Both")
head(both_taxa[, c("Study", "Taxon", "MeanAbundanceControl", "MeanAbundanceRapamycin")])

# Which taxa are present in both groups in more than one study?
both_across_studies <- both_taxa %>%
  group_by(Taxon) %>%
  summarize(n_studies = n_distinct(Study)) %>%
  filter(n_studies > 1)
both_across_studies

all_overlap %>%
  filter(Taxon %in% both_across_studies$Taxon) %>%
  arrange(Taxon, Study) %>%
  select(Study, Taxon, MeanAbundanceControl, MeanAbundanceRapamycin)

#### Filtering approach for meta-analysis#####
filter_for_meta_analysis <- function(treatment_results, separation_analysis, max_se = 5, min_se = 0.01) {
  # Create filtered results
  filtered_results <- list()
  
  # Track filtering statistics
  filter_stats <- data.frame(
    study = character(),
    n_original = integer(),
    n_high_se = integer(),
    n_low_se = integer(),
    n_na_pvalue = integer(),
    n_final = integer(),
    stringsAsFactors = FALSE
  )
  
  for(study_name in names(treatment_results)) {
    study_data <- treatment_results[[study_name]]
    n_original <- nrow(study_data)
    
    # Count exclusions
    n_high_se <- sum(study_data$std_error > max_se, na.rm = TRUE)
    n_low_se <- sum(study_data$std_error < min_se, na.rm = TRUE)
    n_na_pvalue <- sum(is.na(study_data$p_value))
    
    # Filter based on standard error AND p-value
    filtered_results[[study_name]] <- study_data %>%
      filter(
        std_error <= max_se,           # Remove high SE (>5)
        std_error >= min_se,           # Remove tiny SE (<0.01)
        !is.na(p_value)                # Remove NA p-values (model failures)
      )
    
    n_final <- nrow(filtered_results[[study_name]])
    
    cat(study_name, ":\n")
    cat("  Original taxa:", n_original, "\n")
    cat("  Excluded - High SE (>", max_se, "):", n_high_se, "\n")
    cat("  Excluded - Low SE (<", min_se, "):", n_low_se, "\n")
    cat("  Excluded - NA p-value:", n_na_pvalue, "\n")
    cat("  Final taxa:", n_final, "\n")
    cat("  Total removed:", n_original - n_final, "\n\n")
    
    # Store stats
    filter_stats <- rbind(filter_stats, data.frame(
      study = study_name,
      n_original = n_original,
      n_high_se = n_high_se,
      n_low_se = n_low_se,
      n_na_pvalue = n_na_pvalue,
      n_final = n_final,
      stringsAsFactors = FALSE
    ))
  }
  
  # Check how many taxa remain for meta-analysis
  all_taxa <- unique(unlist(lapply(filtered_results, function(x) x$taxon)))
  
  # Count how many studies each taxon appears in
  taxa_counts <- table(unlist(lapply(filtered_results, function(x) x$taxon)))
  multi_study_taxa <- names(taxa_counts[taxa_counts >= 2])
  
  cat("========== OVERALL FILTERING SUMMARY ==========\n")
  cat("Total unique taxa after filtering:", length(all_taxa), "\n")
  cat("Taxa present in ≥2 studies:", length(multi_study_taxa), "\n\n")
  
  # Print filtering stats table
  cat("Filtering statistics by study:\n")
  print(filter_stats)
  
  # Check if filtered taxa with separation issues are actually significant
  cat("\n========== SEPARATION PATTERNS FOR FILTERED TAXA ==========\n")
  significant_filtered <- 0
  
  for(study_name in names(treatment_results)) {
    if(study_name %in% names(separation_analysis)) {
      study_data <- treatment_results[[study_name]]
      high_se_taxa <- study_data$taxon[study_data$std_error > max_se]
      
      sig_high_se <- sum(study_data$p_value[study_data$std_error > max_se] < 0.05, na.rm = TRUE)
      significant_filtered <- significant_filtered + sig_high_se
      
      sep_data <- separation_analysis[[study_name]]
      
      # Count separation types for high SE taxa
      if(nrow(sep_data) > 0) {
        sep_counts <- table(sep_data$separation_type)
        cat(study_name, ": Separation types for", length(high_se_taxa), "filtered taxa", 
            "(", sig_high_se, "were significant at p < 0.05)\n")
        print(sep_counts)
      }
    }
  }
  
  cat("\nTotal significant taxa filtered out due to high SE:", significant_filtered, "\n")
  
  return(list(
    filtered_data = filtered_results,
    filter_stats = filter_stats
  ))
}

# Run filtering approach with BOTH min and max SE thresholds
filter_results <- filter_for_meta_analysis(
  treatment_only_results, 
  separation_analysis,
  max_se = 5,
  min_se = 0.01
)

# Extract filtered results
filtered_results <- filter_results$filtered_data
filter_stats <- filter_results$filter_stats

#### Check LOW SE taxa (convergence failures)####
check_low_se_taxa <- function(treatment_results, ps_list, min_se = 0.01) {
  
  low_se_summary <- list()
  
  for(study_name in names(treatment_results)) {
    study_data <- treatment_results[[study_name]]
    
    # Find taxa with suspiciously low SE or NA p-values
    low_se_taxa <- study_data %>%
      filter(std_error < min_se | is.na(p_value))
    
    if(nrow(low_se_taxa) > 0) {
      cat("\n========== LOW SE / CONVERGENCE FAILURES IN", study_name, "==========\n")
      cat("Found", nrow(low_se_taxa), "problematic taxa\n\n")
      
      # Get phyloseq object
      ps <- ps_list[[study_name]]
      otu <- as(otu_table(ps), "matrix")
      if(!taxa_are_rows(ps)) otu <- t(otu)
      meta <- data.frame(sample_data(ps))
      
      # For each problematic taxon, check prevalence
      for(i in 1:nrow(low_se_taxa)) {
        taxon <- low_se_taxa$taxon[i]
        
        if(taxon %in% rownames(otu)) {
          abund <- otu[taxon, ]
          rel_abund <- abund / colSums(otu) * 100
          
          # Calculate prevalence by group
          control_idx <- meta$StandardTreatment == "Control"
          rapa_idx <- meta$StandardTreatment == "Rapamycin"
          
          control_prev <- sum(abund[control_idx] > 0) / sum(control_idx) * 100
          rapa_prev <- sum(abund[rapa_idx] > 0) / sum(rapa_idx) * 100
          
          low_se_summary[[paste0(study_name, "_", taxon)]] <- data.frame(
            Study = study_name,
            Taxon = taxon,
            StdError = low_se_taxa$std_error[i],
            Pvalue = low_se_taxa$p_value[i],
            Effect = low_se_taxa$effect[i],
            ControlPrevalence = control_prev,
            RapaPrevalence = rapa_prev,
            Issue = ifelse(is.na(low_se_taxa$p_value[i]), "NA p-value", "Tiny SE"),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if(length(low_se_summary) > 0) {
    low_se_df <- do.call(rbind, low_se_summary)
    rownames(low_se_df) <- NULL
    
    cat("\n========== SUMMARY OF CONVERGENCE FAILURES ==========\n")
    print(low_se_df)
    
    cat("\n⚠️  These taxa have ultra-low prevalence and model convergence issues\n")
    cat("   They will be excluded from meta-analysis\n")
    
    return(low_se_df)
  } else {
    cat("\nNo convergence failures detected!\n")
    return(NULL)
  }
}

# Check for convergence failures
low_se_taxa <- check_low_se_taxa(treatment_only_results, ps_genus_dose_sanitized)


#create table of EXCLUDED taxa (both high and low SE)
make_excluded_taxa_table <- function(treatment_results, separation_analysis, meta_list, 
                                     max_se = 5, min_se = 0.01) {
  all_rows <- list()
  
  for (study in names(treatment_results)) {
    res <- treatment_results[[study]]
    
    # Taxa with HIGH SE
    res_high_se <- res %>% filter(std_error > max_se)
    
    # Taxa with LOW SE or NA p-value
    res_low_se <- res %>% filter(std_error < min_se | is.na(p_value))
    
    # Process HIGH SE taxa
    if(nrow(res_high_se) > 0 && study %in% names(separation_analysis)) {
      sep <- separation_analysis[[study]]
      
      for (i in seq_len(nrow(res_high_se))) {
        taxon <- res_high_se$taxon[i]
        se_type <- sep$separation_type[match(taxon, sep$taxon)]
        
        # Get mean abundance
        ps <- meta_list[[study]]
        otu <- as(otu_table(ps), "matrix")
        if (!taxa_are_rows(ps)) otu <- t(otu)
        meta <- data.frame(sample_data(ps))
        abund <- otu[taxon, ]
        means <- tapply(abund, meta$StandardTreatment, mean)
        
        row <- data.frame(
          Study = study,
          Taxon = taxon,
          ExclusionReason = "High SE (>5)",
          SeparationType = se_type,
          MeanAbundanceControl = means["Control"],
          MeanAbundanceRapamycin = means["Rapamycin"],
          StdError = res_high_se$std_error[i],
          Pvalue = res_high_se$p_value[i],
          stringsAsFactors = FALSE
        )
        all_rows[[length(all_rows) + 1]] <- row
      }
    }
    
    # Process LOW SE taxa
    if(nrow(res_low_se) > 0) {
      for (i in seq_len(nrow(res_low_se))) {
        taxon <- res_low_se$taxon[i]
        
        # Get mean abundance
        ps <- meta_list[[study]]
        otu <- as(otu_table(ps), "matrix")
        if (!taxa_are_rows(ps)) otu <- t(otu)
        meta <- data.frame(sample_data(ps))
        abund <- otu[taxon, ]
        means <- tapply(abund, meta$StandardTreatment, mean)
        
        reason <- ifelse(is.na(res_low_se$p_value[i]), 
                         "NA p-value (model failure)", 
                         "Low SE (<0.01, convergence issue)")
        
        row <- data.frame(
          Study = study,
          Taxon = taxon,
          ExclusionReason = reason,
          SeparationType = "Ultra-rare (near-perfect separation)",
          MeanAbundanceControl = means["Control"],
          MeanAbundanceRapamycin = means["Rapamycin"],
          StdError = res_low_se$std_error[i],
          Pvalue = res_low_se$p_value[i],
          stringsAsFactors = FALSE
        )
        all_rows[[length(all_rows) + 1]] <- row
      }
    }
  }
  
  if(length(all_rows) > 0) {
    df <- do.call(rbind, all_rows)
    df <- df %>% arrange(ExclusionReason, Study, Taxon)
    return(df)
  } else {
    return(data.frame())
  }
}

# Create comprehensive excluded taxa table
excluded_taxa_table <- make_excluded_taxa_table(
  treatment_only_results, 
  separation_analysis, 
  ps_genus_dose_sanitized,
  max_se = 5,
  min_se = 0.01
)

# Show table
table.excluded_taxa <- excluded_taxa_table %>%
  mutate(
    MeanAbundanceControl = round(MeanAbundanceControl, 2),
    MeanAbundanceRapamycin = round(MeanAbundanceRapamycin, 2),
    StdError = round(StdError, 4)
  )

kable(
  table.excluded_taxa,
  caption = "Taxa excluded from meta-analysis due to statistical issues"
)


# Print table to console
print(table.excluded_taxa)

# Summary statistics
cat("\n========== EXCLUSION SUMMARY ==========\n")
cat("Total taxa excluded:", nrow(table.excluded_taxa), "\n\n")

cat("By exclusion reason:\n")
print(table(excluded_taxa_table$ExclusionReason))

cat("\nBy study:\n")
print(table(excluded_taxa_table$Study))

cat("\nBy separation type:\n")
print(table(excluded_taxa_table$SeparationType))

# Summary by study and reason
cat("\nBreakdown by study and exclusion reason:\n")
print(table(excluded_taxa_table$Study, excluded_taxa_table$ExclusionReason))

# Which taxa are excluded in multiple studies?
taxa_excluded_multiple_studies <- excluded_taxa_table %>%
  group_by(Taxon, ExclusionReason) %>%
  summarise(
    n_studies = n_distinct(Study),
    studies = paste(sort(unique(Study)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_studies > 1) %>%
  arrange(desc(n_studies), Taxon)

cat("========== TAXA EXCLUDED IN MULTIPLE STUDIES ==========\n")
print(taxa_excluded_multiple_studies)

# Specifically check Akkermansia
cat("\n========== AKKERMANSIA ACROSS ALL STUDIES ==========\n")
excluded_taxa_table %>%
  filter(Taxon == "Akkermansia") %>%
  select(Study, ExclusionReason, SeparationType, MeanAbundanceControl, MeanAbundanceRapamycin)
#### Table of all filter Taxa ####
supp_table_excluded_detailed <- excluded_taxa_table %>%
  mutate(
    Study = case_when(
      Study == "Study2_Food" ~ "Bitto (2016) - Diet",
      Study == "Study2_IP" ~ "Bitto (2016) - IP",
      Study == "Study5" ~ "Yang (2022)",
      Study == "Study7" ~ "Han (2021)",
      TRUE ~ Study
    ),
    Genus = paste0("*", Taxon, "*"),
    `Exclusion Criterion` = case_when(
      grepl("High SE", ExclusionReason) ~ "SE > 5",
      grepl("Low SE", ExclusionReason) ~ "SE < 0.01",
      grepl("NA p-value", ExclusionReason) ~ "NA p-value"
    ),
    Pattern = case_when(
      grepl("Complete: Control Only", SeparationType) ~ "Control only",
      grepl("Control Sometimes", SeparationType) ~ "Mostly control",
      grepl("Rapa Sometimes", SeparationType) ~ "Mostly rapamycin",
      TRUE ~ "Ultra-rare"
    ),
    `Mean % Control` = sprintf("%.1f", MeanAbundanceControl),
    `Mean % Rapamycin` = sprintf("%.1f", MeanAbundanceRapamycin),
    SE = ifelse(StdError > 100, ">100", sprintf("%.2f", StdError))
  ) %>%
  select(Study, Genus, `Exclusion Criterion`, Pattern, 
         `Mean % Control`, `Mean % Rapamycin`, SE) %>%
  arrange(Study, Genus)


## Compare taxa counts before vs after SE filtering

## BEFORE: raw per-study results
taxa_counts_pre  <- table(unlist(lapply(treatment_only_results, `[[`, "taxon")))
multi_pre_2      <- names(taxa_counts_pre[taxa_counts_pre >= 2])
cat("Before filtering: taxa in ≥2 studies =", length(multi_pre_2), "\n")

## AFTER: SE/p-value filtered results
filtered_results <- filter_results$filtered_data
taxa_counts_post <- table(unlist(lapply(filtered_results, `[[`, "taxon")))
multi_post_2     <- names(taxa_counts_post[taxa_counts_post >= 2])
multi_post_3     <- names(taxa_counts_post[taxa_counts_post >= 3])

cat("After filtering: taxa in ≥2 studies  =", length(multi_post_2), "\n")
cat("After filtering: taxa in ≥3 studies  =", length(multi_post_3), "\n")

## Which taxa are in ≥2 but drop out when requiring ≥3?
setdiff(multi_post_2, multi_post_3)


#### META-ANALYSIS meta.taxa for filtered data ######
## 1. Use the properly filtered results
filtered_treatment_SE <- filtered_results  # ← Use the already-filtered data!

## 2. find taxa that appear in ≥ 3 studies  
keep <- unlist(lapply(filtered_treatment_SE, `[[`, "taxon")) |>
  table() |>
  {\(x) names(x[x >= 3])}()

## 3. drop taxa that don't meet that prevalence rule
filtered_treatment_SE <- lapply(
  filtered_treatment_SE,
  \(x) dplyr::filter(x, taxon %in% keep))

cat("Taxa retained for meta-analysis:\n")
cat("- Total unique taxa:", length(keep), "\n")
sapply(filtered_treatment_SE, nrow) |> print()

## helper that converts ONE study’s data frame 
long_one <- \(df, study)
df |>
  transmute(
    id  = taxon,
    study,
    Estimate.StandardTreatmentRapamycin        = effect,
    Std..Error.StandardTreatmentRapamycin      = std_error,
    `Pr(>|t|).StandardTreatmentRapamycin` = p_value,
    showvar = "StandardTreatmentRapamycin",
    l6      = taxon                             # optional genus label
  )

long_tab <- bind_rows(
  long_one(filtered_treatment_SE$Study2_Food, "Study2_Food"),
  long_one(filtered_treatment_SE$Study2_IP,   "Study2_IP"),
  long_one(filtered_treatment_SE$Study5,      "Study5"),
  long_one(filtered_treatment_SE$Study7,      "Study7")
)

long_tab_renamed <- long_tab |>
  dplyr::mutate(
    id = paste0("g__", id),
    study = case_when(
      study == "Study2_Food" ~ "Bitto (2016) – Food",
      study == "Study2_IP"   ~ "Bitto (2016) – IP",
      study == "Study5"      ~ "Yang (2022)",
      study == "Study7"      ~ "Han (2021)",
      TRUE ~ as.character(study)
    ),
    pop = "All"
  )

meta_treatment_only <- meta.taxa(
  taxcomdat        = long_tab_renamed,
  estimate.pattern = "Estimate\\.",          # literal  “Estimate.”
  se.pattern       = "Std\\.\\.Error\\.",    # literal  “Std..Error.”
  summary.measure  = "OR",                   # BEZI = log rate-ratio
  pool.var         = "id",
  studylab         = "study",
  backtransform    = TRUE,                   # ratios instead of logs
  percent.meta     = 0.67,                   # ≥ 2 / 3 studies
  p.adjust.method  = "fdr"
)

## -------- 3·3  inspect pooled output 
meta_treatment_only
table.pooled_treatment_only <- meta_treatment_only$random$StandardTreatmentRapamycin
cat("Number of taxa pooled:", nrow(table.pooled_treatment_only), "\n")
print(head(table.pooled_treatment_only[order(table.pooled_treatment_only$p.adjust), ]))

table.pooled_treatment_only <- table.pooled_treatment_only |>
  mutate(OR = exp(estimate),
         OR_ll = exp(ll),
         OR_ul = exp(ul))
table.pooled_treatment_only %>% arrange(p.adjust) %>%
  select(id, OR, OR_ll, OR_ul, p.adjust) %>% head()

# id        OR     OR_ll     OR_ul   p.adjust
# 1                 g__Muribaculum 0.5519312 0.3808801 0.7998004 0.03550347
# 2                 g__Bacteroides 1.5886137 1.1790790 2.1403940 0.03550347
# 3                g__Unclassified 0.7211150 0.5801086 0.8963956 0.03550347
# 4                g__Ruminococcus 0.6011158 0.4079689 0.8857052 0.08300883
# 5 g__[Eubacterium]_nodatum_group 1.4767044 1.0701332 2.0377427 0.11660641
# 6           g__Lachnoclostridium 1.3270689 1.0154749 1.7342743 0.21023923

## =====  Publication-ready meta table 

# --- 1.  Build a clean data frame 
tempdata <- table.pooled_treatment_only                                   # keep original safe
tempdata <- tempdata |>
  mutate(
    Taxon      = gsub(" ", "", id),     # drop spaces to match other tables
    OddsRatio = OR,                    # <<< CHANGED  (was RateRatio)
    CI_lower  = OR_ll,                 # <<< CHANGED
    CI_upper  = OR_ul,                 # <<< CHANGED
    q_value   = p.adjust               # FDR-adjusted P
  ) |>
  select(Taxon, OddsRatio, CI_lower, CI_upper, q_value) %>%
  arrange(q_value) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))    # nicer digits
table.Meta.Pooled <- tempdata

table.Meta.Pooled
#                     Taxon         RateRatio CI_lower CI_upper q_value
# 1                              A2    1.2201   0.6837   2.1774  0.7375
# 2                       Alistipes    0.9236   0.6861   1.2433  0.8363
# 3                    Anaeroplasma    0.6770   0.3735   1.2272  0.5130
# 4                     Bacteroides    1.5873   1.1780   2.1387  0.0419
# 5                 Bifidobacterium    0.7713   0.2079   2.8613  0.8363
# 6                         Blautia    0.7743   0.5247   1.1425  0.5130
# 7        Candidatus_Saccharimonas    0.9374   0.6610   1.3294  0.8363
# 8              Clostridia_UCG_014    0.7865   0.5350   1.1562  0.5175
# 9      Clostridia_vadinBB60_group    0.8949   0.5207   1.5382  0.8363
# 10               Colidextribacter    0.8290   0.5999   1.1456  0.5566
# 11                  Desulfovibrio    1.2307   0.7907   1.9157  0.6590
# 12                  Enterorhabdus    1.1669   0.7471   1.8225  0.7375
# 13       Eubacterium_brachy_group    0.9495   0.5619   1.6043  0.9555
# 14      Eubacterium_nodatum_group    1.4762   1.0697   2.0372  0.1246
# 15 Eubacterium_xylanophilum_group    0.6962   0.3578   1.3545  0.5566
# 16            Gastranaerophilales    0.7223   0.1877   2.7790  0.8363
# 17                  GCA_900066575    1.4051   0.5160   3.8259  0.7375
# 18                   Helicobacter    1.6293   0.4738   5.6029  0.7309
# 19                 Incertae_Sedis    1.2416   0.9354   1.6481  0.4621
# 20              Lachnoclostridium    1.3266   1.0152   1.7335  0.1922
# 21  Lachnospiraceae_NK4A136_group    0.6927   0.4256   1.1273  0.4621
# 22        Lachnospiraceae_UCG_006    0.8525   0.6400   1.1356  0.5566
# 23                  Lactobacillus    1.0216   0.4446   2.3473  0.9947
# 24                 Muribaculaceae    1.0007   0.8154   1.2280  0.9947
# 25                    Muribaculum    0.5514   0.3803   0.7994  0.0419
# 26                  NK4A214_group    1.1906   0.5312   2.6684  0.8363
# 27                  Oscillibacter    0.4706   0.1466   1.5104  0.5130
# 28                   Ruminococcus    0.6006   0.4071   0.8861  0.0892
# 29                     Tuzzerella    1.0209   0.7713   1.3513  0.9678
# 30                        UCG_009    0.6964   0.4280   1.1332  0.4621
# 31                        UCG_010    0.8575   0.5960   1.2337  0.7131
# 32                 Unclassified_1    1.5955   1.1545   2.2049  0.0543
# 33                 Unclassified_2    0.6623   0.4019   1.0913  0.4621
# 34                 Unclassified_3    0.9796   0.1810   5.3037  0.9947
# 35                   uncultured_5    0.6588   0.4459   0.9732  0.1922

# Create publication-ready table from sig_tab
# Use your existing table.pooled_treatment_only instead
# Fix the arrange/select order
table.sigtax_publication <- table.pooled_treatment_only %>%
  filter(p.adjust < 0.10) %>%  # Filter for significant results
  arrange(p.adjust) %>%        # Sort BEFORE selecting columns
  mutate(
    # Clean genus names
    Genus = gsub("g__", "", id),
    
    # Calculate rate ratios
    OR = OR,
    OR_lower = OR_ll, 
    OR_upper = OR_ul,
    
    # Format columns nicely
    `Log Odds Ratio` = sprintf("%.3f", estimate),
    `Log OR 95% CI` = sprintf("(%.3f, %.3f)", ll, ul),
    `Odds Ratio` = sprintf("%.2f", OR),
    `OR 95% CI` = sprintf("(%.2f, %.2f)", OR_lower, OR_upper),
    `FDR p-value` = sprintf("%.3f", p.adjust)
  ) %>%
  select(
    Genus, `Log Odds Ratio`, `Log OR 95% CI`,
    `Odds Ratio`, `OR 95% CI`, `FDR p-value`
  )
table.treatment_GAMLSS_result <- knitr::kable(table.sigtax_publication, 
             caption = "Table 2. Meta-analysis of rapamycin effects on gut microbiota composition")
table.treatment_GAMLSS_result


# 
# unsplitted
# |Genus        |Log Odds Ratio |Log OR 95% CI    |Odds Ratio |OR 95% CI    |FDR p-value |
#   |:------------|:--------------|:----------------|:----------|:------------|:-----------|
#   |Muribaculum  |-0.594         |(-0.965, -0.223) |0.55       |(0.38, 0.80) |0.036       |
#   |Bacteroides  |0.463          |(0.165, 0.761)   |1.59       |(1.18, 2.14) |0.036       |
#   |Unclassified |-0.327         |(-0.545, -0.109) |0.72       |(0.58, 0.90) |0.036       |
#   |Ruminococcus |-0.509         |(-0.897, -0.121) |0.60       |(0.41, 0.89) |0.083       |
#   > 
# splitted 
# |Genus                  |Log Odds Ratio |Log OR 95% CI    |Odds Ratio |OR 95% CI    |FDR p-value |
# |:----------------------|:--------------|:----------------|:----------|:------------|:-----------|
# |Candidatus_Arthromitus |0.719          |(0.652, 0.786)   |2.05       |(1.92, 2.19) |0.000       |
# |Muribaculum            |-0.655         |(-0.982, -0.328) |0.52       |(0.37, 0.72) |0.002       |
# |Bacteroides            |0.480          |(0.188, 0.772)   |1.62       |(1.21, 2.16) |0.024       |
# |Unclassified           |-0.296         |(-0.483, -0.108) |0.74       |(0.62, 0.90) |0.029       |
# |Ruminococcus           |-0.535         |(-0.921, -0.148) |0.59       |(0.40, 0.86) |0.075       |

#### Calculate Q² Heterogeneity Statistics ===== ####
#Citation for Q https://www.jstor.org/stable/3001666?origin=crossref
#Citation for τ^2 https://www.sciencedirect.com/science/article/pii/0197245686900462?via%3Dihub
#Citation for I^2 https://onlinelibrary.wiley.com/doi/10.1002/sim.1186
   
# Function to calculate Q statistics for a single taxon
calculate_Q_for_taxon <- function(taxon_name, study_results_list) {
  
  # Extract effects and standard errors for this taxon across all studies
  study_effects <- c()
  study_se <- c()
  study_names <- c()
  
  for(study_name in names(study_results_list)) {
    study_data <- study_results_list[[study_name]]
    
    # Find this taxon in the study results
    taxon_row <- which(study_data$taxon == taxon_name)
    
    if(length(taxon_row) > 0) {
      study_effects <- c(study_effects, study_data$effect[taxon_row])
      study_se <- c(study_se, study_data$std_error[taxon_row])
      study_names <- c(study_names, study_name)
    }
  }
  
  # Need at least 2 studies for heterogeneity calculation
  if(length(study_effects) < 2) {
    return(data.frame(
      taxon = taxon_name,
      n_studies = length(study_effects),
      Q = NA,
      Q_df = NA,
      Q_pval = NA,
      I2 = NA,
      tau2 = NA
    ))
  }
  
  # Calculate weights (inverse variance)
  weights <- 1 / (study_se^2)
  
  # Calculate pooled effect (weighted average)
  pooled_effect <- sum(weights * study_effects) / sum(weights)
  
  # Calculate Q statistic (Cochran's Q)
  Q <- sum(weights * (study_effects - pooled_effect)^2)
  
  # Degrees of freedom
  df <- length(study_effects) - 1
  
  # P-value for Q test
  Q_pval <- 1 - pchisq(Q, df)
  
  # I² statistic (percentage of variation due to heterogeneity)
  I2 <- max(0, (Q - df) / Q) * 100
  
  # Tau² (between-study variance)
  if(df > 0) {
    tau2 <- max(0, (Q - df) / (sum(weights) - sum(weights^2)/sum(weights)))
  } else {
    tau2 <- 0
  }
  
  return(data.frame(
    taxon = taxon_name,
    n_studies = length(study_effects),
    Q = round(Q, 3),
    Q_df = df,
    Q_pval = round(Q_pval, 4),
    I2 = round(I2, 1),
    tau2 = round(tau2, 4)
  ))
}

# Calculate Q statistics for all taxa in meta-analysis
meta_taxa <- unique(table.pooled_treatment_only$id)
meta_taxa <- gsub("g__", "", meta_taxa)  # Remove g__ prefix to match

cat("Calculating Q statistics for", length(meta_taxa), "taxa...\n")

# Calculate Q statistics for each taxon
heterogeneity_results <- map_dfr(meta_taxa, ~calculate_Q_for_taxon(.x, treatment_only_results))

# Display results
cat("\n===== HETEROGENEITY STATISTICS =====\n")
print(heterogeneity_results)

#### Enhanced Meta-Analysis Table with Heterogeneity ===== ####

# Add heterogeneity statistics to meta table (FIXED dplyr syntax)
table.Meta.WithHeterogeneity <- table.pooled_treatment_only %>%
  mutate(Taxon_clean = gsub("g__", "", id)) %>%
  left_join(heterogeneity_results, by = c("Taxon_clean" = "taxon")) %>%
  mutate(
    Log_OR     = round(estimate, 3),             # <<< CHANGED
    OddsRatio  = round(OR, 4),                   # <<< CHANGED
    Pct_Change = round(100 * (OR - 1), 2),       # <<< CHANGED
    Pct_CI_lower = round(100 * (OR_ll - 1), 2),  # <<< CHANGED
    Pct_CI_upper = round(100 * (OR_ul - 1), 2),  # <<< CHANGED
    
    # Statistical significance
    FDR_p = round(p.adjust, 4),
    
    # Significance indicators
    Sig = case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01 ~ "**", 
      p.adjust < 0.05 ~ "*",
      p.adjust < 0.1 ~ "†",
      TRUE ~ "ns"
    ),
    
    # Heterogeneity interpretation
    Heterogeneity_Level = case_when(
      is.na(I2) ~ "Cannot calculate",
      I2 < 25 ~ "Low",
      I2 < 50 ~ "Moderate", 
      I2 < 75 ~ "High",
      TRUE ~ "Very High"
    )
  ) %>%
  select(Taxon_clean, Log_OR, OddsRatio, Pct_Change, Pct_CI_lower, Pct_CI_upper,
         n_studies, Q, Q_df, Q_pval, I2, tau2, Heterogeneity_Level, FDR_p, Sig) %>%
  arrange(FDR_p)

# Display enhanced table
cat("\n===== META-ANALYSIS RESULTS WITH HETEROGENEITY =====\n")
print(table.Meta.WithHeterogeneity)

#Heterogeneity Summary Statistics 

heterogeneity_summary <- heterogeneity_results %>%
  filter(!is.na(I2)) %>%
  summarise(
    n_taxa = n(),
    median_I2 = median(I2, na.rm = TRUE),
    mean_I2 = mean(I2, na.rm = TRUE),
    high_heterogeneity = sum(I2 > 50, na.rm = TRUE),
    very_high_heterogeneity = sum(I2 > 75, na.rm = TRUE),
    significant_Q = sum(Q_pval < 0.05, na.rm = TRUE)
  )

cat("\n===== HETEROGENEITY SUMMARY =====\n")
print(heterogeneity_summary)

#Publication-Ready Table for Significant Results 

table.sigtax_with_heterogeneity <- table.pooled_treatment_only %>%
  filter(p.adjust < 0.10) %>%
  arrange(p.adjust) %>%
  # Join with heterogeneity results
  mutate(Taxon_clean = gsub("g__", "", id)) %>%
  left_join(heterogeneity_results, by = c("Taxon_clean" = "taxon")) %>%
  mutate(
    # Clean genus names
    Genus = gsub("g__", "", id),
    
    # Calculate rate ratios
    OR         = OR,
    OR_lower   = OR_ll,
    OR_upper   = OR_ul,
    Pct_Change = 100 * (OR - 1),
    # Format columns for publication
    `Log Odds Ratio` = sprintf("%.3f", estimate),
    `Log OR 95% CI` = sprintf("(%.3f, %.3f)", ll, ul),
    `Odds Ratio` = sprintf("%.2f", OR),
    `OR 95% CI` = sprintf("(%.2f, %.2f)", OR_lower, OR_upper),
    `% Change` = sprintf("%.1f%%", Pct_Change),
    `I² (%)` = ifelse(is.na(I2), "N/A", sprintf("%.1f%%", I2)),
    `Q p-value` = ifelse(is.na(Q_pval), "N/A", sprintf("%.3f", Q_pval)),
    `FDR p-value` = sprintf("%.3f", p.adjust),
    
    # Heterogeneity interpretation
    `Heterogeneity` = case_when(
      is.na(I2) ~ "N/A",
      I2 < 25 ~ "Low",
      I2 < 50 ~ "Moderate",
      I2 < 75 ~ "High", 
      TRUE ~ "Very High"
    )
  ) %>%
  select(Genus, `Log Odds Ratio`, `Log OR 95% CI`,
         `Odds Ratio`, `OR 95% CI`, `% Change`,
         `I² (%)`, Heterogeneity, `Q p-value`, `FDR p-value`)

cat("\n===== PUBLICATION TABLE WITH HETEROGENEITY =====\n")
table.Metaanalysis_rapamycineffects_heterogeneity <- knitr::kable(
  table.sigtax_with_heterogeneity, 
  caption = "Table: Meta-analysis of rapamycin effects on gut microbiome with heterogeneity assessment"
)

print(table.Metaanalysis_rapamycineffects_heterogeneity)

#Key Findings Summary

cat("\n===== KEY HETEROGENEITY FINDINGS =====\n")
cat("• Total taxa analyzed:", nrow(heterogeneity_results), "\n")
cat("• Mean I² across all taxa:", round(heterogeneity_summary$mean_I2, 1), "%\n")
cat("• Taxa with high heterogeneity (I² > 50%):", heterogeneity_summary$high_heterogeneity, "\n")
cat("• Taxa with very high heterogeneity (I² > 75%):", heterogeneity_summary$very_high_heterogeneity, "\n")
cat("• Taxa with significant Q-test (p < 0.05):", heterogeneity_summary$significant_Q, "\n")
cat("• Significant taxa after FDR correction:", sum(table.pooled_treatment_only$p.adjust < 0.10), "\n")

#### META-ANALYSIS: IP Studies Only (excluding dietary intervention) ######

## Create long format table with ONLY IP studies
long_tab_IP_only <- bind_rows(
  long_one(filtered_treatment_SE$Study2_IP,   "Study2_IP"),
  long_one(filtered_treatment_SE$Study5,      "Study5"),
  long_one(filtered_treatment_SE$Study7,      "Study7")
  # Note: Study2_Food is excluded
)

long_tab_IP_renamed <- long_tab_IP_only |>
  dplyr::mutate(
    id = paste0("g__", id),
    study = case_when(
      study == "Study2_IP"   ~ "Bitto (2016) – IP",
      study == "Study5"      ~ "Yang (2022)",
      study == "Study7"      ~ "Han (2021)",
      TRUE ~ as.character(study)
    ),
    pop = "All"
  )

## Run meta-analysis on IP studies only
meta_treatment_IP_only <- meta.taxa(
  taxcomdat        = long_tab_IP_renamed,
  estimate.pattern = "Estimate\\.",
  se.pattern       = "Std\\.\\.Error\\.",
  summary.measure  = "OR",
  pool.var         = "id",
  studylab         = "study",
  backtransform    = TRUE,
  percent.meta     = 0.66,  # ≥ 2/3 studies (so ≥2 out of 3)
  p.adjust.method  = "fdr"
)

## Inspect IP-only results
meta_treatment_IP_only
table.pooled_IP_only <- meta_treatment_IP_only$random$StandardTreatmentRapamycin
cat("Number of taxa pooled (IP only):", nrow(table.pooled_IP_only), "\n")

## Calculate odds ratios
table.pooled_IP_only <- table.pooled_IP_only |>
  mutate(OR = exp(estimate),
         OR_ll = exp(ll),
         OR_ul = exp(ul))

## View top results
table.pooled_IP_only %>% 
  arrange(p.adjust) %>%
  select(id, OR, OR_ll, OR_ul, p.adjust) %>% 
  head(10)

## Create publication-ready table for IP-only analysis
table.sigtax_IP_only <- table.pooled_IP_only %>%
  filter(p.adjust < 0.10) %>%
  arrange(p.adjust) %>%
  mutate(
    Genus = gsub("g__", "", id),
    OR = OR,
    OR_lower = OR_ll, 
    OR_upper = OR_ul,
    `Log Odds Ratio` = sprintf("%.3f", estimate),
    `Log OR 95% CI` = sprintf("(%.3f, %.3f)", ll, ul),
    `Odds Ratio` = sprintf("%.2f", OR),
    `OR 95% CI` = sprintf("(%.2f, %.2f)", OR_lower, OR_upper),
    `FDR p-value` = sprintf("%.3f", p.adjust)
  ) %>%
  select(
    Genus, `Log Odds Ratio`, `Log OR 95% CI`,
    `Odds Ratio`, `OR 95% CI`, `FDR p-value`
  )

table.IP_only_result <- knitr::kable(
  table.sigtax_IP_only, 
  caption = "Table 3. Meta-analysis of rapamycin effects (IP administration only)"
)
table.IP_only_result

## Calculate heterogeneity for IP-only results
heterogeneity_IP_only <- map_dfr(
  unique(table.pooled_IP_only$id),
  ~calculate_Q_for_taxon(gsub("g__", "", .x), 
                         list(Study2_IP = filtered_results$Study2_IP,
                              Study5 = filtered_results$Study5,
                              Study7 = filtered_results$Study7))
)

## Add to table
table.sigtax_IP_heterogeneity <- table.pooled_IP_only %>%
  filter(p.adjust < 0.10) %>%
  mutate(Taxon_clean = gsub("g__", "", id)) %>%
  left_join(heterogeneity_IP_only, by = c("Taxon_clean" = "taxon")) %>%
  mutate(
    Genus = Taxon_clean,
    `Log Odds Ratio` = sprintf("%.3f", estimate),
    `Log OR 95% CI` = sprintf("(%.3f, %.3f)", ll, ul),
    `Odds Ratio` = sprintf("%.2f", OR),
    `OR 95% CI` = sprintf("(%.2f, %.2f)", OR_ll, OR_ul),
    `I² (%)` = sprintf("%.1f%%", I2),
    `Q p-value` = sprintf("%.3f", Q_pval),
    `FDR p-value` = sprintf("%.3f", p.adjust)
  ) %>%
  select(Genus, `Log Odds Ratio`, `Log OR 95% CI`,
         `Odds Ratio`, `OR 95% CI`, `I² (%)`, `Q p-value`, `FDR p-value`)

kable(table.sigtax_IP_heterogeneity, 
      caption = "IP-only meta-analysis with heterogeneity")

# Table: IP-only meta-analysis with heterogeneity
# 
# |Genus                  |Log Odds Ratio |Log OR 95% CI    |Odds Ratio |OR 95% CI    |I² (%) |Q p-value |FDR p-value |
# |:----------------------|:--------------|:----------------|:----------|:------------|:------|:---------|:-----------|
# |Candidatus_Arthromitus |0.718          |(0.651, 0.785)   |2.05       |(1.92, 2.19) |0.0%   |0.561     |0.000       |
# |Muribaculum            |-0.530         |(-0.840, -0.221) |0.59       |(0.43, 0.80) |47.1%  |0.151     |0.024       |

#### Visualisation####

# Before renaming, check what the study column is called and what values it has
cat("Column names in long_tab:\n")
print(names(long_tab_renamed))

cat("\nUnique study values:\n")
print(unique(long_tab_renamed$study))

## [A] relabel split studies (just before plotting)
##     Study2_Food -> "Bitto 2016 - Diet", Study2_IP -> "Bitto 2016 - IP"
long_tab_renamed <- long_tab_renamed %>%
  mutate(study = case_when(
    study == "Bitto (2016) – Food" ~ "Bitto (2016) - Diet",
    study == "Bitto (2016) – IP" ~ "Bitto (2016) - IP",
    TRUE ~ study
  ))

# Verify the change
cat("Updated study values:\n")
print(unique(long_tab_renamed$study))

#1 · add the prefix --
# Keep the list structure but filter the internal data frame
meta_random_pref <- meta_treatment_only$random
meta_random_pref$StandardTreatmentRapamycin <- meta_random_pref$StandardTreatmentRapamycin %>%
  filter(!grepl("Unclassified", id))

# Remove Unclassified_1 from long_tab data
long_tab_pref <- long_tab_renamed %>%
  filter(!grepl("Unclassified", id))

# or "IP"/"Food" if you have that info
head(long_tab_pref)

## 2 · now the plots finally work 
Heatmap.metanalysis_treatment_only <- metatab.show(
  metatab        = meta_random_pref,
  com.pooled.tab = long_tab_pref,
  
  ## tell it exactly how the column names look
  estimate.pattern = "Estimate\\.",
  se.pattern       = "Std\\.\\.Error\\.",      # TWO dots, escaped
  p.pattern        = "Pr\\(>\\|t\\|\\)",       # escape (), | characters
  
  showvar        = "StandardTreatmentRapamycin",
  tax.lev        = "l6",
  p.cutoff.type  = "p.adjust",
  p.cutoff       = 0.05,
  display        = "plot",
  plot           = "heatmap",
  fill.value     = "log(OR)"
) 
# render + add tiny pooled-model label (keeps everything else unchanged)
grid::grid.newpage()
grid::grid.draw(Heatmap.metanalysis_treatment_only)
grid::grid.text("Pooled (RE)",
                x = grid::unit(0.985, "npc"), y = grid::unit(0.985, "npc"),
                just = c("right","top"),
                gp = grid::gpar(fontface = "bold"))

Heatmap.metanalysis_treatment_only 


sig_tab <- metatab.show(
  metatab = meta_random_pref,
  com.pooled.tab = long_tab_pref,
  sumvar = "taxa",               # Explicitly specify taxa
  highest.lev = "g",             # Genus is highest level
  tax.lev = "l6",                # We're using genus level (l6)
  estimate.pattern = "Estimate\\.",
  se.pattern = "Std\\.\\.Error\\.",
  p.pattern = "Pr\\(>\\|t\\|\\)",
  showvar = "StandardTreatmentRapamycin",
  p.cutoff.type = "p.adjust",
  p.cutoff = 0.05,
  readjust.p = FALSE,
  display = "table",
  digit = 3,
  p.digit = 4
)
knitr::kable(sig_tab, caption = "Genera with FDR < 0.05")

# undplitted
# |   |id              | estimate|     ll|     ul|      p| p.adjust|
#   |:--|:---------------|--------:|------:|------:|------:|--------:|
#   |3  |g__Unclassified |   -0.327| -0.545| -0.109| 0.0032|   0.0355|
#   |1  |g__Muribaculum  |   -0.594| -0.965| -0.223| 0.0017|   0.0355|
#   |2  |g__Bacteroides  |    0.463|  0.165|  0.761| 0.0023|   0.0355|
#   |4  |g__Ruminococcus |   -0.509| -0.897| -0.121| 0.0101|   0.0830|

# splitted 
# Table: Genera with FDR < 0.05
# 
# |id                        | estimate|     ll|     ul|      p| p.adjust|
  # |:-------------------------|--------:|------:|------:|------:|--------:|
  # |g__Candidatus_Arthromitus |    0.719|  0.652|  0.786| 0.0000|   0.0000|
  # |g__Muribaculum            |   -0.655| -0.982| -0.328| 0.0001|   0.0025|
  # |g__Bacteroides            |    0.480|  0.188|  0.772| 0.0013|   0.0235|

###Forest Plot ###

# 1. First, make sure you have the significant taxa from sig_tab
taxa_of_interest <- gsub("g__", "", unique(sig_tab$id))

# Create a function to build the forest plot using sig_tab and long_tab_pref
create_forest_from_sigtab <- function(sig_tab_data, long_data, taxa_list, heterogeneity_data = NULL) {
  # Create a data frame for all plot data
  all_plot_data <- data.frame()
  table_data <- data.frame()
  
  # For each taxon, add meta-analysis and study-specific data
  for(taxon in taxa_list) {
    # Get meta-analysis result from sig_tab
    meta_result <- sig_tab_data %>%
      filter(gsub("g__", "", id) == taxon) %>%
      mutate(
        taxon_name = taxon,
        yi = estimate,
        ci.lb = ll,
        ci.ub = ul
      )
    
    # Get study-specific data from long_tab_pref
    study_data <- long_data %>%
      filter(grepl(paste0(taxon, "$"), id)) %>%
      mutate(
        taxon_name = taxon,
        yi = Estimate.StandardTreatmentRapamycin,
        sei = Std..Error.StandardTreatmentRapamycin,
        ci.lb = yi - 1.96 * sei,
        ci.ub = yi + 1.96 * sei
      )
    
    # Create data frames for plotting
    meta_df <- data.frame(
      study = "Pooled (RE)",
      yi = meta_result$yi,
      ci.lb = meta_result$ci.lb,
      ci.ub = meta_result$ci.ub,
      sei        = NA_real_,
      weight     = NA_real_, 
      is_summary = TRUE,
      taxon_name = taxon
    )
    
    study_df <- data.frame(
      study = study_data$study,
      yi = study_data$yi,
      ci.lb = study_data$ci.lb,
      ci.ub = study_data$ci.ub,
      sei        = study_data$sei,
      weight     = 1 / study_data$sei^2,
      is_summary = FALSE,
      taxon_name = taxon
    )
    
    # Add to all plot data
    taxon_data <- rbind(meta_df, study_df)
    all_plot_data <- rbind(all_plot_data, taxon_data)
    
    # Create table data
    meta_row <- data.frame(
      Taxon = taxon,
      Study = "Pooled (RE)",
      LogOR = meta_result$yi,
      CI_lower = meta_result$ci.lb,
      CI_upper = meta_result$ci.ub,
      OR = exp(meta_result$yi),
      OR_lower = exp(meta_result$ci.lb),
      OR_upper = exp(meta_result$ci.ub),
      P_value = meta_result$p,
      Weight    = NA_real_
    )
    
    study_rows <- data.frame(
      Taxon = taxon,
      Study = study_data$study,
      LogOR = study_data$yi,
      CI_lower = study_data$ci.lb,
      CI_upper = study_data$ci.ub,
      OR = exp(study_data$yi),
      OR_lower = exp(study_data$ci.lb),
      OR_upper = exp(study_data$ci.ub),
      P_value = 2 * (1 - pnorm(abs(study_data$yi / study_data$sei))),
      Weight    = 1 / study_data$sei^2
    )
    
    # Add to table data
    table_data <- rbind(table_data, meta_row, study_rows)
  }
  
  # RECODE study names to desired format BEFORE setting factor levels
  all_plot_data$study <- recode(all_plot_data$study,
                                "Bitto (2016) - Diet" = "Bitto 2016 (Diet)",
                                "Bitto (2016) - IP" = "Bitto 2016 (IP)",
                                "Yang (2022)" = "Yang 2022",
                                "Han (2021)" = "Han 2021")
  
  table_data$Study <- recode(table_data$Study,
                             "Bitto (2016) - Diet" = "Bitto 2016 (Diet)",
                             "Bitto (2016) - IP" = "Bitto 2016 (IP)",
                             "Yang (2022)" = "Yang 2022",
                             "Han (2021)" = "Han 2021")
  
  # NOW set factor levels with the RECODED names (reversed for plotting)
  all_plot_data$study <- factor(all_plot_data$study, 
                                levels = rev(c("Pooled (RE)",
                                               "Bitto 2016 (Diet)", 
                                               "Bitto 2016 (IP)", 
                                               "Yang 2022", 
                                               "Han 2021")))
  
   ## Global x-range (add room on the right for text columns)
  all_ci_bounds <- c(all_plot_data$ci.lb, all_plot_data$ci.ub)
  x_min <- floor(min(all_ci_bounds, na.rm = TRUE) * 2) / 2
  x_max <- ceiling(max(all_ci_bounds, na.rm = TRUE) * 2) / 2
  x_min <- max(x_min, -1.5); x_max <- min(x_max,  1.5)

  ci_col_x   <- x_max + 0.10    # FIX 2: place text just outside panel range
  wgt_col_x  <- x_max + 0.55

  # Create plots for each taxon
  forest_plots <- list()
  
  for(i in seq_along(taxa_list)) {
    taxon <- taxa_list[i]
    is_last <- (i == length(taxa_list))
    
    plot_data <- all_plot_data %>% filter(taxon_name == taxon)
    
    # Create CI text for display on right side
    plot_data$ci_text <- sprintf("%.2f [%.2f; %.2f]", 
                                 plot_data$yi, 
                                 plot_data$ci.lb, 
                                 plot_data$ci.ub)
    
    # Calculate weights as percentages for display
    if(sum(!plot_data$is_summary) > 0) {
      study_weights <- plot_data$weight[!plot_data$is_summary]
      total_weight <- sum(study_weights, na.rm = TRUE)
      weight_pct <- round(100 * study_weights / total_weight, 1)
      
      plot_data$weight_text <- ifelse(plot_data$is_summary, 
                                      "—", 
                                      paste0(weight_pct, "%"))
    } else {
      plot_data$weight_text <- ""
    }
    
    # GET HETEROGENEITY DATA
    het_text <- ""
    if(!is.null(heterogeneity_data)) {
      het_result <- heterogeneity_data %>% filter(taxon == !!taxon)
      if(nrow(het_result) > 0) {
        i2 <- round(het_result$I2, 1)
        tau2 <- round(het_result$tau2, 4)
        
        het_text <- paste0("Heterogeneity: I² = ", i2, "%, τ² = ", tau2)
      }
    }
    
    p <- ggplot(plot_data, aes(x = yi, y = study, xmin = ci.lb, xmax = ci.ub)) +
      geom_vline(xintercept = seq(-1.5, 1.5, by = 0.5), color = "gray90") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
      geom_pointrange(aes(size = is_summary, color = is_summary), linewidth = 0.8) +
      geom_text(aes(label = ci_text), x = 1.9, hjust = 0, size = 3) +
      geom_text(aes(label = weight_text), x = 2.6, hjust = 0, size = 3) +
      scale_x_continuous(
        limits = c(-1.8, 3.4),
        breaks = seq(-1.5, 1.5, by = 0.5),
        minor_breaks = NULL
      ) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x = if(is_last) element_text(size = 10, face = "bold", margin = margin(t = 10)) else element_blank(),
        axis.title.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        plot.margin = margin(t = -10, r = -5, b = -5, l = -5)
      ) +
      scale_color_manual(values = c("FALSE" = "royalblue3", "TRUE" = "darkred")) +
      scale_size_manual(values = c("FALSE" = 2, "TRUE" = 2)) +
      labs(
        title = paste0(taxon),
        x = if(is_last) "log(Odds Ratio)" else NULL 
      ) +
      annotate("text", x = 1.9, y = Inf, vjust = 1.5,
               label = "Effect [95% CI]", hjust = 0, size = 3, fontface = "bold") +
      annotate("text", x = 2.6, y = Inf, vjust = 1.5,
               label = "Weight", hjust = 0, size = 3, fontface = "bold")
    
    # ADD HETEROGENEITY TEXT
    if(het_text != "") {
      p <- p + annotate("text", x = -1.70, y = 0.6, 
                        label = het_text, 
                        hjust = 0, size = 3, color = "gray40", fontface = "italic")
    }
    
    forest_plots[[taxon]] <- p
  }
  
  # Combine plots
  if(require(patchwork) && length(forest_plots) > 0) {
    combined_plot <- forest_plots[[taxa_list[1]]]
    
    if(length(forest_plots) > 1) {
      for(i in 2:length(taxa_list)) {
        combined_plot <- combined_plot / forest_plots[[taxa_list[i]]]
      }
    }
    
    combined_plot <- combined_plot +
      plot_annotation(
        title = "Effect of Rapamycin on Gut Microbiota",
        subtitle = "Meta-analysis of three independent studies (FDR < 0.05)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 12, hjust = 0.5))
      )
  }
  
  # Format table data with RECODED names
  table_data$Study <- factor(table_data$Study, 
                             levels = c("Pooled (RE)",
                                        "Bitto 2016 (Diet)", 
                                        "Bitto 2016 (IP)", 
                                        "Yang 2022", 
                                        "Han 2021"))
  
  table_data <- table_data[order(table_data$Taxon, table_data$Study), ]
  
  formatted_table <- table_data %>%
    mutate(
      `Odds Ratio [95% CI]` = sprintf("%.2f [%.2f-%.2f]", OR, OR_lower, OR_upper),
      `Log OR [95% CI]` = sprintf("%.2f [%.2f-%.2f]", LogOR, CI_lower, CI_upper),
      `P-value` = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value)),
      `Weight (%)` = ifelse(is.na(Weight),
                            "—",
                            sprintf("%.1f%%",
                                    100 * Weight / tapply(Weight, Taxon, sum, na.rm = TRUE)[Taxon])) 
    ) %>%
    select(Taxon, Study, `Log OR [95% CI]`, `Odds Ratio [95% CI]`, `P-value`, `Weight (%)`)
  
  # Return plots and table
  return(list(
    individual = forest_plots,
    combined = combined_plot,
    table = formatted_table
  ))
}

# Run with your desired taxa order
results <- create_forest_from_sigtab(
  sig_tab_data = sig_tab,
  long_data = long_tab_pref,
  taxa_list = c("Muribaculum", "Bacteroides", "Candidatus_Arthromitus"),
  heterogeneity_data = heterogeneity_results
)

# Get the forest plot
plot.MetaForest_combined <- results$combined

# Get the matching table
table.Forestplot_combined <- results$table

# Display
print(plot.MetaForest_combined)
print(table.Forestplot_combined)

head(heterogeneity_results)
head(heterogeneity_IP_only)

###Forest Plot with IP-only Pooled Estimate - CLEAN VERSION ###

create_forest_from_sigtab <- function(sig_tab_data, long_data, taxa_list, 
                                      heterogeneity_data_all = NULL,  # ← CHANGED from heterogeneity_data
                                      heterogeneity_data_ip = NULL, 
                                      ip_pooled_data = NULL) {
  
  all_plot_data <- data.frame()
  table_data <- data.frame()
  
  # For each taxon, add meta-analysis and study-specific data
  for(taxon in taxa_list) {
    # Get meta-analysis result from sig_tab (ALL STUDIES)
    meta_result <- sig_tab_data %>%
      filter(gsub("g__", "", id) == taxon) %>%
      mutate(
        taxon_name = taxon,
        yi = estimate,
        ci.lb = ll,
        ci.ub = ul
      )
    het_all <- data.frame(I2 = NA, tau2 = NA)
    if(!is.null(heterogeneity_data_all)) {
      het_result <- heterogeneity_data_all %>% filter(taxon == !!taxon)
      if(nrow(het_result) > 0) {
        het_all <- het_result
      }
    }
    
    # Get heterogeneity for "IP only" analysis
    het_ip <- data.frame(I2 = NA, tau2 = NA)
    if(!is.null(heterogeneity_data_ip)) {
      het_result_ip <- heterogeneity_data_ip %>% filter(taxon == !!taxon)
      if(nrow(het_result_ip) > 0) {
        het_ip <- het_result_ip
      }
    }
    
    # Get study-specific data from long_tab_pref
    study_data <- long_data %>%
      filter(grepl(paste0(taxon, "$"), id)) %>%
      mutate(
        taxon_name = taxon,
        yi = Estimate.StandardTreatmentRapamycin,
        sei = Std..Error.StandardTreatmentRapamycin,
        ci.lb = yi - 1.96 * sei,
        ci.ub = yi + 1.96 * sei
      )
    
    # Create data frame for ALL STUDIES pooled estimate
    meta_df <- data.frame(
      study = "Pooled (All)",
      yi = meta_result$yi,
      ci.lb = meta_result$ci.lb,
      ci.ub = meta_result$ci.ub,
      sei = NA_real_,
      weight = NA_real_, 
      is_summary = TRUE,
      summary_type = "all",
      taxon_name = taxon
    )
    
    # Create data frame for IP-ONLY pooled estimate
    ip_meta_df <- NULL
    if(!is.null(ip_pooled_data)) {
      ip_result <- ip_pooled_data %>%
        filter(gsub("g__", "", id) == taxon)
      
      if(nrow(ip_result) > 0) {
        ip_meta_df <- data.frame(
          study = "Pooled (IP only)",
          yi = ip_result$estimate,
          ci.lb = ip_result$ll,
          ci.ub = ip_result$ul,
          sei = NA_real_,
          weight = NA_real_,
          is_summary = TRUE,
          summary_type = "ip",
          taxon_name = taxon
        )
      }
    }
    
    # Create data frame for individual studies
    study_df <- data.frame(
      study = study_data$study,
      yi = study_data$yi,
      ci.lb = study_data$ci.lb,
      ci.ub = study_data$ci.ub,
      sei = study_data$sei,
      weight = 1 / study_data$sei^2,
      is_summary = FALSE,
      summary_type = "individual",
      taxon_name = taxon
    )
    
    # Combine for plot data
    taxon_data <- rbind(meta_df, study_df)
    if(!is.null(ip_meta_df)) {
      taxon_data <- rbind(taxon_data, ip_meta_df)
    }
    all_plot_data <- rbind(all_plot_data, taxon_data)
    
    # Create table data
    meta_row <- data.frame(
      Taxon = taxon,
      Study = "Pooled (All)",
      LogOR = meta_result$yi,
      CI_lower = meta_result$ci.lb,
      CI_upper = meta_result$ci.ub,
      OR = exp(meta_result$yi),
      OR_lower = exp(meta_result$ci.lb),
      OR_upper = exp(meta_result$ci.ub),
      P_value = meta_result$p,
      Weight = NA_real_,
      I2 = het_all$I2,      # ← ADD
      tau2 = het_all$tau2   # ← ADD
    )
    
    study_rows <- data.frame(
      Taxon = taxon,
      Study = study_data$study,
      LogOR = study_data$yi,
      CI_lower = study_data$ci.lb,
      CI_upper = study_data$ci.ub,
      OR = exp(study_data$yi),
      OR_lower = exp(study_data$ci.lb),
      OR_upper = exp(study_data$ci.ub),
      P_value = 2 * (1 - pnorm(abs(study_data$yi / study_data$sei))),
      Weight = 1 / study_data$sei^2,
      I2 = NA,     # ← ADD
      tau2 = NA    # ← ADD
    )
    
    # Add IP pooled row to table
    if(!is.null(ip_meta_df)) {
      ip_meta_row <- data.frame(
        Taxon = taxon,
        Study = "Pooled (IP only)",
        LogOR = ip_result$estimate,
        CI_lower = ip_result$ll,
        CI_upper = ip_result$ul,
        OR = exp(ip_result$estimate),
        OR_lower = exp(ip_result$ll),
        OR_upper = exp(ip_result$ul),
        P_value = ip_result$p,
        Weight = NA_real_,
        I2 = het_ip$I2,      # ← ADD
        tau2 = het_ip$tau2   # ← ADD
      )
      table_data <- rbind(table_data, meta_row, study_rows, ip_meta_row)
    } else {
      table_data <- rbind(table_data, meta_row, study_rows)
    }
  }
  
  # RECODE study names - SEPARATE, COMPLETE STATEMENTS
  all_plot_data$study <- dplyr::recode(all_plot_data$study,
                                       "Bitto (2016) - Diet" = "Bitto 2016 (Diet)",
                                       "Bitto (2016) - IP" = "Bitto 2016 (IP)",
                                       "Yang (2022)" = "Yang 2022",
                                       "Han (2021)" = "Han 2021")
  
  table_data$Study <- dplyr::recode(table_data$Study,
                                    "Bitto (2016) - Diet" = "Bitto 2016 (Diet)",
                                    "Bitto (2016) - IP" = "Bitto 2016 (IP)",
                                    "Yang (2022)" = "Yang 2022",
                                    "Han (2021)" = "Han 2021")
  
  # Set factor levels (reversed for plotting)
  all_plot_data$study <- factor(all_plot_data$study, 
                                levels = rev(c("Pooled (All)",
                                               "Bitto 2016 (Diet)",
                                               "Pooled (IP only)",
                                               "Bitto 2016 (IP)", 
                                               "Yang 2022", 
                                               "Han 2021")))
  
  ## Global x-range
  all_ci_bounds <- c(all_plot_data$ci.lb, all_plot_data$ci.ub)
  x_min <- floor(min(all_ci_bounds, na.rm = TRUE) * 2) / 2
  x_max <- ceiling(max(all_ci_bounds, na.rm = TRUE) * 2) / 2
  x_min <- max(x_min, -1.5)
  x_max <- min(x_max, 1.5)
  
  # Create plots for each taxon
  forest_plots <- list()
  
  for(i in seq_along(taxa_list)) {
    taxon <- taxa_list[i]
    is_last <- (i == length(taxa_list))
    
    plot_data <- all_plot_data %>% filter(taxon_name == taxon)
    
    # Create CI text
    plot_data$ci_text <- sprintf("%.2f [%.2f; %.2f]", 
                                 plot_data$yi, 
                                 plot_data$ci.lb, 
                                 plot_data$ci.ub)
    
    # Calculate weights (only for individual studies)
    if(sum(!plot_data$is_summary) > 0) {
      study_weights <- plot_data$weight[!plot_data$is_summary]
      total_weight <- sum(study_weights, na.rm = TRUE)
      weight_pct <- round(100 * study_weights / total_weight, 1)
      
      # Initialize weight_text column
      plot_data$weight_text <- ""
      
      # Assign "—" to pooled estimates
      plot_data$weight_text[plot_data$is_summary] <- "—"
      
      # Assign percentages to individual studies (in correct order!)
      plot_data$weight_text[!plot_data$is_summary] <- paste0(weight_pct, "%")
      
    } else {
      plot_data$weight_text <- ""
    }
    
    # Get heterogeneity data
    het_text <- ""
    if(!is.null(heterogeneity_data_all)) {
      het_result <- heterogeneity_data_all %>% filter(taxon == !!taxon)
      if(nrow(het_result) > 0) {
        i2 <- round(het_result$I2, 1)
        tau2 <- round(het_result$tau2, 4)
        het_text <- paste0("Heterogeneity (All): I² = ", i2, "%, τ² = ", tau2)
      }
    }
    
    p <- ggplot(plot_data, aes(x = yi, y = study, xmin = ci.lb, xmax = ci.ub)) +
      geom_vline(xintercept = seq(-1.5, 1.5, by = 0.5), color = "gray90") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
      geom_pointrange(aes(size = is_summary, color = is_summary), linewidth = 0.8) +
      geom_text(aes(label = ci_text), x = 1.9, hjust = 0, size = 3) +
      geom_text(aes(label = weight_text), x = 2.6, hjust = 0, size = 3) +
      scale_x_continuous(
        limits = c(-1.8, 3.4),
        breaks = seq(-1.5, 1.5, by = 0.5),
        minor_breaks = NULL
      ) +
      scale_y_discrete(expand = expansion(add = c(0.8, 0.7))) +  # ADD THIS LINE - expands y-axis at bottom
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x = if(is_last) element_text(size = 10, face = "bold", margin = margin(t = 10)) else element_blank(),
        axis.title.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        plot.margin = margin(t = -10, r = -5, b = -5, l = -5)
      ) +
      scale_color_manual(values = c("FALSE" = "royalblue3", "TRUE" = "darkred")) +
      scale_size_manual(values = c("FALSE" = 2, "TRUE" = 2)) +
      labs(
        title = paste0(taxon),
        x = if(is_last) "log(Odds Ratio)" else NULL 
      ) +
      annotate("text", x = 1.9, y = Inf, vjust = 1.1,
               label = "Effect [95% CI]", hjust = 0, size = 3, fontface = "bold") +
      annotate("text", x = 2.6, y = Inf, vjust = 1.1,
               label = "Weight", hjust = 0, size = 3, fontface = "bold")
    
    # Add heterogeneity text
    if(het_text != "") {
      p <- p + annotate("text", x = -1.70, y = 0.4, 
                        label = het_text, 
                        hjust = 0, size = 3, color = "gray40", fontface = "italic")
    }
    
    forest_plots[[taxon]] <- p
  }
  
  # Combine plots
  if(require(patchwork) && length(forest_plots) > 0) {
    combined_plot <- forest_plots[[taxa_list[1]]]
    
    if(length(forest_plots) > 1) {
      for(i in 2:length(taxa_list)) {
        combined_plot <- combined_plot / forest_plots[[taxa_list[i]]]
      }
    }
    
    combined_plot <- combined_plot +
      plot_annotation(
        title = "Effect of Rapamycin on Gut Microbiota",
        subtitle = "Meta-analysis of three independent studies (FDR < 0.05)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 12, hjust = 0.5))
      )
  }
  
  # Format table data
  table_data$Study <- factor(table_data$Study, 
                             levels = c("Pooled (All)",
                                        "Bitto 2016 (Diet)",
                                        "Pooled (IP only)",
                                        "Bitto 2016 (IP)", 
                                        "Yang 2022", 
                                        "Han 2021"))
  
  table_data <- table_data[order(table_data$Taxon, table_data$Study), ]
  
  # DEBUG: Print raw weights before formatting
  cat("\n=== DEBUG: Raw Prevotellaceae weights ===\n")
  prev_debug <- table_data %>% 
    filter(Taxon == "Prevotellaceae_UCG-001") %>%
    select(Taxon, Study, Weight)
  print(prev_debug)
  
  cat("\n=== DEBUG: Sum of non-NA weights ===\n")
  print(sum(prev_debug$Weight, na.rm = TRUE))
  
  # Format table
  formatted_table <- table_data %>%
    group_by(Taxon) %>%
    mutate(
      total_weight_for_taxon = sum(Weight[!is.na(Weight)], na.rm = TRUE),
      `Weight (%)` = ifelse(
        is.na(Weight),
        "—",
        sprintf("%.1f%%", 100 * Weight / total_weight_for_taxon)
      ),
      Heterogeneity = case_when(
        !is.na(I2) & !is.na(tau2) ~ sprintf("I²=%.1f%%, τ²=%.4f", I2, tau2),
        TRUE ~ "—"
      )
    ) %>%
    ungroup() %>%
    mutate(
      `Odds Ratio [95% CI]` = sprintf("%.2f [%.2f-%.2f]", OR, OR_lower, OR_upper),
      `Log OR [95% CI]` = sprintf("%.2f [%.2f-%.2f]", LogOR, CI_lower, CI_upper),
      `P-value` = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
    ) %>%
    select(Taxon, Study, `Log OR [95% CI]`, `Odds Ratio [95% CI]`, `P-value`, `Weight (%)`,   Heterogeneity)
  
  return(list(
    individual = forest_plots,
    combined = combined_plot,
    table = formatted_table
  ))
}

# Run
results <- create_forest_from_sigtab(
  sig_tab_data = sig_tab,
  long_data = long_tab_pref,
  taxa_list = c("Bacteroides", "Muribaculum", "Candidatus_Arthromitus"),
  heterogeneity_data_all = heterogeneity_results,  # ← CHANGED to _all
  heterogeneity_data_ip = heterogeneity_IP_only, 
  ip_pooled_data = table.pooled_IP_only
)

plot.MetaForest_combined <- results$combined
table.Forestplot_combined <- results$table

print(plot.MetaForest_combined)
print(table.Forestplot_combined, n = Inf, width = Inf)



#### UnnesesayGAMLSS with Total Dosage ####
# GAMLSS function using total dosage instead of AdminRoute
run_gamlss_with_dose <- function(ps_list) {
  # Create a list to store results
  all_results <- list()
  
  for(study_name in names(ps_list)) {
    cat("\n\n========== ANALYZING", study_name, "==========\n")
    
    # Extract phyloseq object
    ps <- ps_list[[study_name]]
    
    # Get OTU table
    otu <- as(otu_table(ps), "matrix")
    if(!taxa_are_rows(ps)) {
      otu <- t(otu)
    }
    
    # Convert to relative abundance
    rel_abund <- apply(otu, 2, function(x) x / sum(x) * 100)  # margin = 2
    
    # Get metadata
    meta <- data.frame(sample_data(ps))
    
    # Get all taxa
    taxa_list <- rownames(rel_abund)
    total_taxa <- length(taxa_list)
    
    cat("Processing", total_taxa, "taxa in", study_name, "\n")
    
    # Create a list to store results for this study
    study_results <- list()
    
    # Processing counters
    success_count <- 0
    
    # Process each taxon
    for(i in 1:total_taxa) {
      taxon <- taxa_list[i]
      
      # Print progress
      if(i %% 10 == 0) {
        cat("Processing taxon", i, "of", total_taxa, "\n")
      }
      
      # Get abundance values
      abundance_values <- as.numeric(rel_abund[rownames(rel_abund) == taxon, ])
      
      # Skip taxa with no variation
      if(length(unique(abundance_values)) <= 1) {
        next
      }
      
      # Store taxon_data locally in this environment
      data_for_model <- data.frame(
        abundance = abundance_values,
        StandardTreatment = factor(meta$StandardTreatment, levels = c("Control", "Rapamycin")),
        row.names = colnames(rel_abund)
      )
      
      # # Handle zero-inflation
      # epsilon <- 0.01
      # data_for_model$abundance[data_for_model$abundance == 0] <- epsilon
      # data_for_model$abundance[data_for_model$abundance == 100] <- 100 - epsilon
      # 
      # Scale to 0-1
      data_for_model$scaled <- data_for_model$abundance / 100
      
      # Add TotalDose for Study2 if available
      if(study_name == "Study2" && "TotalDose_mgkg" %in% colnames(meta)) {
        data_for_model$TotalDose <- meta$TotalDose_mgkg
        
        # Standardize TotalDose to avoid numerical issues
        data_for_model$TotalDose <- scale(data_for_model$TotalDose, center = TRUE, scale = TRUE)
      }
      
      # Run GAMLSS
      model <- NULL
      tryCatch({
        if(study_name == "Study2" && "TotalDose" %in% colnames(data_for_model)) {
          model <- gamlss(
            scaled ~ StandardTreatment,
            family = BEZI(),
            data = data_for_model,
            trace = FALSE,
            control = gamlss.control(n.cyc = 50)
          )
        } else {
          model <- gamlss(
            scaled ~ StandardTreatment,
            family = BEZI(),
            data = data_for_model,
            trace = FALSE
          )
        }
      }, error = function(e) {
        # Skip silently
      })
      
      # Check if model was successfully fit
      if(!is.null(model)) {
        # Try to extract from summary output
        summ_output <- tryCatch({
          capture.output(summary(model))
        }, error = function(e) {
          return(NULL)
        })
        
        if(!is.null(summ_output)) {
          # Extract treatment effect line
          treatment_line <- grep("StandardTreatmentRapamycin", summ_output, value = TRUE)
          
          if(length(treatment_line) > 0) {
            # Find the line with actual numeric values
            for(line in treatment_line) {
              # Split the line
              values <- strsplit(trimws(line), "\\s+")[[1]]
              
              # Check if it has numeric values
              if(length(values) >= 4 && !is.na(suppressWarnings(as.numeric(values[2])))) {
                # First try the standard format
                effect <- as.numeric(values[2])
                std_error <- as.numeric(values[3])
                
                # Check if t-value and p-value are available
                if(length(values) >= 5) {
                  t_value <- as.numeric(values[4])
                  p_value <- as.numeric(values[5])
                } else {
                  # If not, calculate p-value from t-value
                  t_value <- effect / std_error
                  df <- nrow(data_for_model) - length(coef(model))
                  p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
                }
                
                # Store results
                result <- data.frame(
                  taxon = taxon,
                  effect = effect,
                  std_error = std_error,
                  t_value = t_value,
                  p_value = p_value,
                  stringsAsFactors = FALSE
                )
                
                # Add to results list
                study_results[[length(study_results) + 1]] <- result
                success_count <- success_count + 1
                
                # Break the loop once we found the values
                break
              }
            }
          }
        }
      }
    }
    
    cat("Successfully analyzed", success_count, "taxa out of", total_taxa, "\n")
    
    # Combine results for this study
    if(length(study_results) > 0) {
      study_df <- do.call(rbind, study_results)
      
      # Add FDR correction
      study_df$p_adjust <- p.adjust(study_df$p_value, method = "BH")
      
      # Store in all results
      all_results[[study_name]] <- study_df
      
      # Show summary
      cat("\nResults summary for", study_name, ":\n")
      cat("- Successfully analyzed", nrow(study_df), "taxa\n")
      
      # Report significant taxa
      sig_taxa <- sum(study_df$p_value < 0.05, na.rm = TRUE)
      cat("- Significant taxa (p < 0.05):", sig_taxa, "\n")
      
      # Show top significant taxa
      if(sig_taxa > 0) {
        top_taxa <- head(study_df[order(study_df$p_value), ], 5)
        cat("- Top significant taxa:\n")
        print(top_taxa[, c("taxon", "effect", "t_value", "p_value")])
      }
    } else {
      cat("No successful analyses for", study_name, "\n")
    }
  }
  
  return(all_results)
}
dash <- bind_rows(all_results, .id = "study") |>
  group_by(study) |>
  summarise(
    n_tested     = n(),
    n_raw_sig    = sum(p_value   < 0.05, na.rm = TRUE),
    n_fdr_sig    = sum(p_adjust  < 0.10, na.rm = TRUE),
    med_abs_beta = median(abs(effect), na.rm = TRUE),
    .groups = "drop"
  )

print(dash)
# 
# study  n_tested n_raw_sig n_fdr_sig med_abs_beta
# study  n_tested n_raw_sig n_fdr_sig med_abs_beta
# <chr>     <int>     <int>     <int>        <dbl>
# 1 Study2      145        17         7        0.342
# 2 Study5      126        11         1        0.650
# 3 Study7       86        12        12        0.719

#### TotalDosage####
head(ps_genus_dose_sanitized)

# CREATE dose_lookup FIRST!
dose_lookup <- do.call(rbind, lapply(names(ps_genus_dose_sanitized), function(study_name) {
  ps <- ps_genus_dose_sanitized[[study_name]]
  meta <- data.frame(sample_data(ps))
  unique_dose <- unique(meta$TotalDose_mgkg)
  dose <- unique_dose[unique_dose != 0]
  if(length(dose) == 0) dose <- 0
  data.frame(Study = study_name, TotalDose_mgkg = dose)
}))
print(dose_lookup)


## ANALYSIS 1: ALL 4 STUDIES (including dietary)


   # Add doses to filtered_results
   for (study_name in names(filtered_results)) {
       this_dose <- dose_lookup$TotalDose_mgkg[dose_lookup$Study == study_name]
       filtered_results[[study_name]]$TotalDose_mgkg <- this_dose
     }

   cat("Studies in analysis:\n")
   print(names(filtered_results))

   # Combine all 4 studies
  meta_input_all <- bind_rows(filtered_results, .id = "Study") %>%
  filter(!is.na(effect), !is.na(std_error))
# Meta-regression for all studies
 unique_taxa_all <- unique(meta_input_all$taxon)
 meta_reg_results_all <- list()
 
   for (tx in unique_taxa_all) {
       dat <- meta_input_all %>% filter(taxon == tx)
       if (nrow(dat) >= 3) {
           res <- tryCatch({
               rma(yi = effect, sei = std_error, mods = ~ TotalDose_mgkg, data = dat)
             }, error = function(e) NULL)
           meta_reg_results_all[[tx]] <- res
         }
     }
   # Extract results
   meta_reg_summary_all <- data.frame(
       taxon = names(meta_reg_results_all),
       slope = sapply(meta_reg_results_all, function(fit) {
           if(is.null(fit)) return(NA)
           coef(fit)["TotalDose_mgkg"]
         }),
       se = sapply(meta_reg_results_all, function(fit) {
           if(is.null(fit)) return(NA)
           sqrt(vcov(fit)["TotalDose_mgkg","TotalDose_mgkg"])
         }),
       pval = sapply(meta_reg_results_all, function(fit) {
           if(is.null(fit)) return(NA)
           coef(summary(fit))["TotalDose_mgkg","pval"]
         }),
       I2 = sapply(meta_reg_results_all, function(fit) {
           if(is.null(fit)) return(NA)
           fit$I2
         }),
       Q = sapply(meta_reg_results_all, function(fit) {
           if(is.null(fit)) return(NA)
           fit$QE
         }),
       pval_Q = sapply(meta_reg_results_all, function(fit) {
           if(is.null(fit)) return(NA)
           fit$QEp
         })
     ) %>%
     mutate(
    p_adjust = p.adjust(pval, method="BH"),
         RR_per_mgkg = exp(slope),
         pct_change = 100 * (RR_per_mgkg - 1)
       ) %>%
  arrange(pval)
 
 # Significant taxa (all 4 studies)
   significant_dose_response_all <- meta_reg_summary_all %>%
     filter(p_adjust < 0.1) %>%
     arrange(p_adjust) %>%
     select(taxon, slope, se, pval, p_adjust, I2, Q, pval_Q, RR_per_mgkg, pct_change)
 
   table.dose_all <- significant_dose_response_all %>%
     mutate(
         slope = round(slope, 7),
         se = round(se, 7),
         pval = format(pval, scientific = TRUE, digits = 3),
         p_adjust = format(p_adjust, scientific = TRUE, digits = 3),
         I2 = round(I2, 1),
         Q = round(Q, 3),
         pval_Q = round(pval_Q, 3),
         RR_per_mgkg = round(RR_per_mgkg, 6),
         pct_change = round(pct_change, 4)
       )
 
   kable(table.dose_all, 
            caption = "Dose-Response: All 4 Studies (FDR < 0.10)")

## SENSITIVITY ANALYSIS: Excluding Han 2021 (lowest taxa count)####
 
   # Filter out Study7
   meta_input_minus_han <- meta_input_all %>%
     filter(Study != "Study7")
   
   cat("\n\nStudies excluding Han:\n")
   print(unique(meta_input_minus_han$Study))
   print(meta_input_minus_han %>% distinct(Study, TotalDose_mgkg))
   
   # Meta-regression
   unique_taxa_minus_han <- unique(meta_input_minus_han$taxon)
   meta_reg_results_minus_han <- list()
   
   for (tx in unique_taxa_minus_han) {
     dat <- meta_input_minus_han %>% filter(taxon == tx)
     if (nrow(dat) >= 3) {
       res <- tryCatch({
         rma(yi = effect, sei = std_error, mods = ~ TotalDose_mgkg, data = dat)
       }, error = function(e) NULL)
       meta_reg_results_minus_han[[tx]] <- res
     }
   }
   
   # Extract results with ALL stats
   meta_reg_summary_minus_han <- data.frame(
     taxon = names(meta_reg_results_minus_han),
     slope = sapply(meta_reg_results_minus_han, function(fit) {
       if(is.null(fit)) return(NA)
       coef(fit)["TotalDose_mgkg"]
     }),
     se = sapply(meta_reg_results_minus_han, function(fit) {
       if(is.null(fit)) return(NA)
       sqrt(vcov(fit)["TotalDose_mgkg","TotalDose_mgkg"])
     }),
     pval = sapply(meta_reg_results_minus_han, function(fit) {
       if(is.null(fit)) return(NA)
       coef(summary(fit))["TotalDose_mgkg","pval"]
     }),
     I2 = sapply(meta_reg_results_minus_han, function(fit) {
       if(is.null(fit)) return(NA)
       fit$I2
     }),
     Q = sapply(meta_reg_results_minus_han, function(fit) {
       if(is.null(fit)) return(NA)
       fit$QE
     }),
     pval_Q = sapply(meta_reg_results_minus_han, function(fit) {
       if(is.null(fit)) return(NA)
       fit$QEp
     })
   ) %>%
     mutate(
       p_adjust = p.adjust(pval, method="BH"),
       RR_per_mgkg = exp(slope),
       pct_change = 100 * (RR_per_mgkg - 1)
     ) %>%
     arrange(pval)
   
   # Significant taxa (FDR < 0.10)
   significant_dose_minus_han <- meta_reg_summary_minus_han %>%
     filter(p_adjust < 0.10) %>%
     arrange(p_adjust) %>%
     select(taxon, slope, se, pval, p_adjust, I2, Q, pval_Q, RR_per_mgkg, pct_change)
   
   # Formatted table
   table.dose_minus_han <- significant_dose_minus_han %>%
     mutate(
       slope = round(slope, 7),
       se = round(se, 7),
       pval = format(pval, scientific = TRUE, digits = 3),
       p_adjust = format(p_adjust, scientific = TRUE, digits = 3),
       I2 = round(I2, 1),
       Q = round(Q, 3),
       pval_Q = round(pval_Q, 3),
       RR_per_mgkg = round(RR_per_mgkg, 6),
       pct_change = round(pct_change, 4)
     )
   
   cat("\n=== EXCLUDING HAN 2021 RESULTS ===\n")
   kable(table.dose_minus_han, 
         caption = "Dose-Response: Excluding Han 2021 (FDR < 0.10)")
   
###Sensitivity Analysis Bitto food out #####
   
   # Filter to only IP studies
   meta_input_IP <- meta_input_all %>%
     filter(Study != "Study2_Food")
   
   cat("\n\nIP-only studies:\n")
   print(unique(meta_input_IP$Study))
   
   # Meta-regression for IP only
   unique_taxa_IP <- unique(meta_input_IP$taxon)
   meta_reg_results_IP <- list()
   
   for (tx in unique_taxa_IP) {
     dat <- meta_input_IP %>% filter(taxon == tx)
     if (nrow(dat) >= 2) {  # Need at least 2 studies for meta-regression
       res <- tryCatch({
         rma(yi = effect, sei = std_error, mods = ~ TotalDose_mgkg, data = dat, 
             method = "REML")
       }, error = function(e) NULL)
       meta_reg_results_IP[[tx]] <- res
     }
   }
   
   # Extract results - ROBUST VERSION
   meta_reg_summary_IP <- data.frame(
     taxon = names(meta_reg_results_IP),
     slope = sapply(meta_reg_results_IP, function(fit) {
       if(is.null(fit)) return(NA)
       tryCatch({
         coefs <- coef(fit)
         # Get the dose coefficient (second coefficient, first is intercept)
         if(length(coefs) >= 2) coefs[2] else NA
       }, error = function(e) NA)
     }),
     se = sapply(meta_reg_results_IP, function(fit) {
       if(is.null(fit)) return(NA)
       tryCatch({
         vc <- vcov(fit)
         # Get SE for dose coefficient (position [2,2])
         if(nrow(vc) >= 2 && ncol(vc) >= 2) sqrt(vc[2, 2]) else NA
       }, error = function(e) NA)
     }),
     pval = sapply(meta_reg_results_IP, function(fit) {
       if(is.null(fit)) return(NA)
       tryCatch({
         summ <- summary(fit)
         # Get p-value for dose coefficient (row 2)
         if(nrow(summ$beta) >= 2) summ$pval[2] else NA
       }, error = function(e) NA)
     }),
     I2 = sapply(meta_reg_results_IP, function(fit) {
       if(is.null(fit)) return(NA)
       tryCatch(fit$I2, error = function(e) NA)
     }),
     Q = sapply(meta_reg_results_IP, function(fit) {
       if(is.null(fit)) return(NA)
       tryCatch(fit$QE, error = function(e) NA)
     }),
     pval_Q = sapply(meta_reg_results_IP, function(fit) {
       if(is.null(fit)) return(NA)
       tryCatch(fit$QEp, error = function(e) NA)
     })
   ) %>%
     filter(!is.na(slope)) %>%  # Remove failed fits
     mutate(
       p_adjust = p.adjust(pval, method = "BH"),
       RR_per_mgkg = exp(slope),
       pct_change = 100 * (RR_per_mgkg - 1)
     ) %>%
     arrange(pval)
   
   # Show results
   cat("\n\nIP-only dose-response results (FDR < 0.10):\n")
   print(meta_reg_summary_IP %>% 
           filter(p_adjust < 0.10) %>%
           select(taxon, slope, se, pval, p_adjust, I2, pct_change))

# Check presence of significant taxa in IP studies
sig_taxa_all <- significant_dose_response_all$taxon

for(tx in sig_taxa_all) {
  n_ip <- meta_input_IP %>% 
    filter(taxon == tx) %>% 
    nrow()
  cat(tx, ":", n_ip, "IP studies\n")
}
# Check if all 7 significant taxa from "All 4" are in the IP results
sig_taxa_all <- significant_dose_response_all$taxon

cat("=== CHECKING IP RESULTS FOR ALL 7 TAXA ===\n")
for(tx in sig_taxa_all) {
  # Check if it was tested in IP-only
  ip_result <- meta_reg_summary_IP %>% filter(taxon == tx)
  
  if(nrow(ip_result) > 0) {
    cat(tx, ":\n")
    cat("  - p_adjust:", ip_result$p_adjust, "\n")
    cat("  - slope:", ip_result$slope, "\n")
    cat("  - Significant in IP?", ip_result$p_adjust < 0.10, "\n\n")
  } else {
    cat(tx, ": NOT TESTED in IP-only (missing data)\n\n")
  }
}

# # Show how many taxa total were tested
# cat("Total taxa tested in IP-only:", nrow(meta_reg_summary_IP), "\n")
# cat("Significant in IP-only (FDR < 0.10):", nrow(significant_dose_response_IP), "\n")

# Check presence of all significant taxa (from ALL 4 analysis) in each IP study

sig_taxa_all <- c("Lachnospira", "A2", "Faecalibaculum", "Alistipes", 
                  "Roseburia", "Muribaculum", "Butyricicoccus")

ip_studies <- c("Study2_IP", "Study5", "Study7")

# Create presence/absence matrix
presence_matrix <- data.frame()

for(tx in sig_taxa_all) {
  row_data <- data.frame(Taxon = tx)
  
  for(study in ip_studies) {
    # Check if this taxon is present in this study
    n_obs <- meta_input_IP %>%
      filter(taxon == tx, Study == study) %>%
      nrow()
    
    row_data[[study]] <- ifelse(n_obs > 0, "✓", "✗")
  }
  
  # Add total count
  row_data$Total_Studies <- sum(row_data[ip_studies] == "✓")
  
  # Add if testable (need 2+ studies for meta-regression)
  row_data$Testable <- ifelse(row_data$Total_Studies >= 2, "Yes", "No")
  
  presence_matrix <- rbind(presence_matrix, row_data)
}

# Display the table
cat("\n=== TAXA PRESENCE IN IP STUDIES ===\n")
print(presence_matrix)

# More detailed version with actual sample counts
presence_detailed <- data.frame()

for(tx in sig_taxa_all) {
  row_data <- data.frame(Taxon = tx)
  
  for(study in ip_studies) {
    n_obs <- meta_input_IP %>%
      filter(taxon == tx, Study == study) %>%
      nrow()
    
    row_data[[paste0(study, "_n")]] <- n_obs
  }
  
  row_data$Total_n <- sum(row_data[2:4])
  
  presence_detailed <- rbind(presence_detailed, row_data)
}

cat("\n=== DETAILED SAMPLE COUNTS ===\n")
print(presence_detailed)

# Check which taxa have data in ALL 3 IP studies
taxa_in_all_3 <- presence_matrix %>% 
  filter(Total_Studies == 3) %>%
  pull(Taxon)

cat("\n=== TAXA PRESENT IN ALL 3 IP STUDIES ===\n")
cat(paste(taxa_in_all_3, collapse = ", "), "\n")

# Summary
cat("\n=== SUMMARY ===\n")
cat("Taxa in all 3 IP studies:", length(taxa_in_all_3), "\n")
cat("Taxa in 2 IP studies:", sum(presence_matrix$Total_Studies == 2), "\n")
cat("Taxa in 1 IP study:", sum(presence_matrix$Total_Studies == 1), "\n")
cat("Testable via meta-regression (2+):", sum(presence_matrix$Testable == "Yes"), "\n")

# Check ALL taxa (not just the 7 from ALL analysis) present in IP studies
all_taxa_presence <- data.frame()

unique_taxa_IP <- unique(meta_input_IP$taxon)

for(tx in unique_taxa_IP) {
  row_data <- data.frame(Taxon = tx)
  
  for(study in ip_studies) {
    n_obs <- meta_input_IP %>%
      filter(taxon == tx, Study == study) %>%
      nrow()
    
    row_data[[study]] <- ifelse(n_obs > 0, "✓", "✗")
  }
  
  # Count total studies
  row_data$Total_Studies <- sum(row_data[ip_studies] == "✓")
  row_data$Testable <- ifelse(row_data$Total_Studies >= 3, "Yes (3 studies)", 
                              ifelse(row_data$Total_Studies == 2, "Marginal (2 studies)", "No"))
  
  all_taxa_presence <- rbind(all_taxa_presence, row_data)
}

# Summary
cat("\n=== ALL TAXA IN IP STUDIES ===\n")
cat("Total unique taxa:", nrow(all_taxa_presence), "\n")
cat("In all 3 IP studies:", sum(all_taxa_presence$Total_Studies == 3), "\n")
cat("In 2 IP studies:", sum(all_taxa_presence$Total_Studies == 2), "\n")
cat("In 1 IP study:", sum(all_taxa_presence$Total_Studies == 1), "\n")

# Show taxa in all 3 studies
taxa_3_studies <- all_taxa_presence %>% 
  filter(Total_Studies == 3) %>%
  arrange(Taxon)

cat("\n=== TAXA IN ALL 3 IP STUDIES (fully testable) ===\n")
print(taxa_3_studies)

# Compare with what was actually significant in IP-only
cat("\n=== SIGNIFICANT IN IP-ONLY (q < 0.10) ===\n")
print(meta_reg_summary_IP)

#### Effect size (dose) — minimal-change, uses significant_dose_response ####
library(pROC)
## DEBUG: Check what's happening
cat("Number of studies in X_list:", length(X_list), "\n")
print(names(X_list))

# Check each study
for(st in names(X_list)) {
  cat("\n", st, ":\n")
  cat("  Dimensions:", dim(X_list[[st]]), "\n")
  cat("  Taxa:", rownames(X_list[[st]]), "\n")
}

# Check if genus_mat_adj was created
if(exists("genus_mat_adj")) {
  cat("\ngenus_mat_adj dimensions:", dim(genus_mat_adj), "\n")
  cat("Is it a matrix?", is.matrix(genus_mat_adj), "\n")
} else {
  cat("\ngenus_mat_adj doesn't exist!\n")
}

## DEBUG: Check what's happening
cat("Number of studies in X_list:", length(X_list), "\n")
print(names(X_list))

# Check each study
for(st in names(X_list)) {
  cat("\n", st, ":\n")
  cat("  Dimensions:", dim(X_list[[st]]), "\n")
  cat("  Taxa:", rownames(X_list[[st]]), "\n")
}

# Check if genus_mat_adj was created
if(exists("genus_mat_adj")) {
  cat("\ngenus_mat_adj dimensions:", dim(genus_mat_adj), "\n")
  cat("Is it a matrix?", is.matrix(genus_mat_adj), "\n")
} else {
  cat("\ngenus_mat_adj doesn't exist!\n")
}

stopifnot(exists("significant_dose_response_all"))

## DEBUG dose_meta (FIXED)
cat("\nNumber of metadata dataframes:", length(md_list), "\n")
print(names(md_list))

# Manually build dose_meta to debug
dose_meta <- dplyr::bind_rows(md_list)
rownames(dose_meta) <- dose_meta$SampleID

cat("\ndose_meta class:", class(dose_meta), "\n")
cat("dose_meta dimensions:", dim(dose_meta), "\n")
cat("Has SampleID column?", "SampleID" %in% names(dose_meta), "\n")

# Check rownames
cat("First 5 rownames:", head(rownames(dose_meta), 5), "\n")
cat("Number of rownames:", length(rownames(dose_meta)), "\n")

# Check if SampleID matches column names
cat("\nFirst 5 genus_mat_adj colnames:", head(colnames(genus_mat_adj), 5), "\n")
cat("First 5 dose_meta rownames:", head(rownames(dose_meta), 5), "\n")

# Check for match
match_count <- sum(colnames(genus_mat_adj) %in% rownames(dose_meta))
cat("\nMatching samples:", match_count, "out of", ncol(genus_mat_adj), "\n")

# Now try the problematic line
cat("\nAttempting to subset genus_mat_adj...\n")
tryCatch({
  genus_mat_adj <- genus_mat_adj[, rownames(dose_meta), drop = FALSE]
  cat("SUCCESS!\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})

# Check what exists
cat("genus_mat_adj exists?", exists("genus_mat_adj"), "\n")
cat("dose_meta exists?", exists("dose_meta"), "\n")

# If dose_meta doesn't exist, that's the problem!
# The if() block at the top must have been skipped

# Let's check the md_list
cat("md_list exists?", exists("md_list"), "\n")
if(exists("md_list")) {
  cat("md_list length:", length(md_list), "\n")
}
## 0) weights from your dose-response meta (slopes per mg/kg)
w_tbl <- significant_dose_response_all %>%
  dplyr::select(feature = taxon, coef = slope) %>%
  filter(is.finite(coef))
w <- setNames(w_tbl$coef, w_tbl$feature)
stopifnot(length(w) >= 2)

## 1) If needed, build genus_mat_adj (features × samples) + dose_meta (per-sample meta)
if (!exists("genus_mat_adj") || !exists("dose_meta")) {
  stopifnot(exists("ps_genus_dose_sanitized"))
  
  sig_set <- unique(names(w))            # only need weighted genera
  X_list <- list(); md_list <- list()
  
  for (st in names(ps_genus_dose_sanitized)) {
    ps <- ps_genus_dose_sanitized[[st]]
    
    M <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) M <- t(M)    # rows = taxa, cols = samples
    M <- sweep(M, 2, pmax(colSums(M), 1e-12), "/")  # per-sample proportions
    
    rn <- rownames(M)
    rn <- sub("^g__", "", rn)
    rownames(M) <- rn
    
    keep <- intersect(sig_set, rn)
    if (length(keep) == 0) next
    
    M_keep <- M[keep, , drop = FALSE]
    # pad rows to full sig_set so all studies have same taxa order
    M_pad <- matrix(0, nrow = length(sig_set), ncol = ncol(M_keep),
                    dimnames = list(sig_set, colnames(M_keep)))
    M_pad[rownames(M_keep), ] <- M_keep
    colnames(M_pad) <- paste(st, colnames(M_pad), sep = "::")
    X_list[[st]] <- M_pad
    
    md <- as.data.frame(sample_data(ps))
    # match metadata to kept samples (original sample IDs are column names before prefix)
    orig_samp <- colnames(as(otu_table(ps), "matrix"))
    samp_kept <- sub("^.*::", "", colnames(M_pad))
    idx <- match(samp_kept, orig_samp)
    md <- md[idx, , drop = FALSE]
    
    md$SampleID       <- colnames(M_pad)
    md$StudyID_Split  <- st
    md_list[[st]]     <- md
  }
  
  genus_mat_adj <- do.call(cbind, X_list)
  dose_meta     <- dplyr::bind_rows(md_list)
  rownames(dose_meta) <- dose_meta$SampleID
  
  # align sample order
  genus_mat_adj <- genus_mat_adj[, rownames(dose_meta), drop = FALSE]
}

## Ensure required meta columns
if (!("StandardTreatment" %in% names(dose_meta))) {
  stop("dose_meta is missing StandardTreatment.")
}
dose_meta$StandardTreatment <- factor(dose_meta$StandardTreatment, levels = c("Control","Rapamycin"))

# build/repair log-dose
if (!("LogTotalDose" %in% names(dose_meta)) || all(is.na(dose_meta$LogTotalDose))) {
  dose_col <- if ("TotalDose_mgkg" %in% names(dose_meta)) "TotalDose_mgkg" else
    if ("TotalDose"      %in% names(dose_meta)) "TotalDose" else NA
  if (is.na(dose_col)) stop("Need LogTotalDose or TotalDose_mgkg in dose_meta.")
  dose_meta$LogTotalDose <- log10(pmax(as.numeric(dose_meta[[dose_col]]), 1e-8))
}
dose_meta$LogTotalDose <- as.numeric(dose_meta$LogTotalDose)

## 2) DoseScore in Rapamycin samples
rapa_idx <- dose_meta$StandardTreatment == "Rapamycin"
feats    <- intersect(rownames(genus_mat_adj), names(w))
stopifnot(length(feats) >= 2)

Xr  <- t(genus_mat_adj[feats, rapa_idx, drop = FALSE])   # samples × features
mdr <- droplevels(dose_meta[rapa_idx, , drop = FALSE])

score <- as.numeric(Xr %*% w[feats])
mdr$DoseScore        <- as.numeric(scale(score))
mdr$LogTotalDose_std <- as.numeric(scale(mdr$LogTotalDose))

## 3) Global association + partial R^2 + Spearman
fit    <- lm(DoseScore ~ LogTotalDose_std + StudyID_Split, data = mdr)
fit_red<- lm(DoseScore ~ StudyID_Split, data = mdr)
SS_full<- sum(resid(fit)^2); SS_red <- sum(resid(fit_red)^2)
partial_R2 <- (SS_red - SS_full) / SS_red
spe <- suppressWarnings(cor.test(mdr$DoseScore, mdr$LogTotalDose, method = "spearman"))

cat("\n=== GLOBAL DOSE → DOSE-SCORE ASSOCIATION ===\n")
cat(sprintf("Linear (adj. StudyID_Split): slope=%.3f, p=%.3g\n",
            coef(fit)["LogTotalDose_std"], summary(fit)$coefficients["LogTotalDose_std","Pr(>|t|)"]))
cat(sprintf("Partial R^2 for dose (over Study): %.3f\n", partial_R2))
cat(sprintf("Spearman rho: %.3f (p=%.3g)\n", spe$estimate, spe$p.value))

## 4) Within-study slopes
per_study <- lapply(split(mdr, mdr$StudyID_Split), function(df){
  if (nrow(df) < 3 || sd(df$LogTotalDose) == 0) return(NULL)
  fm <- lm(DoseScore ~ scale(LogTotalDose), data = df)
  data.frame(
    StudyID_Split = unique(df$StudyID_Split),
    slope = unname(coef(fm)[2]),
    p     = summary(fm)$coefficients[2,4],
    n     = nrow(df),
    stringsAsFactors = FALSE
  )
})
per_study_tbl <- dplyr::bind_rows(Filter(Negate(is.null), per_study)) %>%
  mutate(sig = ifelse(p < 0.05, "*", "ns"))
cat("\n=== WITHIN-STUDY SLOPES (DoseScore ~ log dose) ===\n")
print(per_study_tbl)

## 5) Visualization
p_dose_score <- ggplot(mdr, aes(x = LogTotalDose, y = DoseScore, color = StudyID_Split)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "Dose–Effect Score vs Log10 Total Dose (Rapamycin only)",
       subtitle = "Linear fit shown in black; points colored by StudyID_Split",
       x = "Log10 Total Dose (mg/kg)", y = "Dose–Effect Score (z)") +
  theme_bw()
print(p_dose_score)

## 6) Optional AUC (extremes only)
q <- quantile(mdr$LogTotalDose, c(1/3, 2/3), na.rm = TRUE)
use <- mdr$LogTotalDose <= q[1] | mdr$LogTotalDose >= q[2]
if (sum(use) >= 6) {
  y_bin <- as.integer(mdr$LogTotalDose[use] >= q[2])    # 1=high dose
  auc_meta <- as.numeric(pROC::auc(y_bin, mdr$DoseScore[use]))
  cat(sprintf("\nAUC (high vs low dose, extremes only) using DoseScore: %.3f\n", auc_meta))
}

####Heatmap#####
# Suppose filtered_results is a list of data frames (CHANGED: was filtered_results_dosage)
meta_input <- do.call(rbind, lapply(names(filtered_results), function(study) {
  df <- filtered_results[[study]]
  df$Study <- study
  df
}))

heatmap_mat <- meta_input %>%
  select(taxon, Study, effect) %>%
  pivot_wider(names_from = Study, values_from = effect) %>%
  column_to_rownames("taxon")

# Define column renaming and ordering
old_names <- c("Study2_Food", "Study2_IP", "Study5", "Study7")
new_names <- c("Bitto2016_Food", "Bitto2016_IP", "Yang2022", "Han2021")
desired_order <- c("Bitto2016_Food", "Bitto2016_IP", "Yang2022", "Han2021")

# Rename columns
rename_lookup <- setNames(new_names, old_names)
for (i in seq_along(old_names)) {
  if (old_names[i] %in% colnames(heatmap_mat)) {
    colnames(heatmap_mat)[colnames(heatmap_mat) == old_names[i]] <- new_names[i]
  }
}

# Remove unclassified/uncultured taxa
clean_mat <- heatmap_mat[!grepl("unclassified|uncultured", rownames(heatmap_mat), ignore.case = TRUE), ]

# Filter for taxa present in at least 2 studies
min_studies <- 2
row_valid <- rowSums(!is.na(clean_mat)) >= min_studies
mat_filtered <- clean_mat[row_valid, ]

# Ensure columns are in desired order
mat_filtered <- mat_filtered[, desired_order, drop = FALSE]

# Clean up bacterial names (remove "g__" prefix if present)
rownames(mat_filtered) <- gsub("^g__", "", rownames(mat_filtered))

# Filter with stricter threshold (min 3 studies)
min_studies_strict <- 3
row_valid_strict <- rowSums(!is.na(clean_mat)) >= min_studies_strict
mat_strict <- clean_mat[row_valid_strict, ]


# Create the heatmap with improved formatting
library(pheatmap)
plot.heatmap_dose_response <- pheatmap(as.matrix(mat_strict),
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Don't cluster columns to maintain your desired order
         main = "Effect Sizes by Taxon and Study/Dose",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         na_col = "lightgray",
         # Make bacterial names smaller
         fontsize_row = 8,        # Smaller row labels (bacterial names)
         fontsize_col = 10,       # Column labels (study names)
         # Angle the study names
         angle_col = 45,          # Angle column labels
         # Add some spacing and borders
         border_color = "white",
         cellwidth = 25,
         cellheight = 8,
         # Legend and title formatting
         fontsize = 10,
         main_fontsize = 14,
         # You can also customize the clustering method if needed
         clustering_distance_rows = "euclidean",
         clustering_method = "complete")


### Scatter plot ####

## Which Lachnospira labels exist?
unique(meta_input_all$taxon[grepl("lachnospira", meta_input_all$taxon, ignore.case = TRUE)])

## Coverage by arm (treated vs control-level effect estimates per arm)
meta_input_all %>%
  dplyr::filter(grepl("^Lachnospira$", taxon, ignore.case = TRUE)) %>%
  dplyr::count(Study, name = "n_rows")
# how many study arms have an effect estimate for Lachnospira?
lach_arm_coverage <- meta_input_all %>%
  dplyr::filter(taxon == "Lachnospira") %>%
  dplyr::distinct(Study) %>%      # or StudyArm if you have Study+Route
  dplyr::arrange(Study)

lach_arm_coverage
n_arms <- nrow(lach_arm_coverage)   # expect 4 if it's in all arms

# Filter for significant taxa present in at least 3 studies
sig_taxa <- meta_reg_summary_all$taxon[meta_reg_summary_all$p_adjust < 0.05]
scatter_data <- meta_input_all %>%
  filter(taxon %in% sig_taxa) %>%
  group_by(taxon) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  filter(!grepl("unclassified|uncultured", taxon, ignore.case = TRUE))

# Faceted scatter plot
plot.scatter_dose_response <- ggplot(scatter_data, aes(x = TotalDose_mgkg, y = effect, color = Study)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "Dose-Response for Significant Genera (≥3 Studies, FDR < 0.05)",
       x = "Total Dose (mg/kg)",
       y = "Effect size")

plot.scatter_dose_response

###Bubble plot ###
# Define significant genera except Unclassified_3
sig_taxa <- meta_reg_summary_all %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%
  pull(taxon)

# Filter your data for just these genera
plot_data <- meta_input_all %>% filter(taxon %in% sig_taxa)

plot.bubble_dose_response <- ggplot(plot_data, aes(x = TotalDose_mgkg, y = effect, color = taxon)) +
  geom_point(aes(size = 1/std_error), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, aes(group = taxon, color = taxon)) +
  theme_bw() +
  labs(title = "Combined Dose-Response Meta-Regression (Significant Genera, Excl. Unclassified_3)",
       x = "Total Dose (mg/kg)",
       y = "Effect size (meta-analysis)",
       color = "Genus")

plot.bubble_dose_response

###Forest plot#####
# Make a barplot of slopes for significant genera
meta_reg_summary_all %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%
  ggplot(aes(x = reorder(taxon, slope), y = slope)) +
  geom_col(aes(fill = slope > 0)) +
  geom_errorbar(aes(ymin = slope - 1.96*NA, ymax = slope + 1.96*NA), width = 0.2) + # Replace NA with SE if you have it
  coord_flip() +
  theme_bw() +
  labs(title = "Regression slopes by genus", x = "Genus", y = "Slope")

plot.barplot_dose_slopes <- meta_reg_summary_all %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%
  ggplot(aes(x = reorder(taxon, slope), y = slope)) +
  geom_col(aes(fill = slope > 0)) +
  geom_errorbar(aes(ymin = slope - 1.96*se, ymax = slope + 1.96*se), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(title = "Regression slopes by genus", x = "Genus", y = "Slope") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"))

# More selective significance threshold and better design
plot.barplot_dose_slopes_improved <- meta_reg_summary_all %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%  # More stringent + remove unclassified
  mutate(
    Direction = ifelse(slope > 0, "Increased with dose", "Decreased with dose"),
    Taxon_clean = gsub("_", " ", taxon)  # Clean names
  ) %>%
  ggplot(aes(x = reorder(Taxon_clean, slope), y = slope)) +
  geom_col(aes(fill = Direction), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = slope - 1.96*se, ymax = slope + 1.96*se), 
                width = 0.3, size = 0.8, color = "black") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray50"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_manual(
    values = c("Decreased with dose" = "#2166AC", "Increased with dose" = "#D73027"),
    name = "Dose-Response Pattern"
  ) +
  labs(
    title = "Dose-Response Relationships in Gut Microbiota",
    subtitle = "Meta-regression slopes for genera with FDR < 0.05",
    x = "Genus",
    y = "Slope (log Rate Ratio change per mg/kg)",
    caption = "Error bars show 95% confidence intervals"
  )

plot.barplot_dose_slopes_improved

print(meta_reg_summary_all %>%filter(p_adjust < 0.05, taxon != "Unclassified_3") )  # More stringent + remove unclassified)
# Even cleaner dot plot version
plot.dotplot_dose_slopes <- meta_reg_summary_all %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%
  mutate(
    Taxon_clean = gsub("_", " ", taxon),
    Significance = case_when(
      p_adjust < 0.001 ~ "***",
      p_adjust < 0.01 ~ "**",
      p_adjust < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = slope, y = reorder(Taxon_clean, slope))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = slope - 1.96*se, xmax = slope + 1.96*se), 
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = slope > 0), size = 4, alpha = 0.8) +
  geom_text(aes(label = Significance), hjust = -0.8, vjust = -0.2, size = 4, fontface = "bold") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#2166AC")) +
  theme(legend.position = "none") +
  labs(
    title = "Significant Dose-Response Effects",
    subtitle = "Meta-regression slopes for genera with FDR < 0.05",
    x = "Slope (log Odds Ratio change per mg/kg)",
    y = "Genus"
  )

plot.dotplot_dose_slopes

# Check significant taxa
print(meta_reg_summary_all %>% filter(p_adjust < 0.05, taxon != "Unclassified_3"))

# Percentage version of the dot plot
plot.dotplot_dose_slopes_percent <- meta_reg_summary_all %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%
  mutate(
    Taxon_clean = gsub("_", " ", taxon),
    # Calculate percentage confidence intervals
    pct_lower = 100 * (exp(slope - 1.96*se) - 1),
    pct_upper = 100 * (exp(slope + 1.96*se) - 1),
    Significance = case_when(
      p_adjust < 0.001 ~ "***",
      p_adjust < 0.01 ~ "**",
      p_adjust < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = pct_change, y = reorder(Taxon_clean, pct_change))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = pct_lower, xmax = pct_upper), 
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = pct_change > 0), size = 4, alpha = 0.8) +
  geom_text(aes(label = Significance), hjust = -0.5, size = 4, fontface = "bold") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#2166AC")) +
  theme(legend.position = "none") +
  labs(
    title = "Significant Dose-Response Effects",
    subtitle = "Meta-regression effects as % change per mg/kg dose (FDR < 0.05)",
    x = "% Change per mg/kg",
    y = "Genus"
  )

plot.dotplot_dose_slopes_percent

## 1) One-time: compute SD of the moderator (across study arms)
dose_sd <- meta_input %>%
  dplyr::distinct(Study, TotalDose_mgkg) %>%
  dplyr::pull(TotalDose_mgkg) %>%
  sd()

## 2) Build the canonical table with per‑SD columns
dose_meta <- meta_reg_summary_all %>%
  dplyr::rename(
    slope_logOR_per_mgkg = slope,
    se_logOR_per_mgkg    = se,
    p_raw                = pval,
    p_fdr                = p_adjust,
    Q_p                  = pval.Q
  ) %>%
  dplyr::mutate(
    slope_logOR_per_SD = slope_logOR_per_mgkg * dose_sd,
    se_logOR_per_SD    = se_logOR_per_mgkg    * dose_sd,
    ci_lo_logOR_per_SD = (slope_logOR_per_mgkg - 1.96*se_logOR_per_mgkg) * dose_sd,
    ci_hi_logOR_per_SD = (slope_logOR_per_mgkg + 1.96*se_logOR_per_mgkg) * dose_sd,
    pct_change_per_SD  = 100 * (exp(slope_logOR_per_SD) - 1),
    pct_lo_per_SD      = 100 * (exp(ci_lo_logOR_per_SD) - 1),
    pct_hi_per_SD      = 100 * (exp(ci_hi_logOR_per_SD) - 1)
  )
plot.dotplot_dose_slopes <- dose_meta %>%
  dplyr::filter(p_fdr < 0.05, taxon != "Unclassified_3") %>%
  dplyr::mutate(
    Taxon_clean = gsub("_"," ", taxon),
    Sig = dplyr::case_when(p_fdr < 0.001 ~ "***",
                           p_fdr < 0.01  ~ "**",
                           p_fdr < 0.05  ~ "*",
                           TRUE ~ "")
  ) %>%
  ggplot(aes(x = slope_logOR_per_SD, y = reorder(Taxon_clean, slope_logOR_per_SD))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo_logOR_per_SD, xmax = ci_hi_logOR_per_SD),
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = slope_logOR_per_SD > 0), size = 4, alpha = 0.85) +
  geom_text(aes(label = Sig), hjust = -0.5, size = 4, fontface = "bold") +
  scale_color_manual(values = c("TRUE"="#D73027","FALSE"="#2166AC")) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(
    title    = "Significant dose–response effects",
    subtitle = paste0("Slopes on log-OR scale per 1 SD dose (SD = ", round(dose_sd, 0), " mg/kg)"),
    x        = "Slope (log OR per 1 SD dose)",
    y        = "Genus"
  )


print(meta_reg_summary_IP %>%filter(p_adjust < 0.05, taxon != "Unclassified_3") )  # More stringent + remove unclassified)
# Even cleaner dot plot version
plot.dotplot_dose_slopes <- meta_reg_summary_IP %>%
  filter(p_adjust < 0.05, taxon != "Unclassified_3") %>%
  mutate(
    Taxon_clean = gsub("_", " ", taxon),
    Significance = case_when(
      p_adjust < 0.001 ~ "***",
      p_adjust < 0.01 ~ "**",
      p_adjust < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = slope, y = reorder(Taxon_clean, slope))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = slope - 1.96*se, xmax = slope + 1.96*se), 
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = slope > 0), size = 4, alpha = 0.8) +
  geom_text(aes(label = Significance), hjust = -0.5, size = 4, fontface = "bold") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#2166AC")) +
  theme(legend.position = "none") +
  labs(
    title = "Significant Dose-Response Effects",
    subtitle = "Meta-regression slopes for genera with FDR < 0.05",
    x = "Slope (log Odds Ratio change per mg/kg)",
    y = "Genus"
  )

plot.dotplot_dose_slopes

#### PLots ###
plot.combined_plot_dose
plot.p_dose_before
plot.p_dose_after
plot.treatment_effect_plot

Heatmap.metanalysis_treatment_only
plot.MetaForest_combined
plot.heatmap_dose_response
plot.scatter_dose_response 
plot.barplot_dose_slopes
plot.dotplot_dose_slopes_percent

### 2025-03-13 Exported as PDF
pdf("MetaRapamycinAnalysis_20250506_V01.pdf", width = 14, height = 10)
plot.combined_plot_dose
plot.p_dose_before
plot.p_dose_after
plot.treatment_effect_plot
plot.MetaForest_combined
Heatmap.metanalysis_treatment_only
plot.heatmap_dose_response
plot.scatter_dose_response 
plot.barplot_dose_slopes
dev.off()

##### tables (in progress) #####
table.high_se_summary ##from treatmentonly
table.sigtax_publication
table.treatment_GAMLSS_result
table.Forestplot_combined
significant_dose_response_table
##### image #####
# Save the current R workspace
save.image("Rapamycin_Metaanalysis_V04_20250728_latest.RData")

# Load the saved workspace for later use
load("Rapamycin_Metaanalysis_V04_20250728_latest.RData")

# Print the session information
sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.4.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Panama
# tzcode source: internal
# 
# attached base packages:
#   [1] parallel  splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# # other attached packages:
# [1] pheatmap_1.0.12             meta_8.1-0                  metafor_4.8-0               numDeriv_2016.8-1.1        
# [5] metadat_1.4-0               Matrix_1.7-1                metamicrobiomeR_1.2         gamlss_5.4-22              
# [9] nlme_3.1-166                gamlss.dist_6.1-1           gamlss.data_6.0-6           MMUPHin_1.20.0             
# [13] DESeq2_1.46.0               SummarizedExperiment_1.36.0 Biobase_2.66.0              MatrixGenerics_1.18.1      
# [17] matrixStats_1.5.0           GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1             
# [21] S4Vectors_0.44.0            BiocGenerics_0.52.0         splitstackshape_1.4.8       ape_5.8-1                  
# [25] vegan_2.6-10                lattice_0.22-6              permute_0.9-7               microbiome_1.28.0          
# [29] phyloseq_1.50.0             knitr_1.50                  viridis_0.6.5               viridisLite_0.4.2          
# [33] patchwork_1.3.0             ggpubr_0.6.0                lubridate_1.9.4             forcats_1.0.0              
# [37] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.4                 readr_2.1.5                
# [41] tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.2               tidyverse_2.0.0                          
# 
