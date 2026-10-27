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

# load required packages 
library(dplyr) # required for '%>%' in 'RelativeAbundancePlots'
library(tidyr) # data processing
library(vegan) # for 'adonis()'
library(ape) # for read.tree()
library(DESeq2) # required for 'DESeq', 'results'
library(ggplot2) # plots
library(ggpubr) # improved ggplots, error bars, wilcoxon/kruskal etc
library(patchwork) # arranging ggplots
library(microbiome) # core microbiome calculations
library(viridis) # for 'scale_fill_viridis()'
library(splitstackshape) # for core microbiome matrix
library(tidyverse) # for core microbiome matrix
library(knitr) # GLM tables
library(phyloseq)
library(MMUPHin)      # For batch correction
library(metamicrobiomeR) # For differential abundance testing with zero-inflation
library(metafor)      # For meta-analysis support
library(reshape2)
library(ggVennDiagram)

packageVersion("MMUPHin")

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


# Check genera in filtered_ps_list (after probiotics removal)
cat("=== Genus Counts After Filtering (Before Standardization) ===\n")
sapply(filtered_ps_list, function(ps) length(unique(tax_table(ps)[,"Genus"])))

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

# Use StandardTreatment as your primary variable of interest to compare Control vs. Rapamycin across all studies
# Include AdminRoute as a covariate to account for different administration methods

#### 4. STREAMLINED TAXON HARMONISATION ####

# 4a) Quick rank counts (unchanged)
rank_vec <- c("Phylum","Class","Order","Family","Genus","Species")
taxa_by_rank <- lapply(ps_list, function(ps) {
  tt <- as.data.frame(tax_table(ps))
  sapply(rank_vec, function(r) {
    if (r %in% names(tt)) length(unique(tt[[r]])) else NA
  })
})
taxa_df <- as.data.frame(do.call(rbind, taxa_by_rank))
print(taxa_df)

# 4b) Spot-check raw Genus names in Study2
suspect   <- "Study2"
raw_genera <- sort(unique(tax_table(standardized_ps_list[[suspect]])[,"Genus"]))
head(raw_genera, 30)

# 4c) Basic whitespace trim on Genus across all studies
standardized_ps_list <- lapply(standardized_ps_list, function(ps) {
  tt <- as.data.frame(tax_table(ps))
  if ("Genus" %in% colnames(tt)) {
    tt$Genus <- trimws(tt$Genus)
    tax_table(ps) <- tax_table(as.matrix(tt))
  }
  ps
})
# sanity check
sapply(standardized_ps_list, function(ps) length(unique(tax_table(ps)[,"Genus"])))

# Check genera in standardized_ps_list (what goes into genus harmonization)
cat("=== Genus Counts After Standardization (Before Harmonization) ===\n") 
sapply(standardized_ps_list, function(ps) length(unique(tax_table(ps)[,"Genus"])))

#### 5.VennDiagramm ####
# Function to clean genus names for Venn diagram (minimal cleaning)
clean_genus_for_venn <- function(genus_vector) {
  cleaned <- genus_vector %>%
    # Remove prefixes
    gsub("^[kpcofgs]__", "", ., ignore.case = TRUE) %>%
    # Trim whitespace
    trimws() %>%
    # Basic collapsing - only the most obvious cases
    gsub("^Uncultured.*$", "Uncultured", ., ignore.case = TRUE) %>%
    gsub("^Unclassified.*$", "Unclassified", ., ignore.case = TRUE) %>%
    # Your specific cases
    gsub("Tuzzerella", "Tyzzerella", .) %>%
    gsub("Erysipelatoclostridiaceae", "Erysipelatoclostridium", .)
  
  # Remove NA, empty, and "NA" strings
  cleaned <- cleaned[!is.na(cleaned) & cleaned != "" & cleaned != "NA"]
  
  # Return unique values
  unique(cleaned)
}
# Extract and clean genus names using the function
study2_genera <- clean_genus_for_venn(tax_table(standardized_ps_list$Study2)[,"Genus"])
study5_genera <- clean_genus_for_venn(tax_table(standardized_ps_list$Study5)[,"Genus"])
study7_genera <- clean_genus_for_venn(tax_table(standardized_ps_list$Study7)[,"Genus"])

# Remove NA and empty values
study2_genera <- study2_genera[!is.na(study2_genera) & study2_genera != "" & study2_genera != "NA"]
study5_genera <- study5_genera[!is.na(study5_genera) & study5_genera != "" & study5_genera != "NA"]
study7_genera <- study7_genera[!is.na(study7_genera) & study7_genera != "" & study7_genera != "NA"]

# Check what we got
cat("=== GENUS COUNTS BY STUDY ===\n")
cat("Study2:", length(study2_genera), "genera\n")
cat("Study5:", length(study5_genera), "genera\n") 
cat("Study7:", length(study7_genera), "genera\n")

# Show some examples of what got collapsed
cat("\n=== EXAMPLES OF COLLAPSED GENERA ===\n")
all_raw_genera <- c(
  tax_table(standardized_ps_list$Study2)[,"Genus"],
  tax_table(standardized_ps_list$Study5)[,"Genus"],
  tax_table(standardized_ps_list$Study7)[,"Genus"]
)
uncultured_variants <- unique(all_raw_genera[grepl("Uncultured", all_raw_genera, ignore.case = TRUE)])
unclassified_variants <- unique(all_raw_genera[grepl("Unclassified", all_raw_genera, ignore.case = TRUE)])

cat("Uncultured variants found:", length(uncultured_variants), "\n")
if(length(uncultured_variants) > 0) cat("Examples:", paste(head(uncultured_variants, 5), collapse = ", "), "\n")
cat("Unclassified variants found:", length(unclassified_variants), "\n")
if(length(unclassified_variants) > 0) cat("Examples:", paste(head(unclassified_variants, 5), collapse = ", "), "\n")

study2_genera <- sort(unique(study2_genera))
study5_genera <- sort(unique(study5_genera))
study7_genera <- sort(unique(study7_genera))

cat("\n== Study2 genera (cleaned) ==\n", paste(study2_genera, collapse=", "), "\n")
cat("\n== Study5 genera (cleaned) ==\n", paste(study5_genera, collapse=", "), "\n")
cat("\n== Study7 genera (cleaned) ==\n", paste(study7_genera, collapse=", "), "\n")

# Create Venn diagram

# Replace scale_fill_brewer() with scale_fill_viridis_c() for continuous values
plot.Venndiagram_allstudys_v3 <- ggVennDiagram(
  list(
    "Bitto (2016)" = study2_genera,
    "Yang (2022)" = study5_genera, 
    "Han (2021)" = study7_genera
  ),
  label = "count",
  label_percent_digit = 0,
  label_size = 5,
  label_alpha = 1,
  edge_lty = "solid",
  edge_size = 0.5,
  set_size = 4.5
) +
  scale_fill_gradient2(
    low = "#FF8080",      # Stronger pink for low values
    mid = "#80B3FF",      # Stronger blue for mid values  
    high = "#80FF80",     # Stronger green for high values
    midpoint = 30,        # Adjust based on your data range
    name = "Count"
  ) +
  labs(
    title = "Genus Overlap Between Rapamycin Studies",
    subtitle = "Total unique genera identified across three independent studies",
    caption = paste("Total genera: Study2 =", length(study2_genera), 
                    "| Study5 =", length(study5_genera), 
                    "| Study7 =", length(study7_genera))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 10, hjust = 0.5, margin = margin(t = 15)),
    legend.position = "bottom",
    plot.margin = margin(25, 50, 25, 50)
  )

print(plot.Venndiagram_allstudys_v3)

# Show overlap statistics
all_genera <- unique(c(study2_genera, study5_genera, study7_genera))
shared_all <- Reduce(intersect, list(study2_genera, study5_genera, study7_genera))

cat("\n=== OVERLAP STATISTICS ===\n")
cat("Total unique genera across all studies:", length(all_genera), "\n")
cat("Genera present in all 3 studies:", length(shared_all), "\n")
cat("Percentage of genera in all studies:", round(100 * length(shared_all) / length(all_genera), 1), "%\n")

cat("\nFirst 10 shared genera:", paste(head(shared_all, 10), collapse = ", "), "\n")

#### 6.Gene Collapse ####
# 4d) True genus‐level collapse, clean & rename OTUs, then merge

# helper to strip prefixes and trim
clean_genus <- function(x) {
  x <- gsub("^[kpcofgs]__", "", x, ignore.case=TRUE)
  trimws(x)
}

# 1) Collapse each study to genus, clean names, overwrite taxa_names()
ps_genus <- lapply(standardized_ps_list, function(ps) {
  ps <- subset_taxa(ps, !Genus %in% c("Chloroplast", "Mitochondria"))
  ps_g <- tax_glom(ps, taxrank="Genus", NArm=FALSE)
  tt   <- as.data.frame(tax_table(ps_g))
  
  if ("Genus" %in% colnames(tt)) {
    # Clean and standardize names
    tt$Genus <- clean_genus(tt$Genus) %>%
      gsub("Tuzzerella", "Tyzzerella", .) %>%  # Tuzzerella -> Tyzzerella
      gsub("Erysipelatoclostridiaceae", "Erysipelatoclostridium", .) %>%  # -> Erysipelatoclostridium
      gsub("^Uncultured.*$", "Uncultured", ., ignore.case = TRUE) %>%
      gsub("^Unclassified.*$", "Unclassified", ., ignore.case = TRUE)
    
    
    tax_table(ps_g) <- tax_table(as.matrix(tt))
    
    # Use make.unique BUT then we'll fix it after merging
    taxa_names(ps_g) <- make.unique(as.character(tt$Genus))
  }
  ps_g
})

# 2) Drop trees
ps_no_tree <- lapply(ps_genus, function(ps) {
  phyloseq(otu_table(ps), tax_table(ps), sample_data(ps))
})

# 3) Merge all studies into one genus‐level phyloseq
merged_genus <- do.call(merge_phyloseq, ps_no_tree)
cat("Before cleaning:", ntaxa(merged_genus), "genera\n")
print(taxa_names(merged_genus))

# 4) DIRECT APPROACH: Aggregate OTU counts by genus name manually
# Get OTU table and tax table
otu_mat <- as.matrix(otu_table(merged_genus))
tax_df <- as.data.frame(tax_table(merged_genus))

# Clean the genus names in tax_table (removing .1, .2, etc. suffixes)
tax_df$Genus_clean <- tax_df$Genus %>%
  gsub("^Uncultured(\\.\\d+)?$", "Uncultured", .) %>%
  gsub("^Unclassified(\\.\\d+)?$", "Unclassified", .) %>%
  gsub("^Tuzzerella(\\.\\d+)?$", "Tyzzerella", .) %>%
  gsub("^Tyzzerella(\\.\\d+)?$", "Tyzzerella", .) %>%
  gsub("^Erysipelatoclostridiaceae(\\.\\d+)?$", "Erysipelatoclostridium", .) %>%
  gsub("^Erysipelatoclostridium(\\.\\d+)?$", "Erysipelatoclostridium", .)
# Convert to data frame with genus info
otu_df <- as.data.frame(otu_mat)
otu_df$Genus_clean <- tax_df$Genus_clean

# Aggregate by genus
otu_collapsed <- otu_df %>%
  group_by(Genus_clean) %>%
  summarise(across(everything(), sum), .groups = 'drop')

# Convert back to matrix format
otu_collapsed_mat <- as.matrix(otu_collapsed[, -1])  # Remove genus column
rownames(otu_collapsed_mat) <- otu_collapsed$Genus_clean

# Create new tax_table with unique genera
unique_genera <- unique(tax_df$Genus_clean)
tax_collapsed <- data.frame(
  Kingdom = "Bacteria",  # or whatever your kingdom is
  Phylum = NA,
  Class = NA,
  Order = NA,
  Family = NA,
  Genus = unique_genera,
  row.names = unique_genera
)

# Create new phyloseq object
merged_genus_collapsed <- phyloseq(
  otu_table(otu_collapsed_mat, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_collapsed)),
  sample_data(merged_genus)
)

cat("After collapsing:", ntaxa(merged_genus_collapsed), "genera\n")

# 4) Normalization 
#see https://bioconductor.statistik.tu-dortmund.de/packages/3.10/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html
merged_genus_normal <- merged_genus_collapsed %>%                     # genus-level object
  transform_sample_counts(function(x) x / sum(x)) %>%  # TSS
  filter_taxa(function(x) sum(x > 1e-5) >= 2,          # ≥ 1 × 10⁻⁵ in ≥ 2 samples
              prune = TRUE)

cat("Merged genera:", ntaxa(merged_genus_normal),
    "across", nsamples(merged_genus_normal), "samples\n")
#Merged genera: 156 across 54 samples

# Print the taxa names (genus names) from merged_genus
cat("=== GENERA IN MERGED_GENUS ===\n")
print(taxa_names(merged_genus_normal))


#### 7.Heatmap of all Studies ####

# Get abundance matrix (already TSS normalized and filtered!)
abundance_matrix <- as.matrix(otu_table(merged_genus_normal))
if(!taxa_are_rows(merged_genus_normal)) {
  abundance_matrix <- t(abundance_matrix)
}

# No additional filtering needed here
cat("Matrix dimensions:", nrow(abundance_matrix), "genera x", ncol(abundance_matrix), "samples\n")

# Get metadata 
meta <- data.frame(sample_data(merged_genus_normal))

# Create group labels
meta$Group <- paste0(meta$Study, "_", meta$StandardTreatment, "_", meta$AdminRoute)
cat("Groups created:", paste(unique(meta$Group), collapse = ", "), "\n")

# Calculate mean abundance per group
group_means_matrix <- sapply(unique(meta$Group), function(group) {
  samples_in_group <- which(meta$Group == group)
  if(length(samples_in_group) == 1) {
    abundance_matrix[, samples_in_group]
  } else {
    rowMeans(abundance_matrix[, samples_in_group, drop = FALSE])
  }
})

# Check dimensions
cat("Group means matrix:", nrow(group_means_matrix), "genera x", ncol(group_means_matrix), "groups\n")

# Convert to long format
heatmap_data <- group_means_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "StudyTreatment", values_to = "Abundance")

# Optional: Filter out very low abundance genera for cleaner visualization
# Uncomment next lines if you want to show only more abundant genera
abundance_threshold <- 0.002  # 0.1% minimum mean abundance
abundant_genera <- heatmap_data %>%
group_by(Genus) %>%
summarise(max_abundance = max(Abundance)) %>%
filter(max_abundance >= abundance_threshold) %>%
pull(Genus)
heatmap_data <- heatmap_data %>% filter(Genus %in% abundant_genera)

heatmap_data_filtered <- heatmap_data %>%
  filter(Genus != "Muribaculaceae") %>%
  filter(!grepl("Unclassified|unclassified|Uncultured|uncultured", Genus, ignore.case = TRUE))

# Change "Food" to "Diet" in StudyTreatment names
heatmap_data_filtered$StudyTreatment <- gsub("_Food", "_Diet", heatmap_data_filtered$StudyTreatment)


# Calculate max excluding Unclassified
max_abundance <- heatmap_data_filtered %>%
  filter(!grepl("Unclassified|unclassified", Genus, ignore.case = TRUE)) %>%
  pull(Abundance) %>%
  max(na.rm = TRUE)

cat("New max abundance (excluding Unclassified):", round(max_abundance * 100, 1), "%\n")

# Reorder StudyTreatment factor
heatmap_data_filtered$StudyTreatment <- factor(
  heatmap_data_filtered$StudyTreatment,
  levels = c("Study7_Rapamycin_IP",
             "Study7_Control_IP",
             "Study5_Rapamycin_IP",
             "Study5_Control_IP",
             "Study2_Rapamycin_IP",
             "Study2_Control_IP",      
             "Study2_Rapamycin_Diet",    # Changed from Food
             "Study2_Control_Diet")      # Changed from Food
)

# Create heatmap with filtered data
plot.heatmap_normalized <- ggplot(heatmap_data_filtered, aes(x = Genus, y = StudyTreatment)) +
  geom_tile(aes(fill = Abundance), colour = NA, linewidth = 0) +
  xlab('') + 
  ylab('') + 
  labs(
    title = paste("Genus Abundance (", length(unique(heatmap_data_filtered$Genus)), "genera)"),
    subtitle = "Normalized and filtered data - grouped by study and treatment",
    fill = "Relative\nAbundance"
  ) +
  scale_fill_gradientn(
    name = "Relative\nAbundance",
    colors = c("white",      # 0%
               "#C6DBEF",    # Light blue
               "#6BAED6",    # Medium blue  
               "#2171B5",    # Dark blue
               "#08306B",    # Very dark blue
               "#FEE0D2",    # Very light pink
               "#FC9272",    # Light-medium red
               "#DE2D26",    # Red
               "#A50F15"),   # Dark red
    values = c(0, 0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.20, 1),  # Explicit breaks!
    limits = c(0, max_abundance),
    breaks = seq(0, max_abundance, by = 0.05),  # Legend breaks every 5%
    labels = function(x) paste0(round(x * 100, 1), "%")
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
    axis.text.y = element_text(hjust = 1, size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(plot.heatmap_normalized)

# Convert to binary presence/absence (exactly like your template)
tempdata <- ifelse(group_means_matrix > 0, 1, 0)

tempdata <- tempdata[!grepl("Unclassified|unclassified|Uncultured|uncultured", rownames(tempdata), ignore.case = TRUE), ]

# Clean row names (copy from your template)
rownames(tempdata) <- gsub(".*ales_", "", rownames(tempdata))
rownames(tempdata) <- gsub("aceae_", "aceae", rownames(tempdata))
rownames(tempdata) <- gsub("aceaeunc", "aceae_unc", rownames(tempdata))
rownames(tempdata) <- gsub("\\[", "", rownames(tempdata))
rownames(tempdata) <- gsub("\\]", "", rownames(tempdata))
rownames(tempdata) <- gsub("_", " ", rownames(tempdata))

# Convert to long format (exactly like your template)
tempdata2 <- tempdata %>% 
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-c(Genus), names_to = "application_method", values_to = "value")

# Change "Food" to "Diet" in application_method names
tempdata2$application_method <- gsub("_Food", "_Diet", tempdata2$application_method)

# EXACT SAME GGPLOT CODE AS YOUR TEMPLATE - just different title
tempplot <- ggplot(tempdata2, aes(Genus, application_method)) +
  geom_tile(aes(fill = value), colour = "black") +
  coord_equal() + 
  xlab('') + 
  ylab('Study Treatment') + 
  labs(
    title = "Genus Presence/Absence Across Studies",
    fill = "Presence"
  ) +
  scale_fill_gradient(low = "white", high = "gray") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        axis.text.y = element_text(hjust = 0),
        legend.position = "none")

print(tempplot)


# Updated plot - remove background, adjust proportions, tilt labels
# Updated plot - with reversed genus order
tempplot <- ggplot(tempdata2, aes(x = application_method, y = Genus)) +
  geom_tile(aes(fill = value), colour = "black") +
  xlab('') + 
  ylab('') + 
  scale_fill_gradient(low = "white", high = "gray") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(position = "right", expand = c(0, 0), 
                   limits = rev(levels(factor(tempdata2$Genus)))) +  # REVERSED
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0, size = 8),
    axis.text.y = element_text(hjust = 0, size = 8),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )
print(tempplot)

# Save with dimensions in CENTIMETERS
ggsave("heatmap_presence_absence.png", 
       plot = tempplot,
       width = 15,      # Width in cm
       height = 50,     # Height in cm
       dpi = 300,       # 300 dots per inch (publication quality)
       units = "cm")    # CENTIMETERS

### 7.5 Split study 2 up #####

meta <- data.frame(sample_data(merged_genus_normal))
meta$StudyID_Split <- ifelse(
  meta$StudyID == "Study2" & meta$AdminRoute == "Food", "Study2_Food",
  ifelse(
    meta$StudyID == "Study2" & meta$AdminRoute == "IP", "Study2_IP",
    meta$StudyID
  )
)
meta$StudyID_Split <- factor(meta$StudyID_Split)

# Update the phyloseq object with split StudyID
sample_data(merged_genus_normal)$StudyID_Split <- meta$StudyID_Split

print(table(meta$StudyID_Split))
print(table(meta$StudyID))

#### 8.BATCH CORRECTION  (MMUPHin_Combat) ####

# 1) Make sure you’re using the merged genus‐level phyloseq
#    (from your harmonise_taxa step)
rapa_genus_ps <- merged_genus_normal

# 2) Extract the genus × sample matrix
genus_mat <- as(otu_table(rapa_genus_ps), "matrix")
if (!taxa_are_rows(rapa_genus_ps)) genus_mat <- t(genus_mat)

# 3) Rebuild your metadata
rapa_meta <- as(sample_data(rapa_genus_ps), "data.frame")
rapa_meta$StudyID_Split          <- factor(trimws(rapa_meta$StudyID_Split))
rapa_meta$SampleID         <- rownames(rapa_meta)
rapa_meta$StandardTreatment <- factor(trimws(rapa_meta$StandardTreatment))
rapa_meta$AdminRoute <- factor(trimws(rapa_meta$AdminRoute))

# 3) Define your three covariates
covar_cols <- c("AdminRoute")

# 4) Run adjust_batch on the genus matrix
fit_genus <- MMUPHin::adjust_batch(
  feature_abd = genus_mat,
  batch       = "StudyID_Split",           
  covariates  = NULL, 
  data        = rapa_meta,
  control     = list(method="combat", zero_correct=TRUE, verbose=FALSE, mean.only=TRUE)
)

# 5) Grab the corrected genus counts
genus_mat_adj <- fit_genus$feature_abd_adj

# 6) Replace the OTU table in your phyloseq
otu_table(rapa_genus_ps) <- otu_table(genus_mat_adj, taxa_are_rows = TRUE)
batch_corrected_ps <- rapa_genus_ps

cat("Batch correction on genus-level complete: ",
    ntaxa(batch_corrected_ps), "genera × ",
    nsamples(batch_corrected_ps), "samples\n")


#### 9.Permanova Test ######

# 1) Get the matrices 
genus_before <- genus_mat  
genus_after  <- fit_genus$feature_abd_adj

# 2) Calculate Bray-Curtis distances 
D_before <- vegdist(t(genus_before), method="bray")
D_after  <- vegdist(t(genus_after), method="bray")

# 3) Make sure we're using the correct metadata
# Use rapa_meta instead of meta
meta_df <- rapa_meta  

# 4) VERIFY SAMPLE ORDER MATCHES METADATA
sample_names_matrix <- colnames(genus_before)
sample_names_meta <- rownames(rapa_meta)

if(!all(sample_names_matrix == sample_names_meta)) {
  warning("Sample order mismatch between matrix and metadata!")
  cat("Matrix samples:", head(sample_names_matrix), "...\n")
  cat("Metadata samples:", head(sample_names_meta), "...\n")
} else {
  cat("Sample order matches between matrix and metadata\n")
}

# Test study effect
set.seed(123)
study_before <- adonis2(D_before ~ StudyID_Split, data = rapa_meta,permutations = 9999)
print(study_before)
set.seed(123)
study_after <- adonis2(D_after ~ StudyID_Split, data = rapa_meta,permutations = 9999)
print(study_after)

# Test treatment effect 
set.seed(123)
treatment_before <- adonis2(D_before ~ StandardTreatment, data = rapa_meta,permutations = 9999)
set.seed(123)
treatment_after <- adonis2(D_after ~ StandardTreatment, data = rapa_meta,permutations = 9999)

# Test combined effect
set.seed(123)
combined_before <- adonis2(D_before ~ StudyID_Split + StandardTreatment, data = rapa_meta,permutations = 9999,by = "terms")
set.seed(123)
combined_after <- adonis2(D_after ~ StudyID_Split + StandardTreatment, data = rapa_meta,permutations = 9999,by = "terms")


extract_stats <- function(adonis_result, row = 1) {
  # Check if row exists and has valid data
  if(row > nrow(adonis_result) || is.na(adonis_result$F[row])) {
    return(data.frame(
      Df = NA, 
      R2 = NA, 
      R2_percent = "ERROR",  # Same type as normal case
      F_stat = NA, 
      P_val = NA, 
      Sig = "ERROR"
    ))
  }
  
  r2 <- adonis_result$R2[row]
  p <- adonis_result$`Pr(>F)`[row]
  
  data.frame(
    Df = adonis_result$Df[row],
    R2 = round(r2, 4),
    R2_percent = paste0(round(r2*100, 1), "%"),
    F_stat = round(adonis_result$F[row], 2),
    P_val = round(p, 5),
    Sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "†", TRUE ~ "ns")
  )
}

# Then check your combined models first:
cat("Combined before rows:", nrow(combined_before), "\n")
cat("Combined after rows:", nrow(combined_after), "\n")
print(combined_before)
print(combined_after)


# Then run the results table:
results_table <- rbind(
  cbind(Effect = "Study", Timing = "Before", extract_stats(study_before, 1)),
  cbind(Effect = "Study", Timing = "After", extract_stats(study_after, 1)),
  cbind(Effect = "Treatment", Timing = "Before", extract_stats(treatment_before, 1)),
  cbind(Effect = "Treatment", Timing = "After", extract_stats(treatment_after, 1)),
  cbind(Effect = "Treatment (adj)", Timing = "Before", extract_stats(combined_before, 2)),
  cbind(Effect = "Treatment (adj)", Timing = "After", extract_stats(combined_after, 2))
)


# 7) Display results
cat("\n=== PERMANOVA RESULTS ===\n")
print(results_table)

table.results_Permanova <- results_table

#### 9B. PERMDISP TESTING (Dispersion Homogeneity) ####

# Function to extract PERMDISP results
get_permdisp <- function(betadisper_result) {
  perm_test <- vegan::permutest(betadisper_result, permutations = 9999)
  data.frame(
    F = round(perm_test$tab$F[1], 3),
    p = round(perm_test$tab$`Pr(>F)`[1], 4),
    sig = ifelse(perm_test$tab$`Pr(>F)`[1] < 0.05, "⚠️", "✓")
  )
}

cat("\n=== PERMDISP (DISPERSION HOMOGENEITY) TESTING ===\n")

# 1. Test STUDY dispersion
set.seed(123)
disp_study_before <- betadisper(D_before, rapa_meta$StudyID_Split)
set.seed(123)
disp_study_after <- betadisper(D_after, rapa_meta$StudyID_Split)

cat("\nStudy dispersion - BEFORE:\n")
print(vegan::permutest(disp_study_before, permutations = 9999))
cat("\nStudy dispersion - AFTER:\n")
print(vegan::permutest(disp_study_after, permutations = 9999))

# 2. Test TREATMENT dispersion (all samples)
set.seed(123)
disp_treat_before <- betadisper(D_before, rapa_meta$StandardTreatment)
set.seed(123)
disp_treat_after <- betadisper(D_after, rapa_meta$StandardTreatment)

cat("\nTreatment dispersion - BEFORE:\n")
print(vegan::permutest(disp_treat_before, permutations = 9999))
cat("\nTreatment dispersion - AFTER:\n")
print(vegan::permutest(disp_treat_after, permutations = 9999))

# 3. Create combined PERMANOVA + PERMDISP table
combined_results <- data.frame(
  Effect = c("Study", "Study", "Treatment", "Treatment"),
  Timing = c("Before", "After", "Before", "After"),
  
  # PERMANOVA
  PERM_R2 = c(
    study_before$R2[1],
    study_after$R2[1],
    treatment_before$R2[1],
    treatment_after$R2[1]
  ),
  PERM_p = c(
    study_before$`Pr(>F)`[1],
    study_after$`Pr(>F)`[1],
    treatment_before$`Pr(>F)`[1],
    treatment_after$`Pr(>F)`[1]
  ),
  
  # PERMDISP
  DISP_F = c(
    vegan::permutest(disp_study_before, permutations = 9999)$tab$F[1],
    vegan::permutest(disp_study_after, permutations = 9999)$tab$F[1],
    vegan::permutest(disp_treat_before, permutations = 9999)$tab$F[1],
    vegan::permutest(disp_treat_after, permutations = 9999)$tab$F[1]
  ),
  DISP_p = c(
    vegan::permutest(disp_study_before, permutations = 9999)$tab$`Pr(>F)`[1],
    vegan::permutest(disp_study_after, permutations = 9999)$tab$`Pr(>F)`[1],
    vegan::permutest(disp_treat_before, permutations = 9999)$tab$`Pr(>F)`[1],
    vegan::permutest(disp_treat_after, permutations = 9999)$tab$`Pr(>F)`[1]
  )
)

# Add interpretation column
combined_results$Interpretation <- with(combined_results, 
                                        ifelse(PERM_p >= 0.05, "No difference",
                                               ifelse(DISP_p >= 0.05, "✓ Location effect",
                                                      "⚠️ Location + Dispersion")))

cat("\n=== COMBINED PERMANOVA + PERMDISP RESULTS ===\n")
print(combined_results)

table.combined_permanova_permdisp <- combined_results


#### 9C. WITHIN-STUDY ANALYSES ####

cat("\n=== WITHIN-STUDY TREATMENT EFFECTS ===\n")

# Get unique studies (4-way split)
studies <- unique(rapa_meta$StudyID_Split)
studies <- studies[!is.na(studies)]

within_study_results <- data.frame()

for (study in studies) {
  cat("\n--- Analyzing:", study, "---\n")
  
  # Subset data for this study
  study_idx <- which(rapa_meta$StudyID_Split == study)
  study_meta <- rapa_meta[study_idx, ]
  
  # Create distance matrices for this study only
  study_D_before <- as.dist(as.matrix(D_before)[study_idx, study_idx])
  study_D_after <- as.dist(as.matrix(D_after)[study_idx, study_idx])
  
  # Check if we have both treatments
  treat_table <- table(study_meta$StandardTreatment)
  if (length(treat_table) < 2) {
    cat("  Skipping - only one treatment group\n")
    next
  }
  
  # Sample sizes
  n_ctrl <- sum(study_meta$StandardTreatment == "Control")
  n_rapa <- sum(study_meta$StandardTreatment == "Rapamycin")
  balanced <- ifelse(n_ctrl == n_rapa, "Yes", "No")
  
  cat("  Sample sizes: Control =", n_ctrl, ", Rapamycin =", n_rapa, "\n")
  
  # PERMANOVA within this study - BEFORE
  set.seed(123)
  perm_before <- adonis2(study_D_before ~ StandardTreatment, 
                         data = study_meta, 
                         permutations = 9999)
  
  # PERMANOVA within this study - AFTER
  set.seed(123)
  perm_after <- adonis2(study_D_after ~ StandardTreatment, 
                        data = study_meta, 
                        permutations = 9999)
  
  # PERMDISP within this study - BEFORE
  set.seed(123)
  disp_before <- betadisper(study_D_before, study_meta$StandardTreatment)
  disp_test_before <- permutest(disp_before, permutations = 9999)
  
  # PERMDISP within this study - AFTER
  set.seed(123)
  disp_after <- betadisper(study_D_after, study_meta$StandardTreatment)
  disp_test_after <- permutest(disp_after, permutations = 9999)
  
  # Store results
  within_study_results <- rbind(within_study_results, 
                                data.frame(
                                  Study = study,
                                  Timing = "Before",
                                  n_control = n_ctrl,
                                  n_rapamycin = n_rapa,
                                  Balanced = balanced,
                                  PERM_R2 = round(perm_before$R2[1], 4),
                                  PERM_p = round(perm_before$`Pr(>F)`[1], 4),
                                  DISP_F = round(disp_test_before$tab$F[1], 3),
                                  DISP_p = round(disp_test_before$tab$`Pr(>F)`[1], 4),
                                  Interpretation = ifelse(perm_before$`Pr(>F)`[1] >= 0.05, "NS",
                                                          ifelse(disp_test_before$tab$`Pr(>F)`[1] >= 0.05, "✓",
                                                                 "⚠️"))
                                ))
  
  within_study_results <- rbind(within_study_results, 
                                data.frame(
                                  Study = study,
                                  Timing = "After",
                                  n_control = n_ctrl,
                                  n_rapamycin = n_rapa,
                                  Balanced = balanced,
                                  PERM_R2 = round(perm_after$R2[1], 4),
                                  PERM_p = round(perm_after$`Pr(>F)`[1], 4),
                                  DISP_F = round(disp_test_after$tab$F[1], 3),
                                  DISP_p = round(disp_test_after$tab$`Pr(>F)`[1], 4),
                                  Interpretation = ifelse(perm_after$`Pr(>F)`[1] >= 0.05, "NS",
                                                          ifelse(disp_test_after$tab$`Pr(>F)`[1] >= 0.05, "✓",
                                                                 "⚠️"))
                                ))
}

cat("\n=== WITHIN-STUDY RESULTS SUMMARY ===\n")
print(within_study_results)

table.within_study_results <- within_study_results

# Count significant studies
sig_before <- sum(within_study_results$Timing == "Before" & 
                    within_study_results$PERM_p < 0.05)
sig_after <- sum(within_study_results$Timing == "After" & 
                   within_study_results$PERM_p < 0.05)

cat("\n=== SUMMARY ===")
cat("\nSignificant treatment effects BEFORE correction:", sig_before, "out of", 
    length(studies), "studies\n")
cat("Significant treatment effects AFTER correction:", sig_after, "out of", 
    length(studies), "studies\n")


#### 9E. SAVE ALL RESULTS ####

# Create a comprehensive results list
permanova_results_full <- list(
  combined_permanova_permdisp = combined_results,
  within_study_results = within_study_results,
  sample_sizes = data.frame(
    Study = names(table(rapa_meta$StudyID_Split)),
    Total = as.numeric(table(rapa_meta$StudyID_Split)),
    Control = as.numeric(table(rapa_meta$StudyID_Split[rapa_meta$StandardTreatment == "Control"])),
    Rapamycin = as.numeric(table(rapa_meta$StudyID_Split[rapa_meta$StandardTreatment == "Rapamycin"]))
  )
)

# Print sample size table
cat("\n=== SAMPLE SIZES BY STUDY ===\n")
print(permanova_results_full$sample_sizes)

# Save for manuscript
saveRDS(permanova_results_full, "permanova_full_results.rds")

#### 10.PCOA ####

# 2) Calculate Bray-Curtis distances 
D_before <- vegdist(t(genus_before), method="bray")
D_after  <- vegdist(t(genus_after), method="bray")

ord_before <- ordinate(merged_genus_normal, "PCoA", D_before)
ord_after <- ordinate(batch_corrected_ps, "PCoA", D_after)

# 2) Extract ordination data with sample metadata
data_before <- plot_ordination(merged_genus, ord_before)$data
data_after <- plot_ordination(batch_corrected_ps, ord_after)$data

data_before$StudyID_Split <- ifelse(
  data_before$StudyID == "Study2" & data_before$AdminRoute == "Food", "Study2_Diet",  # Changed
  ifelse(data_before$StudyID == "Study2" & data_before$AdminRoute == "IP", "Study2_IP",
         data_before$StudyID))

data_after$StudyID_Split <- ifelse(
  data_after$StudyID == "Study2" & data_after$AdminRoute == "Food", "Study2_Diet",  # Changed
  ifelse(data_after$StudyID == "Study2" & data_after$AdminRoute == "IP", "Study2_IP",
         data_after$StudyID))

# Change the label assignments
data_before$StudyID_Split[data_before$StudyID_Split == "Study2_Diet"] <- "Bitto 2016-Diet"
data_before$StudyID_Split[data_before$StudyID_Split == "Study2_IP"] <- "Bitto 2016-IP"
data_before$StudyID_Split[data_before$StudyID_Split == "Study5"] <- "Yang 2022"
data_before$StudyID_Split[data_before$StudyID_Split == "Study7"] <- "Han 2021"

data_after$StudyID_Split[data_after$StudyID_Split == "Study2_Diet"] <- "Bitto 2016-Diet"
data_after$StudyID_Split[data_after$StudyID_Split == "Study2_IP"] <- "Bitto 2016-IP"
data_after$StudyID_Split[data_after$StudyID_Split == "Study5"] <- "Yang 2022"
data_after$StudyID_Split[data_after$StudyID_Split == "Study7"] <- "Han 2021"

# Update factor levels
lvl <- c("Bitto 2016-Diet", "Bitto 2016-IP", "Yang 2022", "Han 2021")
data_before$StudyID_Split <- factor(data_before$StudyID_Split, levels = lvl)
data_after$StudyID_Split <- factor(data_after$StudyID_Split, levels = lvl)
# Lines 319-326 - Keep them separate:
data_before$StudyID_Split[data_before$StudyID_Split == "Study2_Diet"] <- "Bitto 2016 (Diet)"  # Changed
data_before$StudyID_Split[data_before$StudyID_Split == "Study2_IP"] <- "Bitto 2016 (IP)"


# 3) Create separate plotting functions for before and after
create_pcoa_plot_before <- function(data, title, ord_obj) {
  # Create combined grouping variable for ellipses
  data$Group <- paste(data$StudyID_Split, data$StandardTreatment, sep="_")
  
  # Create the plot - BEFORE (shows study differences)
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    # Points with study color and treatment shape
    geom_point(aes(color = StudyID_Split, shape = StandardTreatment), size = 3) +
    # Ellipses with study color but different alpha for treatments
    stat_ellipse(aes(group = Group, color = StudyID_Split, linetype = StandardTreatment),
                 size = 1) +
    # Filled ellipses with low alpha
    stat_ellipse(aes(group = Group, fill = StudyID_Split), 
                 geom = "polygon", alpha = 0.1, show.legend = FALSE) +
    # Custom shapes for treatments
    scale_shape_manual(values = c(16, 17)) +  # Circle for Control, Triangle for Rapamycin
    # Different line types for treatments
    scale_linetype_manual(values = c("solid", "dashed")) +
    # Add panel label
    # Add titles and labels
    labs(
      title = "a",
      color = "Study",
      shape = "Treatment",
      linetype = "Treatment"
    ) +
    theme_bw() +
    # Add variance explained to axes
    xlab(paste0("PCoA1 [", round(100 * ord_obj$values$Relative_eig[1], 1), "%]")) +
    ylab(paste0("PCoA2 [", round(100 * ord_obj$values$Relative_eig[2], 1), "%]"))
  
  return(p)
}
create_pcoa_plot_after <- function(data, title, ord_obj) {
  # Create the plot - AFTER (focuses on treatment differences only)
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    # Points with study color and treatment shape
    geom_point(aes(color = StudyID_Split, shape = StandardTreatment), size = 3) +
    # Ellipses with study color but different alpha for treatments
    # NO ellipses - just focus on the point patterns
    # Custom shapes for treatments
    scale_shape_manual(values = c(16, 17)) +  # Circle for Control, Triangle for Rapamycin
    # Add panel label
    # Add titles and labels
    labs(
      title = "b",
      shape = "Treatment"
    ) +
    theme_bw() +
    theme(legend.position = "none") + 
    # Add variance explained to axes
    xlab(paste0("PCoA1 [", round(100 * ord_obj$values$Relative_eig[1], 1), "%]")) +
    ylab(paste0("PCoA2 [", round(100 * ord_obj$values$Relative_eig[2], 1), "%]"))
  
  return(p)
}

# 4) Create the plots using different functions
plot.p_before <- create_pcoa_plot_before(data_before, "Before Batch Correction", ord_before)
plot.p_after <- create_pcoa_plot_after(data_after, "After Batch Correction", ord_after)

# 5) Combine the plots
library(patchwork)
plot.combined_plot <- plot.p_before + plot.p_after + 
  plot_layout(guides = "collect") 


# 6) Print the combined plot
print(plot.combined_plot)

# Check Han dispersion with 4-study split
han_before <- data_before[data_before$StudyID_Split == "Han 2021", ]
han_after <- data_after[data_after$StudyID_Split == "Han 2021", ]

cat("=== HAN 2021 DISPERSION (4-study split) ===\n")
cat("BEFORE batch correction:\n")
print(apply(han_before[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(han_after[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(han_before[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(han_after[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))

# Check Yang dispersion with 4-study split
Yang_before <- data_before[data_before$StudyID_Split == "Yang 2022", ]
Yang_after <- data_after[data_after$StudyID_Split == "Yang 2022", ]

cat("=== YANG 2022 DISPERSION (4-study split) ===\n")  # Fixed label
cat("BEFORE batch correction:\n")
print(apply(Yang_before[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(Yang_after[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(Yang_before[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(Yang_after[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))

##### VERSION: UNSPLIT (3 studies) - CORRECTED ####
#### 8.BATCH CORRECTION  (MMUPHin_Combat) 

# 1) Make sure you're using the merged genus‐level phyloseq
#    (from your harmonise_taxa step)
rapa_genus_ps_3study <- merged_genus_normal  # ← Added suffix

# 2) Extract the genus × sample matrix
genus_mat_3study <- as(otu_table(rapa_genus_ps_3study), "matrix")  # ← Added suffix
if (!taxa_are_rows(rapa_genus_ps_3study)) genus_mat_3study <- t(genus_mat_3study)

# 3) Rebuild your metadata
rapa_meta_3study <- as(sample_data(rapa_genus_ps_3study), "data.frame")  # ← Added suffix
rapa_meta_3study$StudyID          <- factor(trimws(rapa_meta_3study$StudyID))
rapa_meta_3study$SampleID         <- rownames(rapa_meta_3study)
rapa_meta_3study$StandardTreatment <- factor(trimws(rapa_meta_3study$StandardTreatment))
rapa_meta_3study$AdminRoute <- factor(trimws(rapa_meta_3study$AdminRoute))

# 3) Define your three covariates
covar_cols <- c("AdminRoute")

# 4) Run adjust_batch on the genus matrix
fit_genus_3study <- MMUPHin::adjust_batch(  # ← Added suffix
  feature_abd = genus_mat_3study,  # ← Added suffix
  batch       = "StudyID",         # ← CHANGED: was "StudyID_Split"
  covariates  = covar_cols,        # ← CHANGED: was NULL
  data        = rapa_meta_3study,  # ← Added suffix
  control     = list(method="combat", zero_correct=TRUE, verbose=FALSE, mean.only=TRUE)
)

# 5) Grab the corrected genus counts
genus_mat_adj_3study <- fit_genus_3study$feature_abd_adj  # ← Added suffix

# 6) Replace the OTU table in your phyloseq
otu_table(rapa_genus_ps_3study) <- otu_table(genus_mat_adj_3study, taxa_are_rows = TRUE)
batch_corrected_ps_3study <- rapa_genus_ps_3study  # ← Added suffix

cat("Batch correction on genus-level complete (3-study): ",
    ntaxa(batch_corrected_ps_3study), "genera × ",
    nsamples(batch_corrected_ps_3study), "samples\n")


#### 9.Permanova Test #

# 1) Get the matrices 
genus_before_3study <- genus_mat_3study
genus_after_3study  <- fit_genus_3study$feature_abd_adj

# 2) Calculate Bray-Curtis distances 
D_before_3study <- vegdist(t(genus_before_3study), method="bray")
D_after_3study  <- vegdist(t(genus_after_3study), method="bray")

# 3) Use the correct metadata
meta_df_3study <- rapa_meta_3study

# 4) VERIFY SAMPLE ORDER MATCHES METADATA
sample_names_matrix <- colnames(genus_before_3study)
sample_names_meta <- rownames(rapa_meta_3study)

if(!all(sample_names_matrix == sample_names_meta)) {
  warning("Sample order mismatch between matrix and metadata!")
} else {
  cat("Sample order matches between matrix and metadata\n")
}

# Test study effect - USE STUDYID (3 levels)
set.seed(123)
study_before_3study <- adonis2(D_before_3study ~ StudyID, data = rapa_meta_3study, permutations = 9999)
print(study_before_3study)
set.seed(123)
study_after_3study <- adonis2(D_after_3study ~ StudyID, data = rapa_meta_3study, permutations = 9999)
print(study_after_3study)

# Test treatment effect 
set.seed(123)
treatment_before_3study <- adonis2(D_before_3study ~ StandardTreatment, data = rapa_meta_3study, permutations = 9999)
set.seed(123)
treatment_after_3study <- adonis2(D_after_3study ~ StandardTreatment, data = rapa_meta_3study, permutations = 9999)

# Test combined effect
set.seed(123)
combined_before_3study <- adonis2(D_before_3study ~ StudyID + StandardTreatment, 
                                  data = rapa_meta_3study, permutations = 9999, by = "terms")
set.seed(123)
combined_after_3study <- adonis2(D_after_3study ~ StudyID + StandardTreatment, 
                                 data = rapa_meta_3study, permutations = 9999, by = "terms")


extract_stats <- function(adonis_result, row = 1) {
  if(row > nrow(adonis_result) || is.na(adonis_result$F[row])) {
    return(data.frame(
      Df = NA, 
      R2 = NA, 
      R2_percent = "ERROR",
      F_stat = NA, 
      P_val = NA, 
      Sig = "ERROR"
    ))
  }
  
  r2 <- adonis_result$R2[row]
  p <- adonis_result$`Pr(>F)`[row]
  
  data.frame(
    Df = adonis_result$Df[row],
    R2 = round(r2, 4),
    R2_percent = paste0(round(r2*100, 1), "%"),
    F_stat = round(adonis_result$F[row], 2),
    P_val = round(p, 5),
    Sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "†", TRUE ~ "ns")
  )
}

# Check combined models
cat("Combined before rows:", nrow(combined_before_3study), "\n")
cat("Combined after rows:", nrow(combined_after_3study), "\n")
print(combined_before_3study)
print(combined_after_3study)

# Results table
results_table_3study <- rbind(
  cbind(Effect = "Study", Timing = "Before", extract_stats(study_before_3study, 1)),
  cbind(Effect = "Study", Timing = "After", extract_stats(study_after_3study, 1)),
  cbind(Effect = "Treatment", Timing = "Before", extract_stats(treatment_before_3study, 1)),
  cbind(Effect = "Treatment", Timing = "After", extract_stats(treatment_after_3study, 1)),
  cbind(Effect = "Treatment (adj)", Timing = "Before", extract_stats(combined_before_3study, 2)),
  cbind(Effect = "Treatment (adj)", Timing = "After", extract_stats(combined_after_3study, 2))
)

cat("\n=== PERMANOVA RESULTS (3-STUDY) ===\n")
print(results_table_3study)

table.results_Permanova_3study <- results_table_3study

#### 10.PCOA - VISUALIZE AS 3 STUDIES ##

ord_before_3study <- ordinate(merged_genus_normal, "PCoA", D_before_3study)
ord_after_3study <- ordinate(batch_corrected_ps_3study, "PCoA", D_after_3study)

# Extract ordination data
data_before_3study <- plot_ordination(merged_genus_normal, ord_before_3study)$data
data_after_3study <- plot_ordination(batch_corrected_ps_3study, ord_after_3study)$data

# Create proper study names (3 levels, NOT 4!)
data_before_3study$Study_Name <- ifelse(
  data_before_3study$StudyID == "Study2", "Bitto 2016",
  ifelse(data_before_3study$StudyID == "Study5", "Yang 2022", "Han 2021"))

data_after_3study$Study_Name <- ifelse(
  data_after_3study$StudyID == "Study2", "Bitto 2016",
  ifelse(data_after_3study$StudyID == "Study5", "Yang 2022", "Han 2021"))

# Factor with 3 levels
lvl <- c("Bitto 2016", "Yang 2022", "Han 2021")
data_before_3study$Study_Name <- factor(data_before_3study$Study_Name, levels = lvl)
data_after_3study$Study_Name <- factor(data_after_3study$Study_Name, levels = lvl)

# Plotting functions - 3 STUDIES
create_pcoa_plot_before_3study <- function(data, title, ord_obj) {
  # Create combined grouping variable for ellipses
  data$Group <- paste(data$Study_Name, data$StandardTreatment, sep="_")
  
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = Study_Name, shape = StandardTreatment), size = 3) +
    stat_ellipse(aes(group = Group, color = Study_Name, linetype = StandardTreatment), size = 1) +
    stat_ellipse(aes(group = Group, fill = Study_Name), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(
      title = "A",
      color = "Study",
      shape = "Treatment",
      linetype = "Treatment"
    ) +
    theme_bw() +
    xlab(paste0("PCoA1 [", round(100 * ord_obj$values$Relative_eig[1], 1), "%]")) +
    ylab(paste0("PCoA2 [", round(100 * ord_obj$values$Relative_eig[2], 1), "%]"))
  
  return(p)
}

create_pcoa_plot_after_3study <- function(data, title, ord_obj) {
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = Study_Name, shape = StandardTreatment), size = 3) +
    scale_shape_manual(values = c(16, 17)) +
    labs(
      title = "B",
      shape = "Treatment"
    ) +
    theme_bw() +
    theme(legend.position = "none") + 
    xlab(paste0("PCoA1 [", round(100 * ord_obj$values$Relative_eig[1], 1), "%]")) +
    ylab(paste0("PCoA2 [", round(100 * ord_obj$values$Relative_eig[2], 1), "%]"))
  
  return(p)
}

# Create plots
plot.p_before_3study <- create_pcoa_plot_before_3study(data_before_3study, "Before", ord_before_3study)
plot.p_after_3study <- create_pcoa_plot_after_3study(data_after_3study, "After", ord_after_3study)

# Combine
library(patchwork)
plot.combined_plot_3study <- plot.p_before_3study + plot.p_after_3study + 
  plot_layout(guides = "collect")

print(plot.combined_plot_3study)

#### Check dispersion - USE STUDY_NAME (3 studies) #

han_before_3study <- data_before_3study[data_before_3study$Study_Name == "Han 2021", ]
han_after_3study <- data_after_3study[data_after_3study$Study_Name == "Han 2021", ]

cat("=== HAN 2021 DISPERSION (3-study) ===\n")
cat("BEFORE batch correction:\n")
print(apply(han_before_3study[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(han_after_3study[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(han_before_3study[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(han_after_3study[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))

# Yang dispersion
yang_before_3study <- data_before_3study[data_before_3study$Study_Name == "Yang 2022", ]
yang_after_3study <- data_after_3study[data_after_3study$Study_Name == "Yang 2022", ]

cat("=== YANG 2022 DISPERSION (3-study) ===\n")
cat("BEFORE batch correction:\n")
print(apply(yang_before_3study[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(yang_after_3study[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(yang_before_3study[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(yang_after_3study[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))

# Bitto dispersion (all samples together)
bitto_before_3study <- data_before_3study[data_before_3study$Study_Name == "Bitto 2016", ]
bitto_after_3study <- data_after_3study[data_after_3study$Study_Name == "Bitto 2016", ]

cat("=== BITTO 2016 DISPERSION (3-study) ===\n")
cat("BEFORE batch correction:\n")
print(apply(bitto_before_3study[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(bitto_after_3study[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(bitto_before_3study[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(bitto_after_3study[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))


#### VERSION: IP-ONLY (3 IP studies, Study2_Food filtered out) ####
#### 8.BATCH CORRECTION  (MMUPHin_Combat) ##

# 1) Start with merged genus-level phyloseq
rapa_genus_ps_IPonly <- merged_genus_normal

# 2) Extract matrix
genus_mat_IPonly <- as(otu_table(rapa_genus_ps_IPonly), "matrix")
if (!taxa_are_rows(rapa_genus_ps_IPonly)) genus_mat_IPonly <- t(genus_mat_IPonly)

# 3) Build metadata with StudyID_Split
rapa_meta_IPonly <- as(sample_data(rapa_genus_ps_IPonly), "data.frame")
rapa_meta_IPonly$StudyID          <- factor(trimws(rapa_meta_IPonly$StudyID))
rapa_meta_IPonly$SampleID         <- rownames(rapa_meta_IPonly)
rapa_meta_IPonly$StandardTreatment <- factor(trimws(rapa_meta_IPonly$StandardTreatment))
rapa_meta_IPonly$AdminRoute <- factor(trimws(rapa_meta_IPonly$AdminRoute))

# Create StudyID_Split
rapa_meta_IPonly$StudyID_Split <- ifelse(
  rapa_meta_IPonly$StudyID == "Study2" & rapa_meta_IPonly$AdminRoute == "Food", "Study2_Food",
  ifelse(rapa_meta_IPonly$StudyID == "Study2" & rapa_meta_IPonly$AdminRoute == "IP", "Study2_IP",
         rapa_meta_IPonly$StudyID))
rapa_meta_IPonly$StudyID_Split <- factor(rapa_meta_IPonly$StudyID_Split)

# **KEY STEP: FILTER OUT Study2_Food BEFORE batch correction**
ip_samples <- rapa_meta_IPonly$AdminRoute == "IP"
rapa_meta_IPonly <- rapa_meta_IPonly[ip_samples, ]
genus_mat_IPonly <- genus_mat_IPonly[, ip_samples]

# Verify we have 3 study levels now
cat("Number of study arms after filtering:", length(unique(rapa_meta_IPonly$StudyID_Split)), "\n")
cat("Study levels:", unique(rapa_meta_IPonly$StudyID_Split), "\n")
# Should show: Study2_IP, Study5, Study7 (3 levels)

# Re-factor to drop unused level (Study2_Food)
rapa_meta_IPonly$StudyID_Split <- droplevels(rapa_meta_IPonly$StudyID_Split)

# 4) Run adjust_batch on 3 IP studies
fit_genus_IPonly <- MMUPHin::adjust_batch(
  feature_abd = genus_mat_IPonly,
  batch       = "StudyID_Split",   # 3 levels: Study2_IP, Study5, Study7
  covariates  = NULL,              # No AdminRoute (all IP now)
  data        = rapa_meta_IPonly,
  control     = list(method="combat", zero_correct=TRUE, verbose=FALSE, mean.only=TRUE)
)

# 5) Grab corrected counts
genus_mat_adj_IPonly <- fit_genus_IPonly$feature_abd_adj

# 6) Create phyloseq with IP-only batch corrected data
# First subset the phyloseq to IP samples
sample_data(rapa_genus_ps_IPonly) <- sample_data(rapa_meta_IPonly)
rapa_genus_ps_IPonly <- prune_samples(rownames(rapa_meta_IPonly), rapa_genus_ps_IPonly)
otu_table(rapa_genus_ps_IPonly) <- otu_table(genus_mat_adj_IPonly, taxa_are_rows = TRUE)
batch_corrected_ps_IPonly <- rapa_genus_ps_IPonly

cat("Batch correction (IP-only, 3 studies) complete: ",
    ntaxa(batch_corrected_ps_IPonly), "genera × ",
    nsamples(batch_corrected_ps_IPonly), "samples\n")

#### 9.Permanova Test - IP-ONLY ##

# 1) Get matrices
genus_before_IPonly <- genus_mat_IPonly
genus_after_IPonly  <- fit_genus_IPonly$feature_abd_adj

# 2) Calculate Bray-Curtis distances
D_before_IPonly <- vegdist(t(genus_before_IPonly), method="bray")
D_after_IPonly  <- vegdist(t(genus_after_IPonly), method="bray")

# 3) Metadata
meta_df_IPonly <- rapa_meta_IPonly

# 4) Verify sample order
sample_names_matrix <- colnames(genus_before_IPonly)
sample_names_meta <- rownames(rapa_meta_IPonly)

if(!all(sample_names_matrix == sample_names_meta)) {
  warning("Sample order mismatch!")
} else {
  cat("Sample order matches\n")
}

# Test study effect
set.seed(123)
study_before_IPonly <- adonis2(D_before_IPonly ~ StudyID_Split, data = rapa_meta_IPonly, permutations = 9999)
print(study_before_IPonly)
set.seed(123)
study_after_IPonly <- adonis2(D_after_IPonly ~ StudyID_Split, data = rapa_meta_IPonly, permutations = 9999)
print(study_after_IPonly)

# Test treatment effect
set.seed(123)
treatment_before_IPonly <- adonis2(D_before_IPonly ~ StandardTreatment, data = rapa_meta_IPonly, permutations = 9999)
set.seed(123)
treatment_after_IPonly <- adonis2(D_after_IPonly ~ StandardTreatment, data = rapa_meta_IPonly, permutations = 9999)

# Test combined
set.seed(123)
combined_before_IPonly <- adonis2(D_before_IPonly ~ StudyID_Split + StandardTreatment, 
                                  data = rapa_meta_IPonly, permutations = 9999, by = "terms")
set.seed(123)
combined_after_IPonly <- adonis2(D_after_IPonly ~ StudyID_Split + StandardTreatment, 
                                 data = rapa_meta_IPonly, permutations = 9999, by = "terms")

extract_stats <- function(adonis_result, row = 1) {
  if(row > nrow(adonis_result) || is.na(adonis_result$F[row])) {
    return(data.frame(
      Df = NA, R2 = NA, R2_percent = "ERROR",
      F_stat = NA, P_val = NA, Sig = "ERROR"
    ))
  }
  r2 <- adonis_result$R2[row]
  p <- adonis_result$`Pr(>F)`[row]
  data.frame(
    Df = adonis_result$Df[row],
    R2 = round(r2, 4),
    R2_percent = paste0(round(r2*100, 1), "%"),
    F_stat = round(adonis_result$F[row], 2),
    P_val = round(p, 5),
    Sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "†", TRUE ~ "ns")
  )
}

# Results table
results_table_IPonly <- rbind(
  cbind(Effect = "Study", Timing = "Before", extract_stats(study_before_IPonly, 1)),
  cbind(Effect = "Study", Timing = "After", extract_stats(study_after_IPonly, 1)),
  cbind(Effect = "Treatment", Timing = "Before", extract_stats(treatment_before_IPonly, 1)),
  cbind(Effect = "Treatment", Timing = "After", extract_stats(treatment_after_IPonly, 1)),
  cbind(Effect = "Treatment (adj)", Timing = "Before", extract_stats(combined_before_IPonly, 2)),
  cbind(Effect = "Treatment (adj)", Timing = "After", extract_stats(combined_after_IPonly, 2))
)

cat("\n=== PERMANOVA RESULTS (IP-ONLY) ===\n")
print(results_table_IPonly)

table.results_Permanova_IPonly <- results_table_IPonly

#### 10.PCOA - IP-ONLY ##

# Create phyloseq for BEFORE (IP samples only)
merged_genus_IPonly_before <- prune_samples(rownames(rapa_meta_IPonly), merged_genus_normal)

ord_before_IPonly <- ordinate(merged_genus_IPonly_before, "PCoA", D_before_IPonly)
ord_after_IPonly <- ordinate(batch_corrected_ps_IPonly, "PCoA", D_after_IPonly)

# Extract ordination data
data_before_IPonly <- plot_ordination(merged_genus_IPonly_before, ord_before_IPonly)$data
data_after_IPonly <- plot_ordination(batch_corrected_ps_IPonly, ord_after_IPonly)$data

# Create study names for plotting (3 IP studies)
data_before_IPonly$StudyID_Split <- ifelse(
  data_before_IPonly$StudyID == "Study2", "Bitto 2016 (IP)",
  ifelse(data_before_IPonly$StudyID == "Study5", "Yang 2022", "Han 2021"))

data_after_IPonly$StudyID_Split <- ifelse(
  data_after_IPonly$StudyID == "Study2", "Bitto 2016 (IP)",
  ifelse(data_after_IPonly$StudyID == "Study5", "Yang 2022", "Han 2021"))

# Factor with 3 levels
lvl <- c("Bitto 2016 (IP)", "Yang 2022", "Han 2021")
data_before_IPonly$StudyID_Split <- factor(data_before_IPonly$StudyID_Split, levels = lvl)
data_after_IPonly$StudyID_Split <- factor(data_after_IPonly$StudyID_Split, levels = lvl)

# Plotting functions
create_pcoa_plot_before_IPonly <- function(data, title, ord_obj) {
  data$Group <- paste(data$StudyID_Split, data$StandardTreatment, sep="_")
  
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = StudyID_Split, shape = StandardTreatment), size = 3) +
    stat_ellipse(aes(group = Group, color = StudyID_Split, linetype = StandardTreatment), size = 1) +
    stat_ellipse(aes(group = Group, fill = StudyID_Split), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(title = "A", color = "Study", shape = "Treatment", linetype = "Treatment") +
    theme_bw() +
    xlab(paste0("PCoA1 [", round(100 * ord_obj$values$Relative_eig[1], 1), "%]")) +
    ylab(paste0("PCoA2 [", round(100 * ord_obj$values$Relative_eig[2], 1), "%]"))
  return(p)
}

create_pcoa_plot_after_IPonly <- function(data, title, ord_obj) {
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = StudyID_Split, shape = StandardTreatment), size = 3) +
    scale_shape_manual(values = c(16, 17)) +
    labs(title = "B", shape = "Treatment") +
    theme_bw() +
    theme(legend.position = "none") + 
    xlab(paste0("PCoA1 [", round(100 * ord_obj$values$Relative_eig[1], 1), "%]")) +
    ylab(paste0("PCoA2 [", round(100 * ord_obj$values$Relative_eig[2], 1), "%]"))
  return(p)
}

# Create plots
plot.p_before_IPonly <- create_pcoa_plot_before_IPonly(data_before_IPonly, "Before", ord_before_IPonly)
plot.p_after_IPonly <- create_pcoa_plot_after_IPonly(data_after_IPonly, "After", ord_after_IPonly)

# Combine
library(patchwork)
plot.combined_plot_IPonly <- plot.p_before_IPonly + plot.p_after_IPonly + 
  plot_layout(guides = "collect")

print(plot.combined_plot_IPonly)

#### Check dispersion - IP-ONLY ##

han_before_IPonly <- data_before_IPonly[data_before_IPonly$StudyID_Split == "Han 2021", ]
han_after_IPonly <- data_after_IPonly[data_after_IPonly$StudyID_Split == "Han 2021", ]

cat("\n=== HAN 2021 DISPERSION (IP-only) ===\n")
cat("BEFORE batch correction:\n")
print(apply(han_before_IPonly[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(han_after_IPonly[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(han_before_IPonly[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(han_after_IPonly[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))

# Yang
yang_before_IPonly <- data_before_IPonly[data_before_IPonly$StudyID_Split == "Yang 2022", ]
yang_after_IPonly <- data_after_IPonly[data_after_IPonly$StudyID_Split == "Yang 2022", ]

cat("\n=== YANG 2022 DISPERSION (IP-only) ===\n")
cat("BEFORE batch correction:\n")
print(apply(yang_before_IPonly[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(yang_after_IPonly[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(yang_before_IPonly[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(yang_after_IPonly[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))

# Bitto IP
bitto_before_IPonly <- data_before_IPonly[data_before_IPonly$StudyID_Split == "Bitto 2016 (IP)", ]
bitto_after_IPonly <- data_after_IPonly[data_after_IPonly$StudyID_Split == "Bitto 2016 (IP)", ]

cat("\n=== BITTO 2016 IP DISPERSION (IP-only) ===\n")
cat("BEFORE batch correction:\n")
print(apply(bitto_before_IPonly[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nAFTER batch correction:\n")
print(apply(bitto_after_IPonly[, c("Axis.1", "Axis.2")], 2, sd))

cat("\nPercent change:\n")
before_sd <- apply(bitto_before_IPonly[, c("Axis.1", "Axis.2")], 2, sd)
after_sd <- apply(bitto_after_IPonly[, c("Axis.1", "Axis.2")], 2, sd)
percent_change <- ((after_sd - before_sd) / before_sd) * 100
print(round(percent_change, 1))


#### 11.lm_meta:run ##########

# 1) Extract the batch-corrected abundance data from your phyloseq object
genus_mat_adj <- as(otu_table(batch_corrected_ps), "matrix")
if (!taxa_are_rows(batch_corrected_ps)) genus_mat_adj <- t(genus_mat_adj)

# 2) Make sure your metadata is correctly formatted
rapa_meta <- as(sample_data(batch_corrected_ps), "data.frame")
rapa_meta$StudyID_Split <- factor(trimws(rapa_meta$StudyID_Split))
rapa_meta$StandardTreatment <- factor(trimws(rapa_meta$StandardTreatment))
rapa_meta$AdminRoute <- factor(trimws(rapa_meta$AdminRoute))
# 3) Define covariates (excluding the exposure variable)
covar_cols <- c("AdminRoute")  # Excluding StandardTreatment since it's our exposure

# 4) Verify that all samples match between abundance matrix and metadata
if(!all(colnames(genus_mat_adj) %in% rownames(rapa_meta))) {
  warning("Some samples in abundance matrix are missing from metadata")
  # Keep only matching samples
  genus_mat_adj <- genus_mat_adj[, colnames(genus_mat_adj) %in% rownames(rapa_meta)]
}
if(!all(rownames(rapa_meta) %in% colnames(genus_mat_adj))) {
  warning("Some samples in metadata are missing from abundance matrix")
  # Keep only matching samples
  rapa_meta <- rapa_meta[rownames(rapa_meta) %in% colnames(genus_mat_adj), ]
}

# Ensure same sample order
rapa_meta <- rapa_meta[colnames(genus_mat_adj), ]

# 1) Verify exact match and order
if (!identical(rownames(rapa_meta), colnames(genus_mat_adj))) {
  stop("Sample names/ordering STILL mismatch between metadata and abundance matrix!")
} else {
  message("ample names and ordering match perfectly.")
}

# 2) (Optional) Print a quick head of both to manually eyeball
print(head(rownames(rapa_meta)))
print(head(colnames(genus_mat_adj)))

# 1)running lm_meta 
fit_lm_meta <- MMUPHin::lm_meta(
  feature_abd  = genus_mat_adj,
  batch        = "StudyID_Split",
  exposure     = "StandardTreatment",
  covariates   = NULL,
  data         = rapa_meta,
  control      = list(
    normalization   = "NONE",
    transform       = "AST",
    analysis_method = "LM",      # or "CPLM" if you’re modeling raw counts
    rma_method      = "REML",    # best balance of bias & variance
    rma_conv        = 1e-6,      # tighter convergence criterion
    rma_maxit       = 5000,      # more iterations if needed
    verbose         = TRUE
  )
)

# 1a) Pull out the full results
all_results <- fit_lm_meta$meta_fits %>%
  mutate(qval.fdr = p.adjust(pval, method = "fdr"))

# 1b) Keep only genera present in ≥2 studies
results_k2 <- all_results %>%
  filter(k >= 3)

# Check your results
print(head(fit_lm_meta$meta_fits))
cat("Total genera in ≥3 studies:", nrow(results_k2), "\n")

sig_results <- results_k2 %>%
  filter(qval.fdr < 0.05,
         !grepl("uncultured|unclassified|Unclassified|Other|Mitochondria", feature, ignore.case=TRUE)
  )%>%
arrange(qval.fdr)

# 5) Pick the top 5 by FDR
top_genera <- head(sig_results, 10)

print(top_genera)
print(all)

## LM:meta IP-ONLY STUDIES ####

# 1) Filter to IP-only samples
rapa_meta_IP <- rapa_meta %>%
  filter(StudyID_Split != "Study2_Food") %>%
  droplevels()  # CRITICAL: Drop unused factor levels

# Keep only IP samples in abundance matrix
genus_mat_IP <- genus_mat_adj[, rownames(rapa_meta_IP)]

cat("\n=== IP-ONLY ANALYSIS ===\n")
cat("Studies included:", unique(as.character(rapa_meta_IP$StudyID_Split)), "\n")
cat("Total samples:", ncol(genus_mat_IP), "\n")

# Check treatment levels in each study
cat("\nTreatment distribution by study:\n")
print(table(rapa_meta_IP$StudyID_Split, rapa_meta_IP$StandardTreatment))

# Try running with more lenient filtering
fit_lm_meta_IP <- MMUPHin::lm_meta(
  feature_abd  = genus_mat_IP,
  batch        = "StudyID_Split",
  exposure     = "StandardTreatment",
  covariates   = NULL,
  data         = rapa_meta_IP,
  control      = list(
    normalization   = "NONE",
    transform       = "AST",
    analysis_method = "LM",
    rma_method      = "REML",
    rma_conv        = 1e-6,      # tighter convergence criterion
    rma_maxit       = 5000,      # more iterations if needed
    verbose         = TRUE
  )
)

# 3) Extract significant results (IP-ONLY)
all_results_IP <- fit_lm_meta_IP$meta_fits %>%
  mutate(qval.fdr = p.adjust(pval, method = "fdr"))

results_k2_IP <- all_results_IP %>%
  filter(k >= 2)

sig_results_IP <- results_k2_IP %>%
  filter(qval.fdr < 0.05,
         !grepl("uncultured|unclassified|Unclassified|Other|Mitochondria", 
                feature, ignore.case=TRUE)) %>%
  arrange(qval.fdr)

cat("\n=== IP-ONLY RESULTS ===\n")

print(sig_results_IP %>% select(feature, coef, stderr, pval, qval.fdr, k, I2,tau2))
print(all_results_IP[all_results_IP$feature == "Bacteroides", c("feature", "coef", "stderr","pval", "qval.fdr","k", "I2", "tau2")])
print(all_results_IP[all_results_IP$feature == "Muribaculum", c("feature", "coef", "stderr", "pval", "qval.fdr","k", "I2", "tau2")])
print(all_results_IP[all_results_IP$feature == "Candidatus_Arthromitus", c("feature", "coef", "pval", "qval.fdr","stderr", "k", "I2", "tau2")])
print(all_results_IP[all_results_IP$feature == "ASF356", c("feature", "coef","pval", "qval.fdr", "stderr", "k", "I2", "tau2")])

# Check ALL significant taxa
for(taxon in sig_results$feature) {
  all_res <- sig_results[sig_results$feature == taxon, c("coef", "stderr", "k", "I2", "tau2")]
  ip_res <- all_results_IP[all_results_IP$feature == taxon, c("coef", "stderr", "k", "I2", "tau2")]
  
  cat("\n===", taxon, "===\n")
  cat("ALL:", "k=", all_res$k, ", I²=", round(all_res$I2,1), "%, tau²=", round(all_res$tau2,4), "\n")
  cat("IP: ", "k=", ip_res$k, ", I²=", round(ip_res$I2,1), "%, tau²=", round(ip_res$tau2,4), "\n")
  cat("tau² changed:", abs(all_res$tau2 - ip_res$tau2) > 0.001, "\n")
}

# Extract your 4 main taxa from IP-only k=3 results
main_taxa_IP <- results_k2_IP %>%
  filter(feature %in% c("Bacteroides", "Muribaculum", 
                        "Candidatus_Arthromitus", "ASF356")) %>%
  select(feature, coef, stderr, pval, qval.fdr, k, I2, tau2) %>%
  arrange(qval.fdr)

cat("\n=== YOUR 4 MAIN TAXA: IP-ONLY (k=3) ===\n")
print(main_taxa_IP)

# Compare to ALL-studies results
main_taxa_ALL <- sig_results %>%
  filter(feature %in% c("Bacteroides", "Muribaculum", 
                        "Candidatus_Arthromitus", "ASF356")) %>%
  select(feature, coef, stderr, pval, qval.fdr, k, I2, tau2) %>%
  arrange(qval.fdr)

cat("\n=== YOUR 4 MAIN TAXA: ALL STUDIES (k=4 or 3) ===\n")
print(main_taxa_ALL)

# Side-by-side comparison
cat("\n=== COMPARISON TABLE ===\n")
comparison <- data.frame(
  Taxon = c("Bacteroides", "Muribaculum", "Candidatus_Arthromitus", "ASF356"),
  
  ALL_k = main_taxa_ALL$k,
  ALL_pval = sprintf("%.4f", main_taxa_ALL$pval),
  ALL_FDR = sprintf("%.4f", main_taxa_ALL$qval.fdr),
  ALL_Sig = ifelse(main_taxa_ALL$qval.fdr < 0.05, "✓", "✗"),
  
  IP_k = c(main_taxa_IP$k[main_taxa_IP$feature == "Bacteroides"],
           main_taxa_IP$k[main_taxa_IP$feature == "Muribaculum"],
           main_taxa_IP$k[main_taxa_IP$feature == "Candidatus_Arthromitus"],
           NA),  # ASF356 not in k=3
  IP_pval = c(sprintf("%.4f", main_taxa_IP$pval[main_taxa_IP$feature == "Bacteroides"]),
              sprintf("%.4f", main_taxa_IP$pval[main_taxa_IP$feature == "Muribaculum"]),
              sprintf("%.4f", main_taxa_IP$pval[main_taxa_IP$feature == "Candidatus_Arthromitus"]),
              "NA (k=2)"),
  IP_FDR = c(sprintf("%.4f", main_taxa_IP$qval.fdr[main_taxa_IP$feature == "Bacteroides"]),
             sprintf("%.4f", main_taxa_IP$qval.fdr[main_taxa_IP$feature == "Muribaculum"]),
             sprintf("%.4f", main_taxa_IP$qval.fdr[main_taxa_IP$feature == "Candidatus_Arthromitus"]),
             "NA (k=2)"),
  IP_Sig = c("✗", "✗", "✗", "✗")
)

print(comparison)

###Comparison###
# For Bacteroides specifically
bact_result <- sig_results %>% filter(feature == "Bacteroides")

cat("=== BACTEROIDES DIAGNOSTICS ===\n")
cat("Pooled: coef =", bact_result$coef, ", SE =", bact_result$stderr, 
    ", CI width =", bact_result$stderr * 1.96 * 2, "\n")
cat("I² =", bact_result$I2, "%, τ² =", bact_result$tau2, "\n\n")

# Individual studies
for(study in names(fit_lm_meta$maaslin_fits)) {
  study_res <- fit_lm_meta$maaslin_fits[[study]]
  bact_study <- study_res[study_res$feature == "Bacteroides", ]
  if(nrow(bact_study) > 0) {
    ci_width <- bact_study$stderr * 1.96 * 2
    cat(study, ": coef =", bact_study$coef, ", SE =", bact_study$stderr,
        ", CI width =", ci_width, "\n")
  }
}



# Look at taxa significant at p < 0.05 (before FDR)
pre_fdr_sig <- fit_lm_meta_IP$meta_fits %>%
  filter(pval < 0.05) %>%
  arrange(pval) %>%
  select(feature, coef, stderr, pval, qval.fdr, k, I2)

cat("Taxa significant at p < 0.05 (before FDR):\n")
print(pre_fdr_sig)

cat("\nHow many lost significance after FDR?\n")
cat("p < 0.05:", nrow(pre_fdr_sig), "taxa\n")
cat("FDR < 0.05:", sum(fit_lm_meta_IP$meta_fits$qval.fdr < 0.05, na.rm = TRUE), "taxa\n")

# Count NAs per taxon across studies
na_counts <- data.frame(
  feature = fit_lm_meta_IP$meta_fits$feature,
  k = fit_lm_meta_IP$meta_fits$k,
  stringsAsFactors = FALSE
)

# Check how many taxa are in all 3 IP studies vs fewer
cat("\nTaxa by number of studies (k):\n")
print(table(na_counts$k))

# Find taxa present in all 3 studies
taxa_all_3 <- na_counts %>% filter(k == 3)
cat("\nTaxa present in all 3 IP studies:", nrow(taxa_all_3), "\n")

# Find taxa in only 1-2 studies
taxa_incomplete <- na_counts %>% filter(k < 3)
cat("Taxa in <3 studies:", nrow(taxa_incomplete), "\n")


# Create comparison table
comparison_summary <- data.frame(
  Analysis = c("All 4 studies", "IP-only (3 studies)"),
  Total_Tested = c(nrow(results_k2), nrow(results_k2_IP)),
  Significant = c(nrow(sig_results), nrow(sig_results_IP)),
  Studies_Included = c("Study2_Food, Study2_IP, Study5, Study7", 
                       "Study2_IP, Study5, Study7")
)
print(comparison_summary)

# Check overlap
cat("\n--- Overlap in significant genera ---\n")
all_sig_names <- sig_results$feature
ip_sig_names <- sig_results_IP$feature

cat("Only in ALL studies:", paste(setdiff(all_sig_names, ip_sig_names), collapse=", "), "\n")
cat("Only in IP-only:", paste(setdiff(ip_sig_names, all_sig_names), collapse=", "), "\n")
cat("In both:", paste(intersect(all_sig_names, ip_sig_names), collapse=", "), "\n")

print(all_results_IP[all_results_IP$feature == "Bacteroides", c("I2", "tau2")])

###Diagnostics ####
# TEST 1: Check k values distribution
cat("--- Distribution of k values (IP-only) ---\n")
k_distribution <- table(sig_results_IP$k)
print(k_distribution)

# TEST 2: For k=2 genera, check which study is missing
if(any(sig_results_IP$k == 2)) {
  cat("\n--- Analyzing k=2 genera (present in only 2/3 IP studies) ---\n")
  k2_genera <- sig_results_IP %>% filter(k == 2)
  
  for(genus in k2_genera$feature) {
    cat("\n● ", genus, "\n", sep="")
    
    for(study in c("Study2_IP", "Study5", "Study7")) {
      study_results <- fit_lm_meta_IP$maaslin_fits[[study]]
      study_samples <- rownames(rapa_meta_IP)[rapa_meta_IP$StudyID_Split == study]
      abundances <- genus_mat_IP[genus, study_samples]
      
      if(genus %in% study_results$feature) {
        effect <- study_results[study_results$feature == genus, ]
        cat(sprintf("  %s: ✓ Present (coef=%.3f, p=%.4f, mean=%.4f)\n", 
                    study, effect$coef, effect$pval, mean(abundances)))
      } else {
        prevalence <- 100 * sum(abundances > 0) / length(abundances)
        cat(sprintf("  %s: ✗ Absent (mean=%.4f, prev=%.1f%%)\n", 
                    study, mean(abundances), prevalence))
      }
    }
  }
} else {
  cat("No k=2 genera found.\n")
}

# TEST 3: Check Study2_Food effects on genera significant in ALL studies
if(nrow(sig_results) > 0) {
  cat("\n--- Study2_Food effects on ALL-studies significant genera ---\n")
  
  food_study_data <- fit_lm_meta$maaslin_fits[["Study2_Food"]]
  
  for(genus in head(sig_results$feature, 5)) {  # Check top 5
    if(genus %in% food_study_data$feature) {
      effect <- food_study_data[food_study_data$feature == genus, ]
      cat(sprintf("%s: coef=%.3f, p=%.4f\n", 
                  genus, effect$coef, effect$pval))
    } else {
      cat(sprintf("%s: Not present in Study2_Food\n", genus))
    }
  }
}

# TEST 4: Sample size per study
cat("\n--- Sample sizes per study ---\n")
print(table(rapa_meta$StudyID_Split, rapa_meta$StandardTreatment))

####Viszalizations####
# Create visualizations for the top genera
for (i in 1:nrow(top_genera)) {
  genus <- top_genera$feature[i]
  
  # Extract study-specific effects
  genus_effects <- data.frame()
  
  for (study in names(fit_lm_meta$maaslin_fits)) {
    study_results <- fit_lm_meta$maaslin_fits[[study]]
    
    if (!is.null(study_results)) {
      # Find the genus in this study's results
      genus_row <- study_results[study_results$feature == genus, ]
      
      if (nrow(genus_row) > 0) {
        # Add to our data frame
        genus_effects <- rbind(genus_effects, data.frame(
          Study = study,
          Feature = genus,
          Coefficient = genus_row$coef,
          StdError = genus_row$stderr,
          P_value = genus_row$pval
        ))
      }
    }
  }
  
  # Print study-specific effects
  cat("\nStudy-specific effects for genus:", genus, "\n")
  print(genus_effects)
  
  # Create visualization data
  if (nrow(genus_effects) > 0) {
    # Data for boxplot
    genus_data <- data.frame(
      Abundance = genus_mat_adj[genus, ],
      SampleID = colnames(genus_mat_adj)
    )
    
    # Merge with metadata
    genus_data <- merge(genus_data, rapa_meta, by.x = "SampleID", by.y = 0)
    
    p <- ggplot(genus_data, aes(x = StandardTreatment, y = Abundance, fill = StandardTreatment)) +
      geom_boxplot() +
      facet_wrap(~StudyID_Split, scales = "free_y") +
      labs(
        title = paste('Abundance of', genus, 'by Treatment and Study'),
        subtitle = paste("Meta-analysis p =", round(top_genera$pval[i], 4)),
        y = 'Batch-corrected abundance',
        x = 'Treatment'
      ) +
      theme_bw()
    
    print(p)
    
    # Also create a forest plot
    if (nrow(genus_effects) > 1) {
      # Add meta-analysis result
      meta_row <- top_genera[i, ]
      genus_effects <- rbind(genus_effects, data.frame(
        Study = "Meta-analysis",
        Feature = genus,
        Coefficient = meta_row$coef,
        StdError = meta_row$stderr,
        P_value = meta_row$pval
      ))
      
      # Calculate confidence intervals
      genus_effects$Lower <- genus_effects$Coefficient - 1.96 * genus_effects$StdError
      genus_effects$Upper <- genus_effects$Coefficient + 1.96 * genus_effects$StdError
      
      # Forest plot
      forest_plot <- ggplot(genus_effects, aes(x = Coefficient, y = Study)) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
        labs(
          title = paste("Forest Plot for", genus),
          subtitle = paste("Meta-analysis p =", round(meta_row$pval, 4)),
          x = "Effect Size (Coefficient)",
          y = ""
        ) +
        theme_bw()
      
      print(forest_plot)
    }
  }
}

results_k2 <- results_k2 %>%
  mutate(significant = qval.fdr < 0.05)

# Remove unwanted taxa only (not filtering for significance)
results_plot <- results_k2 %>%
  filter(!grepl("uncultured|unclassified|other|mitochondria", feature, ignore.case=TRUE)) %>%
  arrange(coef)

# Optional: mark which ones are significant for the plot
results_plot <- results_plot %>%
  mutate(significant = qval.fdr < 0.05)

print(sig_results$feature)


####second version

# Check what we're getting for each genus
for(taxon in c("Bacteroides", "Muribaculum", "Candidatus_Arthromitus", "ASF356")) {
  cat("\n=== CHECKING", taxon, "===\n")
  
  # Check sig_tab data
  meta_result <- sig_results %>% filter(feature == taxon)
  cat("Meta-analysis from sig_tab:\n")
  cat("  Coefficient:", meta_result$coef, "\n")
  cat("  Std Error:", meta_result$stderr, "\n")
  cat("  I2:", meta_result$I2, "%\n")
  cat("  tau2:", meta_result$tau2, "\n")
  cat("  Weights: Study2_Food =", meta_result$weight_Study2_Food, "%, Study2_IP =", meta_result$weight_Study2_IP, "%, Study5 =", meta_result$weight_Study5, "%, Study7 =", meta_result$weight_Study7, "%\n")
  # Check individual study data
  
  cat("Individual studies from maaslin_fits:\n")
  for (study in names(fit_lm_meta$maaslin_fits)) {
    study_results <- fit_lm_meta$maaslin_fits[[study]]
    if (!is.null(study_results)) {
      genus_row <- study_results[study_results$feature == taxon, ]
      if (nrow(genus_row) > 0) {
        cat("  ", study, "- coef:", genus_row$coef, ", stderr:", genus_row$stderr, ", pval:", genus_row$pval,"\n")
      }
    }
  }
}

# Create forest plot function for MMUPHin results (using correct data sources)
create_forest_from_mmuphin <- function(mmuphin_fit, sig_tab_data, taxa_list) {
  
  all_plot_data <- data.frame()
  
  # For each taxon, add meta-analysis and study-specific data
  for(taxon in taxa_list) {
    # Get meta-analysis result from sig_tab (filtered, significant results)
    meta_result <- sig_tab_data %>%
      filter(feature == taxon)
    
    # Get study-specific data from MMUPHin maaslin_fits (your exact method)
    study_data <- data.frame()
    
    for (study in names(mmuphin_fit$maaslin_fits)) {
      study_results <- mmuphin_fit$maaslin_fits[[study]]
      
      if (!is.null(study_results)) {
        genus_row <- study_results[study_results$feature == taxon, ]
        
        if (nrow(genus_row) > 0) {
          study_data <- rbind(study_data, data.frame(
            study = study,
            yi = genus_row$coef,
            sei = genus_row$stderr,
            ci.lb = genus_row$coef - 1.96 * genus_row$stderr,
            ci.ub = genus_row$coef + 1.96 * genus_row$stderr,
            pval = genus_row$pval
          ))
        }
      }
    }
    
    # Create data frames for plotting
    meta_df <- data.frame(
      study = "Pooled (RE)",
      yi = meta_result$coef,
      ci.lb = meta_result$coef - 1.96 * meta_result$stderr,
      ci.ub = meta_result$coef + 1.96 * meta_result$stderr,
      is_summary = TRUE,
      taxon_name = taxon,
      sei = meta_result$stderr
    )
    
    study_df <- data.frame(
      study = study_data$study,
      yi = study_data$yi,
      ci.lb = study_data$ci.lb,
      ci.ub = study_data$ci.ub,
      is_summary = FALSE,
      taxon_name = taxon,
      sei = study_data$sei
    )
    
    # Add to all plot data
    taxon_data <- rbind(meta_df, study_df)
    all_plot_data <- rbind(all_plot_data, taxon_data)
  }
  
  # Set factor levels for ordering with proper study names
  all_plot_data$study <- recode(all_plot_data$study,
                                "Study2_Food" = "Bitto 2016 (Diet)",
                                "Study2_IP" = "Bitto 2016 (IP)",
                                "Study5" = "Yang 2022", 
                                "Study7" = "Han 2021"
  )
  
  all_plot_data$study <- factor(all_plot_data$study, 
                                levels = rev(c("Pooled (RE)", "Bitto 2016 (Diet)", "Bitto 2016 (IP)", "Yang 2022", "Han 2021")))
  
  # Calculate GLOBAL x-axis limits across ALL taxa (so all plots have same scale)
  global_x_min <- min(all_plot_data$ci.lb, na.rm = TRUE)
  global_x_max <- max(all_plot_data$ci.ub, na.rm = TRUE)
  global_x_range <- global_x_max - global_x_min
  global_x_limits <- c(global_x_min - 0.05 * global_x_range, 0.26)
  
  # Create breaks for the global range
  break_start <- floor(global_x_min / 0.04) * 0.04
  break_end <- 0.16
  global_x_breaks <- seq(break_start, break_end, by = 0.04)
  
  # Create plots for each taxon
  forest_plots <- list()
  
  # NEW:
  for(i in seq_along(taxa_list)){
    taxon <- taxa_list[i]
    is_last <- (i == length(taxa_list))
    plot_data <- all_plot_data %>% filter(taxon_name == taxon)
    
    # Get meta-analysis result for heterogeneity stats (from sig_tab)
    meta_result <- sig_tab_data %>% filter(feature == taxon)
    
    # Extract heterogeneity stats from sig_tab
    i2 <- paste0(round(meta_result$I2, 0), "%")
    tau2 <- round(meta_result$tau2, 4)
  
    
    # Add effect size text
    plot_data$effect_text <- sprintf("%.2f [%.2f; %.2f]", 
                                     plot_data$yi, 
                                     plot_data$ci.lb, 
                                     plot_data$ci.ub)
    
    # Get correct weights from sig_results (not calculated)
    study2_food_weight <- ifelse(is.na(meta_result$weight_Study2_Food), "—", paste0(round(meta_result$weight_Study2_Food, 1), "%"))
    study2_ip_weight <- ifelse(is.na(meta_result$weight_Study2_IP), "—", paste0(round(meta_result$weight_Study2_IP, 1), "%"))
    study5_weight <- ifelse(is.na(meta_result$weight_Study5), "—", paste0(round(meta_result$weight_Study5, 1), "%"))
    study7_weight <- ifelse(is.na(meta_result$weight_Study7), "—", paste0(round(meta_result$weight_Study7, 1), "%"))
    
    plot_data$weight_text <- case_when(
      plot_data$is_summary ~ "—",
      plot_data$study == "Bitto 2016 (Diet)" ~ study2_food_weight,
      plot_data$study == "Bitto 2016 (IP)" ~ study2_ip_weight,
      plot_data$study == "Yang 2022" ~ study5_weight,
      plot_data$study == "Han 2021" ~ study7_weight,
      TRUE ~ "")
    
    p <- ggplot(plot_data, aes(x = yi, y = study, xmin = ci.lb, xmax = ci.ub)) +
      geom_vline(xintercept = global_x_breaks, color = "gray90") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
      geom_errorbarh(height = 0.08, linewidth = 0.8, aes(color = is_summary)) +
      geom_point(aes(size = is_summary, color = is_summary)) +
      geom_text(aes(label = effect_text), x = 0.18, hjust = 0, size = 3) +
      geom_text(aes(label = weight_text), x = 0.24, hjust = 0, size = 3) +
      scale_x_continuous(
        limits = global_x_limits,
        breaks = global_x_breaks,
        minor_breaks = NULL
      ) +
      scale_y_discrete(expand = expansion(add = c(0.22, 1.1))) +
      theme_minimal() +
      coord_cartesian(clip = "off") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        axis.title.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        plot.margin = margin(6, 28, 12, 70)
        ) +
      scale_color_manual(values = c("FALSE" = "royalblue3", "TRUE" = "darkred")) +
      scale_size_manual(values = c("FALSE" = 5, "TRUE" = 5)) +
      labs(
        title = paste0(taxon),
        x = if(is_last)"Treatment Effect (Linear Coefficient)" else NULL
      ) +
    annotate("text", x = 0.18, y = length(levels(plot_data$study)) + 0.9,
             label = "Effect [95% CI]", hjust = 0, size = 3, fontface = "bold") +
      annotate("text", x = 0.24, y = length(levels(plot_data$study)) + 0.9,
               label = "Weight",          hjust = 0, size = 3, fontface = "bold") +
      annotate("text",  x = -Inf, y = 0.5,
               label = paste("Heterogeneity: I² =", i2, ", τ² =", tau2), 
               hjust = 0.03, size = 3.5, color = "gray40")
    
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
                      plot.subtitle = element_text(size = 12, hjust = 0.5),
                      plot.caption = element_text(size = 10, hjust = 0.5),
                      panel.spacing = unit(0.18, "lines"))
      )
  }
  
  # Return plots
  return(list(
    individual = forest_plots,
    combined = combined_plot
  ))
}

# Run with sig_results for meta-analysis results (3 significant genera)
results <- create_forest_from_mmuphin(
  mmuphin_fit = fit_lm_meta,
  sig_tab_data = sig_results,
  taxa_list = c("Bacteroides", "Muribaculum", "Candidatus_Arthromitus", "ASF356")
)

# Get the COMBINED forest plot 
plot.MetaForest_combined <- results$combined

# Display the combined plot
print(plot.MetaForest_combined)

#### Forest Plot IP ####

for(taxon in sig_results$feature) {
  cat("\n=== CHECKING", taxon, "===\n")
  
  # Check ALL studies meta-analysis
  meta_result_all <- sig_results %>% filter(feature == taxon)
  cat("Meta-analysis ALL (sig_results):\n")
  cat("  Coefficient:", meta_result_all$coef, "\n")
  cat("  Std Error:", meta_result_all$stderr, "\n")
  cat("  I2:", meta_result_all$I2, "%\n")
  cat("  P-value:", meta_result_all$pval, "\n")
  cat("  FDR:", meta_result_all$qval.fdr, "\n")
  
  # Check IP-only meta-analysis from ALL_RESULTS_IP (not just significant)
  meta_result_ip <- all_results_IP %>% filter(feature == taxon)
  if(nrow(meta_result_ip) > 0) {
    cat("Meta-analysis IP-only (all_results_IP):\n")
    cat("  Coefficient:", meta_result_ip$coef, "\n")
    cat("  Std Error:", meta_result_ip$stderr, "\n")
    cat("  I2:", meta_result_ip$I2, "%\n")
    cat("  P-value:", meta_result_ip$pval, "\n")
    cat("  FDR:", meta_result_ip$qval.fdr, "\n")
    if(meta_result_ip$qval.fdr < 0.05) {
      cat("  ✓ SIGNIFICANT in IP-only\n")
    } else {
      cat("  ✗ NOT significant in IP-only (FDR >= 0.05)\n")
    }
  } else {
    cat("Meta-analysis IP-only: NOT FOUND\n")
  }
  
  # Check individual study data
  cat("Individual studies from maaslin_fits:\n")
  for (study in names(fit_lm_meta$maaslin_fits)) {
    study_results <- fit_lm_meta$maaslin_fits[[study]]
    if (!is.null(study_results)) {
      genus_row <- study_results[study_results$feature == taxon, ]
      if (nrow(genus_row) > 0) {
        cat("  ", study, "- coef:", genus_row$coef, ", stderr:", genus_row$stderr, 
            ", pval:", genus_row$pval, "\n")
      }
    }
  }
}

create_forest_with_both_pooled <- function(mmuphin_fit_all, 
                                           mmuphin_fit_ip,
                                           sig_tab_all, 
                                           all_results_ip_data,
                                           taxa_list) {
  
  all_plot_data <- data.frame()
  
  # For each taxon, compile all data
  for(taxon in taxa_list) {
    
    # 1. Get POOLED ALL meta-analysis result
    meta_result_all <- sig_tab_all %>% filter(feature == taxon)
    
    # 2. Get POOLED IP from ALL results (not just significant)
    meta_result_ip <- all_results_ip_data %>% 
      filter(feature == taxon)
    
    # 3. Get individual study data
    study_data <- data.frame()
    
    for (study in names(mmuphin_fit_all$maaslin_fits)) {
      study_results <- mmuphin_fit_all$maaslin_fits[[study]]
      
      if (!is.null(study_results)) {
        genus_row <- study_results[study_results$feature == taxon, ]
        
        if (nrow(genus_row) > 0) {
          study_data <- rbind(study_data, data.frame(
            study = study,
            yi = genus_row$coef,
            sei = genus_row$stderr,
            ci.lb = genus_row$coef - 1.96 * genus_row$stderr,
            ci.ub = genus_row$coef + 1.96 * genus_row$stderr,
            pval = genus_row$pval
          ))
        }
      }
    }
    
    # Create data frame for POOLED ALL
    pooled_all_df <- data.frame(
      study = "Pooled (All)",  # Changed from "Pooled ALL"
      yi = meta_result_all$coef,
      ci.lb = meta_result_all$coef - 1.96 * meta_result_all$stderr,
      ci.ub = meta_result_all$coef + 1.96 * meta_result_all$stderr,
      is_summary = TRUE,
      summary_type = "all",
      taxon_name = taxon,
      sei = meta_result_all$stderr,
      pval = meta_result_all$pval,
      fdr = meta_result_all$qval.fdr,
      i2 = meta_result_all$I2,
      tau2 = meta_result_all$tau2,
      weight_study2_food = meta_result_all$weight_Study2_Food,
      weight_study2_ip = meta_result_all$weight_Study2_IP,
      weight_study5 = meta_result_all$weight_Study5,
      weight_study7 = meta_result_all$weight_Study7
    )
    
    
    # Create data frame for POOLED IP (if exists)
    if(nrow(meta_result_ip) > 0) {
      pooled_ip_df <- data.frame(
        study = "Pooled (IP only)",  # Changed from "Pooled IP"
        yi = meta_result_ip$coef,
        ci.lb = meta_result_ip$coef - 1.96 * meta_result_ip$stderr,
        ci.ub = meta_result_ip$coef + 1.96 * meta_result_ip$stderr,
        is_summary = TRUE,
        summary_type = "ip",
        taxon_name = taxon,
        sei = meta_result_ip$stderr,
        pval = meta_result_ip$pval,
        fdr = meta_result_ip$qval.fdr,
        i2 = meta_result_ip$I2,
        tau2 = meta_result_ip$tau2,
        weight_study2_food = NA,
        weight_study2_ip = NA,
        weight_study5 = NA,
        weight_study7 = NA
      )
    } else {
      pooled_ip_df <- data.frame()
    }
    
    # Create data frame for individual studies
    study_df <- data.frame(
      study = study_data$study,
      yi = study_data$yi,
      ci.lb = study_data$ci.lb,
      ci.ub = study_data$ci.ub,
      is_summary = FALSE,
      summary_type = "individual",
      taxon_name = taxon,
      sei = study_data$sei,
      pval = study_data$pval,               # Already there
      fdr = NA,                             # ADD
      i2 = NA,
      tau2 = NA,
      weight_study2_food = NA,              # ADD
      weight_study2_ip = NA,                # ADD
      weight_study5 = NA,                   # ADD
      weight_study7 = NA                    # ADD
    )
    
    # Combine all data for this taxon
    taxon_data <- rbind(pooled_all_df, pooled_ip_df, study_df)
    all_plot_data <- rbind(all_plot_data, taxon_data)
  }
  
  # Recode study names
  all_plot_data$study <- recode(all_plot_data$study,
                                "Study2_Food" = "Bitto 2016 (Diet)",
                                "Study2_IP" = "Bitto 2016 (IP)",
                                "Study5" = "Yang 2022", 
                                "Study7" = "Han 2021"
  )
  
  # Set factor levels in the ORDER you want (reversed for plotting)
  all_plot_data$study <- factor(all_plot_data$study, 
                                levels = rev(c("Pooled (All)",      # Changed
                                               "Bitto 2016 (Diet)", 
                                               "Pooled (IP only)",  # Changed
                                               "Bitto 2016 (IP)", 
                                               "Yang 2022", 
                                               "Han 2021")))
  
  # Add significance markers
  all_plot_data$sig_marker <- ""
  all_plot_data$sig_marker[all_plot_data$summary_type == "all" & all_plot_data$fdr < 0.05] <- "*"
  all_plot_data$sig_marker[all_plot_data$summary_type == "ip" & all_plot_data$fdr < 0.05] <- "*"

  
  # Calculate GLOBAL x-axis limits
  global_x_min <- min(all_plot_data$ci.lb, na.rm = TRUE)
  global_x_max <- max(all_plot_data$ci.ub, na.rm = TRUE)
  global_x_range <- global_x_max - global_x_min
  global_x_limits <- c(global_x_min - 0.05 * global_x_range, 0.26)
  
  # Create breaks
  break_start <- floor(global_x_min / 0.04) * 0.04
  break_end <- 0.16
  global_x_breaks <- seq(break_start, break_end, by = 0.04)
  
  # Create plots for each taxon
  forest_plots <- list()
  
  for(i in seq_along(taxa_list)) {
    taxon <- taxa_list[i]
    is_last <- (i == length(taxa_list))
    plot_data <- all_plot_data %>% filter(taxon_name == taxon)
    
    # Get heterogeneity stats from POOLED ALL
    meta_result_all <- sig_tab_all %>% filter(feature == taxon)
    i2_all <- paste0(round(meta_result_all$I2, 0), "%")
    tau2_all <- round(meta_result_all$tau2, 4)
    
    # Get IP result for significance check
    meta_result_ip <- all_results_ip_data %>% filter(feature == taxon)
    
    # Add effect size text
    plot_data$effect_text <- sprintf("%.2f [%.2f; %.2f]%s",
                                     plot_data$yi, 
                                     plot_data$ci.lb, 
                                     plot_data$ci.ub,
                                     plot_data$sig_marker)
    
    # Get weights from sig_results for POOLED ALL
    study2_food_weight <- ifelse(is.na(meta_result_all$weight_Study2_Food), 
                                 "—", 
                                 paste0(round(meta_result_all$weight_Study2_Food, 1), "%"))
    study2_ip_weight <- ifelse(is.na(meta_result_all$weight_Study2_IP), 
                               "—", 
                               paste0(round(meta_result_all$weight_Study2_IP, 1), "%"))
    study5_weight <- ifelse(is.na(meta_result_all$weight_Study5), 
                            "—", 
                            paste0(round(meta_result_all$weight_Study5, 1), "%"))
    study7_weight <- ifelse(is.na(meta_result_all$weight_Study7), 
                            "—", 
                            paste0(round(meta_result_all$weight_Study7, 1), "%"))
    
    # Assign weights
    plot_data$weight_text <- case_when(
      plot_data$study == "Pooled (All)" ~ "—",      # Changed
      plot_data$study == "Pooled (IP only)" ~ "—",  # Changed
      plot_data$study == "Bitto 2016 (Diet)" ~ study2_food_weight,
      plot_data$study == "Bitto 2016 (IP)" ~ study2_ip_weight,
      plot_data$study == "Yang 2022" ~ study5_weight,
      plot_data$study == "Han 2021" ~ study7_weight,
      TRUE ~ "")
    
    # Check if IP pooled is significant for this taxon
    ip_is_sig <- if(nrow(meta_result_ip) > 0) meta_result_ip$qval.fdr < 0.05 else FALSE
    
    # Create plot
    p <- ggplot(plot_data, aes(x = yi, y = study, xmin = ci.lb, xmax = ci.ub)) +
      geom_vline(xintercept = global_x_breaks, color = "gray90") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
      geom_errorbarh(height = 0.08, linewidth = 0.8, 
                     aes(color = summary_type)) +
      geom_point(aes(size = is_summary, color = summary_type, shape = summary_type)) +
      geom_text(aes(label = effect_text), x = 0.18, hjust = 0, size = 3) +
      geom_text(aes(label = weight_text), x = 0.24, hjust = 0, size = 3) +
      scale_x_continuous(
        limits = global_x_limits,
        breaks = global_x_breaks,
        minor_breaks = NULL
      ) +
      scale_y_discrete(expand = expansion(add = c(0.22, 1.1))) +
      theme_minimal() +
      coord_cartesian(clip = "off") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        axis.title.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        plot.margin = margin(6, 28, 12, 70)
      ) +
      scale_color_manual(values = c("all" = "darkred", 
                                    "ip" = "darkred", 
                                    "individual" = "royalblue3")) +
      scale_size_manual(values = c("FALSE" = 5, "TRUE" = 5)) +
      scale_shape_manual(values = c("all" = 16,      # Diamond for ALL
                                    "ip" = 16,       # Diamond for IP
                                    "individual" = 16)) +  # Circle for studies
      labs(
        title = paste0(taxon),
        x = if(is_last)"Treatment Effect (Linear Coefficient)" else NULL
      ) +
      annotate("text", x = 0.18, y = length(levels(plot_data$study)) + 0.9,
               label = "Effect [95% CI]", hjust = 0, size = 3, fontface = "bold") +
      annotate("text", x = 0.24, y = length(levels(plot_data$study)) + 0.9,
               label = "Weight", hjust = 0, size = 3, fontface = "bold") +
      annotate("text", x = -Inf, y = 0.3,
               label = paste("Heterogeneity (ALL): I² =", i2_all, ", τ² =", tau2_all), 
               hjust = 0.03, size = 3.5, color = "gray40") 
    
    forest_plots[[taxon]] <- p
  }
  
  # Combine plots
  if(length(forest_plots) > 0) {
    combined_plot <- forest_plots[[taxa_list[1]]]
    
    if(length(forest_plots) > 1) {
      for(i in 2:length(taxa_list)) {
        combined_plot <- combined_plot / forest_plots[[taxa_list[i]]]
      }
    }
    
    combined_plot <- combined_plot +
      plot_annotation(
        title = "Effect of Rapamycin on Gut Microbiota",
        subtitle = "Meta-analysis showing pooled estimates from all 4 studies and IP-only 3 studies",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 11, hjust = 0.5),
                      panel.spacing = unit(0.18, "lines"))
      )
  }
  
  # Return plots
 return(list(
  individual = forest_plots,
  combined = combined_plot,
  plot_data = all_plot_data          # ADD THIS
))
}

# Get list of significant taxa (from ALL studies analysis)
significant_taxa <- sig_results$feature

# Create forest plots with BOTH pooled estimates
# KEY CHANGE: Using all_results_IP instead of sig_results_IP
results <- create_forest_with_both_pooled(
  mmuphin_fit_all = fit_lm_meta,
  mmuphin_fit_ip = fit_lm_meta_IP,
  sig_tab_all = sig_results,
  all_results_ip_data = all_results_IP,  # ← CHANGED FROM sig_results_IP
  taxa_list = significant_taxa
)

# Get the COMBINED forest plot 
plot.MetaForest_combined <- results$combined

# Display the combined plot
print(plot.MetaForest_combined)

cat("\n✓ Forest plot created successfully!\n")
cat("  - Pooled ALL (red diamond) = Meta-analysis of all 4 studies\n")
cat("  - Pooled IP (orange diamond) = Meta-analysis of 3 IP studies\n")
cat("      (Shows IP pooled estimate even if not significant)\n")
cat("  - Individual studies (blue circles) = Study-specific effects\n\n")

# Show which taxa have significant IP pooled
cat("Pooled IP significance status:\n")
for(taxon in significant_taxa) {
  ip_res <- all_results_IP %>% filter(feature == taxon)
  if(nrow(ip_res) > 0) {
    if(ip_res$qval.fdr < 0.05) {
      cat("  ✓", taxon, "- IP pooled IS significant (FDR =", 
          sprintf("%.4f", ip_res$qval.fdr), ")\n")
    } else {
      cat("  ○", taxon, "- IP pooled NOT significant (FDR =", 
          sprintf("%.4f", ip_res$qval.fdr), ")\n")
    }
  }
}
cat("\n")

create_results_table <- function(plot_data_full) {
    table_data <- plot_data_full %>%
      # ADD: Group by taxon to extract pooled weights
      group_by(taxon_name) %>%
      mutate(
        # ADD: Get weights from Pooled ALL row
        pooled_w2f = weight_study2_food[study == "Pooled ALL"][1],
        pooled_w2ip = weight_study2_ip[study == "Pooled ALL"][1],
        pooled_w5 = weight_study5[study == "Pooled ALL"][1],
        pooled_w7 = weight_study7[study == "Pooled ALL"][1]
      ) %>%
      ungroup() %>%  # ADD: Ungroup before continuing
      mutate(
        Effect_CI = sprintf("%.2f [%.2f; %.2f]%s", 
                            yi, ci.lb, ci.ub, sig_marker),
        P_value = case_when(
          is.na(pval) ~ "—",
          pval < 0.001 ~ "<0.001",
          TRUE ~ sprintf("%.3f", pval)
        ),
        # CHANGED: Use pooled_* instead of weight_*
        Weight = case_when(
          study == "Pooled ALL" ~ "—",
          study == "Pooled IP" ~ "—",
          study == "Bitto 2016 (Diet)" & !is.na(pooled_w2f) ~ paste0(round(pooled_w2f, 1), "%"),
          study == "Bitto 2016 (IP)" & !is.na(pooled_w2ip) ~ paste0(round(pooled_w2ip, 1), "%"),
          study == "Yang 2022" & !is.na(pooled_w5) ~ paste0(round(pooled_w5, 1), "%"),
          study == "Han 2021" & !is.na(pooled_w7) ~ paste0(round(pooled_w7, 1), "%"),
          TRUE ~ "—"
        ),
        # ADD: Heterogeneity column
        Heterogeneity = case_when(
          study == "Pooled ALL" ~ sprintf("I²=%d%%, τ²=%.4f", round(i2), tau2),
          study == "Pooled IP" ~ sprintf("I²=%d%%, τ²=%.4f", round(i2), tau2),
          TRUE ~ "—"
        )
      ) %>%
      # CHANGED: Add Heterogeneity to select
      select(Taxon = taxon_name, Study = study, 
             `Effect [95% CI]` = Effect_CI, `P-value` = P_value, 
             Weight = Weight, Heterogeneity = Heterogeneity)
    return(table_data)
}

results_table <- create_results_table(results$plot_data)
print.data.frame(results_table)  # don't pass n=Inf here

### Only Sig, IP only ###
create_forest_ip_only <- function(mmuphin_fit_ip,
                                  sig_tab_ip,
                                  taxa_list) {
  
  all_plot_data <- data.frame()
  
  # For each taxon, compile data
  for(taxon in taxa_list) {
    
    # 1. Get POOLED IP meta-analysis result
    meta_result_ip <- sig_tab_ip %>% filter(feature == taxon)
    
    # 2. Get individual IP study data
    study_data <- data.frame()
    
    for (study in names(mmuphin_fit_ip$maaslin_fits)) {
      study_results <- mmuphin_fit_ip$maaslin_fits[[study]]
      
      if (!is.null(study_results)) {
        genus_row <- study_results[study_results$feature == taxon, ]
        
        if (nrow(genus_row) > 0) {
          study_data <- rbind(study_data, data.frame(
            study = study,
            yi = genus_row$coef,
            sei = genus_row$stderr,
            ci.lb = genus_row$coef - 1.96 * genus_row$stderr,
            ci.ub = genus_row$coef + 1.96 * genus_row$stderr,
            pval = genus_row$pval
          ))
        }
      }
    }
    
    # Create data frame for POOLED IP
    pooled_ip_df <- data.frame(
      study = "Pooled IP",
      yi = meta_result_ip$coef,
      ci.lb = meta_result_ip$coef - 1.96 * meta_result_ip$stderr,
      ci.ub = meta_result_ip$coef + 1.96 * meta_result_ip$stderr,
      is_summary = TRUE,
      taxon_name = taxon,
      sei = meta_result_ip$stderr,
      pval = meta_result_ip$pval,
      fdr = meta_result_ip$qval.fdr,
      i2 = meta_result_ip$I2,
      tau2 = meta_result_ip$tau2
    )
    
    # Create data frame for individual studies
    study_df <- data.frame(
      study = study_data$study,
      yi = study_data$yi,
      ci.lb = study_data$ci.lb,
      ci.ub = study_data$ci.ub,
      is_summary = FALSE,
      taxon_name = taxon,
      sei = study_data$sei,
      pval = study_data$pval,
      fdr = NA,
      i2 = NA,
      tau2 = NA
    )
    
    # Combine
    taxon_data <- rbind(pooled_ip_df, study_df)
    
    
    # Calculate weights (inverse variance)
    taxon_data$weight <- 1 / (taxon_data$sei^2)
    # Normalize weights to percentages within each taxon
    taxon_data$weight_pct <- (taxon_data$weight / sum(taxon_data$weight, na.rm = TRUE)) * 100
    
    
    all_plot_data <- rbind(all_plot_data, taxon_data)
  }
  
  # Recode study names
  all_plot_data$study <- recode(all_plot_data$study,
                                "Study2_IP" = "Bitto 2016",
                                "Study5" = "Yang 2022", 
                                "Study7" = "Han 2021"
  )
  
  # Set factor levels (reversed for plotting)
  all_plot_data$study <- factor(all_plot_data$study, 
                                levels = rev(c("Pooled IP",
                                               "Bitto 2016", 
                                               "Yang 2022", 
                                               "Han 2021")))
  
  # Add significance markers
  all_plot_data$sig_marker <- ""
  all_plot_data$sig_marker[all_plot_data$is_summary & all_plot_data$fdr < 0.05] <- "*"
  
  # Calculate GLOBAL x-axis limits
  global_x_min <- min(all_plot_data$ci.lb, na.rm = TRUE)
  global_x_max <- max(all_plot_data$ci.ub, na.rm = TRUE)
  global_x_range <- global_x_max - global_x_min
  global_x_limits <- c(global_x_min - 0.05 * global_x_range, 
                       global_x_max + 0.15 * global_x_range)
  
  # Create breaks
  break_start <- floor(global_x_min / 0.04) * 0.04
  break_end <- ceiling(global_x_max / 0.04) * 0.04
  global_x_breaks <- seq(break_start, break_end, by = 0.04)
  
  # Create plots for each taxon
  forest_plots <- list()
  
  for(i in seq_along(taxa_list)) {
    taxon <- taxa_list[i]
    is_first <- (i == 1)
    is_last <- (i == length(taxa_list))
    plot_data <- all_plot_data %>% filter(taxon_name == taxon)
    
    # Get heterogeneity stats
    meta_result_ip <- sig_tab_ip %>% filter(feature == taxon)
    i2_ip <- paste0(round(meta_result_ip$I2, 0), "%")
    tau2_ip <- round(meta_result_ip$tau2, 4)
    k_ip <- meta_result_ip$k
    
    # Add effect size text
    plot_data$effect_text <- sprintf("%.2f [%.2f; %.2f]%s",
                                     plot_data$yi, 
                                     plot_data$ci.lb, 
                                     plot_data$ci.ub,
                                     plot_data$sig_marker)
    plot_data$weight_text <- ifelse(is.na(plot_data$weight_pct), 
                                    "-", 
                                    sprintf("%.1f%%", plot_data$weight_pct))
    
    # Create plot
    p <- ggplot(plot_data, aes(x = yi, y = study, xmin = ci.lb, xmax = ci.ub)) +
      geom_vline(xintercept = global_x_breaks, color = "gray90") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
      geom_errorbarh(height = 0.08, linewidth = 0.8, 
                     aes(color = is_summary)) +
      geom_point(aes(size = is_summary, color = is_summary, shape = is_summary)) +
      geom_text(aes(label = effect_text), x = Inf, hjust = 0.7, size = 3) +
      geom_text(aes(label = weight_text), x = Inf, hjust = -1, size = 3) +  
      scale_x_continuous(
        limits = global_x_limits,
        breaks = global_x_breaks,
        minor_breaks = NULL
      ) +
      scale_y_discrete(expand = expansion(add = c(0.9, 1.1))) +  # Increased from 0.22 to 0.5
      theme_minimal() +
      coord_cartesian(clip = "off") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        axis.title.x = if(is_last) element_text(size = 10, face = "bold") else element_blank(),
        plot.margin = margin(6, 100, 20, 70)
      ) +
      scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "royalblue3")) +
      scale_size_manual(values = c("FALSE" = 5, "TRUE" = 5)) +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 16)) +
      labs(
        title = paste0(taxon),
        x = if(is_last) "Effect Size" else NULL
      ) +
      annotate("text", x = -Inf, y = -0.1,
               label = paste0("Heterogeneity: k=", k_ip, ", I²=", i2_ip, ", τ²=", tau2_ip), 
               hjust = 0.05, size = 3.5, color = "gray40")
    # Add column headers only to first plot
    if(is_first) {
      p <- p +
        annotate("text", x = Inf, y = Inf, 
                 label = "Effect [95% CI]", 
                 hjust = 0.7, vjust = -0.5, size = 3.5, fontface = "bold") +
        annotate("text", x = Inf, y = Inf, 
                 label = "Weight", 
                 hjust = -0.9, vjust = -0.5, size = 3.5, fontface = "bold")
    }
    
    forest_plots[[taxon]] <- p
  }
  
  # Combine plots
  if(length(forest_plots) > 0) {
    combined_plot <- forest_plots[[taxa_list[1]]]
    
    if(length(forest_plots) > 1) {
      for(i in 2:length(taxa_list)) {
        combined_plot <- combined_plot / forest_plots[[taxa_list[i]]]
      }
    }
    
    combined_plot <- combined_plot +
      plot_annotation(
        title = "IP-Route Significant Taxa",
        subtitle = "Meta-analysis of 3 IP-administration studies (* FDR < 0.05)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 11, hjust = 0.5),
                      panel.spacing = unit(0.18, "lines"))
      )
  }
  
  return(list(
    individual = forest_plots,
    combined = combined_plot,
    plot_data = all_plot_data
  ))
}

# RUN IT


# Get IP-significant taxa
ip_significant_taxa <- sig_results_IP$feature

# Create forest plot
results_ip <- create_forest_ip_only(
  mmuphin_fit_ip = fit_lm_meta_IP,
  sig_tab_ip = sig_results_IP,
  taxa_list = ip_significant_taxa
)

# Display
plot.MetaForest_IP <- results_ip$combined
print(plot.MetaForest_IP)

#### no weights###
create_sign_forest <- function(mmuphin_fit, sig_tab_data, taxa_list,
                               study_name_map = c("Study2_Food" = "Bitto 2016 (Diet)",
                                                  "Study2_IP" = "Bitto 2016 (IP)",
                                                  "Study5" = "Yang 2022",
                                                  "Study7" = "Han 2021")) {
  
  # helper to coerce possible Rle/list scalars to numeric
  to_num <- function(x) {
    if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
    as.numeric(x)
  }
  
  # coerce sig_tab_data to plain data.frame/tibble
  sig_df <- as.data.frame(sig_tab_data)
  
  all_plot_data <- data.frame()
  
  for (taxon in taxa_list) {
    # ---- meta row
    meta_rows <- sig_df[sig_df$feature == taxon, , drop = FALSE]
    if (nrow(meta_rows) > 0) {
      # take the first matching row; pull coef
      mcoef <- to_num(meta_rows$coef)[1]
      meta_row <- data.frame(
        study = "Pooled (RE)",
        yi = mcoef,
        is_summary = TRUE,
        taxon_name = taxon,
        stringsAsFactors = FALSE
      )
    } else {
      # no meta; still allow plot with studies
      meta_row <- data.frame(
        study = "Meta-analysis",
        yi = NA_real_,
        is_summary = TRUE,
        taxon_name = taxon,
        stringsAsFactors = FALSE
      )
    }
    
    # ---- study rows keep study 2 seperate
    study_data <- data.frame()
    for (study in names(mmuphin_fit$maaslin_fits)) {
      st_res <- mmuphin_fit$maaslin_fits[[study]]
      if (is.null(st_res)) next
      st_df <- as.data.frame(st_res)
      g_row <- st_df[st_df$feature == taxon, , drop = FALSE]
      if (nrow(g_row) > 0) {
        study_data <- rbind(
          study_data,
          data.frame(
            study = study,
            yi = to_num(g_row$coef),
            is_summary = FALSE,
            taxon_name = taxon,
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    taxon_df <- rbind(meta_row, study_data)
    all_plot_data <- rbind(all_plot_data, taxon_df)
  }
  
  # pretty study labels
  all_plot_data$study <- dplyr::recode(all_plot_data$study, !!!study_name_map, .default = all_plot_data$study)
  all_plot_data$study <- factor(all_plot_data$study,
                                levels = rev(c("Pooled (RE)","Bitto 2016 (Food)","Bitto 2016 (IP)","Yang 2022","Han 2021")))
  
  # direction-only labels
  all_plot_data <- all_plot_data %>%
    mutate(direction_label = case_when(
      is.na(yi) ~ "—",
      yi > 0 ~ "↑ Increase",
      yi < 0 ~ "↓ Decrease",
      TRUE ~ "0 / No change"
    ))
  all_plot_data$direction_label <- factor(all_plot_data$direction_label,
                                          levels = c("↓ Decrease","0 / No change","↑ Increase","—"))
  
  # build plots
  plots <- list()
  for (taxon in taxa_list) {
    plot_data <- all_plot_data %>% filter(taxon_name == taxon)
    p <- ggplot(plot_data, aes(x = direction_label, y = study)) +
      geom_point(aes(shape = is_summary), size = 3) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) + # circle vs diamond
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_blank()
      ) +
      labs(title = taxon, x = "Direction (no effect size)")
    plots[[taxon]] <- p
  }
  
  return(plots)
}

# Run
sign_plots <- create_sign_forest(
  mmuphin_fit = fit_lm_meta,
  sig_tab_data = sig_results,
  taxa_list = c("Bacteroides", "Muribaculum", "Candidatus_Arthromitus", "ASF356")
)
#### Investigation ####
# Extract exact numbers for all taxa
for(taxon in c("Bacteroides", "Muribaculum", "Candidatus_Arthromitus", "ASF356")) {
  cat("\n==================", taxon, "==================\n")
  
  # Get pooled meta-analysis result
  meta_result <- sig_results %>% filter(feature == taxon)
  pooled_se <- meta_result$stderr
  pooled_coef <- meta_result$coef
  pooled_ci_lower <- pooled_coef - 1.96 * pooled_se
  pooled_ci_upper <- pooled_coef + 1.96 * pooled_se
  pooled_ci_width <- pooled_ci_upper - pooled_ci_lower
  
  cat("POOLED (All studies):\n")
  cat("  Coefficient:", round(pooled_coef, 4), "\n")
  cat("  SE:", round(pooled_se, 4), "\n")
  cat("  95% CI: [", round(pooled_ci_lower, 4), ",", round(pooled_ci_upper, 4), "]\n")
  cat("  CI width:", round(pooled_ci_width, 4), "\n")
  cat("  I²:", round(meta_result$I2, 1), "%, τ²:", round(meta_result$tau2, 6), "\n\n")
  
  # Get individual study results
  cat("INDIVIDUAL STUDIES:\n")
  all_ses <- c()
  
  for (study in names(fit_lm_meta$maaslin_fits)) {
    study_results <- fit_lm_meta$maaslin_fits[[study]]
    if (!is.null(study_results)) {
      genus_row <- study_results[study_results$feature == taxon, ]
      if (nrow(genus_row) > 0) {
        study_se <- genus_row$stderr
        study_coef <- genus_row$coef
        study_ci_lower <- study_coef - 1.96 * study_se
        study_ci_upper <- study_coef + 1.96 * study_se
        study_ci_width <- study_ci_upper - study_ci_lower
        
        all_ses <- c(all_ses, study_se)
        
        cat("  ", study, ":\n")
        cat("    Coefficient:", round(study_coef, 4), "\n")
        cat("    SE:", round(study_se, 4), "\n")
        cat("    95% CI: [", round(study_ci_lower, 4), ",", round(study_ci_upper, 4), "]\n")
        cat("    CI width:", round(study_ci_width, 4), "\n")
      }
    }
  }
  
  # Compare pooled vs individual
  cat("\nCOMPARISON:\n")
  cat("  Pooled SE:", round(pooled_se, 4), "\n")
  cat("  Individual SEs:", paste(round(all_ses, 4), collapse=", "), "\n")
  cat("  Min individual SE:", round(min(all_ses), 4), "\n")
  cat("  Pooled SE < Min individual SE?", pooled_se < min(all_ses), "\n")
  
  # Manual calculation of what SE should be with I²=0%
  if(meta_result$I2 < 5) {  # Essentially zero heterogeneity
    expected_se <- 1/sqrt(sum(1/all_ses^2))
    cat("  Expected SE (fixed-effects):", round(expected_se, 4), "\n")
    cat("  Actual vs Expected:", round(pooled_se, 4), "vs", round(expected_se, 4), "\n")
  }
}

# Check if taxa exist in the abundance matrix
problem_taxa <- sig_results_IP$feature

for(taxon in problem_taxa) {
  cat("\n═══", taxon, "═══\n")
  
  # Check MaAsLin2 results per study
  cat("\nMaAsLin2 Results (included in differential abundance):\n")
  for(study in names(fit_lm_meta_IP$maaslin_fits)) {
    study_res <- fit_lm_meta_IP$maaslin_fits[[study]]
    if(!is.null(study_res)) {
      taxon_row <- study_res[study_res$feature == taxon, ]
      if(nrow(taxon_row) > 0) {
        cat("  ✓", study, ": coef =", round(taxon_row$coef, 4), 
            ", pval =", sprintf("%.2e", taxon_row$pval), "\n")
      } else {
        cat("  ✗", study, ": FILTERED OUT by MaAsLin2\n")
      }
    }
  }
  
  # Check raw abundances
  if(exists("genus_mat_adj") && taxon %in% rownames(genus_mat_adj)) {
    cat("\nRaw Abundance (batch-corrected genus_mat_adj):\n")
    
    for(study in c("Study2_IP", "Study5", "Study7")) {
      # Get samples for this study
      study_samples <- rapa_meta %>% 
        filter(StudyID_Split == study) %>% 
        pull(SampleID)
      
      study_samples <- study_samples[study_samples %in% colnames(genus_mat_adj)]
      
      if(length(study_samples) > 0) {
        abundances <- genus_mat_adj[taxon, study_samples]
        
        n_total <- length(abundances)
        n_detected <- sum(abundances > 0)
        prevalence <- (n_detected / n_total) * 100
        mean_all <- mean(abundances)
        mean_nonzero <- ifelse(n_detected > 0, mean(abundances[abundances > 0]), 0)
        
        cat(sprintf("  %s: %d/%d samples (%.1f%% prevalence)\n", 
                    study, n_detected, n_total, prevalence))
        cat(sprintf("           Mean (all) = %.6f, Mean (nonzero) = %.6f\n", 
                    mean_all, mean_nonzero))
      }
    }
    
    # Compare by treatment
    cat("\nAbundance by Treatment:\n")
    for(study in c("Study2_IP", "Study5", "Study7")) {
      study_meta <- rapa_meta %>% filter(StudyID_Split == study)
      
      for(trt in c("Control", "Rapamycin")) {
        trt_samples <- study_meta %>% 
          filter(StandardTreatment == trt) %>% 
          pull(SampleID)
        
        trt_samples <- trt_samples[trt_samples %in% colnames(genus_mat_adj)]
        
        if(length(trt_samples) > 0) {
          abundances <- genus_mat_adj[taxon, trt_samples]
          n_detected <- sum(abundances > 0)
          prevalence <- (n_detected / length(abundances)) * 100
          mean_val <- mean(abundances)
          
          cat(sprintf("  %s %s: %.1f%% prevalence, mean = %.6f\n", 
                      study, trt, prevalence, mean_val))
        }
      }
    }
    
  } else {
    cat("\n⚠️  NOT FOUND in genus_mat_adj!\n")
  }
}


cat("\n\n═══ PREVALENCE SUMMARY TABLE ═══\n")

prevalence_table <- data.frame()

for(taxon in problem_taxa) {
  if(taxon %in% rownames(genus_mat_adj)) {
    for(study in c("Study2_IP", "Study5", "Study7")) {
      study_samples <- rapa_meta %>% 
        filter(StudyID_Split == study) %>% 
        pull(SampleID)
      
      study_samples <- study_samples[study_samples %in% colnames(genus_mat_adj)]
      
      if(length(study_samples) > 0) {
        abundances <- genus_mat_adj[taxon, study_samples]
        prev <- (sum(abundances > 0) / length(abundances)) * 100
        mean_abund <- mean(abundances)
        
        # Check if in MaAsLin results
        in_maaslin <- taxon %in% fit_lm_meta_IP$maaslin_fits[[study]]$feature
        
        prevalence_table <- rbind(prevalence_table, data.frame(
          Taxon = taxon,
          Study = study,
          Prevalence = round(prev, 1),
          Mean_Abundance = sprintf("%.6f", mean_abund),
          In_MaAsLin = in_maaslin
        ))
      }
    }
  }
}

print(prevalence_table)

cat("\n═══ KEY FINDINGS ═══\n")
cat("MaAsLin2 default filters:\n")
cat("  - Minimum prevalence: 10% (taxa must be in ≥10% of samples)\n")
cat("  - Minimum abundance: typically 0.0001 (relative abundance)\n\n")
cat("If a taxon has <10% prevalence in a study, MaAsLin2 excludes it.\n")
cat("This is why some taxa appear in Yang/Han but not Bitto.\n\n")



#### 11.5 Predictivity #####
#### PREDICTIVE ANALYSIS USING META-ANALYSIS RESULTS 
#### 12. ADD DOSAGE INFORMATION AND BATCH CORRECTION ####

# Create a copy of the merged genus phyloseq object
merged_genus_dose <- merged_genus_normal

# Extract the metadata with a new name
dose_meta <- as(sample_data(merged_genus_dose), "data.frame")
dose_meta$StudyID <- factor(trimws(dose_meta$StudyID))
dose_meta$StudyID_Split  <- factor(dose_meta$StudyID_Split)
dose_meta$SampleID <- rownames(dose_meta)
dose_meta$StandardTreatment <- factor(trimws(dose_meta$StandardTreatment))
dose_meta$AdminRoute <- factor(trimws(dose_meta$AdminRoute))


# Add dosage information
# 1) Initialize dose columns
dose_meta$DailyDose_mgkg <- 0
dose_meta$TotalDose_mgkg <- 0
dose_meta$DurationDays <- 0

# 2) Set duration values
dose_meta$DurationDays[dose_meta$StudyID == "Study2"] <- 90
dose_meta$DurationDays[dose_meta$StudyID == "Study5"] <- 14
dose_meta$DurationDays[dose_meta$StudyID == "Study7"] <- 84

# 3) Set daily and total doses based on study and administration route
# For Study2 - IP group
dose_meta$DailyDose_mgkg[dose_meta$StudyID == "Study2" & dose_meta$AdminRoute == "IP"] <- 8.0
dose_meta$TotalDose_mgkg[dose_meta$StudyID == "Study2" & dose_meta$AdminRoute == "IP"] <- 720

# For Study2 - feeding group
dose_meta$DailyDose_mgkg[dose_meta$StudyID == "Study2" & dose_meta$AdminRoute == "Food"] <- 11.0
dose_meta$TotalDose_mgkg[dose_meta$StudyID == "Study2" & dose_meta$AdminRoute == "Food"] <- 990

# For Study5
dose_meta$DailyDose_mgkg[dose_meta$StudyID == "Study5"] <- 2.5
dose_meta$TotalDose_mgkg[dose_meta$StudyID == "Study5"] <- 35

# For Study7
dose_meta$DailyDose_mgkg[dose_meta$StudyID == "Study7"] <- 0.5
dose_meta$TotalDose_mgkg[dose_meta$StudyID == "Study7"] <- 42

# 4) Set doses to 0 for control samples
dose_meta$DailyDose_mgkg[dose_meta$StandardTreatment == "Control"] <- 0
dose_meta$TotalDose_mgkg[dose_meta$StandardTreatment == "Control"] <- 0

# 5) Create log-transformed versions for analysis
dose_meta$LogDailyDose <- log10(dose_meta$DailyDose_mgkg + 1)  # Add 1 to avoid log(0)
dose_meta$LogTotalDose <- log10(dose_meta$TotalDose_mgkg + 1)  # Add 1 to avoid log(0)

dose_meta$StandardTreatment <- factor(dose_meta$StandardTreatment, levels=c("Control", "Rapamycin"))
dose_meta$AdminRoute <- factor(dose_meta$AdminRoute, levels=c("Food", "IP"))

table(dose_meta$DailyDose_mgkg, dose_meta$StandardTreatment)
table(dose_meta$TotalDose_mgkg, dose_meta$StandardTreatment)

#1) Make the rownames of dose_meta match the phyloseq sample IDs
rownames(dose_meta) <- dose_meta$SampleID

# 6) Update the sample data in the phyloseq object
sample_data(merged_genus_dose) <- sample_data(dose_meta)

# 7) Print summary of dosage information
cat("Dosage information by study and administration route:\n")
print(aggregate(cbind(DailyDose_mgkg, TotalDose_mgkg) ~ StudyID + AdminRoute, 
                data = dose_meta[dose_meta$StandardTreatment == "Rapamycin", ], 
                FUN = function(x) round(mean(x, na.rm = TRUE), 1)))

# StudyID AdminRoute DailyDose_mgkg TotalDose_mgkg
# 1  Study2       Food           11.0            990
# 2  Study2         IP            8.0            720
# 3  Study5         IP           17.5            245
# 4  Study7         IP            0.5            42

# Confirm metadata and genus matrix alignment
stopifnot(all(rownames(dose_meta) == sample_names(merged_genus_dose)))


#### 13.BATCH CORRECTION WITH DOSE INFORMATION ####

# 1) Extract genus × sample matrix
genus_mat_dose <- as(otu_table(merged_genus_dose), "matrix")
if (!taxa_are_rows(merged_genus_dose)) genus_mat_dose <- t(genus_mat_dose)

stopifnot(
  !anyNA(genus_mat_dose),
  all(genus_mat_dose >= 0)
)

# 2) Define your covariates - now including dosage information
covar_cols_dose <- c("AdminRoute") #run without a logtotaldose
#"StandardTreatment"

# 3) Run adjust_batch on the FULL genus matrix (not just core genera)
fit_genus_dose <- MMUPHin::adjust_batch(
  feature_abd = genus_mat_dose,  # Use full matrix with all genera
  batch = "StudyID_Split",
  covariates = NULL,
  data = dose_meta,
  control = list(method = "combat", zero_correct = TRUE, verbose = FALSE, mean.only= TRUE)
)

genus_mat_adj <- fit_genus_dose$feature_abd_adj


# 4) Check if batch correction actually modified the data
# Calculate the sum of absolute differences between original and corrected matrices
diff_sum <- sum(abs(genus_mat_dose - fit_genus_dose$feature_abd_adj))
cat("Sum of absolute differences between original and corrected matrices:", diff_sum, "\n")
#419591 

# If difference is very small, batch correction might not be working
if (diff_sum < 1e-10) {
  warning("Batch correction may not be working - very small difference between original and corrected matrices")
}

# 5) Create a NEW phyloseq object for the batch-corrected data (don't overwrite the original)
batch_corrected_ps_dose <- merged_genus_dose  # Start with a copy
otu_table(batch_corrected_ps_dose) <- otu_table(fit_genus_dose$feature_abd_adj, taxa_are_rows = TRUE)

cat("✓ Batch correction with dose adjustment complete: ",
    ntaxa(batch_corrected_ps_dose), "genera × ",
    nsamples(batch_corrected_ps_dose), "samples\n")
#✓ Batch correction with dose adjustment complete:  112 genera ×  54 samples

# 6) Now you can directly compare the original and batch-corrected data
original_D <- distance(merged_genus_dose, method="bray")
corrected_D <- distance(batch_corrected_ps_dose, method="bray")

# Check if distance matrices are identical
dist_diff <- sum(abs(as.matrix(original_D) - as.matrix(corrected_D)))
cat("Difference between original and batch-corrected distance matrices:", dist_diff, "\n")
#310.0229  

# Check after adjustment
cat("Summary of adjusted genus counts (should be non-negative, no NAs):\n")
summary(as.vector(genus_mat_adj))

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000000 0.000000 0.000000 0.008924 0.003039 0.713924 

# fit_lm_meta$meta_fits %>% 
#   filter(feature == "Muribaculum")


#### 14. PCOA WITH DOSE-ADJUSTED DATA ####
genus_before <- genus_mat_dose  
genus_after <- as(otu_table(batch_corrected_ps_dose), "matrix")
if (!taxa_are_rows(batch_corrected_ps_dose))
    genus_after <- t(genus_after)
# 1) Compute distances & ordinations
D_dose_before <- vegdist(t(genus_before), method="bray")
D_dose_after <- vegdist(t(genus_after), method="bray")

ord_dose_before <- ordinate(merged_genus_dose, "PCoA", D_dose_before)
ord_dose_after <- ordinate(batch_corrected_ps_dose, "PCoA", D_dose_after)

# 2) Extract ordination data with sample metadata
data_dose_before <- plot_ordination(merged_genus_dose, ord_dose_before)$data
data_dose_after <- plot_ordination(batch_corrected_ps_dose, ord_dose_after)$data
# Calculate the range of each axis
x_range <- range(c(data_dose_before$Axis.1, data_dose_after$Axis.1))
y_range <- range(c(data_dose_before$Axis.2, data_dose_after$Axis.2))

# Add 20% padding on each side to accommodate ellipses
x_padding <- diff(x_range) * 0.1
y_padding <- diff(y_range) * 0.1

# Create expanded limits
xlims <- c(x_range[1] - x_padding, x_range[2] + x_padding)
ylims <- c(y_range[1] - y_padding, y_range[2] + y_padding)

# 3) Create a plotting function to ensure consistent styling
create_pcoa_plot <- function(data, ord, title, xlims, ylims) { 
  # Use the correct ordination object for axis labels
  var_explained1 <- round(100 * ord$values$Relative_eig[1], 1)
  var_explained2 <- round(100 * ord$values$Relative_eig[2], 1)
  
  # Create combined grouping variable for ellipses
  data$Group <- paste(data$StudyID_Split, data$StandardTreatment, sep="_")
  
  # Create the plot
  p <- ggplot(data, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = StudyID_Split, shape = StandardTreatment), size = 3) +
    stat_ellipse(aes(group = Group, color = StudyID_Split, linetype = StandardTreatment), size = 1) +
    stat_ellipse(aes(group = Group, fill = StudyID_Split), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(
      title = title,
      subtitle = "Bray-Curtis Dissimilarity",
      color = "Study",
      shape = "Treatment",
      linetype = "Treatment"
    ) +
    theme_bw() +
    coord_equal(xlim = xlims, ylim = ylims)+
    xlab(paste0("PCoA1 [", var_explained1, "%]")) +
    ylab(paste0("PCoA2 [", var_explained2, "%]"))
  
  return(p)
}

# 4) Create the plots with correct ordination objects
p_dose_before <- create_pcoa_plot(data_dose_before, ord_dose_before, "Before Batch Correction (with Dose)",xlims, ylims)
p_dose_after <- create_pcoa_plot(data_dose_after, ord_dose_after, "After Batch Correction (with Dose)", xlims, ylims)

# 5) Combine the plots
library(patchwork)
combined_plot_dose <- p_dose_before + p_dose_after + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Bray-Curtis PCoA Before & After Batch Correction (with Dose Adjustment)")

# 6) Display the combined plot
print(combined_plot_dose)

plot.combined_plot_dose <- combined_plot_dose
plot.p_dose_before <- p_dose_before
plot.p_dose_after <- p_dose_after

#### 15. PERMANOVA Test with Dose-Adjusted Batch Correction ######
# 1) Get the matrices 
genus_before <- genus_mat_dose  
genus_after <- fit_genus_dose$feature_abd_adj

# 2) Calculate Bray-Curtis distances using vegan
D_before <- vegdist(t(genus_before), method="bray")
D_after <- vegdist(t(genus_after), method="bray")

# 3) Use dose_meta which has all the dose information
meta_df <- dose_meta  

# 4) VERIFY SAMPLE ORDER MATCHES METADATA
sample_names_matrix <- colnames(genus_before)
sample_names_meta <- rownames(dose_meta)
if(!all(sample_names_matrix == sample_names_meta)) {
  warning("Sample order mismatch between matrix and metadata!")
  cat("Matrix samples:", head(sample_names_matrix), "...\n")
  cat("Metadata samples:", head(sample_names_meta), "...\n")
} else {
  cat("✓ Sample order matches between matrix and metadata\n")
}

# Test study effect
set.seed(123)
study_before <- adonis2(D_before ~ StudyID_Split, data = dose_meta)
print(study_before)
set.seed(123)
study_after <- adonis2(D_after ~ StudyID_Split, data = dose_meta)
print(study_after)

# Test treatment effect 
set.seed(123)
treatment_before <- adonis2(D_before ~ StandardTreatment, data = dose_meta)
print(treatment_before)
set.seed(123)
treatment_after <- adonis2(D_after ~ StandardTreatment, data = dose_meta)
print(treatment_after)

# Test combined effect
set.seed(123)
combined_before <- adonis2(D_before ~ StudyID_Split + StandardTreatment, data = dose_meta, by= "margin")
print(combined_before)
set.seed(123)
combined_after <- adonis2(D_after ~ StudyID_Split + StandardTreatment, data = dose_meta, by = "margin")
print(combined_after)

# Dose effect (only in treated samples)
dose_samples <- dose_meta$StandardTreatment == "Rapamycin"
dose_meta_treated <- dose_meta[dose_samples, ]
dose_D_before <- as.dist(as.matrix(D_before)[dose_samples, dose_samples])
dose_D_after <- as.dist(as.matrix(D_after)[dose_samples, dose_samples])

cat("\nDose effect before batch correction (treated samples only):\n")
set.seed(123)
dose_before <- adonis2(dose_D_before ~ LogTotalDose, data = dose_meta_treated)
print(dose_before)

cat("\nDose effect after batch correction (treated samples only):\n")
set.seed(123)
dose_after <- adonis2(dose_D_after ~ LogTotalDose, data = dose_meta_treated)
print(dose_after)

# Use same extract_stats function
extract_stats <- function(adonis_result, row = 1) {
  if(row > nrow(adonis_result) || is.na(adonis_result$F[row])) {
    return(data.frame(
      Df = NA, 
      R2 = NA, 
      R2_percent = "ERROR",
      F_stat = NA, 
      P_val = NA, 
      Sig = "ERROR"
    ))
  }
  
  r2 <- adonis_result$R2[row]
  p <- adonis_result$`Pr(>F)`[row]
  
  data.frame(
    Df = adonis_result$Df[row],
    R2 = round(r2, 3),
    R2_percent = paste0(round(r2*100, 1), "%"),
    F_stat = round(adonis_result$F[row], 2),
    P_val = round(p, 4),
    Sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "†", TRUE ~ "ns")
  )
}

# Create results table
results_table_dose <- rbind(
  cbind(Effect = "Study", Timing = "Before", extract_stats(study_before, 1)),
  cbind(Effect = "Study", Timing = "After", extract_stats(study_after, 1)),
  cbind(Effect = "Treatment", Timing = "Before", extract_stats(treatment_before, 1)),
  cbind(Effect = "Treatment", Timing = "After", extract_stats(treatment_after, 1)),
  cbind(Effect = "Treatment (adj)",  Timing = "Before",
        extract_stats(combined_before, which(rownames(combined_before) == "StandardTreatment"))),
  cbind(Effect = "Treatment (adj)",  Timing = "After",
        extract_stats(combined_after,  which(rownames(combined_after)  == "StandardTreatment"))),
  cbind(Effect = "Dose (treated)", Timing = "Before", extract_stats(dose_before, 1)),
  cbind(Effect = "Dose (treated)", Timing = "After", extract_stats(dose_after, 1))
)

# Display results
cat("\n=== PERMANOVA RESULTS (Dose-Adjusted Batch Correction) ===\n")
print(results_table_dose)
table.results_dose_Permanova <- results_table_dose

set.seed(123)
adonis2(dose_D_after ~ StudyID_Split + LogTotalDose,
        data = dose_meta_treated, by = "margin")


#### 16. PERMDISPS####
# Simple function to extract PERMDISP stats
get_permdisp <- function(disp_obj) {
  result <- vegan::permutest(disp_obj, permutations = 999)
  f_stat <- as.numeric(result$statistic)
  p_val <- result$tab$`Pr(>F)`[1]  # P-value is here, not in result$p.value!
  
  data.frame(F_stat = round(f_stat, 4), 
             P_val = round(p_val, 5))
  
}

## 1. CONTROL vs RAPAMYCIN (all samples) ───────────────────────────────────────
set.seed(123)
disp_before  <- betadisper(D_before,  dose_meta$StandardTreatment)
set.seed(123)
disp_after   <- betadisper(D_after,   dose_meta$StandardTreatment)

cat("\nPERMDISP  |  StandardTreatment  |  BEFORE correction\n")
set.seed(123)
print( vegan::permutest(disp_before, permutations = 9999) )

cat("\nPERMDISP  |  StandardTreatment  |  AFTER  correction\n")
set.seed(123)
print( vegan::permutest(disp_after,  permutations = 9999) )


### 2. RAPAMYCIN-ONLY DOSE GRADIENT ─────────────────────────────────────────────
## betadisper needs a *factor*, so bin LogTotalDose into tertiles (“Low/Med/High”);
## change ‘breaks = 3’ to 2 or 4 if you prefer coarser/finer bins.
dose_bins <- cut(dose_meta_treated$LogTotalDose,
                 breaks = 3,
                 labels = c("Low", "Medium", "High"))
set.seed(123)
disp_dose_before <- betadisper(dose_D_before, dose_bins)
set.seed(123)
disp_dose_after  <- betadisper(dose_D_after,  dose_bins)

cat("\nPERMDISP  |  Dose tertiles (Rapamycin only)  |  BEFORE correction\n")
set.seed(123)
print( vegan::permutest(disp_dose_before, permutations = 9999) )

cat("\nPERMDISP  |  Dose tertiles (Rapamycin only)  |  AFTER  correction\n")
set.seed(123)
print( vegan::permutest(disp_dose_after,  permutations = 9999) )

# Add the missing Study objects:
set.seed(123)
disp_study_before <- betadisper(D_before, dose_meta$StudyID_Split)
set.seed(123)
disp_study_after  <- betadisper(D_after,  dose_meta$StudyID_Split)

# Create simple table
permdisp_table <- rbind(
  cbind(Effect = "Treatment", Timing = "Before", get_permdisp(disp_before)),
  cbind(Effect = "Treatment", Timing = "After",  get_permdisp(disp_after)),
  cbind(Effect = "Study", Timing = "Before", get_permdisp(disp_study_before)),
  cbind(Effect = "Study", Timing = "After",  get_permdisp(disp_study_after)),
  cbind(Effect = "Dose", Timing = "Before", get_permdisp(disp_dose_before)),
  cbind(Effect = "Dose", Timing = "After",  get_permdisp(disp_dose_after))
)
cat("\n=== PERMDISP RESULTS ===\n")
print(permdisp_table)
table.permidsp_table <- permdisp_table

#### 17.Dose-Response Relationship #########

# Create dose-response visualization
# make sure 'meta_after' exists
meta_after <- as(sample_data(batch_corrected_ps_dose), "data.frame")

# ensure split study is available
stopifnot("StudyID_Split" %in% colnames(meta_after))

# identify groups
control_samples  <- meta_after$StandardTreatment == "Control"
rapa_samples     <- meta_after$StandardTreatment == "Rapamycin"

# matrix: samples x features (to match row indexing by samples below)
X_after <- t(genus_mat_adj)

# Calculate centroid for control samples (mean abundance of each genus)
control_centroid <- colMeans(t(genus_mat_adj)[control_samples,drop=FALSE ])

# Calculate Bray-Curtis distance from each Rapamycin sample to control centroid
dist_to_control <- sapply(which(rapa_samples), function(i) {
  sample_vec <- t(genus_mat_adj)[i, ]
  vegdist(rbind(sample_vec, control_centroid), method = "bray")[1]
})

# Create data frame for plotting
plot_data <- data.frame(
  SampleID = rownames(meta_after)[rapa_samples],
  LogTotalDose = meta_after$LogTotalDose[rapa_samples],
  TotalDose = meta_after$TotalDose_mgkg[rapa_samples],
  StudyID = meta_after$StudyID[rapa_samples],
  AdminRoute = meta_after$AdminRoute[rapa_samples],
  Distance = dist_to_control
)

# Clean up StudyID labels (remove any trailing characters like "a")
plot_data$StudyID <- gsub("[a-z]$", "", plot_data$StudyID)
plot_data$StudyID <- trimws(plot_data$StudyID)

# Create detailed study labels that differentiate dose-route combinations
plot_data <- plot_data %>%
  mutate(
    detailed_label = case_when(
      # For studies with multiple doses/routes, show dose and route
      StudyID == "Study2" ~ paste0(StudyID, " (", round(TotalDose), "mg/kg, ", AdminRoute, ")"),
      # For studies with single dose, just show study and dose
      TRUE ~ paste0(StudyID, " (", round(TotalDose), "mg/kg)")
    )
  )

# Calculate correlation statistics
cor_result <- cor.test(plot_data$LogTotalDose, plot_data$Distance, method = "spearman")
r_value <- cor_result$estimate
p_value <- cor_result$p.value

# Create dose-response plot
p_dose_response <- ggplot(plot_data, aes(x = LogTotalDose, y = Distance, color = detailed_label, shape = AdminRoute)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  # Add text labels with actual dose values
  labs(
    title = "Dose-Response Relationship Between Rapamycin and Microbiome Composition",
    subtitle = paste0("\np-value = ", round(p_value, 3)),
    x = "Log10 Total Dose (mg/kg)",
    y = "Bray-Curtis Distance to Control Centroid",
    color = "Study (Dose Range)",
    shape = "Admin Route"
  ) +
  theme_bw()

# Add correlation info to the plot
p_dose_response <- p_dose_response +
  annotate("text", x = min(plot_data$LogTotalDose) + 0.1, y = max(plot_data$Distance) - 0.02,
           label = paste0("Spearman's rho = ", round(r_value, 3)),
           hjust = 0, vjust = 1, size = 3.5)

# Display the plot
print(p_dose_response)

# Analysis: Dose-Response Relationship
# Method: Spearman Correlation  
# Variable_X: Log10 Total Dose (mg/kg)
# Variable_Y: Bray-Curtis Distance to Control Centroid
# Sample_Size: 27
# Spearmans_Rho: 0.200
# P_Value: 0.318
# Significance: ns
# Interpretation: No significant dose-response

#### 18. LM META WITH DOSE-ADJUSTED DATA ####
# 1) Extract the batch-corrected abundance data from your phyloseq object
genus_mat_adj <- as(otu_table(batch_corrected_ps_dose), "matrix")
if (!taxa_are_rows(batch_corrected_ps_dose)) genus_mat_adj <- t(genus_mat_adj)

# 2) Make sure your metadata is correctly formatted
rapa_meta <- as(sample_data(batch_corrected_ps_dose), "data.frame")
rapa_meta$StudyID <- factor(trimws(rapa_meta$StudyID))
rapa_meta$StandardTreatment <- factor(trimws(rapa_meta$StandardTreatment))
rapa_meta$AdminRoute <- factor(trimws(rapa_meta$AdminRoute))

# 3) Verify that all samples match between abundance matrix and metadata
if(!all(colnames(genus_mat_adj) %in% rownames(rapa_meta))) {
  warning("Some samples in abundance matrix are missing from metadata")
  # Keep only matching samples
  genus_mat_adj <- genus_mat_adj[, colnames(genus_mat_adj) %in% rownames(rapa_meta)]
}
if(!all(rownames(rapa_meta) %in% colnames(genus_mat_adj))) {
  warning("Some samples in metadata are missing from abundance matrix")
  # Keep only matching samples
  rapa_meta <- rapa_meta[rownames(rapa_meta) %in% colnames(genus_mat_adj), ]
}

#### 19.TREATMENT EFFECT ANALYSIS (not usable due to missing factors)######
# Try with default settings and CPLM

# Run lm_meta with more robust settings
fit_treatment_meta <- MMUPHin::lm_meta(
  feature_abd = genus_mat_adj,     # Your abundance matrix (genera in rows, samples in columns)
  batch = "StudyID",               # Study identifier column
  exposure = "LogTotalDose",  # Treatment variable 
  covariates = "AdminRoute",               # No additional covariates
  data = rapa_meta,                # Metadata dataframe
  control = list(
    normalization = "NONE",         # Total Sum Scaling normalization
    transform = "AST",             # Arcsin Square Root transformation
    analysis_method = "LM",      # Compound Poisson Linear Model
    random_effects = TRUE,         # Use random effects model
    verbose = TRUE                 # Print progress
  )
)
  
# Extract significant results and filter out unclassified genera
# sig_features <- fit_treatment_meta$meta_fits %>% 
#   filter(pval < 0.1) %>%             # Filter to genera with p < 0.1
#   filter(!grepl("unclassified", feature, ignore.case = TRUE)) %>%  # Remove unclassified genera
#   filter(!grepl("uncultured", feature, ignore.case = TRUE)) %>%    # Remove uncultured genera
#   filter(!grepl("unknown", feature, ignore.case = TRUE)) %>%       # Remove unknown genera  
#   arrange(coef)                      # Sort by coefficient value

sig_features <- fit_treatment_meta$meta_fits %>%
  mutate(qval.fdr = p.adjust(pval, method = "fdr"))%>% 
filter(qval.fdr < 0.05,
         !grepl("uncultured|unclassified|Unclassified|Other|Mitochondria", feature, ignore.case=TRUE)
  )%>%
 arrange(coef)  # Sort by effect size instead

print(sig_features)

# Add this section after you create sig_features but before the plot:
sig_features$p_label <- case_when(
  sig_features$pval < 0.001 ~ "***",
  sig_features$pval < 0.01 ~ "**", 
  sig_features$pval < 0.05 ~ "*",
  TRUE ~ ""
)

# Make feature names cleaner for plotting
sig_features$clean_name <- gsub("s__", "", sig_features$feature)
sig_features$clean_name <- gsub("_", " ", sig_features$clean_name)

# Make feature a factor to preserve sorting in plot
sig_features$clean_name <- factor(sig_features$clean_name, levels = sig_features$clean_name)


# Create treatment effect plot
treatment_effect_plot <- ggplot(sig_features, aes(x = clean_name, y = coef)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  geom_errorbar(aes(ymin = coef - 1.96 * stderr, ymax = coef + 1.96 * stderr), width = 0.2) +
  geom_text(aes(label = p_label), hjust = -0.2, size = 4) +  # <- ADD THIS LINE
  coord_flip() +
  labs(
    title = "Effect of Rapamycin on Microbial Genera",
    subtitle = "Meta-analysis across three studies (n=27 rapamycin samples, FDR < 0.05)",
    x = "",
    y = "Effect Size (Coefficient)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)



# Display and save plot
print(treatment_effect_plot)

plot.treatment_effect_plot <- treatment_effect_plot

#### 20. Analyze dose-response relationship (only in rapamycin samples) #####

# Check the complete structure of dose_meta
cat("=== CHECKING DATA STRUCTURE ===\n")
cat("Total samples in dose_meta:", nrow(dose_meta), "\n")
cat("Unique StudyID_Split values:\n")
print(table(dose_meta$StudyID_Split))

cat("\nSample distribution by StudyID_Split and Treatment:\n")
print(table(dose_meta$StudyID_Split, dose_meta$StandardTreatment))

cat("\nDose values by StudyID_Split (Rapamycin only):\n")
rapa_dose_summary <- dose_meta %>%
  filter(StandardTreatment == "Rapamycin") %>%
  group_by(StudyID_Split) %>%
  summarise(
    n = n(),
    DailyDose = unique(DailyDose_mgkg),
    TotalDose = unique(TotalDose_mgkg),
    LogTotalDose = round(unique(LogTotalDose), 2)
  )
print(as.data.frame(rapa_dose_summary))

# Check if genus_mat_adj exists
if (!exists("genus_mat_adj")) {
  stop("genus_mat_adj not found. Please run section 13 (batch correction) first.")
}

cat("\ngenus_mat_adj dimensions:", dim(genus_mat_adj), "\n")
cat("Sample alignment check:", 
    all(colnames(genus_mat_adj) == rownames(dose_meta)), "\n\n")

# Function to run linear model for each genus with SPLIT study as covariate
run_dose_model <- function(genus_name) {
  # Extract abundances from batch-corrected matrix
  abund <- genus_mat_adj[genus_name, ]
  
  # Subset to Rapamycin samples only
  rapa_idx <- dose_meta$StandardTreatment == "Rapamycin"
  
  if (sum(rapa_idx) == 0) return(NULL)
  
  # Extract data for Rapamycin samples
  rapa_abund <- abund[rapa_idx]
  rapa_meta_subset <- dose_meta[rapa_idx, , drop = FALSE]
  
  # Remove NA values
  keep <- !is.na(rapa_abund) & !is.na(rapa_meta_subset$LogTotalDose)
  
  if (sum(keep) < 5) return(NULL)  # Need sufficient samples
  
  rapa_abund_clean <- rapa_abund[keep]
  rapa_meta_clean <- rapa_meta_subset[keep, , drop = FALSE]
  
  # Check variance
  if (var(rapa_abund_clean, na.rm = TRUE) < 1e-10) return(NULL)
  
  # Standardize LogTotalDose
  rapa_meta_clean$LogTotalDose_std <- as.numeric(scale(rapa_meta_clean$LogTotalDose))
  
  # Drop unused levels
  rapa_meta_clean$StudyID_Split <- droplevels(rapa_meta_clean$StudyID_Split)
  
  # Fit model with StudyID_Split as covariate
  # This treats Study2_IP, Study2_Food, Study5, Study7 as separate groups
  if (length(unique(rapa_meta_clean$StudyID_Split)) > 1) {
    model <- lm(rapa_abund_clean ~ LogTotalDose_std + StudyID_Split, data = rapa_meta_clean)
  } else {
    model <- lm(rapa_abund_clean ~ LogTotalDose_std, data = rapa_meta_clean)
  }
  
  # Extract results
  summ <- summary(model)
  
  if (!"LogTotalDose_std" %in% rownames(summ$coefficients)) return(NULL)
  
  return(data.frame(
    feature = genus_name,
    coef = summ$coefficients["LogTotalDose_std", "Estimate"],
    stderr = summ$coefficients["LogTotalDose_std", "Std. Error"],
    pval = summ$coefficients["LogTotalDose_std", "Pr(>|t|)"],
    n_samples = sum(keep),
    stringsAsFactors = FALSE
  ))
}

# Run analysis for all genera
cat("=== RUNNING DOSE-RESPONSE ANALYSIS ===\n")
valid_genera <- rownames(genus_mat_adj)
cat("Testing", length(valid_genera), "genera...\n")

dose_results_list <- lapply(valid_genera, function(g) {
  tryCatch(run_dose_model(g), error = function(e) NULL)
})

# Combine results
non_null <- Filter(Negate(is.null), dose_results_list)

if (length(non_null) == 0) {
  stop("No valid models. Check data alignment.")
}

cat("Successfully fit models for", length(non_null), "genera\n\n")

dose_results <- dplyr::bind_rows(non_null)

# Apply FDR correction
dose_results$qval <- p.adjust(dose_results$pval, method = "BH")

# Show summary
cat("P-value distribution:\n")
print(summary(dose_results$pval))
cat("\nQ-value distribution:\n")
print(summary(dose_results$qval))

# Get significant results
sig_dose <- dose_results %>% 
  filter(qval < 0.05) %>%
  filter(!grepl("unclassified|uncultured|unknown|Muribaculaceae", feature, ignore.case = TRUE)) %>%
  arrange(coef)

if (nrow(sig_dose) > 0) {
  # Add clean names
  sig_dose$clean_name <- gsub("s__|_", " ", sig_dose$feature)
  sig_dose$clean_name <- factor(sig_dose$clean_name, levels = sig_dose$clean_name)
  
  cat("\n=== SIGNIFICANT DOSE-RESPONSE EFFECTS (FDR < 0.05) ===\n")
  cat("Found", nrow(sig_dose), "significant genera:\n\n")
  print(sig_dose)
  
} else {
  cat("\n=== NO SIGNIFICANT EFFECTS AT FDR < 0.05 ===\n")
  
  # Try relaxed threshold
  sig_dose <- dose_results %>% 
    filter(qval < 0.1) %>%
    filter(!grepl("unclassified|uncultured|unknown|Muribaculaceae", feature, ignore.case = TRUE)) %>%
    arrange(coef)
  
  if (nrow(sig_dose) > 0) {
    sig_dose$clean_name <- gsub("s__|_", " ", sig_dose$feature)
    sig_dose$clean_name <- factor(sig_dose$clean_name, levels = sig_dose$clean_name)
    
    cat("\nGenera with FDR < 0.1:\n")
    print(sig_dose)
  }
}

# If you have significant results, continue with visualization
if(nrow(sig_dose) > 0) {
  # Your visualization code here (the plots from your original code)
  # ...
}

#### Dose Analysis IP only 

# Filter to IP-only samples
dose_meta_IP <- dose_meta %>%
  filter(StudyID_Split != "Study2_Food") %>%
  droplevels()

cat("=== IP-ONLY DATA STRUCTURE ===\n")
cat("Total samples in dose_meta_IP:", nrow(dose_meta_IP), "\n")
cat("Studies included:", paste(unique(dose_meta_IP$StudyID_Split), collapse=", "), "\n\n")

cat("Sample distribution by StudyID_Split and Treatment:\n")
print(table(dose_meta_IP$StudyID_Split, dose_meta_IP$StandardTreatment))

# Check dose-study collinearity
dose_by_study_IP <- dose_meta_IP %>%
  filter(StandardTreatment == "Rapamycin") %>%
  group_by(StudyID_Split) %>%
  summarise(
    sd_log_dose = sd(LogTotalDose, na.rm = TRUE),
    n_unique_doses = n_distinct(LogTotalDose)
  )

# Determine if we can use study covariate
within_study_variation <- any(dose_by_study_IP$sd_log_dose > 0.01, na.rm = TRUE)
use_study_covariate_IP <- within_study_variation

cat("\nDose values by StudyID_Split (Rapamycin only):\n")
rapa_dose_summary_IP <- dose_meta_IP %>%
  filter(StandardTreatment == "Rapamycin") %>%
  group_by(StudyID_Split) %>%
  summarise(
    n = n(),
    DailyDose = unique(DailyDose_mgkg),
    TotalDose = unique(TotalDose_mgkg),
    LogTotalDose = round(unique(LogTotalDose), 2)
  )
print(as.data.frame(rapa_dose_summary_IP))

# Align with genus_mat_adj
dose_meta_IP <- dose_meta_IP[rownames(dose_meta_IP) %in% colnames(genus_mat_adj), ]
cat("\nSample alignment check:", 
    all(rownames(dose_meta_IP) %in% colnames(genus_mat_adj)), "\n\n")

# Function to run linear model (IP-only)
run_dose_model_IP <- function(genus_name) {
  abund <- genus_mat_adj[genus_name, rownames(dose_meta_IP)]
  rapa_idx <- dose_meta_IP$StandardTreatment == "Rapamycin"
  if (sum(rapa_idx) == 0) return(NULL)
  
  rapa_abund <- abund[rapa_idx]
  rapa_meta_subset <- dose_meta_IP[rapa_idx, , drop = FALSE]
  
  keep <- !is.na(rapa_abund) & !is.na(rapa_meta_subset$LogTotalDose)
  if (sum(keep) < 5) return(NULL)
  
  rapa_abund_clean <- rapa_abund[keep]
  rapa_meta_clean <- rapa_meta_subset[keep, , drop = FALSE]
  
  if (var(rapa_abund_clean, na.rm = TRUE) < 1e-10) return(NULL)
  
  rapa_meta_clean$LogTotalDose_std <- as.numeric(scale(rapa_meta_clean$LogTotalDose))
  rapa_meta_clean$StudyID_Split <- droplevels(rapa_meta_clean$StudyID_Split)
  
  # Choose model based on collinearity check
  if (use_study_covariate_IP && length(unique(rapa_meta_clean$StudyID_Split)) > 1) {
    model <- lm(rapa_abund_clean ~ LogTotalDose_std + StudyID_Split, data = rapa_meta_clean)
    model_type <- "with_study"
    
    # Calculate VIF
    dose_on_study <- lm(LogTotalDose_std ~ StudyID_Split, data = rapa_meta_clean)
    vif_value <- 1 / (1 - summary(dose_on_study)$r.squared)
  } else {
    model <- lm(rapa_abund_clean ~ LogTotalDose_std, data = rapa_meta_clean)
    model_type <- "no_study"
    vif_value <- NA
  }
  
  summ <- summary(model)
  if (!"LogTotalDose_std" %in% rownames(summ$coefficients)) return(NULL)
  
  return(data.frame(
    feature = genus_name,
    coef = summ$coefficients["LogTotalDose_std", "Estimate"],
    stderr = summ$coefficients["LogTotalDose_std", "Std. Error"],
    pval = summ$coefficients["LogTotalDose_std", "Pr(>|t|)"],
    n_samples = sum(keep),
    model_type = model_type,
    vif = vif_value,
    stringsAsFactors = FALSE
  ))
}

# Run analysis for all genera
cat("=== RUNNING IP-ONLY DOSE-RESPONSE ANALYSIS ===\n")
valid_genera <- rownames(genus_mat_adj)
cat("Testing", length(valid_genera), "genera...\n")

dose_results_IP_list <- lapply(valid_genera, function(g) {
  tryCatch(run_dose_model_IP(g), error = function(e) NULL)
})

# Combine results
non_null_IP <- Filter(Negate(is.null), dose_results_IP_list)

if (length(non_null_IP) == 0) {
  stop("No valid models for IP-only. Check data alignment.")
}

cat("Successfully fit models for", length(non_null_IP), "genera\n\n")

dose_results_IP <- dplyr::bind_rows(non_null_IP)

# Apply FDR correction
dose_results_IP$qval <- p.adjust(dose_results_IP$pval, method = "BH")

# Show summary
cat("=== IP-ONLY RESULTS ===\n")
cat("P-value distribution:\n")
print(summary(dose_results_IP$pval))
cat("\nQ-value distribution:\n")
print(summary(dose_results_IP$qval))

# Get significant results
sig_dose_IP <- dose_results_IP %>% 
  filter(qval < 0.05) %>%
  filter(!grepl("unclassified|uncultured|unknown|Muribaculaceae", feature, ignore.case = TRUE)) %>%
  arrange(coef)

if (nrow(sig_dose_IP) > 0) {
  sig_dose_IP$clean_name <- gsub("s__|_", " ", sig_dose_IP$feature)
  sig_dose_IP$clean_name <- factor(sig_dose_IP$clean_name, levels = sig_dose_IP$clean_name)
  
  cat("\n=== SIGNIFICANT IP-ONLY DOSE-RESPONSE EFFECTS (FDR < 0.05) ===\n")
  cat("Found", nrow(sig_dose_IP), "significant genera:\n\n")
  print(sig_dose_IP)
  
} else {
  cat("\n=== NO SIGNIFICANT IP-ONLY EFFECTS AT FDR < 0.05 ===\n")
  
  # Try relaxed threshold
  sig_dose_IP <- dose_results_IP %>% 
    filter(qval < 0.1) %>%
    filter(!grepl("unclassified|uncultured|unknown|Muribaculaceae", feature, ignore.case = TRUE)) %>%
    arrange(coef)
  
  if (nrow(sig_dose_IP) > 0) {
    sig_dose_IP$clean_name <- gsub("s__|_", " ", sig_dose_IP$feature)
    sig_dose_IP$clean_name <- factor(sig_dose_IP$clean_name, levels = sig_dose_IP$clean_name)
    
    cat("\nIP-only genera with FDR < 0.1:\n")
    print(sig_dose_IP)
  }
}

# COMPARISON: ALL vs IP-ONLY

# Check overlap
if(nrow(sig_dose) > 0 & nrow(sig_dose_IP) > 0) {
  all_sig_names <- sig_dose$feature
  ip_sig_names <- sig_dose_IP$feature
  
  cat("Overlap in significant genera:\n")
  cat("  Only in ALL:", paste(setdiff(all_sig_names, ip_sig_names), collapse=", "), "\n")
  cat("  Only in IP-only:", paste(setdiff(ip_sig_names, all_sig_names), collapse=", "), "\n")
  cat("  In both:", paste(intersect(all_sig_names, ip_sig_names), collapse=", "), "\n")
}

#### VERIFICATION OF DOSE-RESPONSE ANALYSIS 

cat("=== COMPREHENSIVE VERIFICATION ===\n\n")

# 1. Check dose distribution across studies
cat("1. DOSE DISTRIBUTION:\n")
dose_summary <- dose_meta %>%
  filter(StandardTreatment == "Rapamycin") %>%
  group_by(StudyID_Split) %>%
  summarise(
    n = n(),
    TotalDose = unique(TotalDose_mgkg),
    LogTotalDose = round(unique(LogTotalDose), 3)
  ) %>%
  arrange(TotalDose)
print(dose_summary)

cat("\nDose variation check:\n")
cat("- Dose range:", range(dose_summary$TotalDose), "mg/kg\n")
cat("- Fold difference:", max(dose_summary$TotalDose)/min(dose_summary$TotalDose), "\n")
cat("- This should provide good power for dose-response analysis\n\n")

# 2. Manually verify the top significant genus
cat("2. MANUAL VERIFICATION OF TOP HIT:\n")
test_genus <- "Rikenellaceae_RC9_gut_group"

# Get the data
abund <- genus_mat_adj[test_genus, ]
rapa_idx <- dose_meta$StandardTreatment == "Rapamycin"
test_data <- data.frame(
  abundance = abund[rapa_idx],
  dose = dose_meta$TotalDose_mgkg[rapa_idx],
  log_dose = dose_meta$LogTotalDose[rapa_idx],
  study = dose_meta$StudyID_Split[rapa_idx]
)

# Check correlation
cat("Raw correlation (abundance vs log dose):", 
    cor(test_data$abundance, test_data$log_dose), "\n\n")


#### AUC tests####
suppressPackageStartupMessages({ library(dplyr); library(pROC) })

# treated only
treated <- dose_meta$StandardTreatment == "Rapamycin"
X  <- t(genus_mat_adj[, treated, drop=FALSE])   # samples x genera
md <- droplevels(dose_meta[treated, , drop=FALSE])
md$study <- droplevels(factor(md$StudyID_Split))

# make a binary label: top tertile vs bottom tertile of LogTotalDose
q <- quantile(md$LogTotalDose, probs = c(1/3, 2/3), na.rm = TRUE)
use <- md$LogTotalDose <= q[1] | md$LogTotalDose >= q[2]
y   <- as.integer(md$LogTotalDose[use] >= q[2])  # 1 = high dose

X2  <- X[use, , drop=FALSE]
res_auc <- lapply(colnames(X2), function(g){
  pred <- X2[, g]
  if (all(pred == pred[1]) || anyNA(pred)) return(NULL)
  auc  <- as.numeric(pROC::auc(y, pred))
  data.frame(feature=g, auc=auc)
})
res_auc <- dplyr::bind_rows(Filter(Negate(is.null), res_auc)) %>%
  arrange(desc(auc))

head(res_auc, 15)


# pick significant or top-|coef| genera as weights
w_tbl <- dose_results %>%
  filter(!is.na(coef)) %>%
  arrange(dplyr::desc(abs(coef))) %>%
  slice_head(n = 50)    # tweak K if you like
w <- setNames(w_tbl$coef, w_tbl$feature)

# keep only weighted features
feats <- intersect(colnames(X2), names(w))
if (length(feats) >= 2) {
  score <- as.numeric(X2[, feats, drop=FALSE] %*% w[feats])
  meta_auc <- as.numeric(pROC::auc(y, score))
  cat(sprintf("Meta-score AUC (high vs low dose): %.3f\n", meta_auc))
}

partial_r2_for_genus <- function(genus_name){
  abund <- genus_mat_adj[genus_name, ]
  idx   <- dose_meta$StandardTreatment == "Rapamycin"
  yv    <- abund[idx]
  df    <- droplevels(dose_meta[idx, , drop=FALSE])
  keep  <- is.finite(yv) & is.finite(df$LogTotalDose)
  yv    <- yv[keep]; df <- df[keep, , drop=FALSE]
  if (length(yv) < 5) return(NULL)
  
  df$LogTotalDose_std <- as.numeric(scale(df$LogTotalDose))
  full <- lm(yv ~ LogTotalDose_std + StudyID_Split, data=df)
  red  <- lm(yv ~ StudyID_Split, data=df)
  SSfull <- sum(resid(full)^2); SSred <- sum(resid(red)^2)
  R2partial <- (SSred - SSfull) / SSred
  data.frame(feature=genus_name, partial_R2=R2partial)
}

r2_tbl <- dplyr::bind_rows(lapply(rownames(genus_mat_adj), partial_r2_for_genus)) %>%
  arrange(dplyr::desc(partial_R2))
head(r2_tbl, 15)



###Visualization####
# Create the dose effect plot
if(nrow(sig_dose) > 0) {
  dose_effect_plot <- ggplot(sig_dose, aes(x = clean_name, y = coef)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_errorbar(aes(ymin = coef - 1.96 * stderr, ymax = coef + 1.96 * stderr), width = 0.2) +
    coord_flip() +
    labs(
      title = "Dose Effect of Rapamycin on Microbial Genera",
      subtitle = "Linear model across studies (FDR < 0.05)",
      x = "",
      y = "Effect Size (Coefficient)"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)
  
  # Display and save the plot
  print(dose_effect_plot)
}

##Forests plot
# Make a barplot of slopes for significant genera
plot.barplot_dose_coefficients <- sig_dose %>%
  ggplot(aes(x = reorder(clean_name, coef), y = coef)) +
  geom_col(aes(fill = coef > 0)) +
  geom_errorbar(aes(ymin = coef - 1.96*stderr, ymax = coef + 1.96*stderr), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(title = "Dose-Response Coefficients by Genus", x = "Genus", y = "Coefficient") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"))

# Display the plot
print(plot.barplot_dose_coefficients)


# Even cleaner dot plot version
plot.dotplot_dose_coefficient <- sig_dose %>%
  mutate(
    clean_name = str_remove_all(clean_name, "\\[|\\]"),  # Remove brackets
    Significance = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01 ~ "**",
      qval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = coef, y = reorder(clean_name,coef))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = coef - 1.96*stderr, xmax = coef + 1.96*stderr), 
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = coef > 0), size = 4, alpha = 0.8) +
  geom_text(aes(label = Significance), hjust = -0.8, vjust = -0.1, size = 4, fontface = "bold") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#2166AC")) +
  theme(legend.position = "none") +
  labs(
    title = "Significant Dose-Response Effects",
    subtitle = "Linear model coefficients for genera with FDR < 0.05",
    x = bquote("Coefficient (proportion per SD " * log[10] * " dose)"),
    y = "Genus"
  )


plot.dotplot_dose_coefficient

plot.dotplot_dose_coefficients_percent <- sig_dose %>%
  mutate(
    # Calculate percentage change and confidence intervals
    # Note: Since your coefficients are already in transformed space, 
    # we'll convert them to approximate percentage changes
    pct_change = 100 * coef,  # Simple conversion for transformed data
    pct_lower = 100 * (coef - 1.96*stderr),
    pct_upper = 100 * (coef + 1.96*stderr),
    
    # Create significance stars
    Significance = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01 ~ "**",
      qval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = pct_change, y = reorder(clean_name, pct_change))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = pct_lower, xmax = pct_upper), 
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = pct_change > 0), size = 4, alpha = 0.8) +
  geom_text(aes(label = Significance), hjust = -0.5, size = 4, fontface = "bold") +
  scale_x_continuous(limits = c(-0.6, 0.2)) +  # Add this line
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#2166AC")) +
  theme(legend.position = "none") +
  labs(
    title = "Significant Dose-Response Effects",
    subtitle = "Linear model effects as % change per standardized log dose (FDR < 0.05)",
    x = "% Change per standardized log dose",
    y = "Genus"
  )

# Display the plot
print(plot.dotplot_dose_coefficients_percent)

####IP only ###

# Dot plot for IP-only dose-response results
# Dot plot for IP-only dose-response results
plot.dotplot_dose_coefficient_IP <- sig_dose_IP %>%
  mutate(
    clean_name = str_remove_all(clean_name, "\\[|\\]"),  # Remove brackets
    Significance = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01 ~ "**",
      qval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = coef, y = reorder(clean_name, coef))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = coef - 1.96*stderr, xmax = coef + 1.96*stderr), 
                 height = 0.2, size = 0.8) +
  geom_point(aes(color = coef > 0), size = 4, alpha = 0.8) +
  geom_text(aes(label = Significance),hjust = -0.7, vjust = -0.01, size = 4, fontface = "bold") +
  scale_x_continuous(limits = c(-0.003, 0.001)) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#2166AC")) +
  theme(legend.position = "none") +
  labs(
    title = "Dose-Associated Effects (IP-Route Only)",
    subtitle = "Linear model coefficients for genera with FDR < 0.05",
    x = bquote("Coefficient (proportion per SD " * log[10] * " dose)"),  # CHANGED to "proportion"
    y = "Genus"
  )

plot.dotplot_dose_coefficient_IP

#### plotted together 
#### plotted together 
combined <- (plot.dotplot_dose_coefficient + 
               labs(title = "All studies", subtitle = NULL, x = NULL)) +  # Add panel title
  (plot.dotplot_dose_coefficient_IP + 
     labs(title = "IP only", subtitle = NULL, x = NULL, y = NULL)) +  # Add panel title
  plot_annotation(
    tag_levels = 'a',
    caption = bquote("Coefficient (proportion per SD " * log[10] * " dose)")
  ) +
  plot_layout(widths = c(1, 1.3)) &
  theme(plot.tag = element_text(size = 14, face = 'bold'),
        plot.caption = element_text(hjust = 0.5, size = 12))

# Display
print(combined)

# Save with high resolution
ggsave("figure5_dose_response.png", combined,
       width = 11, height = 4, dpi = 300)
### Clipped

lim <- c(-0.004, 0.002)

sig_dose_IP %>%
  dplyr::mutate(
    ci_lo = coef - 1.96*stderr,
    ci_hi = coef + 1.96*stderr
  ) %>%
  dplyr::filter(
    !is.finite(ci_lo) | !is.finite(ci_hi) |
      ci_lo < lim[1] | ci_hi > lim[2] |
      coef  < lim[1] | coef  > lim[2]
  ) %>%
  dplyr::select(feature, clean_name, coef, stderr, ci_lo, ci_hi) %>%
  dplyr::arrange(coef)

sum(dose_meta_IP$StandardTreatment == "Rapamycin")


#### Effect size ###
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(pROC) })

## 0) pick weights (prefer your significant dose-response taxa)
if (exists("sig_dose") && nrow(sig_dose) > 0) {
  w_tbl <- sig_dose %>% dplyr::select(feature, coef)
} else if (exists("dose_results") && nrow(dose_results) > 0) {
  w_tbl <- dose_results %>% 
    filter(!is.na(coef)) %>% 
    arrange(dplyr::desc(abs(coef))) %>% 
    slice_head(n = 50) %>%          # keep top-|coef| (tweak K if you want)
    dplyr::select(feature, coef)
} else {
  stop("Need sig_dose or dose_results to build a dose-effect score.")
}
w <- setNames(w_tbl$coef, w_tbl$feature)

## 1) build a sample-level “dose-effect” score in Rapamycin samples
rapa_idx <- dose_meta$StandardTreatment == "Rapamycin"
feats    <- intersect(rownames(genus_mat_adj), names(w))
stopifnot(length(feats) >= 2)

Xr  <- t(genus_mat_adj[feats, rapa_idx, drop = FALSE])   # samples × features
mdr <- droplevels(dose_meta[rapa_idx, , drop = FALSE])

score <- as.numeric(Xr %*% w[feats])     # linear meta-score
mdr$DoseScore <- as.numeric(scale(score))  # z-score for readability
mdr$LogTotalDose_std <- as.numeric(scale(mdr$LogTotalDose))

## 2) global association: linear (adj. for study) + Spearman
fit <- lm(DoseScore ~ LogTotalDose_std + StudyID_Split, data = mdr)
fit_red <- lm(DoseScore ~ StudyID_Split, data = mdr)

SS_full <- sum(resid(fit)^2); SS_red <- sum(resid(fit_red)^2)
partial_R2 <- (SS_red - SS_full) / SS_red

spe <- suppressWarnings(cor.test(mdr$DoseScore, mdr$LogTotalDose, method = "spearman"))

cat("\n=== GLOBAL DOSE → DOSE-SCORE ASSOCIATION ===\n")
cat(sprintf("Linear (adj. StudyID_Split): slope=%.3f, p=%.3g\n",
            coef(fit)["LogTotalDose_std"], summary(fit)$coefficients["LogTotalDose_std","Pr(>|t|)"]))
cat(sprintf("Partial R^2 for dose (over Study): %.3f\n", partial_R2))
cat(sprintf("Spearman rho: %.3f (p=%.3g)\n", spe$estimate, spe$p.value))

## 3) per-study slopes (sanity: does trend hold inside each split study?)
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

## 4) visualization
p_dose_score <- ggplot(mdr, aes(x = LogTotalDose, y = DoseScore, color = StudyID_Split)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "Dose–Effect Score vs Log10 Total Dose (Rapamycin only)",
       subtitle = "Linear fit shown in black; points colored by StudyID_Split",
       x = "Log10 Total Dose (mg/kg)", y = "Dose–Effect Score (z)") +
  theme_bw()
print(p_dose_score)

## 5) optional: AUC for high- vs low-dose using the same score (extremes only)
q <- quantile(mdr$LogTotalDose, c(1/3, 2/3), na.rm = TRUE)
use <- mdr$LogTotalDose <= q[1] | mdr$LogTotalDose >= q[2]
if (sum(use) >= 6) {
  y_bin <- as.integer(mdr$LogTotalDose[use] >= q[2])    # 1=high dose
  auc_meta <- as.numeric(pROC::auc(y_bin, mdr$DoseScore[use]))
  cat(sprintf("\nAUC (high vs low dose, extremes only) using DoseScore: %.3f\n", auc_meta))
}


#######

plot.combined_plot_dose <- combined_plot_dose
plot.p_dose_before <- p_dose_before
plot.p_dose_after <- p_dose_after
plot.treatment_effect_plot <- treatment_effect_plot
plot.rapamycin_effect_plot


#### plots (in progress) #####

### 2025-03-13 Exported as PDF
pdf("MetaRapamycinAnalysisMMUPHin_20250506_V01.pdf", width = 14, height = 10)
plot.combined_plot_dose <- combined_plot_dose
plot.p_dose_before <- p_dose_before
plot.p_dose_after <- p_dose_after
plot.treatment_effect_plot <- treatment_effect_plot
plot.rapamycin_effect_plot

dev.off()

##### tables (in progress) #####
table.results_table 
table.results_Permanova
table.results_dose_Permanova
table.permidsp_table 

##### image #####
# Save the current R workspace
save.image("MetaAnalysisMMUPHIN_20250605_V01.RData")

# Load the saved workspace for later use
load("MetaAnalysisMMUPHIN_20250605_V01.RData")

# Print the session information
sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] reshape2_1.4.4              metafor_4.8-0               numDeriv_2016.8-1.1        
# [4] metadat_1.4-0               Matrix_1.7-1                metamicrobiomeR_1.2        
# [7] gamlss_5.4-22               nlme_3.1-166                gamlss.dist_6.1-1          
# [10] gamlss.data_6.0-6           MMUPHin_1.20.0              knitr_1.50                 
# [13] lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1              
# [16] purrr_1.0.4                 readr_2.1.5                 tibble_3.2.1               
# [19] tidyverse_2.0.0             splitstackshape_1.4.8       viridis_0.6.5              
# [22] viridisLite_0.4.2           microbiome_1.28.0           phyloseq_1.50.0            
# [25] patchwork_1.3.0             ggpubr_0.6.0                ggplot2_3.5.2              
# [28] DESeq2_1.46.0               SummarizedExperiment_1.36.0 Biobase_2.66.0             
# [31] MatrixGenerics_1.18.1       matrixStats_1.5.0           GenomicRanges_1.58.0       
# [34] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0           
# [37] BiocGenerics_0.52.0         ape_5.8-1                   vegan_2.6-10               
# [40] lattice_0.22-6              permute_0.9-7               tidyr_1.3.1                
# [43] dplyr_1.1.4                




