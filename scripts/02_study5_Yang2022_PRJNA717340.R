### 2024-12-06 Rapamycin reanalysis 2024
#for Dataset 5.PRJNA717340
### 2024-12-06 to reload
# run chunk 1
# load image
# > load("Rapamycin_reanalysis_5.PRJNA717340_V02_20250129_latest.RData") # have to change to new rdata file
# free memory
# > gc()


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
  library(stringr)# for str_split_fixed()
  library(reshape2) # for 'melt()'; unifrac plot custom function 
  library(tidyr) # data processing
  library(phyloseq) # BIOM and sample data processing
  library(vegan) # for 'adonis()'
  library(ape) # for read.tree()
  library(DESeq2) # required for 'DESeq', 'results'
  library(ggplot2) # plots
  library(ggpubr) # improved ggplots, error bars, wilcoxon/kruskal etc
  library(patchwork) # arranging ggplots
  #library("dendextend") # color branches and labels of dendrograms
  library(openxlsx) # for createWorkbook(), addWorksheet(), writeData(), saveWorkbook()
  library(microbiome) # core microbiome calculations
  library(viridis) # for 'scale_fill_viridis()'
  library(splitstackshape) # for core microbiome matrix
  library(tidyverse) # for core microbiome matrix
  library(jsonlite) # for getDisbiome()
  library(rstatix) # for wilcox_test() etc in SCFA calculations
  library(data.table) # # rowid() in clin parameters
  library(randomForest)  # RF classifiers
  library(knitr) # GLM tables

suppressPackageStartupMessages({
  library(dplyr) # required for '%>%' in 'RelativeAbundancePlots'
  library(stringr)# for str_split_fixed()
  library(reshape2) # for 'melt()'; unifrac plot custom function 
  library(tidyr) # data processing
  library(phyloseq) # BIOM and sample data processing
  library(vegan) # for 'adonis()'
  library(ape) # for read.tree()
  library(DESeq2) # required for 'DESeq', 'results'
  library(ggplot2) # plots
  library(ggpubr) # improved ggplots, error bars, wilcoxon/kruskal etc
  library(patchwork) # arranging ggplots
  #library("dendextend") # color branches and labels of dendrograms
  library(openxlsx) # for createWorkbook(), addWorksheet(), writeData(), saveWorkbook()
  library(microbiome) # core microbiome calculations
  library(viridis) # for 'scale_fill_viridis()'
  library(splitstackshape) # for core microbiome matrix
  library(tidyverse) # for core microbiome matrix
  library(jsonlite) # for getDisbiome()
  library(rstatix) # for wilcox_test() etc in SCFA calculations
  library(data.table) # # rowid() in clin parameters
  library(randomForest)  # RF classifiers
  library(knitr) # GLM tables
}) 


##### Constants and Functions #####

### Constants
# FDR cutoff for differential abundances
fdr.threshold <- 0.01

# working dir
setwd(here::here("data/study5_PRJNA717340"))

# load functions
source(here::here("scripts/MicrobiomeFunctions.R"))

#### Load Data #### 
### input file paths

inputfile.metadata <- "filereport_read_run_PRJNA717340.tsv"
inputfile.otu <- "feature-table.tsv"
inputfile.taxonomy <- "taxonomy.tsv"
inputfile.tree <- "tree.nwk"


##### 16S phyloseq object (ok) #####
file_path <- "filereport_read_run_PRJNA717340.tsv"
metadata <- read.delim(file_path, header = TRUE, sep = "\t")

# Example: drop columns where *all* entries are NA
metadata <- metadata[, colSums(is.na(metadata)) < nrow(metadata)]

write.table(metadata, file = "shifted_sample_metadata.tsv",sep = "\t",row.names = FALSE,quote = FALSE)

### mapping
SampleData <- import_qiime_sample_data("shifted_sample_metadata.tsv")


### 2) FEATURE (OTU) TABLE

otumat <- read.csv(inputfile.otu, sep = "\t", header = TRUE, stringsAsFactors=FALSE)

# Rename the "X.OTU.ID" column to "OTU.ID"
colnames(otumat)[colnames(otumat) == "X.OTU.ID"] <- "OTU.ID"

# Move the first column (OTU.ID) to row names
rownames(otumat) <- otumat$OTU.ID
otumat$OTU.ID  <- NULL

# Now convert to a phyloseq 'otu_table', specifying rows = taxa (OTUs)
OTU <- otu_table(as.matrix(otumat), taxa_are_rows = TRUE)

# Check that the column names of otumat match your sample IDs:
colnames(otumat)


### 3) TAXONOMY TABLE

taxonomy <- read.table(
  inputfile.taxonomy,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = ""
)

# The taxonomy file presumably has these columns: OTU_ID, Taxon, Confidence
# We only need the OTU_ID and Taxon columns to build the tax_table.
# Let's move OTU_ID into row names so it matches the row names in OTU.

colnames(taxonomy)[colnames(taxonomy) == "Feature.ID"] <- "OTU.ID" # Rename the "Feature.ID" column to "OTU.ID"

rownames(taxonomy) <- taxonomy$OTU.ID
taxonomy$OTU.ID <- NULL

# Now split the "Taxon" column by semicolons
tax_levels_list <- strsplit(taxonomy$Taxon, ";", fixed = TRUE)

# For each row, pad or trim to 7 ranks:
tax_levels_fixed <- lapply(tax_levels_list, function(x) {
  # If more than 7, keep only first 7
  if (length(x) > 7) {
    x <- x[1:7]
  }
  
  # If fewer than 7, fill the rest with "Unclassified"
  needed <- 7 - length(x)
  if (needed > 0) {
    x <- c(x, rep("Unclassified", needed))
  }
  
  return(x)
})

# Turn this into a data frame
taxmat <- as.data.frame(do.call(rbind, tax_levels_fixed), stringsAsFactors = FALSE)

# Give columns sensible names
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Make sure row names line up with your OTUs
# (they already do, because we set rownames(taxonomy) above)
rownames(taxmat) <- rownames(taxonomy)

# Convert to phyloseq tax_table
TAX <- tax_table(as.matrix(taxmat))

### tree data
TreeData <- read.tree(inputfile.tree)
# transforming all multichotomies into a series of dichotomies with one (or several) branch(es)
# of length zero, to avoid UniFrac warning messages about unrooted trees
# see https://github.com/joey711/phyloseq/issues/936
TreeData <- ape::multi2di(TreeData)

### phyloseq object
qiime2.ps <- phyloseq(OTU, TAX, TreeData)
qiime2.ps <- merge_phyloseq(qiime2.ps, SampleData)

# sample data for visualization
tempsampledata2 <- data.frame(sample_data(qiime2.ps))

qiime2.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1004 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 40 sample variables ]
# tax_table()   Taxonomy Table:    [ 1004 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1004 tips and 1003 internal nodes ]

### clean phyloseq object
temp.ps <- tidy_phyloseq(qiime2.ps)
# further clean mitochondria and chloroplast
# latest SILVA classifies them at the genus level
temp.ps <-  subset_taxa(temp.ps, Genus != "Chloroplast")
temp.ps <-  subset_taxa(temp.ps, Genus != "Mitochondria")
temp.ps <-  subset_taxa(temp.ps, Phylum != "")
# also removing the "aqua" negative control
#temp.ps <- subset_samples(temp.ps, X.SampleID != "negativ")

# 1) Save original sample names in a vector
original_samples <- sample_names(temp.ps)

# 2) Prune low-coverage samples
temp.ps <- prune_samples(sample_sums(temp.ps) >= 2000, temp.ps)

# 3) Check which samples were removed
removed_samples <- setdiff(original_samples, sample_names(temp.ps))

# 4) Print the result
if (length(removed_samples) == 0) {
  cat("No samples removed.\n")
} else {
  cat("Removed the following samples:\n")
  print(removed_samples)
}
# Removed the following samples:
# No samples removed.
allsamples.ps <- temp.ps #for correlation analysis

temp.ps <-subset_samples(temp.ps, Group %in% c("WT", "WT-RAPA"))
sample_data(temp.ps)$Group %>% unique()

# Add this at the beginning of your analysis, right after creating temp.ps
# This changes the group names in the actual data
sample_data(temp.ps)$Group <- plyr::mapvalues(
  sample_data(temp.ps)$Group,
  from = c("WT", "WT-RAPA"),
  to = c("Control", "Rapamycin (IP)")
)

# Update grouporder
grouporder <- c("Control", "Rapamycin (IP)")

# THEN generate comparisons using the renamed groups
groups <- levels(factor(sample_data(temp.ps)$Group))
my_comparisons <- combn(groups, 2, simplify = FALSE)

temp.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1004 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 40 sample variables ]
# tax_table()   Taxonomy Table:    [ 1004 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1004 tips and 1003 internal nodes ]

# ## outlier  # In phyloseq:
# alpha_df <- estimate_richness(temp.ps, measures = "Shannon")
# alpha_df$Sample <- rownames(alpha_df)
# 
# alpha_df %>%
#   arrange(Shannon) %>%
#   head(10)  # see the "lowest Shannon" samples
# sample_comp <- psmelt(temp.ps) %>%
#   filter(Sample == "ERR1428197") %>%
#   group_by(Taxonomy = Genus) %>%
#   summarize(Reads = sum(Abundance)) %>%
#   arrange(desc(Reads))
#head(sample_comp, 20)

summary_results <- microbiome::summarize_phyloseq(temp.ps)
print(summary_results)
# Compositional = NO2
# 1] Min. number of reads = 170862] Max. number of reads = 283553] Total number of reads = 2185864] Average number of reads = 21858.65] Median number of reads = 19579.57] Sparsity = 0.7259960159362556] Any OTU sum to 1 or less? YES8] Number of singletons = 2749] Percent of OTUs that are singletons 
# (i.e. exactly one read detected across all samples)0.19920318725099610] Number of sample variables are: 40Sample.idsecondary_study_accessionsample_accessionsecondary_sample_accessionexperiment_accessionstudy_accessionsubmission_accessiontax_idscientific_nameinstrument_platforminstrument_modellibrary_namelibrary_layoutlibrary_strategylibrary_sourcelibrary_selectionread_countbase_countcenter_namefirst_publiclast_updatedexperiment_titlestudy_titlestudy_aliasexperiment_aliasrun_aliasfastq_bytesfastq_md5fastq_ftpfastq_asperafastq_galaxysra_bytessra_md5sra_ftpsra_asperasra_galaxysample_aliasGroupsample_titlefirst_created2
# > print(summary_results)
# [[1]]
# [1] "1] Min. number of reads = 17086"
# 
# [[2]]
# [1] "2] Max. number of reads = 28355"
# 
# [[3]]
# [1] "3] Total number of reads = 218586"
# 
# [[4]]
# [1] "4] Average number of reads = 21858.6"
# 
# [[5]]
# [1] "5] Median number of reads = 19579.5"
# 
# [[6]]
# [1] "7] Sparsity = 0.725996015936255"
# 
# [[7]]
# [1] "6] Any OTU sum to 1 or less? YES"
# 
# [[8]]
# [1] "8] Number of singletons = 274"
# 
# [[9]]
# [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)0.199203187250996"
# 
# [[10]]
# [1] "10] Number of sample variables are: 40"
# 
# [[11]]
# [1] "Sample.id"                  "secondary_study_accession"  "sample_accession"          
# [4] "secondary_sample_accession" "experiment_accession"       "study_accession"           
# [7] "submission_accession"       "tax_id"                     "scientific_name"           
# [10] "instrument_platform"        "instrument_model"           "library_name"              
# [13] "library_layout"             "library_strategy"           "library_source"            
# [16] "library_selection"          "read_count"                 "base_count"                
# [19] "center_name"                "first_public"               "last_updated"              
# [22] "experiment_title"           "study_title"                "study_alias"               
# [25] "experiment_alias"           "run_alias"                  "fastq_bytes"               
# [28] "fastq_md5"                  "fastq_ftp"                  "fastq_aspera"              
# [31] "fastq_galaxy"               "sra_bytes"                  "sra_md5"                   
# [34] "sra_ftp"                    "sra_aspera"                 "sra_galaxy"                
# [37] "sample_alias"               "Group"                      "sample_title"              
# [40] "first_created"             


####Save Phyloseq object externally####
## What phyloseq objects are in memory?
ls()                     # lists all objects
class(temp.ps)           # should print "phyloseq" if it exists


# 1.  Make sure the clean object is called `clean_ps`
clean_ps <- temp.ps            # or qiime2.ps / ps_final …

# 2.  Tag samples with the study ID  ← change to your study label
sample_data(clean_ps)$StudyID <- "PRJNA717340"


saveRDS(
  clean_ps,
  file = here::here("data/processed/ps_PRJNA717340_clean.rds"),
  compress = "xz"
)

message("✓ Saved cleaned phyloseq for PRJNA717340")


---------------------------------------------
  
  #----------- Phylum-level Extraction -----------#
  # 1) Extract the OTU table
  otumat <- as.data.frame(otu_table(temp.ps))  # convert to data.frame

# 2) If taxa_are_rows == TRUE, rownames are taxa (ASVs/OTUs), colnames are samples.
#    Make that explicit:
otumat$Taxon <- rownames(otumat)

# 3) Pivot to a "long" format: Sample, Taxon, Count
df_otu <- otumat %>%
  pivot_longer(
    cols = -Taxon,
    names_to = "Sample",
    values_to = "Count"
  )

# 4) Convert tax_table to a data frame
df_tax <- as.data.frame(tax_table(temp.ps))

# 5) Keep row names as "Taxon" so we can join
df_tax$Taxon <- rownames(df_tax)

# 6) Join with df_otu by Taxon
df_otu_tax <- left_join(df_otu, df_tax, by = "Taxon")

# 7) Group by Sample and Phylum, then sum all counts
df_phylum <- df_otu_tax %>%
  group_by(Sample, Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop")

# 8) Calculate relative abundance
df_phylum <- df_phylum %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Count / sum(Count)) %>%
  ungroup()

# 9) Write phylum-level table to file
write.table(df_phylum, 
            file = here::here("results/study5_phyloseq_phylum_table.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#----------- Genus-level Extraction -----------#
# Follow the same approach as above but group by Genus

# We can reuse df_otu_tax from above since it already has the joined taxonomic information

# 10) Group by Sample and Genus, then sum all counts
df_genus <- df_otu_tax %>%
  group_by(Sample, Genus) %>%
  summarise(Count = sum(Count), .groups = "drop")

# 11) Calculate relative abundance
df_genus <- df_genus %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Count / sum(Count)) %>%
  ungroup()

# 12) Write genus-level table to file
write.table(df_genus, 
            file = here::here("results/study5_phyloseq_genus_table.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

-----------------------------------------------------

  # cleanup
  rm(OTU, TAX, taxonomy, taxmat)

##### Phylum/Genera compositions (ok) #####

### levels for ordered plots


grouporder <- c("Control", "Rapamycin (IP)")

### Phylum by comparison group
pseq <- microbiome::transform(temp.ps, "compositional")
pseq <- aggregate_rare(pseq, level = "Phylum", detection = 1/1000, prevalence = 50/100)
mysubtitle <- "Detection 0.1%, Prevalence 50%, Kruskal-Wallis p-value "
tempplot <- plot_composition(pseq, average_by = "Group") +
  scale_fill_brewer("Phyla", palette = "Paired") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) +
  labs(x = "Group ", 
       y = "Relative abundance",
       title = "Phylum composition") 
tempplot$data$Sample <- as.character(tempplot$data$Sample)
tempplot$data$Sample <- factor(x = tempplot$data$Sample, levels = grouporder)
ps.kruskal <- kruskal.test(Abundance ~ Sample, data=tempplot$data)
tempplot <- tempplot + labs(subtitle = paste(mysubtitle, round(ps.kruskal$p.value, digits = 4)))
plot.PhylumComp <- tempplot

print(plot.PhylumComp)

bubble_table <- psmelt(pseq) %>%
  arrange(Phylum)
bubbleTable <- bubble_table[, c("Sample","Phylum", "Abundance", "Group")] # columns: Sample, Phylum, Abundance

# Remove rows with Abundance == 0 or NA
bubbleTable <- bubbleTable %>%
  filter(!is.na(Abundance) & Abundance > 0)
# Reorder factor levels here:
bubbleTable$Group <- factor(
  bubbleTable$Group,
  levels = c("Control", "Rapamycin (IP)", "RhoB", "RhoB-RAPA")
)
myBubblePlot <- ggplot(bubbleTable, aes(x = Group, y = Phylum)) + 
  geom_point(aes(color = Phylum, size = Abundance), alpha = 1) +
  xlab("Group") +
  scale_size_continuous(limits = c(0.005, 0.6), range = c(1,25)) + 
  theme_bw() + # remove background color
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  # Adjust x-axis label alignment
  labs(
    title = "Phylum Bubble Plot",
    y = "Phylum",
    x = "Group"
  )

myBubblePlot 
plot.PhylumBubble <- myBubblePlot

# to retrieve plot data
#tempdata <- tempplot$data

### Genus by comparison group
pseq <- microbiome::transform(temp.ps, "compositional")
# fix error by assigning empty to "unclassified"
# src https://github.com/microbiome/microbiome/issues/113
# Step 1: Convert to character and replace empty values
tax_table(pseq)[, "Genus"] <- as.character(tax_table(pseq)[, "Genus"])
tax_table(pseq)[tax_table(pseq)[, "Genus"] == "", "Genus"] <- " Unclassified"
tax_table(pseq)[tax_table(pseq)[, "Genus"] == "Other", "Genus"] <- " Unclassified"
tax_table(pseq)[grepl("uncultured|unclassified", tax_table(pseq)[, "Genus"], ignore.case=TRUE), "Genus"] <- " Unclassified"

pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/1000, prevalence = 50/100)
mysubtitle <- "Detection 1.25%, Prevalence 50%, Kruskal-Wallis p-value "
tempplot <- plot_composition(pseq, average_by = "Group") +
  scale_fill_viridis_d("Genera") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) +
  labs(x = "Group", 
       y = "Relative abundance",
       title = "Genus composition") 

# ADD THESE TWO LINES HERE:
tempplot$data$Tax[tempplot$data$Tax == "Other"] <- " Unclassified"
tempplot$data <- aggregate(Abundance ~ Sample + Tax, data = tempplot$data, FUN = sum)

tempplot$data$Sample <- as.character(tempplot$data$Sample)
tempplot$data$Sample <- factor(x = tempplot$data$Sample, levels = grouporder)
ps.kruskal <- kruskal.test(Abundance ~ Sample, data=tempplot$data)
tempplot <- tempplot + labs(subtitle = paste(mysubtitle, round(ps.kruskal$p.value, digits = 4)))
plot.GenusComp <- tempplot

print(plot.GenusComp)

#phyloseq object to relative abundances (Genus level)

bubble_table <- psmelt(pseq) %>%
  arrange(Genus)

bubble_table$Genus <- ifelse(
  # Use case-insensitive regex:
  grepl("(?i)(uncultured|unclassified)", bubble_table$Genus),
  "Uncharacterized",
  bubble_table$Genus
)

bubbleTable <- bubble_table[, c("Sample","Genus", "Abundance", "Group")] # columns: Sample, Phylum, Abundance
# Remove rows with Abundance == 0 or NA
bubbleTable <- bubbleTable %>%
  filter(!is.na(Abundance) & Abundance > 0)

bubbleTable <- bubbleTable %>%
  mutate(
    Genus = str_trim(Genus),                          # trim hidden spaces
    Genus = case_when(                                # collapse variants
      Genus %in% c("uncultured", "Uncultured","Uncharacterized",
                   "Other", "Unclassified") ~ "Unclassified",
      TRUE                                  ~ Genus
    ),
    Genus = factor(Genus, 
                     levels = c(setdiff(sort(unique(Genus)), "Unclassified"), 
                                "Unclassified"))
  ) %>%
  droplevels()                                         # clean out unused levels

# Reorder factor levels here:
bubbleTable$Group<- factor(
  bubbleTable$Group,
  levels = c("Control", "Rapamycin (IP)")
)

myBubblePlot <- ggplot(bubbleTable, aes(x = Group, y = Genus)) + 
  geom_point(aes(color = Genus, size = Abundance), alpha = 1) +
  xlab("Group") +
  scale_size_continuous(limits = c(0.001, 0.6), range = c(1,25)) + 
  theme_bw() + # remove background color
  theme(axis.text.x = element_text(angle= 45, hjust= 1)) +
  labs(x = "Group", 
       y = "Relative abundance",
       title = "Genus Bubble Plot") 

myBubblePlot 
plot.GenusBubble <- myBubblePlot

# to retrieve plot data
# tempdata <- tempplot$data


### cleanup
rm(pseq, tempplot, myBubblePlot, bubbleTable, bubble_table)


### Abundance tables for publication
tempdata <- plot.PhylumComp$data
tempdata$xlabel <- NULL
table.Phylum.Abundance<- reshape2::dcast(tempdata, Tax ~ Sample, value.var="Abundance", fun.aggregate = sum)
table.Phylum.Abundance$Tax <- gsub(" ", "", table.Phylum.Abundance$Tax)

tempdata <- plot.GenusComp$data
tempdata$xlabel <- NULL
table.Genus.Abundance <- reshape2::dcast(tempdata, Tax ~ Sample, value.var="Abundance", fun.aggregate = sum)
table.Genus.Abundance$Tax <- gsub(" ", "", table.Genus.Abundance$Tax)

rm(tempdata)

###combined 
# Phylum plot with relabeled x-axis
# Phylum plot with correct relabeling
plot.PhylumComp <- plot.PhylumComp + 
  scale_x_discrete(labels = c("Control" = "Control-IP", 
                              "Rapamycin (IP)" = "Rapamycin-IP"),
                   expand = c(0, 0)) +  # Remove gap
  scale_y_continuous(expand = c(0, 0)) +  # Remove gap on y-axis too) +
  labs(tag = "a", title = NULL, subtitle = NULL,x = NULL) +
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"))+
theme(plot.tag.position = c(-0.01, 1),
      plot.tag = element_text(size = 14, face = "plain"),
      panel.background = element_blank(),
      panel.grid = element_blank())



# Genus plot with correct relabeling  
plot.GenusComp <- plot.GenusComp + 
  scale_x_discrete(labels = c("Control" = "Control-IP", 
                              "Rapamycin (IP)" = "Rapamycin-IP"),
                   expand = c(0, 0)) +  # Remove gap
  scale_y_continuous(expand = c(0, 0)) +  # Remove gap on y-axis too)) +
  labs(tag = "b", title = NULL, subtitle = NULL,x = NULL) +
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"))+
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        panel.background = element_blank(),
        panel.grid = element_blank())


# Combine side by side
library(patchwork)
plot.Composition.combined <- plot.PhylumComp | plot.GenusComp
plot.Composition.combined


##### FB ratios (ok) #####

pseq <- microbiome::transform(temp.ps, "compositional")
pseq <- aggregate_rare(pseq, level = "Phylum", detection = 1/1000, prevalence = 50/100)
tempplot <- plot_composition(pseq) 
tempdata <- tempplot$data
tempdata <- tempdata[, c(-4)]

temp.list <- paste(c("Firmicutes", "Bacteroidota"), collapse = '|') 
tempdata <- tempdata %>% filter(str_detect(Tax, temp.list))

tempdata <- tempdata %>% 
  #select(-Group) %>% 
  group_by(Sample) %>%
  spread(Tax, Abundance) 

colnames(tempdata) <- gsub(" ", "", x = colnames(tempdata))
tempdata$FBratio <- tempdata$Firmicutes / tempdata$Bacteroidota
# Before merging, let's see what we have
print("Structure of tempdata before merge:")
print(str(tempdata))

tempdata <- merge(tempdata, tempsampledata2, by.x=1, by.y=1)
tempdata$logFB <- log(tempdata$FBratio) # using logs because otherwise range way too big in plot 

# After the merge but before plotting, rename the groups
tempdata$Group <- plyr::mapvalues(
  tempdata$Group,
  from = c("WT", "WT-RAPA"),
  to = c("Control", "Rapamycin (IP)")
)

# Then use your existing my_comparisons list in the plot

plot.FBratio <- ggboxplot(tempdata, x = "Group", y = "logFB", color = "Group", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons) + 
  #ylim(-1, 10) + # Add pairwise comparisons p-value
  ggtitle("Firmicutes - Bacteroidetes ratios") +
  xlab("Group") +
  ylab("log(FB ratio)")+
  # Remove redundant color legend
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

plot.FBratio
table.FBratio <- tempdata

### mean values
mean(table.FBratio[table.FBratio$Group == 'Control', 'FBratio'])
mean(table.FBratio[table.FBratio$Group == 'Rapamycin (IP)', 'FBratio'])

# mean(table.FBratio[table.FBratio$Group == 'WT', 'FBratio'])
# [1] 1.113575
# > mean(table.FBratio[table.FBratio$Group == 'WT-RAPA', 'FBratio'])
# [1] 1.163256

rm(pseq, tempplot, temp.list, tempdata)


##### core microbiomes (ok)  #####

### dataset cleanup
# Create dataset without low coverage samples (<150 reads)
psHighCov <- prune_samples(sample_sums(temp.ps) >= 150, temp.ps)

# Prune unidentified Genuses; 
# Filter out Eukaryotes, chloroplasts and mitochondria
# ref http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
psHighCov <- subset_taxa(psHighCov,!grepl("(?i)(uncultured|unclassified)", Genus))
psHighCov <-  subset_taxa(psHighCov, Genus != "NA")
psHighCov <-  subset_taxa(psHighCov, Genus != "Unclassified")
psHighCov <-  subset_taxa(psHighCov, Genus != "unclassified")
psHighCov <-  subset_taxa(psHighCov, Genus != "uncultured")
psHighCov <-  subset_taxa(psHighCov, Genus != "-1")
psHighCov <-  subset_taxa(psHighCov, Genus != "Mitochondria")
psHighCov <-  subset_taxa(psHighCov, Genus != "Chloroplast")


### WT core
pseq <- subset_samples(psHighCov, Group == "Control") %>%
  microbiome::aggregate_taxa("Genus") %>% # needed to get genera names
  microbiome::transform("compositional")
core_taxa_standard <- core_members(pseq, detection = 0.001, prevalence = 0.8)

# cleanup taxa names from the core
core_taxa_standard <- gsub(".*ales_", "", core_taxa_standard)
core_taxa_standard <- gsub("___", "", core_taxa_standard)
core_taxa_standard <- gsub("aceae_", "aceae", core_taxa_standard)
core_taxa_standard <- gsub("aceaeunc", "aceae_unc", core_taxa_standard)
core_taxa_standard <- gsub(".*eria_", "", core_taxa_standard)

# core heatmap options
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
gray <- gray(seq(0,1,length=5))

myplot <- plot_core(pseq,
                    plot.type = "heatmap", 
                    colours = gray,
                    prevalences = prevalences,
                    detections = detections, 
                    min.prevalence = 0.3) + 
  labs(x = "Detection Threshold (Relative Abundance (%))")

# edit long names and resort according to prevalence
myplot$data$Taxa <- gsub(".*ales_", "", myplot$data$Taxa)
myplot$data$Taxa <- gsub("___", "", myplot$data$Taxa)
myplot$data$Taxa <- gsub("aceae_", "aceae", myplot$data$Taxa)
myplot$data$Taxa <- gsub("aceaeunc", "aceae_unc", myplot$data$Taxa)
myplot$data <- myplot$data %>% mutate(Taxa = fct_reorder(Taxa, Prevalence))

# 1) Add underscore between "Lachnospiraceae" and "NK4A136_group"
myplot$data$Taxa <- gsub(
  "LachnospiraceaeNK4A136_group", 
  "Lachnospiraceae_NK4A136_group", 
  myplot$data$Taxa)

# 2) Replace "LachnospiraceaeUnclassified" with "Lachnospiraceae Unclassified"
myplot$data$Taxa <- gsub(
  "LachnospiraceaeUnclassified", 
  "Lachnospiraceae Unclassified", 
  myplot$data$Taxa)

# 3) Unify any "Unclassified" text to "uncultured"
#    (using ignore.case = TRUE to catch uppercase "Unclassified" or "UNCLASSIFIED")
myplot$data$Taxa <- gsub(
  "Unclassified", 
  "uncultured", 
  myplot$data$Taxa, 
  ignore.case = TRUE)

# finally reorder or do any final manipulations
myplot$data <- myplot$data %>%
  mutate(Taxa = fct_reorder(Taxa, Prevalence))
# replot
myplot <- myplot  + scale_fill_viridis(option = "E") +
  scale_x_discrete(labels=scaleFUN) +
  labs(title = "Core microbiome: Control", 
       subtitle = "detection 0.1%, prevalence 80%, min. prevalence 0.3")
#myplot

list.Core.WT <- core_taxa_standard
plot.Core.WT <- myplot
plot.Core.WT
length(taxa(pseq))
# [1] 104
length(list.Core.WT)
# [1] 25
# therefore:
# "25 out of 104 genera are shared by 80% of the subjects
# at the minimum detection threshold of 0.1% relative abundance"
# Shetty et al. 2017


### WT-RAPA core
pseq <- subset_samples(psHighCov, Group == "Rapamycin (IP)") %>%
  microbiome::aggregate_taxa("Genus") %>% # needed to get genera names
  microbiome::transform("compositional")
core_taxa_standard <- core_members(pseq, detection = 0.001, prevalence = 0.8)

# cleanup taxa names from the core
core_taxa_standard <- gsub(".*ales_", "", core_taxa_standard)
core_taxa_standard <- gsub("___", "", core_taxa_standard)
core_taxa_standard <- gsub("aceae_", "aceae", core_taxa_standard)
core_taxa_standard <- gsub("aceaeunc", "aceae_unc", core_taxa_standard)
core_taxa_standard <- gsub(".*eria_", "", core_taxa_standard)

# core heatmap options
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
gray <- gray(seq(0,1,length=5))

myplot <- plot_core(pseq,
                    plot.type = "heatmap", 
                    colours = gray,
                    prevalences = prevalences,
                    detections = detections, 
                    min.prevalence = 0.3) + 
  labs(x = "Detection Threshold (Relative Abundance (%))")

# edit long names and resort according to prevalence
myplot$data$Taxa <- gsub(".*ales_", "", myplot$data$Taxa)
myplot$data$Taxa <- gsub("___", "", myplot$data$Taxa)
myplot$data$Taxa <- gsub("aceae_", "aceae", myplot$data$Taxa)
myplot$data$Taxa <- gsub("aceaeunc", "aceae_unc", myplot$data$Taxa)
myplot$data <- myplot$data %>% mutate(Taxa = fct_reorder(Taxa, Prevalence))

# replot
myplot <- myplot  + scale_fill_viridis(option = "E") +
  scale_x_discrete(labels=scaleFUN) +
  labs(title = "Core microbiome: Rapamycin(IP)", 
       subtitle = "detection 0.1%, prevalence 80%, min. prevalence 0.3")
# myplot

list.Core.WT_RAPA <- core_taxa_standard
plot.Core.WT_RAPA <- myplot
length(taxa(pseq))
# [1] 104
length(list.Core.WT_RAPA)
# [1] 22
# therefore:
# "22 out of 104 genera are shared by 80% of the subjects
# at the minimum detection threshold of 0.1% relative abundance"
# Shetty et al. 2017

### cleanup
rm(myplot, core_taxa_standard, prevalences, detections, gray, pseq)

# Better formatting function for percentages - removes trailing zeros
scaleFUN <- function(x) {
  num <- as.numeric(x) * 100
  # Use sprintf with %g which removes trailing zeros automatically
  paste0(sprintf("%g", num), "%")
}

##combined
# Update Control plot - remove title, reduce right margin
plot.Core.Control <- plot.Core.WT + 
  scale_x_discrete(labels = scaleFUN) +
  labs(title = "Control", subtitle = NULL, tag = "a") +  # Changed title to NULL
  theme(plot.tag.position = c(0.1, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        plot.margin = unit(c(5, 2, 5, 5), "pt"),
        legend.position = "none")  # top, right, bottom, left

# Update Rapamycin plot - remove title, reduce left margin  
plot.Core.Rapamycin <- plot.Core.WT_RAPA + 
  scale_x_discrete(labels = scaleFUN) +
  labs(title = "Rapamycin-IP",subtitle = NULL,  tag = "b") +  # Changed title to NULL
  theme(plot.tag.position = c(0.1, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        plot.margin = unit(c(5, 5, 5, 2), "pt"))  # top, right, bottom, left

# Combine side by side
plot.Core.combined <- plot.Core.Control | plot.Core.Rapamycin
plot.Core.combined

# Save at 300 dpi
ggsave(here::here("results/figS9_core_combined.png"), 
       plot.Core.combined,
       width = 13, height = 4, dpi = 300)

##### comparison of core microbiomes (ok) #####

### core vectors to presence/absence matrix
# src https://stackoverflow.com/a/53575310
ls(pattern="list.*Core.*")
# [1] "list.Core.WT"   "list.Core.WT_RAPA"    

temp_vector <- ls(pattern="list.*Core.*") 
tempdata <- t(splitstackshape:::charMat(listOfValues = mget(temp_vector), fill = 0L)) # requires library(splitstackshape)
colnames(tempdata) <- gsub(".*Core.", "", temp_vector)
colnames(tempdata) <- gsub(".list", "", colnames(tempdata))
rownames(tempdata) <- gsub(".*ales_", "", rownames(tempdata))
rownames(tempdata) <- gsub("aceae_", "aceae", rownames(tempdata))
rownames(tempdata) <- gsub("aceaeunc", "aceae_unc", rownames(tempdata))
rownames(tempdata) <- gsub("\\[", "", rownames(tempdata))
rownames(tempdata) <- gsub("\\]", "", rownames(tempdata))
rownames(tempdata) <- gsub("_", " ", rownames(tempdata))
rownames(tempdata) <- gsub("application_method", "gr.", rownames(tempdata))

colnames(tempdata) <- c("WT", "WT-RAPA")
tempdata2 <- tempdata %>% 
  as.data.frame() %>%
  rownames_to_column("Genus") %>% # requires library(tidyverse)
  pivot_longer(-c(Genus), names_to = "Group", values_to = "value")
tempdata2$Group <- factor(tempdata2$Group, levels=c("WT", "WT-RAPA"))

tempplot <- ggplot(tempdata2, aes(Genus, Group)) +
  geom_tile(aes(fill = value), colour = "black") +
  coord_equal() + 
  xlab('') + 
  ylab('Group') + 
  labs(
    title = "Core Microbiome Comparison",
    fill = "Presence"
  ) + # tiles as squares
  # coord_flip() +
  #geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "gray") +
  scale_x_discrete(position = "top") + # x labels on top
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        axis.text.y = element_text(hjust = 0),
        legend.position = "none")   # hide legend

table.Core.comparison.matrix <- tempdata
table.Core.comparison.matrix_for_ggplot2 <- tempdata2
plot.Core.comparison <- tempplot +
  scale_y_discrete(labels = c(
    WT        = "Control",
    `WT-RAPA`  = "Rapamycin (IP)"
  ))

### cleanup
rm(tempdata, tempdata2, tempplot)

##### bugbase () #####

# data acquisition and preprocessing
bugbase_data <- read.csv(inputfile.bugbase, sep = "\t", header = TRUE, row.names = 1)
rownames(tempsampledata) <- tempsampledata$X.SampleID
bugbase_data <- merge(bugbase_data, tempsampledata, by = "row.names")
rownames(bugbase_data) <- bugbase_data$Row.names
bugbase_data <- bugbase_data[, -c(1)]
bugbase_data <- bugbase_data[-c(1), ] # discard negative control
bugbase_data$Gruppe <- factor(bugbase_data$Gruppe, levels = c("Gesund", "Gingivitis", "Parodontitis"))
names(bugbase_data)[names(bugbase_data) == 'Gruppe'] <- 'Group'

# calling adhoc function
plotBugbase2 <- function(inputTable){
  my_comparisons <- my_comparisons
  aerobic.plot <- ggboxplot(bugbase_data, x = "Group", y = "Aerobic", color = "Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    xlab("") +
    ylab("proportion") +
    ggtitle("aerobic") +
    theme(legend.position = "none")
  
  anaerobic.plot <- ggboxplot(bugbase_data, x = "Group", y = "Anaerobic", color = "Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    xlab("") +
    ylab("proportion") +
    ggtitle("anaerobic") +
    theme(legend.position = "none")
  
  facultative.plot <- ggboxplot(bugbase_data, x = "Group", y = "Facultatively_Anaerobic", color = "Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    xlab("") +
    ylab("proportion") +
    ggtitle("facultatively anaerobic") +
    theme(legend.position = "none")
  
  # stresstolerant.plot <- ggboxplot(bugbase_data, x = "Group", y = "Stress_Tolerant", color = "Group", palette = "jco") +
  #   stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  #   xlab("") +
  #   ylab("proportion") +
  #   ggtitle("stress tolerant") +
  #   theme(legend.position = "none")
  
  pathogenic.plot <- ggboxplot(bugbase_data, x = "Group", y = "Potentially_Pathogenic", color = "Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    xlab("") +
    ylab("proportion") +
    ggtitle("potentially pathogenic") +
    theme(legend.position = "none")
  
  biofilms.plot <- ggboxplot(bugbase_data, x = "Group", y = "Forms_Biofilms", color = "Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    xlab("") +
    ylab("proportion") +
    ggtitle("forming biofilms") +
    theme(legend.position = "none")
  
  mobile.plot <- ggboxplot(bugbase_data, x = "Group", y = "Contains_Mobile_Elements", color = "Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    xlab("") +
    ylab("proportion") +
    ggtitle("contains mobile elements") +
    theme(legend.position = "none")
  
  # combined plot
  bugbase.plot1 <- (aerobic.plot + anaerobic.plot + facultative.plot) / 
    (pathogenic.plot + biofilms.plot + mobile.plot)
  # join legends
  bugbase.plot1 <- bugbase.plot1 + 
    plot_layout(guides = "collect") + theme(legend.position = "right") #+
  #plot_annotation(tag_levels = 'A')
  
  return(bugbase.plot1)
}

# all plots
bugbase_data_plots <- plotBugbase2(bugbase_data)
plot.BugBase <- bugbase_data_plots

# backup table
table.bugbase <- bugbase_data

# cleanup
rm(bugbase_data_plots, bugbase_data)

##### shannon diversity (ok) #####

# shannon diversity by Group
season_std_shannon <- plot_richness(temp.ps, x = "Group", color = "Group", measures = c("Shannon")) 
tempdata <- season_std_shannon$data

tempdata$Group <- factor(tempdata$Group, levels = c("Control", "Rapamycin (IP)"))

plot.Shannon.AlphaDiversity <- ggboxplot(tempdata, x = "Group", y = "value",
                                         color = "Group", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  # 2) Overall test label (Kruskal-Wallis), placed at x=1
  stat_compare_means(
    method = "kruskal.test",
    label.x = 1,
    label.y = (1.2 + max(tempdata$value))
  ) +
  ggtitle("Shannon diversity") +
  xlab("Group") +
  ylab("Shannon diversity index")+
  coord_cartesian(ylim = c(4, 5.5))+
  theme(legend.position = "none")


plot.Shannon.AlphaDiversity
# backup table
table.shannon <- tempdata

# cleanup
rm(tempdata)

###combine

# Update F/B ratio plot - consistent font sizes
plot.FBratio <- plot.FBratio + 
  scale_x_discrete(labels = c("Control" = "Control-IP", 
                              "Rapamycin (IP)" = "Rapamycin-IP")) +
  labs(tag = "a", title = NULL, subtitle = NULL, x = NULL, y = "log(F/B ratio)") +
  theme(plot.tag.position = c(-0.0, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        legend.position = "none",
        axis.text.x = element_text(size = 11),    # Add consistent size
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))

# Update Shannon diversity plot - consistent font sizes
plot.Shannon.AlphaDiversity <- plot.Shannon.AlphaDiversity + 
  scale_x_discrete(labels = c("Control" = "Control-IP", 
                              "Rapamycin (IP)" = "Rapamycin-IP")) +
  labs(tag = "b", title = NULL, x = NULL, subtitle = NULL) +
  theme(plot.tag.position = c(-0.0, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        legend.position = "none",
        axis.text.x = element_text(size = 11),    # Add consistent size
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))

# Combine side by side
library(patchwork)
plot.Shannon.FBratio.combined <- plot.FBratio | plot.Shannon.AlphaDiversity 
plot.Shannon.FBratio.combined
##### beta diversity (Bray divergence) (ok) #####

# src https://microbiome.github.io/tutorials/Betadiversity.html
# requires library(reshape)
# ref Salonen et al. ISME J 2014
# http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html

# Using temp.ps for the phyloseq object
pseq <- temp.ps

# Calculate Bray divergence for each group
b.WT <- divergence(subset_samples(pseq, Group == "Control"),
                             apply(abundances(subset_samples(pseq, Group == "Control")), 1, median))

b.WT_RAPA <- divergence(subset_samples(pseq, Group == "Rapamycin (IP)"),
                               apply(abundances(subset_samples(pseq, Group == "Rapamycin (IP)")), 1, median))


# Combine distances into a single table
l <- list(b.WT, b.WT_RAPA)
df <- reshape2::melt(l) # Requires library(reshape)

# Rename groups for clarity
df$L1[df$L1 == '1'] <- 'Control'
df$L1[df$L1 == '2'] <- 'Rapamycin (IP)'

# Set proper factor levels for ordering
df$L1 <- factor(df$L1, levels = c("Control","Rapamycin (IP)"))

# Rename columns for clarity
bray_table <- df
colnames(bray_table) <- c("value", "Group")
# Rename columns for clarity
bray_table <- df
colnames(bray_table) <- c("value", "Group")

# Plot Bray Beta Diversity
plot.Bray.BetaDiversity <- ggboxplot(bray_table, 
                                     x = "Group", 
                                     y = "value",
                                     color = "Group", 
                                     palette = "jco",
                                     ylab = "Divergence") +  # Add y-axis label
  scale_x_discrete(labels = c("Control" = "Control-IP",
                              "Rapamycin (IP)" = "Rapamycin-IP")) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = (0.1 + max(bray_table$value)), 
                     label.x = 1) +
  labs(x = NULL, title = NULL, subtitle = NULL) +  # Remove x-axis label
  theme(legend.position = "none")

plot.Bray.BetaDiversity

# Save the table for backup
table.betadiv <- bray_table


# cleanup
rm(pseq, bray_table, b.WT, b.WT_RAPA, l, df)

##### Beta Diversity by Franzi #####

# Set up phyloseq object for analysis
pseq <- temp.ps

# Ensure the phylogenetic tree is included
phy_tree(pseq)

# Calculate Weighted UniFrac distances
iDist <- phyloseq::distance(pseq, method = "wunifrac")

# Perform PCoA ordination
iMDS <- ordinate(pseq, method = "PCoA", distance = iDist)

# Run PERMANOVA (adonis) to test for significant differences in beta diversity between groups
wunifrac.adonis <- adonis(iDist ~ Group, 
                          data = as(sample_data(pseq), "data.frame"), 
                          permutations = 10000)

# Extract p-value from PERMANOVA results
tempSubtitle.here <- paste0("Adonis PERMANOVA p-value = ", 
                            formatC(wunifrac.adonis$aov.tab$`Pr(>F)`[1], format = "f", digits = 5))

# Plot the PCoA ordination, coloring by application method
plot.unifrac.pcoa <- plot_ordination(pseq, iMDS, color = "Group") +
  labs(title = "PCoA Ordination (Weighted UniFrac)", subtitle = tempSubtitle.here) +
  stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.15) +  # Add ellipses for each group
  theme_bw()

# Display the plot
print(plot.unifrac.pcoa)

dispersion <- betadisper(iDist, sample_data(pseq)$Group)
anova(dispersion)
permutest(dispersion)

# Clean up
rm(iDist, iMDS, wunifrac.adonis, tempSubtitle.here)

#####Unweighted UniFrac#######
# Set up phyloseq object for analysis
pseq <- temp.ps

# Ensure the phylogenetic tree is included
phy_tree(pseq)

# Calculate Unweighted UniFrac distances
iDist <- phyloseq::distance(pseq, method = "unifrac")

# Perform PCoA ordination
iMDS <- ordinate(pseq, method = "PCoA", distance = iDist)

# Run PERMANOVA (adonis) to test for significant differences in beta diversity between groups
unifrac.adonis2 <- adonis2(
  iDist ~ Group,
  data         = as(sample_data(pseq), "data.frame"),
  permutations = 10000,
  by           = "margin"    # tests each term marginally
)

# Correct extraction of the p-value
pval <- unifrac.adonis2["Group", "Pr(>F)"]

# Extract p-value from PERMANOVA results
tempSubtitle.here <- paste0("Adonis PERMANOVA p-value = ", 
                            formatC(pval, format = "f", digits = 5))

# Plot the PCoA ordination, coloring by application method
plot.unweightedunifrac.pcoa <- plot_ordination(pseq, iMDS, color = "Group") +
  labs(
    title    = "PCoA Ordination (Unweighted UniFrac)",
    subtitle = tempSubtitle.here,
    color    = "Group"   # relabel the points/colour legend
  )  +
  stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.15, 
               show.legend = FALSE) +     # drop the fill legend here
  theme_bw()

dispersion <- betadisper(iDist, sample_data(pseq)$Group)
anova(dispersion)
permutest(dispersion)

# Display the plot
print(plot.unweightedunifrac.pcoa)

# Clean up
rm(iDist, iMDS, unifrac.adonis2, tempSubtitle.here)

##combine
# Update Weighted UniFrac plot - remove title/subtitle, change labels, add tag, remove legend
# Update Weighted UniFrac plot
plot.unifrac.pcoa <- plot.unifrac.pcoa + 
  scale_color_discrete(labels = c("Control" = "Control-IP", 
                                  "Rapamycin (IP)" = "Rapamycin-IP")) +
  scale_fill_discrete(labels = c("Control" = "Control-IP", 
                                 "Rapamycin (IP)" = "Rapamycin-IP")) +
  
  labs(tag = "a", title = NULL, subtitle = NULL) +
  theme(plot.tag.position = c(-0.00, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        legend.position = "none")

# Update Unweighted UniFrac plot
plot.unweightedunifrac.pcoa <- plot.unweightedunifrac.pcoa + 
  scale_color_discrete(labels = c("Control" = "Control-IP", 
                                  "Rapamycin (IP)" = "Rapamycin-IP")) +
  scale_fill_discrete(labels = c("Control" = "Control-IP", 
                                 "Rapamycin (IP)" = "Rapamycin-IP")) +
  
  labs(tag = "b", title = NULL, subtitle = NULL) +
  theme(plot.tag.position = c(-0.00, 1),
        plot.tag = element_text(size = 14, face = "plain"))

# Combine side by side
library(patchwork)
plot.unifrac.combined <- plot.unifrac.pcoa | plot.unweightedunifrac.pcoa
plot.unifrac.combined

##### Differential abundances: WT vs WT-RAPA  (ok) #####
# FDR cutoff for differential abundances
fdr.threshold <- 0.01

# Subset to include only the two groups we want to compare
temp.diagdds <- subset_samples(temp.ps, Group %in% c("Control", "Rapamycin (IP)"))

# Create DESeq2 object
diagdds <- phyloseq_to_deseq2(temp.diagdds, ~ Group)

# Set WT as the reference level
diagdds$Group <- relevel(factor(diagdds$Group), ref = "Control")

# Avoid zeros error with geometric means
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local")

#diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
res <- results(diagdds, cooksCutoff = FALSE)
nrow(res)

print(paste("Total number of ASVs tested:", nrow(res)))

# # Filter out low abundance taxa Love (2014), Weiss S, et al. (2017).
# min_count <- 5  # Minimum count in at least one group
# # Get normalized counts
# normalized_counts <- counts(diagdds, normalized=TRUE)
# # Calculate mean counts by group
# mean_counts_control <- rowMeans(normalized_counts[, colData(diagdds)$application_method == "Control_food"])
# mean_counts_rapamycin <- rowMeans(normalized_counts[, colData(diagdds)$application_method == "Rapamycin_food"])
# # Apply abundance filter before statistical filtering
# abundant_taxa <- which(mean_counts_control > min_count | mean_counts_rapamycin > min_count)
# res_filtered <- res[abundant_taxa, ]

# Get total counts per sample
sample_totals <- colSums(counts(diagdds))

# Calculate relative abundance (proportion of reads)
rel_abundance <- sweep(counts(diagdds), 2, sample_totals, "/")

# Set minimum relative abundance threshold (e.g., 0.001 = 0.1%)
min_rel_abundance <- 0.0025  # Adjust this value as needed

# Calculate mean relative abundance by group
mean_rel_abundance_control <- rowMeans(rel_abundance[, colData(diagdds)$Group == "Control"])
mean_rel_abundance_rapamycin <- rowMeans(rel_abundance[, colData(diagdds)$Group == "Rapamycin (IP)"])

# Apply relative abundance filter
abundant_taxa <- which(mean_rel_abundance_control > min_rel_abundance | 
                         mean_rel_abundance_rapamycin > min_rel_abundance)
res_filtered <- res[abundant_taxa, ]

# > head(res)
# log2 fold change (MLE): application method Rapamycin food vs Control food 
# Wald test p-value: application method Rapamycin food vs Control food 
# i.e. log2(WT-RAPA/WT )

# "The text, condition treated vs untreated, tells you that the 
# estimates are of the logarithmic fold change log2(treated/untreated)."
# source https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

print(paste("ASVs after minimum count filtering:", nrow(res_filtered)))

# filter results by FDR (p-adj) < 0.01
# Apply both statistical significance and fold change filters
sigtab <- res_filtered[which(!is.na(res_filtered$padj) & 
                               !is.na(res_filtered$log2FoldChange) & 
                               res_filtered$padj < fdr.threshold & 
                               abs(res_filtered$log2FoldChange) > 3), ]

sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(temp.diagdds)[rownames(sigtab), ], "matrix"))

nrow(sigtab)
length(unique(sigtab$Genus))


# # correcting cell without genera info
# # sigtab[rownames(sigtab) == "f33bc888923d4c7b93d34b684d1cc2c7",]
# sigtab$Genus <- as.character(sigtab$Genus)
# sigtab[rownames(sigtab) == "f33bc888923d4c7b93d34b684d1cc2c7", "Genus"] <- "[Prevotellaceae]"
# sigtab$Genus <- as.factor(sigtab$Genus)

# Merge unclassified and uncultured taxa
sigtab$Genus <- as.character(sigtab$Genus)  # Convert to character
# Replace 'unclassified' and variations
unclassified_patterns <- c("unclassified", " Unclassified","Unclassified", "uncultured", "Uncultured","Mitochondria")
for (pattern in unclassified_patterns) {
  sigtab$Genus[grepl(pattern, sigtab$Genus)] <- "Unclassified"
}


# Phylum order
# x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
# x <- sort(x, TRUE)
# sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Sort by max LFC at genus level
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# plot
tempTitle <- "Differential Abundance by Genera"
tempSubtitle <- "Rapamycin (IP) vs Control; FDR < 0.05"
plot.DiffAbund.WT <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  labs(title = tempTitle, subtitle = tempSubtitle) +
  geom_hline(yintercept=0)
plot.DiffAbund.WT




# diffex table
table.DiffAbund.WT <- sigtab[,c(7:13,1:6)]
# DiffAbund.Gingivitis.Table <- sigtab[,c(7:13,1:6)]
# run in case it's needed on excel
# clip <- pipe("pbcopy", "w")
# write.table(DiffAbund.Gingivitis.Table[,c(2,6,9,13)], file=clip, 
#             col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
# close(clip)

# Take top 50 instead of 100
# Sort by significance (lowest padj)
ordered_sigtab <- sigtab[order(sigtab$padj), ]
top50_sigtab <- head(ordered_sigtab, 50)

# Sort these by fold change for better visualization
x_top50 <- tapply(top50_sigtab$log2FoldChange, top50_sigtab$Genus, function(x) max(x))
x_top50 <- sort(x_top50, TRUE)
top50_sigtab$Genus = factor(as.character(top50_sigtab$Genus), levels=names(x_top50))

# Create plot for top 50
tempTitle_top50 <- "Top 50 Differentially Abundant Genera"
tempSubtitle_top50 <- paste0("Rapamycin vs Control; FDR < ", fdr.threshold, ", |log2FC| > 3")

# Update plot with top 50
plot.Top50.DiffAbund.WT <- ggplot(top50_sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=6) +  # Slightly larger points since fewer data
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=9)) +  # Slightly larger text
  theme(legend.position = "right") +  # Move legend to right
  labs(title = tempTitle_top50, subtitle = tempSubtitle_top50) +
  geom_hline(yintercept=0)

# Display the plot
plot.Top50.DiffAbund.WT

####new ####
# Merge unclassified and uncultured taxa
sigtab$Genus <- as.character(sigtab$Genus)
unclassified_patterns <- c("unclassified", " Unclassified","Unclassified", "uncultured", "Uncultured","Mitochondria")
for (pattern in unclassified_patterns) {
  sigtab$Genus[grepl(pattern, sigtab$Genus)] <- "Unclassified"
}

# Sort by max LFC at genus level
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# Add direction
sigtab$Direction <- ifelse(sigtab$log2FoldChange > 0, "Enriched", "Depleted")

# plot

plot.DiffAbund.WT.clean <- ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Direction)) + 
  geom_vline(xintercept=0, linetype="solid", color="black", linewidth=0.5) +
  geom_point(size=5) + 
  scale_color_manual(values = c("Enriched" = "#00BFC4", "Depleted" = "#F8766D")) +
  scale_x_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by=5)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=10),
    legend.position = "none",
    plot.title = element_text(face="bold", hjust=0, size=12)
  ) +
  labs(title = tempTitle, x = expression(log[2]~"fold change"), y = NULL)

plot.DiffAbund.WT.clean

##### Differential abundances: Pathways () #####

# Remove only pathway-related objects
rm(path_abund, picrust2Abund, abundance.table, wt_samples, wt_conditions, 
   wt_samples_meta, wt_samples_meta, wt_dds, table.DiffPathw.WT, 
   plot.DiffPathw.WT, clean_df)


# table abundances with KO ids
# FDR cutoff for differential abundances
fdr.threshold <- 0.05


# Load PICRUSt2 pathway abundance data with the full file path
path_abund <- read.table(here::here("data/study5_PRJNA717340/picrust2/path_abun_unstrat.tsv.gz"), 
                         sep="\t",
                         header=TRUE,
                         row.names=1, 
                         stringsAsFactors=FALSE,
                         comment.char="")

# Format similar to original picrust1Abund
picrust2Abund <- data.frame(path_abund)

# Print dimension for sanity check
cat("Dimension of the raw PICRUSt2 abundance table:\n")
print(dim(picrust2Abund))

# REPLACE with:
picrust2Abund$path.id <- rownames(picrust2Abund)  # Use the real pathway IDs
picrust2Abund$kegg.desc <- rownames(picrust2Abund)  # Also use for descriptions

# Extract the abundance table
abundance.table <- picrust2Abund[, 1:(ncol(picrust2Abund)-2)]
rownames(abundance.table) <- picrust2Abund$path.id

cat("Dimension of raw abundance.table:\n")
print(dim(abundance.table))
cat("First few column names of abundance.table:\n")
print(head(colnames(abundance.table)))

# Make the row names of SampleData = Sample.id for easy subsetting
rownames(metadata) <- metadata$Sample.id

cat("\nSampleData:\n")
print(head(metadata))

common_samples <- intersect(colnames(abundance.table), rownames(metadata))
length(common_samples)  # how many samples overlap

# 2) Subset abundance.table to only those samples
abundance.table <- abundance.table[, common_samples, drop=FALSE]

# 3) Order the columns to match the row order of SampleData exactly
#    (if you want them in that specific order)
abundance.table <- abundance.table[, rownames(metadata)]

cat("\nDimension of abundance.table after aligning samples:\n")
print(dim(abundance.table))
cat("Colnames now:\n")
print(colnames(abundance.table))

minCount <- 1000
abundance.table <- abundance.table[rowSums(abundance.table) > minCount, ]

abundance.table.rounded <- round(abundance.table)

cat("\nDimension of 'abundance.table' after minCount filter:\n")
print(dim(abundance.table.rounded))

# # Additional per-sample filter as suggested by professor
# min_per_sample <- 10  # Set to 100-300 as suggested 
# 
# # Filter to keep only pathways with at least min_per_sample in each sample
# # This function checks if all values in a row are >= min_per_sample
# keep_pathways <- apply(abundance.table, 1, function(pathway_counts) {
#   all(pathway_counts >= min_per_sample)
# })
# 
# # Apply the per-sample filter
# abundance.table.per_sample <- abundance.table[keep_pathways, ]
# abundance.table.per_sample.rounded <- round(abundance.table.per_sample)
# 
# cat("\nDimension after additional per-sample filter (min", min_per_sample, "per sample):\n")
# print(dim(abundance.table.per_sample.rounded))
# 
# # Replace the existing filtered table with the per-sample filtered one
# abundance.table <- abundance.table.per_sample
# abundance.table.rounded <- abundance.table.per_sample.rounded
# # Suppose your metadata data frame is 'metadata' (with rownames = SRR IDs).
# # 'Group' is the column that has "WT", "WT-RAPA", "RhoB", etc.

# 1) Identify which samples are WT or WT-RAPA
samples_of_interest <- rownames(metadata)[ metadata$Group %in% c("WT", "WT-RAPA") ]

# 2) Subset the abundance table to these samples
wt_samples <- abundance.table.rounded[, samples_of_interest, drop=FALSE ]

# 3) Subset the metadata accordingly
wt_samples_meta <- metadata[samples_of_interest, , drop=FALSE ]

# Check dimension
dim(wt_samples)
# This should be (#pathways) x (number of WT + number of WT-RAPA samples)

# Also confirm the distribution of groups
table(wt_samples_meta$Group)

wt_dds <- DESeqDataSetFromMatrix(
  countData = wt_samples,
  colData   = wt_samples_meta,
  design    = ~ Group
)


# IMPORTANT: Set reference level to "WT" 
wt_dds$Group <- relevel(wt_dds$Group, ref = "WT")

# 2) Run DESeq2 (only once)
wt_dds <- DESeq(wt_dds, test="Wald", fitType="parametric")


# After DESeq2:
res <- results(wt_dds, contrast = c("Group", "WT-RAPA", "WT"))
res_df <- as.data.frame(res)
res_df$feature <- rownames(res_df)  # keep track of feature ID

# Merge in annotations (if applicable)
res_annot <- merge(
  res_df,
  picrust2Abund[, c("path.id", "kegg.desc")],
  by.x = "feature",  # or "row.names" if you didn't copy them
  by.y = "path.id",
  all.x = TRUE
)

# Sort by raw p-value
res_annot <- res_annot[order(res_annot$pvalue), ]

# Look at top 20 unadjusted p-values
head(res_annot, 20)

# 3) Extract results for WT-RAPA vs WT
table.DiffPathw.WT <- results(wt_dds, contrast = c("Group", "WT-RAPA", "WT"))
# Convert to data frame, filter by FDR
table.DiffPathw.WT <- as.data.frame(table.DiffPathw.WT)
table.DiffPathw.WT <- subset(table.DiffPathw.WT, padj < fdr.threshold)

cat("\nNumber of significant pathways at padj < ", fdr.threshold, ":\n", nrow(table.DiffPathw.WT), "\n")

hist(res$pvalue, breaks=50, main="P-value distribution", xlab="pvalue")

res_noIF <- results(
  wt_dds, 
  contrast = c("Group","WT-RAPA","WT"),
  independentFiltering = FALSE
)
# Now see how many are NA
table(is.na(res_noIF$padj))

summary(res)


# 4) Merge with 'picrust2Abund' for descriptions
table.DiffPathw.WT <- merge(
  table.DiffPathw.WT,
  picrust2Abund[, c("path.id", "kegg.desc")],
  by.x = "row.names",
  by.y = "path.id"
)

# Rename for clarity
colnames(table.DiffPathw.WT)[colnames(table.DiffPathw.WT) == "Row.names"] <- "path.id"

cat("\nFinal columns in table.DiffPathw.WT:\n")
print(colnames(table.DiffPathw.WT))

# 5) Plot a bar chart
plot.DiffPathw.WT <- ggplot(
  data = table.DiffPathw.WT,
  aes(
    x   = reorder(kegg.desc, log2FoldChange),
    y   = log2FoldChange,
    fill = log2FoldChange > 0
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  guides(fill = "none") +
  theme_minimal() +
  labs(
    x = "MetaCyc Pathways",
    y = "Log2(Fold Change)",
    title = "Rapamycin (IP) vs Control",
    subtitle = paste0("FDR < ", fdr.threshold)
  )

print(plot.DiffPathw.WT)

# ONLY ADDITION: Simple printout of significant pathways sorted by fold change
cat("\n\nSignificant pathways sorted by fold change magnitude:\n")
pathways_sorted <- table.DiffPathw.WT[order(-abs(table.DiffPathw.WT$log2FoldChange)), ]
result_table <- data.frame(
  Pathway = pathways_sorted$kegg.desc,
  Log2FC = round(pathways_sorted$log2FoldChange, 3),
  Padj = format(pathways_sorted$padj, digits=3, scientific=TRUE)
)
print(result_table)

tempTitle <- "Differential Pathway Abundance"
tempSubtitle <- paste0("Rapamycin_IP vs Vehicle_IP; FDR < ", fdr.threshold)
plot.DiffPathw.IP <- ggplot(table.DiffPathw.WT,
                            aes(x = kegg.desc, y = log2FoldChange,  color = log2FoldChange > 0)) +
  geom_point(size = 6) +
  scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#00BFC4"),
                     breaks = c(FALSE, TRUE),
                     labels = c("Negative", "Positive")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        legend.position = "none") +
  labs(title = tempTitle,
       subtitle = tempSubtitle,
       x = "MetaCyc Pathways",
       y = "Log2(Fold Change)") +
  geom_hline(yintercept = 0)
print(plot.DiffPathw.IP)
#### new ####
##### Pathway merged ####
# Add direction variable
table.DiffPathw.WT$Direction <- ifelse(table.DiffPathw.WT$log2FoldChange > 0, "Enriched", "Depleted")

# Sort by fold change
table.DiffPathw.WT <- table.DiffPathw.WT[order(table.DiffPathw.WT$log2FoldChange, decreasing = FALSE), ]
table.DiffPathw.WT$kegg.desc <- factor(table.DiffPathw.WT$kegg.desc, levels = table.DiffPathw.WT$kegg.desc)

# Create bar plot
plot.DiffPathw.IP.clean <- ggplot(table.DiffPathw.WT,
                                  aes(y = kegg.desc, x = log2FoldChange, fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Enriched" = "#17BECF", "Depleted" = "#FC8D62")) +
  scale_x_continuous(limits = c(-2, 6), breaks = seq(-4, 6, by = 1)) +  # Zero shifted right
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = expression(log[2]~"fold change"), 
       y = NULL)

print(plot.DiffPathw.IP.clean)

##### clinical markers data () #####

### preparing clinical markers' data
# (clinical-marker block below is from an analysis template and is not used in this study)
# and replacing commas for points
table.ClinicalMarkers <- read.delim(pipe("pbpaste"), dec = ",")

# simplified table of clinical markers
table.ClinicalMarkers.simplified <- table.ClinicalMarkers[, -c(2:4)]
# add here leading string 
table.ClinicalMarkers.simplified$Sample_ID <- sub("^", "Zahn", table.ClinicalMarkers.simplified$Sample_ID)
colnames(table.ClinicalMarkers.simplified)[colnames(table.ClinicalMarkers.simplified) == 'STATUS'] <- 'Gruppe'
table.ClinicalMarkers.simplified$Gruppe <- gsub("healthy", "Gesund", table.ClinicalMarkers.simplified$Gruppe)
table.ClinicalMarkers.simplified$Gruppe <- gsub("gingivitis", "Gingivitis", table.ClinicalMarkers.simplified$Gruppe)
table.ClinicalMarkers.simplified$Gruppe <- gsub("periodontitis", "Parodontitis", table.ClinicalMarkers.simplified$Gruppe)
colnames(table.ClinicalMarkers.simplified)[colnames(table.ClinicalMarkers.simplified) == 'IL1β'] <- 'IL1beta'
colnames(table.ClinicalMarkers.simplified)[colnames(table.ClinicalMarkers.simplified) == 'TNFα'] <- 'TNFalpha'
colnames(table.ClinicalMarkers.simplified)[colnames(table.ClinicalMarkers.simplified) == 'IFNγ'] <- 'IFNgamma'
colnames(table.ClinicalMarkers.simplified)[colnames(table.ClinicalMarkers.simplified) == 'TGFβ1'] <- 'TGFbeta1'
table.ClinicalMarkers.simplified <- table.ClinicalMarkers.simplified[complete.cases(table.ClinicalMarkers.simplified),]

# obtaining the means by Gruppe
table.ClinicalMarkers.means <- table.ClinicalMarkers.simplified %>%
  group_by(Gruppe) %>%
  summarise(across(everything(), mean))


##### correlation: Gingi vs clinical markers (ok -corr) #####

# otu table
otu_corr <- aggregate_taxa(temp.ps, 'Genus')
otu_corr <- as(otu_table(otu_corr), "matrix")
if(taxa_are_rows(temp.ps)){otu_corr <- t(otu_corr)}
otu_corr <- as.data.frame(otu_corr)

# Define data sets to cross-correlate
x <- log10(otu_corr) # OTU Log10
#y <- as.matrix(table.clinical.markers.sorted[, c(6:18,20)]) # diff clin markers

### select only those for Gingi
temptable_Gingi <- table.ClinicalMarkers.simplified %>% filter(Gruppe == "Gingivitis")
x <- merge(temptable_Gingi, x, by.x = 1, by.y=0)
rownames(x) <- x[, 1]
x <- as.matrix(x[, c(27:ncol(x))])
y <- as.matrix(temptable_Gingi[, 3:26])
rownames(y) <- temptable_Gingi$Sample_ID

### corr
table.correlation.Gingi <- microbiome::associate(x, y, 
                                                 method = "spearman", 
                                                 mode = "table",
                                                 p.adj.threshold = 0.05)
# > table.correlation.Gingi[table.correlation.Gingi$p.adj < 0.05,]
# X1    X2 Correlation      p.adj
# 3193  Megasphaera CXCL8   0.8475657 0.01601225


##### correlation: Perod genera vs clinical markers (ok -no corr) #####

# otu table
otu_corr <- aggregate_taxa(temp.ps, 'Genus')
otu_corr <- as(otu_table(otu_corr), "matrix")
if(taxa_are_rows(temp.ps)){otu_corr <- t(otu_corr)}
otu_corr <- as.data.frame(otu_corr)

# Define data sets to cross-correlate
x <- log10(otu_corr) # OTU Log10
#y <- as.matrix(table.clinical.markers.sorted[, c(6:18,20)]) # diff clin markers

### select only those for Gingi
temptable_Perod <- table.ClinicalMarkers.simplified %>% filter(Gruppe == "Parodontitis")
x <- merge(temptable_Perod, x, by.x = 1, by.y=0)
rownames(x) <- x[, 1]
x <- as.matrix(x[, c(27:ncol(x))])
y <- as.matrix(temptable_Perod[, 3:26])
rownames(y) <- temptable_Perod$Sample_ID

### corr
table.correlation.Perod <- microbiome::associate(x, y, 
                                                 method = "spearman", 
                                                 mode = "table",
                                                 p.adj.threshold = 0.05)
# No significant correlations 


##### correlation: Healthy vs clinical markers (ok -no corr) #####

# otu table
otu_corr <- aggregate_taxa(temp.ps, 'Genus')
otu_corr <- as(otu_table(otu_corr), "matrix")
if(taxa_are_rows(temp.ps)){otu_corr <- t(otu_corr)}
otu_corr <- as.data.frame(otu_corr)

# Define data sets to cross-correlate
x <- log10(otu_corr) # OTU Log10
#y <- as.matrix(table.clinical.markers.sorted[, c(6:18,20)]) # diff clin markers

### select only those for Gingi
temptable_Gesnd <- table.ClinicalMarkers.simplified %>% filter(Gruppe == "Gesund")
x <- merge(temptable_Gesnd, x, by.x = 1, by.y=0)
rownames(x) <- x[, 1]
x <- as.matrix(x[, c(27:ncol(x))])
y <- as.matrix(temptable_Gesnd[, 3:26])
rownames(y) <- temptable_Gesnd$Sample_ID

### corr
table.correlation.Gesnd <- microbiome::associate(x, y, 
                                                 method = "spearman", 
                                                 mode = "table",
                                                 p.adj.threshold = 0.05)
# No significant correlations

# cleanup
rm(x, y, otu_corr, temptable_Gingi, temptable_Perod, temptable_Gesnd)






##### Supervised ML: microbial biomarkers via RF classifiers for P vs Gingi (ok) #####
# src tempScript07

### whole random forest vs microbiome code here:
# src https://github.com/michberr/randomforest-microbe/blob/master/microbe_rf.md

# > tempdata <- microbiome::summarize_phyloseq(temp.ps)[[1]]
# > tempdata
# [1] "1] Min. number of reads = 26514"
# > rm(tempdata)

# Set minlib
minlib = 25000

# Select unfractionated samples (whole bacterial community)
# Select our two sites; WE2 is nearshore, WE4 is offshore
# Scale reads to minlib library size
sites.scale  <- temp.ps %>%
  subset_samples(Gruppe %in% c("Gingivitis", "Parodontitis")) %>%
  scale_reads(n = minlib, round = "round")  

# Set prunescale 
prunescale = 0.0001 # why this value??

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(sites.scale)/nsamples(sites.scale)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, sites.scale)

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(sites.prune))

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(sites.prune)$Gruppe)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)
# > dim(rf.data)
# [1]  39 919
# hence we have 39 samples and 918 OTUs

# Now we will use the randomForest package to train and test our random forest model 
# using the "out of bag" error to estimate our model error. OOB is a nice feature of 
# random forest models whereby since the training data is bootstrapped, you only use 
# approximately 2/3 of the data at each iteration. The remaining 1/3 or "out of bag" 
# data can be used to validate your model. This removes the need to use another form 
# of cross-validation such as using a separate validation set or k-folds.
# 
# It is important to set a seed for reproducability By default, randomForest uses p/3 
# variables when building a random forest of regression trees and root(p) variables 
# when building a random forest of classification trees. In this case, p/3 = 129

set.seed(2)
dental.classify <- randomForest(response~., data = rf.data, ntree = 100) # library(randomForest)
# > print(dental.classify)
# Call:
#   randomForest(formula = response ~ ., data = rf.data, ntree = 100) 
# Type of random forest: classification
# Number of trees: 100
# No. of variables tried at each split: 30
# 
# OOB estimate of  error rate: 25.64%
# Confusion matrix:
#   Gingivitis Parodontitis class.error
# Gingivitis           16            3   0.1578947
# Parodontitis          7           13   0.3500000

# Make a data frame with predictor names and their importance
imp <- importance(dental.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]

# edited ids
imp.20$predictors <- gsub("^X", "", imp.20$predictors)

### adding taxa information
imp.20 <- merge(imp.20, taxmat[, c(5:7)], by.x=1, by.y=0, all.x=T)

# removing whitespaces
imp.20 <- imp.20 %>% 
  mutate(across(where(is.character), str_trim))

# filter out undefined species names
imp.20 <- imp.20 %>%
  filter(!grepl('_gut', Species)) %>%
  filter(!grepl('uncultured', Species)) %>%
  filter(!grepl('metagenome', Species))

# cleanup
imp.20$Family <- gsub("f__", "", imp.20$Family)
imp.20$Genus <- gsub("g__", "", imp.20$Genus)
imp.20$Species <- gsub("s__", "", imp.20$Species)
imp.20[imp.20==""] <- NA

# lowercase genera and species (uncultured etc) to NA
imp.20$Genus <- ifelse(str_detect(imp.20$Genus,"^[[:lower:]]"),
                       NA,
                       as.character(imp.20$Genus))
imp.20$Species <- ifelse(str_detect(imp.20$Species,"^[[:lower:]]"),
                         NA,
                         as.character(imp.20$Species))

# if species, species; otherwise genus or family info
imp.20$taxid <- ifelse(is.na(imp.20$Species),  
                       ifelse(is.na(imp.20$Genus), 
                              as.character(imp.20$Family), 
                              as.character(imp.20$Genus) ),
                       as.character(imp.20$Species))
imp.20$taxid <- gsub("_", " ", imp.20$taxid)

# add number if genus name repeats with different OTU ids and different LFCs
imp.20 <- imp.20 %>% 
  group_by(taxid) %>% 
  mutate(taxid2 = paste(taxid, row_number(), sep = "_"))


# # ggplot
# ggplot(imp.20, aes(x = taxid, y = MeanDecreaseGini)) +
#   geom_bar(stat = "identity", fill = "indianred") +
#   coord_flip() +
#   ggtitle("Most important OTUs for classifying into PLC_T2 or INT_T2")
# 
# 
# ggbarplot(imp.20, x = "taxid2", y = "MeanDecreaseGini",
#           # fill = "mpg_grp",           # change fill color by mpg_level
#           # color = "white",            # Set bar border colors to white
#           palette = "jco",            # jco journal color palett. see ?ggpar
#           sort.val = "desc",          # Sort the value in descending order
#           sort.by.groups = FALSE,     # Don't sort inside each group
#           x.text.angle = 90,          # Rotate vertically x axis texts
#           ylab = "Mean decrease Gini",
#           xlab = "Taxa",
#           # legend.title = "Most important OTUs for classifying into PLC_T2 or INT_T2",
#           rotate = TRUE,
#           ggtheme = theme_minimal())

# output plot
plot.RF.Parod_Gingi <- ggdotchart(imp.20, x = "taxid2", y = "MeanDecreaseGini",
                                  color = "#00AFBB",                                # Color by groups
                                  # palette = c("#00AFBB"), # Custom color palette
                                  sorting = "descending",                       # Sort value in descending order
                                  rotate = TRUE,                                # Rotate vertically
                                  dot.size = 2,                                 # Large dot size
                                  y.text.col = FALSE,                            # Color y text by groups
                                  ggtheme = theme_pubr())  +                      # ggplot2 theme) 
  theme_cleveland() #+ 
# ggtitle("Microbial biomarkers via supervised ML (RF)") 


# output tables
table.RF.Parod_Gingi <- imp.20
output.RF.Parod_Gingi <- dental.classify

# cleanup
rm(tax.mean, sites.prune, sites.scale, predictors, response, imp, imp.sort, imp.20, dental.classify)








##### Supervised ML: microbial biomarkers via RF classifiers for ALL groups (ok) #####

# Set minlib
minlib = 25000

# Select unfractionated samples (whole bacterial community)
# Select our 3 types
# Scale reads to minlib library size
sites.scale  <- temp.ps %>%
  subset_samples(Gruppe %in% c("Gesund", "Gingivitis", "Parodontitis")) %>%
  scale_reads(n = minlib, round = "round")  

# Set prunescale 
prunescale = 0.0001 # why this value??

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(sites.scale)/nsamples(sites.scale)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, sites.scale)

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(sites.prune))

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(sites.prune)$Gruppe)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)
# > dim(rf.data)
# [1]  59 886
# hence we have 59 samples and 885 OTUs

# train + test RF model 
set.seed(2)
temp.classify <- randomForest(response~., data = rf.data, ntree = 100) # library(randomForest)

# > print(temp.classify)
# 
# Call:
#   randomForest(formula = response ~ ., data = rf.data, ntree = 100) 
# Type of random forest: classification
# Number of trees: 100
# No. of variables tried at each split: 29
# 
# OOB estimate of  error rate: 57.63%
# Confusion matrix:
#   Gesund Gingivitis Parodontitis class.error
# Gesund            8          8            4   0.6000000
# Gingivitis       12          4            3   0.7894737
# Parodontitis      4          3           13   0.3500000


# Make a data frame with predictor names and their importance
imp <- importance(temp.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]

# edited ids
imp.20$predictors <- gsub("^X", "", imp.20$predictors)

### adding taxa information
imp.20 <- merge(imp.20, taxmat[, c(5:7)], by.x=1, by.y=0, all.x=T)

# removing whitespaces
imp.20 <- imp.20 %>% 
  mutate(across(where(is.character), str_trim))

# filter out undefined species names
imp.20 <- imp.20 %>%
  filter(!grepl('_gut', Species)) %>%
  filter(!grepl('uncultured', Species)) %>%
  filter(!grepl('metagenome', Species))

# cleanup
imp.20$Family <- gsub("f__", "", imp.20$Family)
imp.20$Genus <- gsub("g__", "", imp.20$Genus)
imp.20$Species <- gsub("s__", "", imp.20$Species)
imp.20[imp.20==""] <- NA

# lowercase genera and species (uncultured etc) to NA
imp.20$Genus <- ifelse(str_detect(imp.20$Genus,"^[[:lower:]]"),
                       NA,
                       as.character(imp.20$Genus))
imp.20$Species <- ifelse(str_detect(imp.20$Species,"^[[:lower:]]"),
                         NA,
                         as.character(imp.20$Species))

# if species, species; otherwise genus or family info
imp.20$taxid <- ifelse(is.na(imp.20$Species),  
                       ifelse(is.na(imp.20$Genus), 
                              as.character(imp.20$Family), 
                              as.character(imp.20$Genus) ),
                       as.character(imp.20$Species))
imp.20$taxid <- gsub("_", " ", imp.20$taxid)

# add number if genus name repeats with different OTU ids and different LFCs
imp.20 <- imp.20 %>% 
  group_by(taxid) %>% 
  mutate(taxid2 = paste(taxid, row_number(), sep = "_"))


# output plot
plot.RF.all <- ggdotchart(imp.20, x = "taxid2", y = "MeanDecreaseGini",
                          color = "#00AFBB",                                # Color by groups
                          # palette = c("#00AFBB"), # Custom color palette
                          sorting = "descending",                       # Sort value in descending order
                          rotate = TRUE,                                # Rotate vertically
                          dot.size = 2,                                 # Large dot size
                          y.text.col = FALSE,                            # Color y text by groups
                          ggtheme = theme_pubr())  +                      # ggplot2 theme) 
  theme_cleveland() #+ 


# output tables
table.RF.all <- imp.20
output.RF.all <- temp.classify

# cleanup
rm(tax.mean, sites.prune, sites.scale, predictors, response, imp, imp.sort, imp.20, temp.classify)


### Testing partial dependence

# Using highest importance, Tannerella forsythia_1
# > table.RF.all[table.RF.all$Genus == "Tannerella", c("predictors", "taxid2")]
# # A tibble: 3 × 2
# predictors                       taxid2                
# <chr>                            <chr>                 
#   1 53f1dc3d04077a8ba28384d2886e2aab Tannerella forsythia_1
# 2 bee31f158e7de8077bdefd5905066c8f Tannerella forsythia_2
# 3 NA                               NA    
partialPlot(output.RF.all, rf.data, X53f1dc3d04077a8ba28384d2886e2aab, "Parodontitis",
            xlab = "Tannerella forsythia 1",
            main = "Parodontitis: Partial dependence on Tannerella forsythia 1")
partialPlot(output.RF.all, rf.data, X53f1dc3d04077a8ba28384d2886e2aab, "Gingivitis",
            xlab = "Tannerella forsythia 1",
            main = "Gingivitis: Partial dependence on Tannerella forsythia 1")
partialPlot(output.RF.all, rf.data, X53f1dc3d04077a8ba28384d2886e2aab, "Gesund",
            xlab = "Tannerella forsythia 1",
            main = "Gesund: Partial dependence on Tannerella forsythia 1")

# > partialPlot(output.RF.all, rf.data, X53f1dc3d04077a8ba28384d2886e2aab, "Gingivitis")
# > partialPlot(output.RF.all, rf.data, X53f1dc3d04077a8ba28384d2886e2aab, "Gesund")




##### Mediation (GLMs): genera assoc to differences between groups (ok) #####
# src tempScript19.R
# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/differential-abundance-testing-for-univariate-data.html#compare-results-between-parametric-and-non-parametric-tests
# library(microbiome)

### core
tempdata <- as.data.frame(table.Core.comparison.matrix)
temp.core.list <- gsub(" ", "", rownames(tempdata[tempdata$Gesund ==0, ]))
# [1] " Aggregatibacter" " Capnocytophaga"  " Porphyromonas"   " Selenomonas"     " Tannerella"  

### DAGs
tempdata <- rbind(table.DiffAbund.Gingivitis, table.DiffAbund.Parodontitis, table.DiffAbund.ParodGingi)
tempdata$Genus <- as.character(tempdata$Genus)
tempdata <- tempdata %>% mutate(across(where(is.character), str_trim))
temp.DAGs.list <- unique(tempdata$Genus)

### RF
temp.RF.list <- unique(table.RF.Parod_Gingi$Genus)

### all targets
input.GLM.list <- unique(c(temp.core.list, temp.DAGs.list, temp.RF.list))
input.GLM.list <- trimws(gsub('\\[.*', '', input.GLM.list))
input.GLM.list <- trimws(gsub('_', '', input.GLM.list))
input.GLM.list <- trimws(gsub('-', '', input.GLM.list))
input.GLM.list <- trimws(gsub("\\w*[0-9].*", "", input.GLM.list))
input.GLM.list <- input.GLM.list[!is.na(input.GLM.list)]
input.GLM.list <- input.GLM.list[nzchar(input.GLM.list)]

# > input.GLM.list
# [1] "Aggregatibacter"     "Capnocytophaga"      "Porphyromonas"       "Selenomonas"         "Tannerella"         
# [6] "Corynebacterium"     "Streptococcus"       "Fusobacterium"       "Prevotella"          "Rothia"             
# [11] "Veillonella"         "Desulfobulbus"       "Fretibacterium"      "Treponema"           "Filifactor"         
# [16] "Haemophilus"         "Gemella"             "Dialister"           "Eikenella"           "Propionivibrio"     
# [21] "Undibacterium"       "Lachnoanaerobaculum" "Sphingomonas"        "Serratia"            "Saccharimonadales"  
# [26] "Cutibacterium"  

# cleanup
rm(tempdata, temp.core.list, temp.DAGs.list, temp.RF.list)

### (1) genera abundance table by group and sample
test.temp.ps <- microbiome::transform(temp.ps, "log10p")
tax_table(test.temp.ps)[tax_table(test.temp.ps)[, "Genus"] == "", "Genus"] <- "Unclassified"
tax_table(test.temp.ps)[, colnames(tax_table(test.temp.ps))] <- gsub(tax_table(test.temp.ps)[, colnames(tax_table(test.temp.ps))], pattern = "[a-z]__", replacement = "")
test.temp.ps <- aggregate_taxa(test.temp.ps, 'Genus')

test.composition.genus <- data.frame(otu_table(test.temp.ps)) 
test.composition.genus <- as.data.frame(t(test.composition.genus))
colnames(test.composition.genus) <- gsub(" ", "", colnames(test.composition.genus))
test.composition.genus <- test.composition.genus[, input.GLM.list]
test.composition.genus <- merge(test.composition.genus, tempsampledata, by.x=0, by.y=0)
rownames(test.composition.genus) <- test.composition.genus$Row.names
test.composition.genus <- test.composition.genus[, -c(28)]

table.glm.input <- test.composition.genus
rm(test.temp.ps, test.composition.genus)

### (2) preparing an input table for GLM analysis 
colnames(table.glm.input)[colnames(table.glm.input) == 'Row.names'] <- 'Sample'
colnames(table.glm.input)[colnames(table.glm.input) == 'Gruppe'] <- 'Group'
# cols <- c("Group", "Geschlecht", "Genese", "Nikotin", "Alkohol")
cols <- c("Group")
table.glm.input[cols] <- lapply(table.glm.input[cols], factor)  ## as.factor() could also be used
rm(cols)

### (3) Finding genera individual contributions to the differences between groups
# Aggregatibacter significantly explains differences between groups and Group:Parodontitis
res <- glm(Aggregatibacter ~ Group, data = table.glm.input, family = "gaussian")
res <- coefficients(summary(res))
res[res[,1] > 0 & res[,4] < 0.05, ]
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)       5.620556   1.041392 5.397159 1.420481e-06
# GroupParodontitis 3.045520   1.472750 2.067913 4.327941e-02
# kable(res, digits = 3, caption = "Aggregatibacter significantly explains groups (Estimate of 56.20) and for GroupParodontitis")
res.Aggregatibacter <- res[res[,1] > 0 & res[,4] < 0.05, ]
res.Aggregatibacter <- as.data.frame(res.Aggregatibacter)
# res.Aggregatibacter$taxa <- "Aggregatibacter"

# 2024-07-22 loop through all genera
for (i in 2:(ncol(table.glm.input) - 1)) {
  tmp <- as.formula(paste0(colnames(table.glm.input)[i], "~ Group"))
  res <- glm(tmp, family = "gaussian", data = table.glm.input)
  print(colnames(table.glm.input)[i])
  # print(summary(res))
  res <- coefficients(summary(res))
  print(res[res[,1] > 0 & res[,4] < 0.05, ])
}

# [1] "Aggregatibacter"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)       5.620556   1.041392 5.397159 1.420481e-06
# GroupParodontitis 3.045520   1.472750 2.067913 4.327941e-02
# [1] "Capnocytophaga"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.853753e+01 3.262232e+00 5.682468e+00 4.964180e-07 
# [1] "Porphyromonas"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)       8.632507   1.975099 4.370670 5.424089e-05
# GroupParodontitis 7.880497   2.793212 2.821303 6.605915e-03
# [1] "Selenomonas"
# Estimate Std. Error  t value    Pr(>|t|)
# (Intercept)       14.03865   3.857617 3.639203 0.000597419
# GroupParodontitis 13.45559   5.455494 2.466429 0.016733543
# [1] "Tannerella"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)        3.11152   1.044902 2.977810 4.284017e-03
# GroupParodontitis  8.45677   1.477715 5.722871 4.273252e-07
# [1] "Corynebacterium"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.298728e+01 2.004760e+00 6.478223e+00 2.508915e-08 
# [1] "Streptococcus"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.205866e+01 1.606051e+00 1.373471e+01 1.414278e-19 
# [1] "Fusobacterium"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)       15.138731   1.721845 8.792157 3.915583e-12
# GroupParodontitis  7.979207   2.435057 3.276806 1.805768e-03
# [1] "Prevotella"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 3.301060e+01 6.163934e+00 5.355444e+00 1.654642e-06 
# [1] "Rothia"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.600779e+01 1.531331e+00 1.045351e+01 8.941032e-15 
# [1] "Veillonella"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.747777e+01 2.220384e+00 7.871509e+00 1.258380e-10 
# [1] "Desulfobulbus"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.821519e+00 5.721579e-01 4.931364e+00 7.655221e-06 
# [1] "Fretibacterium"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 8.416188e+00 1.631959e+00 5.157106e+00 3.402909e-06 
# [1] "Treponema"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)        8.422188   4.110492 2.048949 0.0451638471
# GroupParodontitis 20.807339   5.813113 3.579380 0.0007200884
# [1] "Filifactor"
# Estimate  Std. Error     t value    Pr(>|t|) 
# 2.218523499 0.656867380 3.377429854 0.001336503 
# [1] "Haemophilus"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.438892e+01 1.780002e+00 8.083653e+00 5.632680e-11 
# [1] "Gemella"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 8.518180e+00 1.046548e+00 8.139308e+00 4.563327e-11 
# [1] "Dialister"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)       3.574263   1.010325 3.537734 0.0008192916
# GroupParodontitis 3.785970   1.428816 2.649726 0.0104510946
# [1] "Eikenella"
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept)       3.068997  0.6649719 4.615228 2.333993e-05
# GroupParodontitis 2.360653  0.9404123 2.510232 1.498118e-02
# [1] "Propionivibrio"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.5707045851 0.4362565106 3.6004152307 0.0006744426 
# [1] "Undibacterium"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.635414e+00 3.297325e-01 7.992583e+00 7.951906e-11 
# [1] "Lachnoanaerobaculum"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.9776675357 0.5511720074 3.5881131644 0.0007007906 
# [1] "Sphingomonas"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 6.352610e+00 7.295527e-01 8.707541e+00 5.374086e-12 
# [1] "Serratia"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.400339e+00 4.754141e-01 5.048944e+00 5.025392e-06 
# [1] "Saccharimonadales"
# Estimate Std. Error t value Pr(>|t|)
# [1] "Cutibacterium"
# Estimate Std. Error  t value   Pr(>|t|)
# (Intercept)     1.065565  0.4063394 2.622352 0.01122670
# GroupGingivitis 1.304343  0.5821628 2.240512 0.02904276



### Summary
# Explain diffs between groups: all except Saccharimonadales
# Explain GroupParodontitis: Aggregatibacter, Porphyromonas, Selenomonas, Tannerella, Fusobacterium, Treponema, Dialister, Eikenella
# Explain GroupGingivitis: Cutibacterium


### cleanup
rm(tmp, res, i)





##### plots (in progress) #####

### 2024-07-22 Exported as PDF
pdf("5.RapamycinAnalysis_20250330_V02.pdf",width = 14, height = 10)
plot.PhylumComp # phylum composition 0,1% detection 50% Prevelance
plot.GenusComp # genus composition
plot.PhylumBubble #Phylum Bubbleplot
plot.GenusBubble #Genus Bubbleplot
plot.FBratio # FB plot
plot.Core.WT # Core microbiome control_food
plot.Core.WT_RAPA # Core microbiome Rapamycin_food
plot.Core.comparison# Core microbiome comparison
#plot.BugBase # BugBase
plot.Shannon.AlphaDiversity # Shannon alpha diversity
plot.Bray.BetaDiversity # Bray divergence (beta div)
plot.unifrac.pcoa #new BEta diversity made after Franzis script (betadiv #2)
plot.unweightedunifrac.pcoa
plot.DiffAbund.WT + geom_segment(aes(xend = Genus, yend = 0)) # Diff Abd G vs H
plot.Top50.DiffAbund.WT
plot.DiffPathw.WT
dev.off()
-------------------------------------------------------------------
dev.off()
plot.RF.Parod_Gingi + 
  labs(title="Microbial biomarkers via supervised ML (RF)", subtitle="Parodontitis and Gingivitis")
plot.RF.all +
  labs(title="Microbial biomarkers via supervised ML (RF)", subtitle="all groups")

dev.off()

##### tables (in progress) #####

table.Phylum.Abundance # phylum composition
table.Genus.Abundance # genus composition
table.FBratio # FB ratios
list.Core.WT # Core microbiome Healthy
list.Core.WT_RAPA # Core microbiome Gingivitis
table.Core.comparison.matrix # Core microbiome comparison
table.Core.comparison.matrix_for_ggplot2 # Core microbiome comparison f ggplot2
#table.bugbase # BugBase
table.shannon # Shannon alpha diversity
table.betadiv # Bray divergence (beta div)
table.DiffAbund.WT# Diff abundance Gingivitis vs Healthy
table.DiffPathw.WT
---------------------------------------------------------------------
table.DiffPathw.Gingi # Diff pathways Gingivitis vs Healthy (EMPTY TABLE)
table.DiffPathw.Parod # Diff pathways Parod vs Healthy
table.DiffPathw.ParodGingi # Diff pathways Parod vs Gingivitis

table.ClinicalMarkers # clinical markers; use for publication
table.ClinicalMarkers.simplified # clinical markers; used for analyses
table.ClinicalMarkers.means # means by group of 'table.ClinicalMarkers.simplified'

table.correlation.Gesnd # correlations genera vs healthy (no signif correlations)
table.correlation.Gingi # correlations genera vs gingivitis 
table.correlation.Perod # correlations genera vs perodontitis (no signif correlations)

table.RF.Parod_Gingi # RF classifiers (P vs Gingi)
table.RF.all # RF classifiers (all groups)

##### excel output (ok) #####

# new excel object
wb <- createWorkbook("dentalmicrobiome_20220722_01.xlsx")

# copy table list then cmd line:
# $ pbpaste | gsed 's/\s\#.*//' | gawk '{gsub ("table.", ""); print "addWorksheet(wb, \x27"$1"\x27)"}' | pbcopy
addWorksheet(wb, 'Phylum.Abundance')
addWorksheet(wb, 'Genus.Abundance')
addWorksheet(wb, 'FBratio')
# addWorksheet(wb, 'list.Core.Gesund')
# addWorksheet(wb, 'list.Core.Gingi')
# addWorksheet(wb, 'list.Core.Parod')
# addWorksheet(wb, 'Core.comparison.matrix')
# addWorksheet(wb, 'Core.comparison.matrix_for_ggplot2')
addWorksheet(wb, 'bugbase')
addWorksheet(wb, 'shannon')
addWorksheet(wb, 'betadiv')
addWorksheet(wb, 'DiffAbund.Gingivitis')
addWorksheet(wb, 'DiffAbund.Parodontitis')
addWorksheet(wb, 'DiffAbund.ParodGingi')
addWorksheet(wb, 'DiffPathw.Gingi')
addWorksheet(wb, 'DiffPathw.Parod')
addWorksheet(wb, 'DiffPathw.ParodGingi')
# addWorksheet(wb, 'ClinicalMarkers')
addWorksheet(wb, 'ClinicalMarkers.simplified')
# addWorksheet(wb, 'ClinicalMarkers.means')
addWorksheet(wb, 'correlation.Gesnd')
addWorksheet(wb, 'correlation.Gingi')
addWorksheet(wb, 'correlation.Perod')
addWorksheet(wb, 'RF.Parod_Gingi')

# $ pbpaste | gsed 's/\s\#.*//' | gawk '{gsub ("table.", ""); print "writeData(wb, \x27"$1"\x27, "$1", rowNames = FALSE)"}' | pbcopy
writeData(wb, 'Phylum.Abundance', table.Phylum.Abundance, rowNames = FALSE)
writeData(wb, 'Genus.Abundance', table.Genus.Abundance, rowNames = FALSE)
writeData(wb, 'FBratio', table.FBratio, rowNames = FALSE)
# writeData(wb, 'list.Core.Gesund', list.Core.Gesund, rowNames = FALSE)
# writeData(wb, 'list.Core.Gingi', list.Core.Gingi, rowNames = FALSE)
# writeData(wb, 'list.Core.Parod', list.Core.Parod, rowNames = FALSE)
# writeData(wb, 'Core.comparison.matrix', table.Core.comparison.matrix, rowNames = FALSE)
# writeData(wb, 'Core.comparison.matrix_for_ggplot2', table.Core.comparison.matrix_for_ggplot2, rowNames = FALSE)
writeData(wb, 'bugbase', table.bugbase, rowNames = FALSE)
writeData(wb, 'shannon', table.shannon, rowNames = FALSE)
writeData(wb, 'betadiv', table.betadiv, rowNames = FALSE)
writeData(wb, 'DiffAbund.Gingivitis', table.DiffAbund.Gingivitis, rowNames = FALSE)
writeData(wb, 'DiffAbund.Parodontitis', table.DiffAbund.Parodontitis, rowNames = FALSE)
writeData(wb, 'DiffAbund.ParodGingi', table.DiffAbund.ParodGingi, rowNames = FALSE)
writeData(wb, 'DiffPathw.Gingi', table.DiffPathw.Gingi, rowNames = FALSE)
writeData(wb, 'DiffPathw.Parod', table.DiffPathw.Parod, rowNames = FALSE)
writeData(wb, 'DiffPathw.ParodGingi', table.DiffPathw.ParodGingi, rowNames = FALSE)
# writeData(wb, 'ClinicalMarkers', table.ClinicalMarkers, rowNames = FALSE)
writeData(wb, 'ClinicalMarkers.simplified', table.ClinicalMarkers.simplified, rowNames = FALSE)
# writeData(wb, 'ClinicalMarkers.means', table.ClinicalMarkers.means, rowNames = FALSE)
writeData(wb, 'correlation.Gesnd', table.correlation.Gesnd, rowNames = FALSE)
writeData(wb, 'correlation.Gingi', table.correlation.Gingi, rowNames = FALSE)
writeData(wb, 'correlation.Perod', table.correlation.Perod, rowNames = FALSE)
writeData(wb, 'RF.Parod_Gingi', table.RF.Parod_Gingi, rowNames = FALSE)

# saving excel
saveWorkbook(wb, file = "dentalmicrobiome_20220722_01.xlsx", overwrite = TRUE)

# cleanup
rm(wb)

##### image #####
# Save the current R workspace
save.image("Rapamycin_reanalysis_5.PRJNA717340_V04_20250502_latest.RData")

# Load the saved workspace for later use
load("Rapamycin_reanalysis_5.PRJNA717340_V04_20250502_latest.RData")

# Print the session information
sessionInfo()

# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.2
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] knitr_1.49                  randomForest_4.7-1.2        data.table_1.16.4          
# [4] rstatix_0.7.2               jsonlite_1.8.9              lubridate_1.9.4            
# [7] forcats_1.0.0               purrr_1.0.2                 readr_2.1.5                
# [10] tibble_3.2.1                tidyverse_2.0.0             splitstackshape_1.4.8      
# [13] viridis_0.6.5               viridisLite_0.4.2           microbiome_1.28.0          
# [16] openxlsx_4.2.7.1            patchwork_1.3.0             ggpubr_0.6.0               
# [19] ggplot2_3.5.1               DESeq2_1.46.0               SummarizedExperiment_1.36.0
# [22] Biobase_2.66.0              MatrixGenerics_1.18.1       matrixStats_1.5.0          
# [25] GenomicRanges_1.58.0        GenomeInfoDb_1.42.1         IRanges_2.40.1             
# [28] S4Vectors_0.44.0            BiocGenerics_0.52.0         ape_5.8-1                  
# [31] vegan_2.6-8                 lattice_0.22-6              permute_0.9-7              
# [34] phyloseq_1.50.0             tidyr_1.3.1                 reshape2_1.4.4             
# [37] stringr_1.5.1               dplyr_1.1.4                
# 
# loaded via a namespace (and not attached):
#   [1] gridExtra_2.3           rlang_1.1.4             magrittr_2.0.3          ade4_1.7-22            
# [5] compiler_4.4.2          mgcv_1.9-1              vctrs_0.6.5             pkgconfig_2.0.3        
# [9] crayon_1.5.3            backports_1.5.0         XVector_0.46.0          labeling_0.4.3         
# [13] utf8_1.2.4              tzdb_0.4.0              UCSC.utils_1.2.0        xfun_0.50              
# [17] zlibbioc_1.52.0         biomformat_1.34.0       rhdf5filters_1.18.0     DelayedArray_0.32.0    
# [21] Rhdf5lib_1.28.0         BiocParallel_1.40.0     broom_1.0.7             parallel_4.4.2         
# [25] cluster_2.1.8           R6_2.5.1                RColorBrewer_1.1-3      stringi_1.8.4          
# [29] car_3.1-3               Rcpp_1.0.14             iterators_1.0.14        timechange_0.3.0       
# [33] Matrix_1.7-1            splines_4.4.2           igraph_2.1.2            tidyselect_1.2.1       
# [37] rstudioapi_0.17.1       abind_1.4-8             codetools_0.2-20        plyr_1.8.9             
# [41] withr_3.0.2             evaluate_1.0.3          Rtsne_0.17              survival_3.8-3         
# [45] zip_2.3.1               Biostrings_2.74.1       pillar_1.10.1           carData_3.0-5          
# [49] foreach_1.5.2           generics_0.1.3          hms_1.1.3               munsell_0.5.1          
# [53] scales_1.3.0            glue_1.8.0              tools_4.4.2             locfit_1.5-9.10        
# [57] ggsignif_0.6.4          rhdf5_2.50.2            grid_4.4.2              colorspace_2.1-1       
# [61] nlme_3.1-166            GenomeInfoDbData_1.2.13 Formula_1.2-5           cli_3.6.3              
# [65] S4Arrays_1.6.0          gtable_0.3.6            ggsci_3.2.0             digest_0.6.37          
# [69] SparseArray_1.6.0       farver_2.1.2            multtest_2.62.0         lifecycle_1.0.4        
# [73] httr_1.4.7              MASS_7.3-64   