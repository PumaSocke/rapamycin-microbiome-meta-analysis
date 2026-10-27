### 2024-12-06 Rapamycin reanalysis 2024

### 2024-12-06 to reload
# run chunk 1
# load image
# > load("dientes_reanalysis_V01_20240620_latest.RData") # have to change to new rdata file
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
setwd(here::here("data/study2_PRJEB13256"))

# load functions
source(here::here("scripts/MicrobiomeFunctions.R"))

#### Load Data #### 
### input file paths

inputfile.metadata <- "sample_metadata.tsv"
inputfile.otu <- "feature-table.tsv"
inputfile.taxonomy <- "taxonomy.tsv"
inputfile.tree <- "tree.nwk"

# comparisons (for ggpubr)
my_comparisons <- list( c("Rapamycin_IP", "Vehicle_IP"),
                        c("Rapamycin_food", "Control_food"))
                       #c("Gingivitis", "Parodontitis") )

##### 16S phyloseq object (ok) #####
file_path <- "sample_metadata.tsv"
metadata <- read.delim(file_path, header = TRUE, sep = "\t")

# Step 1: Switch headers of 'application_method' and 'sex' columns
colnames(metadata)[colnames(metadata) == "application_method"] <- "temp"
colnames(metadata)[colnames(metadata) == "sex"] <- "application_method"
colnames(metadata)[colnames(metadata) == "temp"] <- "sex"

shifted_data <- rbind(names(metadata), metadata) # 2. Insert the original column names as the first row

shifted_data <- cbind(NA, shifted_data) # 3. Add an empty column on the left, but immediately rename it to #SampleID
colnames(shifted_data) <- c("#SampleID", colnames(metadata))


shifted_data <- shifted_data[-1, ] # 4. Remove the row that contained the original header names (now duplicated)
shifted_data$`#SampleID` <- shifted_data$run_accession # 5. Fill the new #SampleID column with the values from 'run_accession'

# 6. Shift all columns one row to the right
shifted_data <- rbind(names(metadata), metadata)  # Add the header as the first row
shifted_data <- cbind(NA, shifted_data)           # Add an empty column to the left
names(shifted_data) <- names(metadata)            # Restore the original column names

shifted_data <- shifted_data[-1, ] # 7. Remove the temporary first row (original header)

shifted_data$`SampleID` <- shifted_data$run_accession # 8. Copy 'run_accession' to '#SampleID' (first column)

shifted_data <- shifted_data[, -ncol(shifted_data)] #9.Remove the last column, no matter what it is

rownames(shifted_data) <- NULL #10.Clear out the old row names so R will display numeric 1,2,3,... 

colnames(shifted_data)

shifted_data <- shifted_data %>%
  mutate(
    application_method = if_else(
      SampleID == "ERR1428190",      # ←  the sample you need to fix
      "Control_food",                # ←  new value
      application_method             # ←  leave all others unchanged
    )
  )

shifted_data %>%
  filter(SampleID == "ERR1428190") %>%
  select(SampleID, application_method, sex)

# 11. Write out the final shifted metadata
write.table(shifted_data, file = "shifted_sample_metadata.tsv",sep = "\t",row.names = FALSE,quote = FALSE)
 head(shifted_data)
 
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
# otu_table()   OTU Table:         [ 4180 taxa and 36 samples ]
# sample_data() Sample Data:       [ 36 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 4180 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4180 tips and 4179 internal nodes ]


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
#[1] "ERR1428194" "ERR1428197"
temp.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4180 taxa and 34 samples ]
# sample_data() Sample Data:       [ 34 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 4180 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4180 tips and 4179 internal nodes ]

### outlier  # In phyloseq:
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
# 
# head(sample_comp, 20)

summary_results <- microbiome::summarize_phyloseq(temp.ps)
print(summary_results)

# > summary_results <- microbiome::summarize_phyloseq(temp.ps)
# Compositional = NO2
# 1] Min. number of reads = 79802] Max. number of reads = 327933] Total number of reads = 8272844] Average number of reads = 24331.88235294125] Median number of reads = 24996.57] Sparsity = 0.9278637770897836] Any OTU sum to 1 or less? YES8] Number of singletons = 229] Percent of OTUs that are singletons 
# (i.e. exactly one read detected across all samples)010] Number of sample variables are: 18SampleIDstudy_accessionsample_accessionexperiment_accessionrun_accessiontax_idscientific_namelibrary_namenominal_lengthlibrary_layoutread_countbase_countexperiment_titlefrpt_typesample_aliassample_titlesexapplication_method2
# > print(summary_results)
# [[1]]
# [1] "1] Min. number of reads = 7980"
# 
# [[2]]
# [1] "2] Max. number of reads = 32793"
# 
# [[3]]
# [1] "3] Total number of reads = 827284"
# 
# [[4]]
# [1] "4] Average number of reads = 24331.8823529412"
# 
# [[5]]
# [1] "5] Median number of reads = 24996.5"
# 
# [[6]]
# [1] "7] Sparsity = 0.927863777089783"
# 
# [[7]]
# [1] "6] Any OTU sum to 1 or less? YES"
# 
# [[8]]
# [1] "8] Number of singletons = 22"
# 
# [[9]]
# [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)0"
# 
# [[10]]
# [1] "10] Number of sample variables are: 18"
# 
# [[11]]
# [1] "SampleID"             "study_accession"      "sample_accession"     "experiment_accession"
# [5] "run_accession"        "tax_id"               "scientific_name"      "library_name"        
# [9] "nominal_length"       "library_layout"       "read_count"           "base_count"          
# [13] "experiment_title"     "frpt_type"            "sample_alias"         "sample_title"        
# [17] "sex"                  "application_method"  
# 

####Save Phyloseq object externally 
## What phyloseq objects are in memory?
ls()                     # lists all objects
class(temp.ps)           # should print "phyloseq" if it exists


# 1.  Make sure the clean object is called `clean_ps`
clean_ps <- temp.ps            # or qiime2.ps / ps_final …

# 2.  Tag samples with the study ID  ← change to your study label
sample_data(clean_ps)$StudyID <- "PRJEB13256"

saveRDS(
  clean_ps,
  file = here::here("data/processed/ps_PRJEB13256_clean.rds"),
  compress = "xz"
)

message("✓ Saved cleaned phyloseq for PRJEB13256")



# cleanup
rm(OTU, TAX, shifted_data, taxonomy, taxmat)

####----------- Phylum-level Extraction -----------####
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
            file = here::here("results/study2_phyloseq_phylum_table.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#####----------- Genus-level Extraction -----------#####
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
            file = here::here("results/study2_phyloseq_genus_table.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

---------------------------------------------


##### Phylum/Genera compositions (ok) #####

set.seed(1)

ps0 <- temp.ps
stopifnot(!is.null(phy_tree(ps0)))  # ensure a tree exists

## Derive route & treatment from application_method
md <- data.frame(sample_data(ps0))
route <- ifelse(grepl("_IP$",   md$application_method), "IP",
                ifelse(grepl("_food$", md$application_method, ignore.case = TRUE), "Food", NA))
treat <- ifelse(grepl("^Rapamycin", md$application_method, ignore.case = TRUE), "Rapamycin", "Control")
sample_data(ps0)$route     <- route
sample_data(ps0)$treatment <- treat

## Helper to include only covariates that exist
rhs_builder <- function(base, covars, ps){
  present <- covars[covars %in% sample_variables(ps)]
  paste(base, paste(present, collapse = " + "), sep = if (length(present)) " + " else "")
}

## Runner (weighted & unweighted UniFrac)
run_perm <- function(ps, rhs){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)  # filter AFTER subsetting
  Dw  <- UniFrac(ps, weighted = TRUE,  normalized = TRUE)
  Duw <- UniFrac(ps, weighted = FALSE, normalized = TRUE)
  md  <- data.frame(sample_data(ps))
  list(
    wunifrac = adonis2(Dw  ~ ., data = md[, all.vars(as.formula(paste("~", rhs)))],
                       permutations = 9999, by = "margin"),
    uunifrac = adonis2(Duw ~ ., data = md[, all.vars(as.formula(paste("~", rhs)))],
                       permutations = 9999, by = "margin")
  )
}

## 1) Injection only: treatment (+ sex if present)
ps_ip   <- subset_samples(ps0, route == "IP")
rhs_ip  <- rhs_builder("treatment", c("sex"), ps_ip)
perm_ip <- run_perm(ps_ip, rhs_ip)

## 2) Feeding only: treatment (+ sex if present)
ps_food   <- subset_samples(ps0, route == "Food")
rhs_food  <- rhs_builder("treatment", c("sex"), ps_food)
perm_food <- run_perm(ps_food, rhs_food)

## 3) Pooled: route + treatment (+ sex if present)
rhs_all  <- rhs_builder("route + treatment", c("sex"), ps0)
perm_all <- run_perm(ps0, rhs_all)

## ---- Summarise key terms (no cat) 
grab <- function(tab, term){
  rn <- rownames(tab); i <- match(term, rn)
  stopifnot(!is.na(i)); data.frame(R2 = tab$R2[i], p = tab[["Pr(>F)"]][i])
}

summ <- rbind(
  cbind(contrast="Injection", metric="Unweighted UniFrac", term="treatment", grab(perm_ip$uunifrac, "treatment")),
  cbind(contrast="Injection", metric="Weighted UniFrac",   term="treatment", grab(perm_ip$wunifrac, "treatment")),
  cbind(contrast="Feeding",   metric="Unweighted UniFrac", term="treatment", grab(perm_food$uunifrac, "treatment")),
  cbind(contrast="Feeding",   metric="Weighted UniFrac",   term="treatment", grab(perm_food$wunifrac, "treatment")),
  cbind(contrast="Pooled",    metric="Unweighted UniFrac", term="route",     grab(perm_all$uunifrac, "route")),
  cbind(contrast="Pooled",    metric="Unweighted UniFrac", term="treatment", grab(perm_all$uunifrac, "treatment")),
  cbind(contrast="Pooled",    metric="Weighted UniFrac",   term="route",     grab(perm_all$wunifrac, "route")),
  cbind(contrast="Pooled",    metric="Weighted UniFrac",   term="treatment", grab(perm_all$wunifrac, "treatment"))
)
print(summ, row.names = FALSE, digits = 3)

## ---- (Optional) dispersion checks to report alongside PERMANOVA 
# anova(betadisper(UniFrac(ps_ip,   weighted=TRUE,  normalized=TRUE), sample_data(ps_ip)$treatment))
# anova(betadisper(UniFrac(ps_food, weighted=TRUE,  normalized=TRUE), sample_data(ps_food)$treatment))
# anova(betadisper(UniFrac(ps0,     weighted=TRUE,  normalized=TRUE), sample_data(ps0)$treatment))


### levels for ordered plots
grouporder <- c("Vehicle_IP", "Rapamycin_IP","Control_food", "Rapamycin_food")

### Phylum by comparison group
pseq <- microbiome::transform(temp.ps, "compositional")
pseq <- aggregate_rare(pseq, level = "Phylum", detection = 1/1000, prevalence = 50/100)
mysubtitle <- "Detection 0.1%, Prevalence 50%, Kruskal-Wallis p-value "
tempplot <- plot_composition(pseq, average_by = "application_method") +
  scale_fill_brewer("Phyla", palette = "Paired") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) +
  labs(x = "Group", 
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
bubbleTable <- bubble_table[, c("Sample","Phylum", "Abundance", "application_method")] # columns: Sample, Phylum, Abundance

# Remove rows with Abundance == 0 or NA
bubbleTable <- bubbleTable %>%
  filter(!is.na(Abundance) & Abundance > 0)
# Reorder factor levels here:
bubbleTable$application_method <- factor(
  bubbleTable$application_method,
  levels = c("Vehicle_IP", "Rapamycin_IP", "Control_food", "Rapamycin_food")
 )
myBubblePlot <- ggplot(bubbleTable, aes(x = application_method, y = Phylum)) + 
  geom_point(aes(color = Phylum, size = Abundance), alpha = 1) +
  xlab("Application method") +
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
tax_table(pseq)[tax_table(pseq)[, "Genus"] == "", "Genus"] <- " Unclassified"

# If you already have a category called "Unclassified" without a space, rename it too
tax_table(pseq)[tax_table(pseq)[, "Genus"] == "Unclassified", "Genus"] <- " Unclassified"

pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/1000, prevalence = 50/100)
mysubtitle <- "Detection 1.25%, Prevalence 50%, Kruskal-Wallis p-value "
tempplot <- plot_composition(pseq, average_by = "application_method") +
  scale_fill_viridis_d("Genera") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) +
  labs(x = "Group", 
       y = "Relative abundance",
       title = "Genus composition") 

# Now directly modify the plot data
tempplot$data$Tax[tempplot$data$Tax %in% c("uncultured"," uncultured", "Other", "Unclassified")] <- " Unclassified"

tempplot$data$Sample <- as.character(tempplot$data$Sample)
tempplot$data$Sample <- factor(x = tempplot$data$Sample, levels = grouporder)
ps.kruskal <- kruskal.test(Abundance ~ Sample, data=tempplot$data)
tempplot <- tempplot + labs(subtitle = paste(mysubtitle, round(ps.kruskal$p.value, digits = 4)))
plot.GenusComp <- tempplot

print(plot.GenusComp)

#phyloseq object to relative abundances (Genus level)

bubble_table <- psmelt(pseq) %>%
  arrange(Genus)

bubbleTable <- bubble_table[, c("Sample","Genus", "Abundance", "application_method")] # columns: Sample, Phylum, Abundance
# Remove rows with Abundance == 0 or NA
bubbleTable <- bubbleTable %>%
  filter(!is.na(Abundance) & Abundance > 0)

bubbleTable <- bubbleTable %>%
  mutate(
    Genus = str_trim(Genus),                          # trim hidden spaces
    Genus = case_when(                                # collapse variants
      Genus %in% c("uncultured", "Uncultured",
                   "Other", "Unclassified") ~ "Unclassified",
      TRUE                                  ~ Genus
    )
  ) %>%
  droplevels()                                         # clean out unused levels

# Reorder factor levels here:
bubbleTable$application_method <- factor(
  bubbleTable$application_method,
  levels = c("Vehicle_IP", "Rapamycin_IP", "Control_food", "Rapamycin_food")
)

myBubblePlot <- ggplot(bubbleTable, aes(x = application_method, y = Genus)) + 
  geom_point(aes(color = Genus, size = Abundance), alpha = 1) +
  xlab("Application method") +
  scale_size_continuous(limits = c(0.001, 0.6), range = c(1,25)) + 
  theme_bw() + # remove background color
  theme(axis.text.x = element_text(angle= 45, hjust= 1)) +
  labs(x = "Group", 
       y = "Relative abundance",
       title = "Genus Bubble Plot"
       ) 

  myBubblePlot 
plot.GenusBubble <- myBubblePlot

# to retrieve plot data
# tempdata <- tempplot$data


### cleanup
rm(pseq, tempplot, myBubblePlot, bubbleTable, bubble_table)

###Combine

# Phylum plot - no gap between bars and axis
plot.PhylumComp <- plot.PhylumComp + 
  scale_x_discrete(labels = c("Vehicle_IP" = "Control-IP", 
                              "Rapamycin_IP" = "Rapamycin-IP",
                              "Control_food" = "Control-Diet",
                              "Rapamycin_food" = "Rapamycin-Diet"),
                   expand = c(0, 0)) +  # Remove gap
  scale_y_continuous(expand = c(0, 0)) +  # Remove gap on y-axis too
  labs(tag = "a", title = NULL, subtitle = NULL, x = NULL) +
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Genus plot - no gap between bars and axis 
plot.GenusComp <- plot.GenusComp + 
  scale_x_discrete(labels = c("Vehicle_IP" = "Control-IP", 
                              "Rapamycin_IP" = "Rapamycin-IP",
                              "Control_food" = "Control-Diet",
                              "Rapamycin_food" = "Rapamycin-Diet"),
                   expand = c(0, 0)) +  # Remove gap
  scale_y_continuous(expand = c(0, 0)) +  # Remove gap on y-axis too
  labs(tag = "b", title = NULL, subtitle = NULL, x = NULL) +
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine side by side
library(patchwork)
plot.Composition.combined <- plot.PhylumComp | plot.GenusComp
plot.Composition.combined

### Abundance tables for publication

## Phylum table (percent, no double-aggregation)
ph_dat <- plot.PhylumComp$data[, c("Tax", "Sample", "Abundance")]
ph_wide <- reshape2::dcast(ph_dat, Tax ~ Sample, value.var = "Abundance")
ph_wide[,-1] <- round(100 * ph_wide[,-1], 2)
table.Phylum.Abundance <- ph_wide

# Collapse duplicates per Tax × Sample (use mean; identical rows will be unchanged)
gn_dat <- plot.GenusComp$data[, c("Tax", "Sample", "Abundance", "xlabel")]
gn_dat <- aggregate(Abundance ~ Tax + Sample, data = gn_dat, FUN = mean)
gn_wide <- reshape2::dcast(gn_dat, Tax ~ Sample, value.var = "Abundance")
gn_wide[,-1] <- round(100 * gn_wide[,-1], 4)
table.Genus.Abundance <- gn_wide

# Aggregate at family level to look for Rikenellaceae and Helicobacteraceae
pseq_family <- microbiome::transform(temp.ps, "compositional")
pseq_family <- aggregate_rare(pseq_family, level = "Family", detection = 1/100000, prevalence = 1/100)

# Get abundance table at family level
family_abundance <- abundances(pseq_family)
family_metadata <- meta(pseq_family)

# Create family abundance table
family_table <- data.frame(
  Tax = rownames(family_abundance),
  Vehicle_IP = rowMeans(family_abundance[, family_metadata$application_method == "Vehicle_IP"]),
  Rapamycin_IP = rowMeans(family_abundance[, family_metadata$application_method == "Rapamycin_IP"]),
  Control_food = rowMeans(family_abundance[, family_metadata$application_method == "Control_food"]),
  Rapamycin_food = rowMeans(family_abundance[, family_metadata$application_method == "Rapamycin_food"])
)

# Look specifically for the two missing families
family_table[grep("Rikenellaceae|Helicobacteraceae", family_table$Tax), ]
tempdata <- plot.PhylumComp$data
tempdata$xlabel <- NULL
table.Phylum.Abundance <- reshape2::dcast(tempdata, Tax ~ Sample, value.var = "Abundance", fun.aggregate = sum)
table.Phylum.Abundance$Tax <- gsub(" ", "", table.Phylum.Abundance$Tax)

tempdata <- plot.GenusComp$data
tempdata$xlabel <- NULL
table.Genus.Abundance <-reshape2::dcast(tempdata, Tax ~ Sample, value.var = "Abundance", fun.aggregate = sum) 
table.Genus.Abundance$Tax <- gsub(" ", "", table.Genus.Abundance$Tax)

rm(tempdata)


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
tempdata <- merge(tempdata, tempsampledata2, by.x=1, by.y=1)
tempdata$logFB <- log(tempdata$FBratio) # using logs because otherwise range way too big in plot 

plot.FBratio <- ggboxplot(tempdata, x = "application_method", y = "logFB", color = "application_method", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons) + 
  #ylim(-1, 10) + # Add pairwise comparisons p-value
  ggtitle("Firmicutes - Bacteroidetes ratios") +
  xlab("Group") +
  ylab("log(FB ratio)") +
  labs(color = "Group")
  

plot.FBratio
table.FBratio <- tempdata

# For the IP comparison
wilcox.test(FBratio ~ application_method, 
            data = subset(tempdata, application_method %in% c("Rapamycin_IP", "Vehicle_IP")))

# For the food comparison
wilcox.test(FBratio ~ application_method, 
            data = subset(tempdata, application_method %in% c("Rapamycin_food", "Control_food")))

### mean values
mean(table.FBratio[table.FBratio$application_method == 'Rapamycin_IP', 'FBratio'])
mean(table.FBratio[table.FBratio$application_method == 'Vehicle_IP', 'FBratio'])
mean(table.FBratio[table.FBratio$application_method == 'Rapamycin_food', 'FBratio'])
mean(table.FBratio[table.FBratio$application_method == 'Control_food', 'FBratio'])
# > mean(table.FBratio[table.FBratio$application_method == 'Rapamycin_IP', 'FBratio'])
# [1] 0.6725011
# > mean(table.FBratio[table.FBratio$application_method == 'Vehicle_IP', 'FBratio'])
# [1] 0.8945744
# > mean(table.FBratio[table.FBratio$application_method == 'Rapamycin_food', 'FBratio'])
# [1] 0.8281735
# > mean(table.FBratio[table.FBratio$application_method == 'Control_food', 'FBratio'])
# [1] 1.354303
# >

rm(pseq, tempplot, temp.list, tempdata)


##### core microbiomes (ok)  #####

### dataset cleanup
# Create dataset without low coverage samples (<150 reads)
psHighCov <- prune_samples(sample_sums(temp.ps) >= 150, temp.ps)

# Prune unidentified Genuses; 
# Filter out Eukaryotes, chloroplasts and mitochondria
# ref http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
psHighCov <-  subset_taxa(psHighCov, Genus != "NA")
psHighCov <-  subset_taxa(psHighCov, Genus != "Unassigned")
psHighCov <-  subset_taxa(psHighCov, Genus != "unidentified")
psHighCov <-  subset_taxa(psHighCov, Genus != "-1")
psHighCov <-  subset_taxa(psHighCov, Genus != "Mitochondria")
psHighCov <-  subset_taxa(psHighCov, Genus != "Chloroplast")


### Control food core
pseq <- subset_samples(psHighCov, application_method == "Control_food") %>%
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
  labs(title = "Core microbiome: Control food", 
       subtitle = "detection 0.1%, prevalence 80%, min. prevalence 0.3")
#myplot

list.Core.Controlfood <- core_taxa_standard
plot.Core.Controlfood <- myplot
plot.Core.Controlfood
length(taxa(pseq))
# [1] 147
length(list.Core.Controlfood)
# [1] 23
# therefore:
# "22 out of 147 genera are shared by 80% of the subjects
# at the minimum detection threshold of 0.1% relative abundance"
# Shetty et al. 2017


### Rapamycin_food core
pseq <- subset_samples(psHighCov, application_method == "Rapamycin_food") %>%
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
  labs(title = "Core microbiome: Rapamycin food", 
       subtitle = "detection 0.1%, prevalence 80%, min. prevalence 0.3")
# myplot

list.Core.Rapamycinfood <- core_taxa_standard
plot.Core.Rapamycinfood <- myplot
length(taxa(pseq))
# [1] 147
 length(list.Core.Rapamycinfood)
# [1] 20
# therefore:
# "25 out of 145 genera are shared by 80% of the subjects
# at the minimum detection threshold of 0.1% relative abundance"
# Shetty et al. 2017


### VehicleIP core
pseq <- subset_samples(psHighCov, application_method == "Vehicle_IP") %>%
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
  labs(title = "Core microbiome: Vehicle IP", 
       subtitle = "detection 0.1%, prevalence 80%, min. prevalence 0.3")
# myplot

list.Core.VehicleIP <- core_taxa_standard
plot.Core.VehicleIP <- myplot

#length(taxa(pseq))
# [1] 147
#length(list.Core.VehicleIP)
# [1] 24
# therefore:
# "24 out of 145 genera are shared by 80% of the subjects
# at the minimum detection threshold of 0.1% relative abundance"
# Shetty et al. 2017


### Rapamycin IP core
pseq <- subset_samples(psHighCov, application_method == "Rapamycin_IP") %>%
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
  labs(title = "Core microbiome: Rapamycin IP", 
       subtitle = "detection 0.1%, prevalence 80%, min. prevalence 0.3")
# myplot

list.Core.RapamycinIP <- core_taxa_standard
plot.Core.RapamycinIP <- myplot

length(taxa(pseq))
# [1] 147
length(list.Core.RapamycinIP)
# [1] 26
# therefore:
# "26 out of 145 genera are shared by 80% of the subjects
# at the minimum detection threshold of 0.1% relative abundance"
# Shetty et al. 2017

### cleanup
rm(myplot, core_taxa_standard, prevalences, detections, gray, pseq)


# Create plots WITHOUT legends for a and c
plot.Core.Control <- plot.Core.VehicleIP + 
  scale_x_discrete(labels = scaleFUN) +
  labs(title = "Control-IP", subtitle = NULL,  # JUST ADD TITLE
       x = "Detection threshold (Relative abundance (%))",
       tag = "a") +
  theme(plot.tag.position = c(0.1, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        plot.title = element_text(size = 12, hjust = 0),  # Left-align title
        legend.position = "none")

plot.Core.Rapamycin <- plot.Core.RapamycinIP + 
  scale_x_discrete(labels = scaleFUN) +
  labs(title = "Rapamycin-IP", subtitle = NULL,  # JUST ADD TITLE
       x = "Detection threshold (Relative abundance (%))",
       tag = "b") +
  theme(plot.tag.position = c(0.05, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        plot.title = element_text(size = 12, hjust = 0))  # Left-align title

plot.Core.ControlDiet <- plot.Core.Controlfood + 
  scale_x_discrete(labels = scaleFUN) +
  labs(title = "Control-Diet", subtitle = NULL,  # JUST ADD TITLE
       x = "Detection threshold (Relative abundance (%))",
       tag = "c") +
  theme(plot.tag.position = c(0.1, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        plot.title = element_text(size = 12, hjust = 0),  # Left-align title
        legend.position = "none")

plot.Core.RapamycinDiet <- plot.Core.Rapamycinfood + 
  scale_x_discrete(labels = scaleFUN) +
  labs(title = "Rapamycin-Diet", subtitle = NULL,  # JUST ADD TITLE
       x = "Detection threshold (Relative abundance (%))",
       tag = "d") +
  theme(plot.tag.position = c(0.05, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        plot.title = element_text(size = 12, hjust = 0))  # Left-align title

# Combine
plot.Core.combined <- (plot.Core.Control | plot.Core.Rapamycin) / 
  (plot.Core.ControlDiet | plot.Core.RapamycinDiet)

plot.Core.combined

# Save
ggsave(here::here("results/fig_s2.png"), 
       plot.Core.combined,
       width = 10, height = 10, dpi = 300)

##### comparison of core microbiomes (ok) #####

### core vectors to presence/absence matrix
# src https://stackoverflow.com/a/53575310
ls(pattern="list.*Core.*")
# [1] "list.Core.Controlfood"   "list.Core.Rapamycinfood" "list.Core.RapamycinIP"  
# [4] "list.Core.VehicleIP"    

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

colnames(tempdata) <- c("Control_food", "Rapamycin_food", "Rapamycin_IP", "Vehicle_IP")
tempdata2 <- tempdata %>% 
  as.data.frame() %>%
  rownames_to_column("Genus") %>% # requires library(tidyverse)
  pivot_longer(-c(Genus), names_to = "application_method", values_to = "value")
tempdata2$application_method <- factor(tempdata2$application_method, levels=c("Control_food", "Rapamycin_food", "Vehicle_IP", "Rapamycin_IP"))

tempplot <- ggplot(tempdata2, aes(Genus, application_method)) +
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
plot.Core.comparison <- tempplot

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
season_std_shannon <- plot_richness(temp.ps, x = "application_method", color = "application_method", measures = c("Shannon")) 
tempdata <- season_std_shannon$data
names(tempdata)[names(tempdata) == 'application_method'] <- 'Group'

tempdata$Group <- factor(tempdata$Group, levels = c("Control_food", "Rapamycin_food", "Vehicle_IP", "Rapamycin_IP"))

plot.Shannon.AlphaDiversity <- ggboxplot(tempdata, x = "Group", y = "value",
                                         color = "Group", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  # 2) Overall test label (Kruskal-Wallis), placed at x=1
  stat_compare_means(
    method = "kruskal.test",
    label.x = 1,
    label.y = (1.2 + max(tempdata$value)))+
  labs(
    title = "Shannon diversity",
    x     = "Group",
    y     = "Shannon diversity index",
    color = NULL            # drops the legend title entirely
  ) +
  theme(legend.position = "none")+
  ggtitle("Shannon diversity") +
  xlab("Group") +
  ylab("Shannon diversity index")+
  coord_cartesian(ylim = c(4, 7))


plot.Shannon.AlphaDiversity
# backup table
table.shannon <- tempdata

# cleanup
rm(tempdata)

##combine

# Update F/B ratio plot - same labels as Shannon, correct y-axis
plot.FBratio <- plot.FBratio + 
  scale_x_discrete(labels = c("Control_food" = "Control-Diet",
                              "Rapamycin_food" = "Rapamycin-Diet",
                              "Vehicle_IP" = "Control-IP", 
                              "Rapamycin_IP" = "Rapamycin-IP")) +
  labs(tag = "a", title = NULL, subtitle = NULL, x = NULL, y = "log(F/B ratio)") +  # Changed y-axis
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        legend.position = "none",
        axis.text.x = element_text())

# Update Shannon diversity plot
plot.Shannon.AlphaDiversity <- plot.Shannon.AlphaDiversity + 
  scale_x_discrete(labels = c("Control_food" = "Control-Diet",
                              "Rapamycin_food" = "Rapamycin-Diet",
                              "Vehicle_IP" = "Control-IP", 
                              "Rapamycin_IP" = "Rapamycin-IP")) +
  labs(tag = "b", title = NULL, subtitle = NULL, x = NULL) +
  theme(plot.tag.position = c(-0.02, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        axis.text.x = element_text())

# Combine
library(patchwork)
plot.Shannon.FBratio.combined <- plot.FBratio | plot.Shannon.AlphaDiversity
plot.Shannon.FBratio.combined
##### beta diversity (Bray divergence) (ok) #####

# src https://microbiome.github.io/tutorials/Betadiversity.html
# requires library(reshape)
# ref Salonen et al. ISME J 2014
# http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html

##### beta diversity (Bray divergence) (ok) #####
pseq <- temp.ps

# Calculate Bray divergence for each group
b.Control_food <- divergence(subset_samples(pseq, application_method == "Control_food"),
                             apply(abundances(subset_samples(pseq, application_method == "Control_food")), 1, median))

b.Rapamycin_food <- divergence(subset_samples(pseq, application_method == "Rapamycin_food"),
                               apply(abundances(subset_samples(pseq, application_method == "Rapamycin_food")), 1, median))

b.Vehicle_IP <- divergence(subset_samples(pseq, application_method == "Vehicle_IP"),
                           apply(abundances(subset_samples(pseq, application_method == "Vehicle_IP")), 1, median))

b.Rapamycin_IP <- divergence(subset_samples(pseq, application_method == "Rapamycin_IP"),
                             apply(abundances(subset_samples(pseq, application_method == "Rapamycin_IP")), 1, median))

# Combine distances
l <- list(b.Control_food, b.Rapamycin_food, b.Vehicle_IP, b.Rapamycin_IP)
df <- reshape2::melt(l)

# Rename groups - keep as Vehicle-IP to match your data
df$L1[df$L1 == '1'] <- 'Control-Diet'
df$L1[df$L1 == '2'] <- 'Rapamycin-Diet'
df$L1[df$L1 == '3'] <- 'Control-IP'
df$L1[df$L1 == '4'] <- 'Rapamycin-IP'

# Set factor levels
df$L1 <- factor(df$L1, levels = c("Control-Diet", "Rapamycin-Diet", "Control-IP", "Rapamycin-IP"))

# Rename columns
bray_table <- df
colnames(bray_table) <- c("value", "Group")

# Define comparisons - NOW MATCHING THE ACTUAL GROUP NAMES
my_comparisons_bray <- list(
  c("Control-Diet", "Rapamycin-Diet"),
  c("Control-IP", "Rapamycin-IP")  # Changed from Control-IP to Vehicle-IP
)

# Plot - ADDED x = "Group"
##inter-individialDivergence
plot.Bray.BetaDiversity <- ggboxplot(bray_table, 
                                     x = "Group",  # ADDED THIS
                                     y = "value",
                                     color = "Group", 
                                     palette = "jco",
                                     ylab = "Divergence") +
  stat_compare_means(comparisons = my_comparisons_bray, method = "wilcox.test") +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(bray_table$value) * 1.15, 
                     label.x = 2) +
  labs(title = NULL, subtitle = NULL, x = NULL) +
  theme(legend.position = "none")

plot.Bray.BetaDiversity

# Save the table
table.betadiv <- bray_table

# cleanup
rm(pseq, bray_table, b.Control_food, b.Rapamycin_food, b.Vehicle_IP, b.Rapamycin_IP, l, df)
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
wunifrac.adonis2 <- adonis2(
  iDist ~ application_method,
  data         = as(sample_data(pseq), "data.frame"),
  permutations = 10000,
  by           = "margin"    # tests each term marginally
)

# Correct extraction of the p-value (first row, column "Pr(>F)")
pval <- wunifrac.adonis2["application_method", "Pr(>F)"]

# Extract p-value from PERMANOVA results
tempSubtitle.here <- paste0("Adonis PERMANOVA p-value = ", 
                            formatC(pval, format = "f", digits = 5))

# Plot the PCoA ordination, coloring by application method
plot.unifrac.pcoa <- plot_ordination(pseq, iMDS, color = "application_method") +
  labs(
    title    = "PCoA Ordination (Weighted UniFrac)",
    subtitle = tempSubtitle.here,
    color    = "Group"   # relabel the points/colour legend
  )  +
  stat_ellipse(geom = "polygon", aes(fill = application_method), alpha = 0.15 , 
  show.legend = FALSE  ) +     # ← drop the fill legend here# Add ellipses for each group
  theme_bw()

dispersion <- betadisper(iDist, sample_data(pseq)$application_method)
anova(dispersion)
permutest(dispersion)


# Display the plot
print(plot.unifrac.pcoa)

# Clean up
rm(iDist, iMDS, wunifrac.adonis2, tempSubtitle.here)

##### Unweighted UniFrac ####
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
  iDist ~ application_method,
  data         = as(sample_data(pseq), "data.frame"),
  permutations = 10000,
  by           = "margin"    # tests each term marginally
)

# Correct extraction of the p-value
pval <- unifrac.adonis2["application_method", "Pr(>F)"]

# Extract p-value from PERMANOVA results
tempSubtitle.here <- paste0("Adonis PERMANOVA p-value = ", 
                            formatC(pval, format = "f", digits = 5))

# Plot the PCoA ordination, coloring by application method
plot.weightedunifrac.pcoa <- plot_ordination(pseq, iMDS, color = "application_method") +
  labs(
    title    = "PCoA Ordination (Unweighted UniFrac)",
    subtitle = tempSubtitle.here,
    color    = "Group"   # relabel the points/colour legend
  )  +
  stat_ellipse(geom = "polygon", aes(fill = application_method), alpha = 0.15, 
               show.legend = FALSE) +     # drop the fill legend here
  theme_bw()

dispersion <- betadisper(iDist, sample_data(pseq)$application_method)
anova(dispersion)
permutest(dispersion)

# Display the plot
print(plot.weightedunifrac.pcoa)

# Clean up
rm(iDist, iMDS, unifrac.adonis2, tempSubtitle.here)
##combine
# Update Weighted UniFrac plot (a) - NO LEGEND, keep original axes
plot.unifrac.pcoa <- plot.unifrac.pcoa + 
  scale_color_discrete(labels = c("Control_food" = "Control-Diet",
                                  "Rapamycin_food" = "Rapamycin-Diet",
                                  "Vehicle_IP" = "Control-IP", 
                                  "Rapamycin_IP" = "Rapamycin-IP")) +
  scale_fill_discrete(labels = c("Control_food" = "Control-Diet",
                                 "Rapamycin_food" = "Rapamycin-Diet",
                                 "Vehicle_IP" = "Control-IP", 
                                 "Rapamycin_IP" = "Rapamycin-IP")) +
  labs(tag = "a", title = NULL, subtitle = NULL) +
  theme(plot.tag.position = c(-0.00, 1),
        plot.tag = element_text(size = 14, face = "plain"),
        legend.position = "none")

# Update Unweighted UniFrac plot (b) - KEEP LEGEND, keep original axes
plot.weightedunifrac.pcoa <- plot.weightedunifrac.pcoa + 
  scale_color_discrete(labels = c("Control_food" = "Control-Diet",
                                  "Rapamycin_food" = "Rapamycin-Diet",
                                  "Vehicle_IP" = "Control-IP", 
                                  "Rapamycin_IP" = "Rapamycin-IP")) +
  scale_fill_discrete(labels = c("Control_food" = "Control-Diet",
                                 "Rapamycin_food" = "Rapamycin-Diet",
                                 "Vehicle_IP" = "Control-IP", 
                                 "Rapamycin_IP" = "Rapamycin-IP")) +
  labs(tag = "b", title = NULL, subtitle = NULL) +
  theme(plot.tag.position = c(0.01, 1),
        plot.tag = element_text(size = 14, face = "plain"))

# Combine side by side
library(patchwork)
plot.unifrac.combined <- plot.unifrac.pcoa | plot.weightedunifrac.pcoa
plot.unifrac.combined

##### Differential abundances: Rapamycin IP vs Vehicle IP  (ok) #####
# FDR cutoff for differential abundances
fdr.threshold <- 0.05

# subsetting 
temp.diagdds <- subset_samples(temp.ps, application_method %in% c("Vehicle_IP", "Rapamycin_IP"))
# deseq object
diagdds <- phyloseq_to_deseq2(temp.diagdds, ~ application_method)

# Set Vehicle_IP as the reference level
diagdds$application_method <- relevel(factor(diagdds$application_method), ref = "Vehicle_IP")


# code added to avoid zeroes error; see functions' section
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local")

#diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
res <- results(diagdds, cooksCutoff = FALSE)
head(res)

print(paste("Total number of ASVs tested:", nrow(res)))

# # Filter out low abundance taxa
# min_count <- 5  # Minimum count in at least one group
# # Get normalized counts
# normalized_counts <- counts(diagdds, normalized = TRUE)
# 
# # Calculate mean counts by group
# mean_counts_vehicle <- rowMeans(normalized_counts[, colData(diagdds)$application_method == "Vehicle_IP"])
# mean_counts_rapamycin <- rowMeans(normalized_counts[, colData(diagdds)$application_method == "Rapamycin_IP"])
# 
# # Apply abundance filter before statistical filtering
# abundant_taxa <- which(mean_counts_vehicle > min_count |
#                          mean_counts_rapamycin > min_count)
# res_filtered <- res[abundant_taxa, ]

# Get total counts per sample
sample_totals <- colSums(counts(diagdds))

# Calculate relative abundance (proportion of reads)
rel_abundance <- sweep(counts(diagdds), 2, sample_totals, "/")

# Set minimum relative abundance threshold (e.g., 0.001 = 0.1%)
min_rel_abundance <- 0.0001  # Adjust this value as needed

# Calculate mean relative abundance by group
mean_rel_abundance_control <- rowMeans(rel_abundance[, colData(diagdds)$application_method == "Vehicle_IP"])
mean_rel_abundance_rapamycin <- rowMeans(rel_abundance[, colData(diagdds)$application_method == "Rapamycin_IP"])

# Apply relative abundance filter
abundant_taxa <- which(mean_rel_abundance_control > min_rel_abundance | 
                         mean_rel_abundance_rapamycin > min_rel_abundance)
res_filtered <- res[abundant_taxa, ]

print(paste("ASVs after minimum count filtering:", nrow(res_filtered)))

# > head(res)
# log2 fold change (MLE): application method Rapamycin food vs Control food 
# Wald test p-value: application method Rapamycin food vs Control food 
# i.e. log2(Rapamycin_IP/Vehicle_IP )

# "The text, condition treated vs untreated, tells you that the 
# estimates are of the logarithmic fold change log2(treated/untreated)."
# source https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

# filter results by FDR (p-adj) < 0.01
# Apply both statistical significance and fold change filters
sigtab <- res_filtered[which(!is.na(res_filtered$padj) & 
                               !is.na(res_filtered$log2FoldChange) & 
                               res_filtered$padj < fdr.threshold & 
                               abs(res_filtered$log2FoldChange) > 1), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(temp.diagdds)[rownames(sigtab), ], "matrix"))

# Merge unclassified and uncultured taxa
sigtab$Genus <- as.character(sigtab$Genus)  # Convert to character
# Replace 'unclassified' and variations
unclassified_patterns <- c("unclassified", "Unclassified", "uncultured", "Uncultured")
for (pattern in unclassified_patterns) {
  sigtab$Genus[grepl(pattern, sigtab$Genus)] <- "Unclassified/Uncultured"
}

# Find rows with unknown genera
unknown_genera <- sigtab[is.na(sigtab$Genus) | 
                           sigtab$Genus == "" | 
                           sigtab$Genus == "Unknown" |
                           sigtab$Genus == "Unclassified/Uncultured", ]

# Print full taxonomy for these ASVs
if(nrow(unknown_genera) > 0) {
  print("Unknown genera with significant differential abundance:")
  print(unknown_genera[, c("Phylum", "Class", "Order", "Family", "Genus", "log2FoldChange", "padj")])
}

# [1] "Unknown genera with significant differential abundance:"
# Phylum       Class           Order           Family                   Genus
# e485bd0337c3e09152d416e1b4ecca5e  Firmicutes  Clostridia  Lachnospirales  Lachnospiraceae Unclassified/Uncultured
# log2FoldChange         padj
# e485bd0337c3e09152d416e1b4ecca5e       22.68019 6.346198e-12

# Get the sequence for the unclassified Lachnospiraceae ASV
target_asv_id <- "e485bd0337c3e09152d416e1b4ecca5e"

# Find the specific ASV and update its genus
sigtab$Genus[rownames(sigtab) == "e485bd0337c3e09152d416e1b4ecca5e"] <- "Lacrimispora sp."

x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# plot
tempTitle <- "Differential Abundance by Genera"
tempSubtitle <- paste0("Rapamycin IP vs Vehicle IP; FDR < ", fdr.threshold)
plot.DiffAbund.IP <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  labs(title = tempTitle, subtitle = tempSubtitle) +
  geom_hline(yintercept=0)
plot.DiffAbund.IP

# Sort by max LFC at genus level - ascending for horizontal layout
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, FALSE)  # FALSE for ascending order (bottom to top)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# Horizontal plot - just swap x and y
plot.DiffAbund.IP_flip <- ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Genus)) + 
  geom_point(size=6) + 
  geom_vline(xintercept=0) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Differential Abundance by Genera",
    subtitle = paste0("Rapamycin IP vs Vehicle IP; FDR < ", fdr.threshold),
    x = expression(log[2]~fold~change),
    y = NULL
  )
plot.DiffAbund.IP_flip

sigtab_IP <- sigtab

# diffex table
table.DiffAbund.IP <- sigtab[,c(7:13,1:6)]

##### Differential abundances: Control_food vs Rapamycin_food (ok) #####
# subsetting 
temp.diagdds <- subset_samples(temp.ps, application_method %in% c("Control_food", "Rapamycin_food"))
# deseq object
diagdds <- phyloseq_to_deseq2(temp.diagdds, ~ application_method)

# Set Control_food as the reference level
diagdds$application_method <- relevel(factor(diagdds$application_method), ref = "Control_food")

# code added to avoid zeroes error; see functions' section
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local")

#diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
res <- results(diagdds, cooksCutoff = FALSE)
head(res)

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
min_rel_abundance <- 0.0001  # Adjust this value as needed

# Calculate mean relative abundance by group
mean_rel_abundance_control <- rowMeans(rel_abundance[, colData(diagdds)$application_method == "Control_food"])
mean_rel_abundance_rapamycin <- rowMeans(rel_abundance[, colData(diagdds)$application_method == "Rapamycin_food"])

# Apply relative abundance filter
abundant_taxa <- which(mean_rel_abundance_control > min_rel_abundance | 
                         mean_rel_abundance_rapamycin > min_rel_abundance)
res_filtered <- res[abundant_taxa, ]

# Statistical significance and fold change filter
sigtab <- res_filtered[which(!is.na(res_filtered$padj) & 
                               !is.na(res_filtered$log2FoldChange) & 
                               res_filtered$padj < fdr.threshold & 
                               abs(res_filtered$log2FoldChange) > 3), ]

print(paste("ASVs after minimum count filtering:", nrow(res_filtered)))

# > head(res)
# log2 fold change (MLE): Gruppe Control_food vs Rapamycin_food
# Wald test p-value: Gruppe Control_food vs Rapamycin_food
# i.e. log2(Control_food/Rapamycin_food)

sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(temp.diagdds)[rownames(sigtab), ], "matrix"))

# Merge unclassified and uncultured taxa
sigtab$Genus <- as.character(sigtab$Genus)  # Convert to character
# Replace 'unclassified' and variations
unclassified_patterns <- c("unclassified", "Unclassified", "uncultured", "Uncultured")
for (pattern in unclassified_patterns) {
  sigtab$Genus[grepl(pattern, sigtab$Genus)] <- "Unclassified/Uncultured"
}

# Find rows with unknown genera
unknown_genera <- sigtab[is.na(sigtab$Genus) | 
                           sigtab$Genus == "" | 
                           sigtab$Genus == "Unknown" |
                           sigtab$Genus == "Unclassified/Uncultured", ]

# correcting cell without genera info
sigtab$Genus <- as.character(sigtab$Genus)
sigtab[rownames(sigtab) == "374aa6542d0867e9dae3170b044e024c", "Genus"] <- "[Prevotellaceae]"
sigtab[rownames(sigtab) == "20a26307116045d45a923bbfbf22de3a", "Genus"] <- "[Prevotellaceae]"
sigtab$Genus[rownames(sigtab) == "fb6fce060e3008b919dad776ec02e05f"] <- "Lachnospiraceae"

# Print full taxonomy for these ASVs
if(nrow(unknown_genera) > 0) {
  print("Unknown genera with significant differential abundance:")
  print(unknown_genera[, c("Phylum", "Class", "Order", "Family", "Genus", "log2FoldChange", "padj")])
}

# Sort by max LFC at genus level (only once)
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# plot
tempTitle <- "Differential Abundance by Genera"
tempSubtitle <- paste0("Rapamycin_food vs Control_food; FDR < ", fdr.threshold)
plot.DiffAbund.food <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  labs(title = tempTitle, subtitle = tempSubtitle) +
  geom_hline(yintercept=0)
plot.DiffAbund.food

x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

plot.DiffAbund.food_flip <- ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Genus)) + 
  geom_point(size=6) + 
  geom_vline(xintercept=0) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = expression(log[2]~fold~change), y = NULL)


# diffex table
table.DiffAbund.food <- sigtab[,c(7:13,1:6)]
#### Differential abundances merged #####

# Get shared x-axis limits
x_limits <- range(c(sigtab_IP$log2FoldChange, sigtab$log2FoldChange))

# Apply same scale to both
plot.DiffAbund.IP_flip <- plot.DiffAbund.IP_flip + scale_x_continuous(limits=x_limits)
plot.DiffAbund.food_flip <- plot.DiffAbund.food_flip + scale_x_continuous(limits=x_limits)

# Remove titles and add tags with positioning
plot.DiffAbund.IP_flip <- plot.DiffAbund.IP_flip + 
  labs(title = NULL, subtitle = NULL, tag = "a") +
  theme(plot.tag.position = c(0.05, 1))

plot.DiffAbund.food_flip <- plot.DiffAbund.food_flip + 
  labs(title = NULL, subtitle = NULL, tag = "b") +
  theme(plot.tag.position = c(0.05, 1))

# Combine
library(patchwork)
plot.DiffAbund.combined <- plot.DiffAbund.IP_flip | plot.DiffAbund.food_flip
plot.DiffAbund.combined

##### Differential abundances: Pathways (ok) #####


# Remove only pathway-related objects
rm(path_abund, picrust2Abund, abundance.table, wt_samples, wt_conditions, 
   wt_samples_meta, wt_samples_meta, wt_dds, table.DiffPathw.WT, 
   plot.DiffPathw.WT, clean_df)


# table abundances with KO ids
# FDR cutoff for differential abundances
fdr.threshold <- 0.05


# Load PICRUSt2 pathway abundance data with the full file path
path_abund <- read.table(here::here("data/study2_PRJEB13256/picrust2/path_abun_unstrat.tsv.gz"), 
                         sep="\t",
                         header=TRUE,
                         row.names=1, 
                         stringsAsFactors=FALSE,
                         comment.char="")

# Format similar to original picrust1Abund
picrust2Abund <- data.frame(path_abund)
# 
# 
# # Create pathway IDs in the same format as original code
# numlimit <- nrow(picrust2Abund)
# preffix <- "PWY"
# suffix <- seq(1:numlimit)
# picrust2Abund$path.id <- paste(preffix, suffix, sep = "")
# picrust2Abund$ <- rownames(picrust2Abund)  # Use pathway IDs as descriptions for now

# Use the real MetaCyc pathway IDs
picrust2Abund$path.id <- rownames(picrust2Abund)  # Real MetaCyc IDs
picrust2Abund$metacyc.desc <- rownames(picrust2Abund)  # Same thing, but better name

# Extract the abundance table
abundance.table <- picrust2Abund[, 1:(ncol(picrust2Abund)-2)]
rownames(abundance.table) <- picrust2Abund$path.id
cat("Dimension of raw abundance.table:\n")

print(dim(abundance.table))
cat("First few column names of abundance.table:\n")
print(head(colnames(abundance.table)))

common_samples <- intersect(colnames(abundance.table), rownames(metadata))
length(common_samples)  # how many samples overlap

minCount <- 1000
abundance.table <- subset(abundance.table, rowSums(abundance.table) > minCount)
abundance.table <- abundance.table[ , order(names(abundance.table))]

colnames(abundance.table)
rownames(metadata)
abundance.table <- round(abundance.table)

# diff abundant PATHWAYS
abundance.ds <- DESeqDataSetFromMatrix(countData = abundance.table, 
                                       colData = metadata, 
                                       design = ~sex)
abundance.ds <- DESeq(abundance.ds, test="Wald", fitType="parametric")

# filter results by FDR (p-adj) < 0.01
# fdr.threshold <- 0.01

rm(table.DiffPathw.Food, table.DiffPathw.IP)

# Contrasts for IP and Food:
# IP: "Rapamycin_IP" vs "Vehicle_IP"
table.DiffPathw.IP <- results(
  abundance.ds, 
  contrast = c("sex", "Rapamycin_IP", "Vehicle_IP")
)

# Food: "Rapamycin_food" vs "Control_food"
table.DiffPathw.Food <- results(
  abundance.ds,
  contrast = c("sex", "Rapamycin_food", "Control_food")
)

# Filter by FDR
table.DiffPathw.IP   <- subset(table.DiffPathw.IP,   padj < fdr.threshold)
table.DiffPathw.Food <- subset(table.DiffPathw.Food, padj < fdr.threshold)

# Convert to data frames
table.DiffPathw.IP   <- as.data.frame(table.DiffPathw.IP)
table.DiffPathw.Food <- as.data.frame(table.DiffPathw.Food)

# Merge IP results with path ID + description
table.DiffPathw.IP <- merge(
  table.DiffPathw.IP,
  picrust2Abund[, c("path.id", "metacyc.desc")],
  by.x = "row.names",   # row names from DESeq2 result
  by.y = "path.id"      # column in your picrust data
)

# Rename "Row.names" to "path.id"
names(table.DiffPathw.IP)[names(table.DiffPathw.IP) == "Row.names"] <- "path.id"

# Same for the food contrast
table.DiffPathw.Food <- merge(
  table.DiffPathw.Food,
  picrust2Abund[, c("path.id", "metacyc.desc")],
  by.x = "row.names",
  by.y = "path.id"
)
names(table.DiffPathw.Food)[names(table.DiffPathw.Food) == "Row.names"] <- "path.id"


# IP bar plot
plot.DiffPathw.IP <- ggplot(
  table.DiffPathw.IP,
  aes(
    x   = reorder(metacyc.desc, log2FoldChange),
    y   = log2FoldChange,
    fill = log2FoldChange > 0
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  guides(fill = "none") +
  theme_minimal() +
  labs(
    x        = "MetaCyc Pathways",
    y        = "Log2(Fold Change)",
    title    = "Rapamycin IP vs Vehicle IP",
    subtitle = paste0("FDR < ", fdr.threshold)
  )

print(plot.DiffPathw.IP)

# Food bar plot
plot.DiffPathw.Food <- ggplot(
  table.DiffPathw.Food,
  aes(
    x   = reorder(metacyc.desc, log2FoldChange),
    y   = log2FoldChange,
    fill = log2FoldChange > 0
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  guides(fill = "none") +
  theme_minimal() +
  labs(
    x        = "MetaCyc Pathways",
    y        = "Log2(Fold Change)",
    title    = "Rapamycin food vs Control food",
    subtitle = paste0("FDR < ", fdr.threshold)
  )

print(plot.DiffPathw.Food)


 (table.plot.DiffPathw.IP | plot.DiffPathw.Food)


tempTitle <- "Differential Pathway Abundance"
tempSubtitle <- paste0("Rapamycin_food vs Control_food; FDR < ", fdr.threshold)


plot.DiffPathw.Food <- ggplot(table.DiffPathw.Food,
                              aes(x = metacyc.desc, y = log2FoldChange, color = log2FoldChange > 0)) +
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
print(plot.DiffPathw.Food)

# IP plot in genus style
tempTitle <- "Differential Pathway Abundance"
tempSubtitle <- paste0("Rapamycin_IP vs Vehicle_IP; FDR < ", fdr.threshold)

plot.DiffPathw.IP <- ggplot(table.DiffPathw.IP,
                            aes(x = metacyc.desc, y = log2FoldChange,  color = log2FoldChange > 0)) +
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

##### Pathway merged ####
# Flip both pathway plots with bars
plot.DiffPathw.IP_flip <- ggplot(table.DiffPathw.IP,
                                 aes(x = log2FoldChange, 
                                     y = reorder(metacyc.desc, log2FoldChange), 
                                     fill = log2FoldChange > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "#F8766D")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9)) +
  labs(x = expression(log[2]~fold~change), y = NULL, tag = "a") +
  theme(plot.tag.position = c(0.05, 1))

plot.DiffPathw.Food_flip <- ggplot(table.DiffPathw.Food,
                                   aes(x = log2FoldChange, 
                                       y = reorder(metacyc.desc, log2FoldChange), 
                                       fill = log2FoldChange > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "#F8766D")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9)) +
  labs(x = expression(log[2]~fold~change), y = NULL, tag = "b") +
  theme(plot.tag.position = c(0.05, 1))

# Get shared scale and combine
x_limits <- range(c(table.DiffPathw.IP$log2FoldChange, table.DiffPathw.Food$log2FoldChange))
plot.DiffPathw.IP_flip <- plot.DiffPathw.IP_flip + scale_x_continuous(limits = x_limits)
plot.DiffPathw.Food_flip <- plot.DiffPathw.Food_flip + scale_x_continuous(limits = x_limits)

plot.DiffPathw.combined <- plot.DiffPathw.IP_flip | plot.DiffPathw.Food_flip
plot.DiffPathw.combined

##### clinical markers data (ok) #####

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

### 2025-03-13 Exported as PDF
pdf("2.RapamycinAnalysis_20250506_V03.pdf", width = 14, height = 10)
plot.PhylumComp # phylum composition 0,1% detection 50% Prevelance
plot.GenusComp # genus composition
plot.PhylumBubble #Phylum Bubbleplot
plot.GenusBubble #Genus Bubbleplot
plot.FBratio # FB plot
plot.Core.Controlfood # Core microbiome control_food
plot.Core.Rapamycinfood # Core microbiome Rapamycin_food
plot.Core.VehicleIP # Core microbiome Vehicle_IP
plot.Core.RapamycinIP #Core microbiome Rapamycin_IP
plot.Core.comparison# Core microbiome comparison
#plot.BugBase # BugBase
plot.Shannon.AlphaDiversity # Shannon alpha diversity
plot.Bray.BetaDiversity # Bray divergence (beta div)
plot.unifrac.pcoa #new BEta diversity made after Franzis script (betadiv #2)
plot.weightedunifrac.pcoa
plot.DiffAbund.IP + geom_segment(aes(xend = Genus, yend = 0)) # Diff Abd G vs H
plot.DiffAbund.food + geom_segment(aes(xend = Genus, yend = 0)) # Diff abd P vs H
plot.DiffPathw.Food
plot.DiffPathw.IP
dev.off ()

plot.DiffPathw.Parod # Diff pathways Parod vs Healthy
plot.DiffPathw.ParodGingi # Diff pathways Parod vs Gingivitis
plot.RF.Parod_Gingi + 
  labs(title="Microbial biomarkers via supervised ML (RF)", subtitle="Parodontitis and Gingivitis")
plot.RF.all +
  labs(title="Microbial biomarkers via supervised ML (RF)", subtitle="all groups")

dev.off()

##### tables (in progress) #####

table.Phylum.Abundance # phylum composition
table.Genus.Abundance # genus composition
table.FBratio # FB ratios
list.Core.Controlfood # Core microbiome Healthy
list.Core.Rapamycinfood # Core microbiome Gingivitis
list.Core.VehicleIP# Core microbiome Parod
list.Core.RapamycinIP
table.Core.comparison.matrix # Core microbiome comparison
table.Core.comparison.matrix_for_ggplot2 # Core microbiome comparison f ggplot2
#table.bugbase # BugBase
table.shannon # Shannon alpha diversity
table.betadiv # Bray divergence (beta div)
table.DiffAbund.IP # Diff abundance Gingivitis vs Healthy
table.DiffAbund.food # Diff abundance Parod vs Healthy
table.DiffPathw.IP
table.DiffPathw.Food
------------------------------------------------------------
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
save.image("Rapamycin_reanalysis_2.PRJEB13256_V03_20250428_latest.RData")

# Load the saved workspace for later use
load("Rapamycin_reanalysis_2.PRJEB13256_V03_20250428_latest.RData")

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
