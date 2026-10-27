##### Functions V2.2 #####
# 2023-07-27 isradelacon@gmail.com

# clean all functions from memory
# detachAllPackages <- function() {
#   basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
#   package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
#   package.list <- setdiff(package.list,basic.packages)
#   if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
# }

# to paste without NAs
# source https://stackoverflow.com/a/40994869
paste_noNA <- function(x,sep=", ") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

# filtering GreenGenes data
# requires library(phyloseq)
filter_gg <- function(psObject){
  
  # Prune OTUs that are not present in any of the samples 
  temp.psObject <- prune_taxa(taxa_sums(psObject) > 0, psObject)
  
  # Filter out Eukaryotes, Archaea, chloroplasts and mitochondria
  # ref http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
  temp.psObject <-  subset_taxa(temp.psObject, 
                                Kingdom == "Bacteria" &
                                  Family  != "mitochondria" &
                                  Class   != "Chloroplast") 
  
  # Prune unidentified Genuses
  temp.psObject <-  subset_taxa(temp.psObject, Genus != "NA")
  temp.psObject <-  subset_taxa(temp.psObject, Genus != "-1")
  
  # output
  return(temp.psObject)
}

# filtering Silva data
# note: filters out Archaea, use with caution
# requires library(phyloseq)
filter_silva <- function(psObject){
  
  # Prune OTUs that are not present in any of the samples 
  temp.psObject <- prune_taxa(taxa_sums(psObject) > 0, psObject)
  
  # Prune non-bacterial OTUs
  # temp.psObject <- subset_taxa(temp.psObject, Rank1 != "Unassigned")
  temp.psObject <- subset_taxa(temp.psObject, Rank1 == "D_0__Bacteria")
  
  # Prune unidentified Genuses ("Rank6" for SILVA)
  temp.psObject <- subset_taxa(temp.psObject, Rank6 != "NA")
  temp.psObject <- subset_taxa(temp.psObject, Rank6 != "-1")
  temp.psObject <- subset_taxa(temp.psObject, Rank7 != "Ambiguous_taxa")
  
  # output
  return(temp.psObject)
}

# silva DESeq2 normalization and counts
# requires library(DESeq2)
normalize_DESeq2 <- function(psObject, silvaflag){
  
  # normalized data
  temp.psObject <- psObject
  temp.psObject <- phyloseq_to_deseq2(temp.psObject, ~ Group)
  
  # code added to avoid zeroes error; see above
  # calculate geometric means prior to estimate size factors
  geoMeans <- apply(counts(temp.psObject), 1, gm_mean)
  temp.psObject <- estimateSizeFactors(temp.psObject, geoMeans = geoMeans)
  
  # normalized counts
  psObject.normalized <- counts(temp.psObject, normalized=TRUE)
  psObject.normalized <- as.data.frame(psObject.normalized)
  psObject.normalized$otuid <- rownames(psObject.normalized)
  
  # raw counts
  psObject.counts <- counts(temp.psObject, normalized=FALSE)
  psObject.counts <- as.data.frame(psObject.counts)
  psObject.counts$otuid <- rownames(psObject.counts)
  
  # id assignment at the GENUS level
  # (agglomerating data at Species generates error since most have no species assignments)
  if(silvaflag == "TRUE"){
    #myTaxRank <- "Rank6"
    specnames <- psObject %>%
      tax_glom(taxrank = "Rank6") %>%                   # agglomerate at GENUS level
      psmelt() %>%                                      # Melt to long format
      arrange(Rank6)                                    # Sort data frame alphabetically by GENUS
  } else {
    specnames <- psObject %>%
      tax_glom(taxrank = "Genus") %>%                 
      psmelt() %>%                                      
      arrange(Genus)                                  
  }
  
  specnames <- as.data.frame(specnames)
  specnames <- unique(specnames[,c(1,16)])
  
  # Changed to "all = TRUE" to report all OTUs with reads, including those not identified until lower taxa level
  psObject.normalized <- merge(psObject.normalized, specnames, by.x = "otuid", by.y = "OTU", all = TRUE)
  psObject.counts <- merge(psObject.counts, specnames, by.x = "otuid", by.y = "OTU", all = TRUE)
  
  result <- list(psObject.counts, psObject.normalized)
  return(result)
}

# gm_mean -- Function to avoid zeroes error in differential abundances:
# "Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#  every gene contains at least one zero, cannot compute log geometric means"
# source https://github.com/joey711/phyloseq/issues/387
# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## Relative Abundances by Phyla
# requires library(phyloseq), library(dplyr)
DoRelativeAbundancePhylum <- function(psObject, myComparisonColumn, myFilter, myPlotTitle){

  # src https://github.com/joey711/phyloseq/issues/1089#issuecomment-668580364  
  # > Can you help me understand a bit more why I need to use transform_sample_counts twice?
  # The merging operation simply sums the abundances. If you want the abundances to sum to 1 prior to 
  # plotting, it is necessary to transform the abundances to proportions after the merging operation 
  # (hence the 2nd call to transform_sample_counts). So then why did I also transform to proportions 
  # before merging samples? I did that because I wanted to end up with an unweighted average of the 
  # proportions across samples within a sample group. You can certainly skip this step; in that case 
  # you would be taking a weighted average of the proportions across samples, where each sample is 
  # weighted by its read depth. Either way is reasonable.
    
  ps <- tax_glom(psObject, "Phylum")      # agglomerate at phylum level
  ps <- transform_sample_counts(ps, function(x) x / sum(x))   # Transform to rel. abundance
  ps <- merge_samples(ps, myComparisonColumn, fun=sum)
  ps <- transform_sample_counts(ps, function(x) x / sum(x))
  
  ps <- ps %>%
    psmelt() %>%                            # Melt to long format
    filter(Abundance > myFilter) %>%        # Filter out low abundance taxa (Phyla < myFilter)
    arrange(Phylum)                         # Sort data frame alphabetically by phylum
  
  # Kruskal-Wallis tests for relative abundance data
  ps$Sample <- as.factor(ps$Sample)
  ps.kruskal <- kruskal.test(Abundance ~ Sample, data=ps)
  # add these results to the subtitle
  # tempSubtitle.here <- "Relative abundances (Phyla > 1%)"
  # tempSubtitle.here <- paste0(tempSubtitle.here, "; Kruskal-Wallis p-value: ", round(ps.kruskal$p.value, digits = 4))
  tempSubtitle.here <- paste0("Kruskal-Wallis p-value: ", round(ps.kruskal$p.value, digits = 4))
  
  # saving relative abundances as a separate table
  myTable <- ps
  
  # Plot 
  tempTitle <- myPlotTitle
  myPlot <- ggplot(ps, aes(x = Sample, y = Abundance, fill = Phylum, label = Phylum)) + 
    geom_bar(stat = "identity", colour="black") +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    #scale_fill_manual("Legend", values = myphylacolors) + 
    labs(title = tempTitle, subtitle = tempSubtitle.here, y = "Relative Abundance (%)") #+
  
  result <- list(myTable, myPlot)
  return(result)
}

## Relative Abundances by Genera
# requires library(phyloseq), library(dplyr)
DoRelativeAbundanceGenus <- function(psObject, myComparisonColumn, myFilter, myPlotTitle){
  
  ps <- tax_glom(psObject, "Genus")      # agglomerate at phylum level
  ps <- transform_sample_counts(ps, function(x) x / sum(x))   # Transform to rel. abundance
  ps <- merge_samples(ps, myComparisonColumn, fun=sum)
  ps <- transform_sample_counts(ps, function(x) x / sum(x))
  
  ps <- ps %>%
    psmelt() %>%                                        # Melt to long format
    filter(Abundance > myFilter) %>%                        # Filter out low abundance taxa (Phyla < 1%)
    arrange(Genus)                                     # Sort data frame alphabetically by phylum
  
  # Kruskal-Wallis tests for relative abundance data
  ps$Sample <- as.factor(ps$Sample)
  ps.kruskal <- kruskal.test(Abundance ~ Sample, data=ps)
  # add these results to the subtitle
  # tempSubtitle.here <- "Relative abundances (Phyla > 1%)"
  # tempSubtitle.here <- paste0(tempSubtitle.here, "; Kruskal-Wallis p-value: ", round(ps.kruskal$p.value, digits = 4))
  tempSubtitle.here <- paste0("Kruskal-Wallis p-value: ", round(ps.kruskal$p.value, digits = 4))
  
  # saving relative abundances as a separate table
  myTable <- ps
  
  # Plot 
  tempTitle <- myPlotTitle
  myPlot <- ggplot(ps, aes(x = Sample, y = Abundance, fill = Genus, label = Genus)) + 
    geom_bar(stat = "identity", colour="black") +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    #scale_fill_manual("Legend", values = myphylacolors) + 
    labs(title = tempTitle, subtitle = tempSubtitle.here, y = "Relative Abundance (%)") #+
  
  result <- list(myTable, myPlot)
  return(result)
}

# UniFrac distances
# modified from https://github.com/joey711/phyloseq/issues/714
# and https://github.com/jeffkimbrel/jakR/blob/master/R/plotDistances.R
# requires library(phyloseq), library(dplyr)
plotDistances <- function(p, m, s, d) {
  
  # p: phyloseq object
  # m: distance method (wunifrac, unifrac, etc)
  # s: sample column
  # d: sample data for groups
  
  # calc distances
  wu = phyloseq::distance(p, m)
  wu.m = melt(as.matrix(wu))
  
  # remove self-comparisons
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor,as.character)
  
  # get sample data (S4 error OK and expected)
  sd = sample_data(p) %>%
    select(s, d) %>%
    mutate_if(is.factor,as.character)
  
  # combined distances with sample data
  colnames(sd) = c("Var1", "Type1")
  wu.sd = left_join(wu.m, sd, by = "Var1")
  
  colnames(sd) = c("Var2", "Type2")
  wu.sd = left_join(wu.sd, sd, by = "Var2")
  
  # plot
  ggplot(wu.sd, aes(x = Type2, y = value)) +
    theme_bw() +
    geom_point() +
    #geom_boxplot() +
    geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
    scale_color_identity() +
    facet_wrap(~ Type1, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggtitle(paste0("Distance Metric: ", m))
}

# UniFrac distances for qiime2 data
# modified from the code above
# requires library(phyloseq), library(dplyr)
plotDistancesQiime2 <- function(p, m, s, d) {
  
  # p: phyloseq object
  # m: distance method (wunifrac, unifrac, etc)
  # s: sample column
  # d: sample data for groups
  
  # calc distances
  wu = phyloseq::distance(p, m)
  wu.m = melt(as.matrix(wu))
  
  # remove self-comparisons
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor,as.character)
  
  # get sample data (S4 error OK and expected)
  sd = sample_data(temp.ps)
  class(sd) <- 'data.frame' # this made the difference!!
  sd = sd %>%
    select(s, d) %>%
    mutate_if(is.factor,as.character)
  
  # combined distances with sample data
  colnames(sd) = c("Var1", "Type1")
  wu.sd = left_join(wu.m, sd, by = "Var1")
  
  colnames(sd) = c("Var2", "Type2")
  wu.sd = left_join(wu.sd, sd, by = "Var2")
  
  # plot
  ggplot(wu.sd, aes(x = Type2, y = value)) +
    theme_bw() +
    geom_point() +
    #geom_boxplot() +
    geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
    scale_color_identity() +
    facet_wrap(~ Type1, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggtitle(paste0("Distance Metric: ", m))
}


# TidyPhyloseq --cleans phyloseq object
# modified from https://github.com/davidelliott/core-microbiome/blob/master/core-microbiome.Rmd
# requires library(phyloseq), library(dplyr)
tidy_phyloseq <- function(my_phyloseq){
  # set ranks. 
  # colnames(tax_table(my_phyloseq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # fix taxa names
  tax_table(my_phyloseq)[,colnames(tax_table(my_phyloseq))] <- gsub(tax_table(my_phyloseq)[,colnames(tax_table(my_phyloseq))],pattern="[a-z]__",replacement="")
  
  # NAs and empty phyla change it to "Unidentified"
  tax_table(my_phyloseq)[is.na(tax_table(my_phyloseq)[,"Phylum"]), "Phylum"] <- "Unidentified"
  tax_table(my_phyloseq)[tax_table(my_phyloseq)[,"Phylum"] == "__", "Phylum"] <- "Unidentified"
  
  # Prune OTUs that are not present in any of the samples 
  my_phyloseq <- prune_taxa(taxa_sums(my_phyloseq) > 0, my_phyloseq)
  
  # Prune unidentified Genuses; 
  # my_phyloseq <-  subset_taxa(my_phyloseq, Genus != "NA")
  # my_phyloseq <-  subset_taxa(my_phyloseq, Genus != "Unassigned")
  # my_phyloseq <-  subset_taxa(my_phyloseq, Genus != "unidentified")
  # my_phyloseq <-  subset_taxa(my_phyloseq, Genus != "-1")
  
  # Filter out chloroplasts and mitochondria
  # ref http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
  my_phyloseq <-  subset_taxa(my_phyloseq, Family != "mitochondria" & Class != "Chloroplast") 
  
  # Filter out eukaryotes and archaea
  # my_phyloseq <-  subset_taxa(my_phyloseq, Kingdom == "Bacteria") 
  
  return(my_phyloseq)
}


# add_global_label: global x and y-axis labels for patchwork plots
# source https://github.com/thomasp85/patchwork/issues/43#issuecomment-692781479
add_global_label <- function(pwobj, Xlab = NULL, Ylab = NULL, Xgap = 0.03, Ygap = 0.03, ...) {
  ylabgrob <- patchwork::plot_spacer()
  if (!is.null(Ylab)) {
    ylabgrob <- ggplot() +
      geom_text(aes(x = .5, y = .5), label = Ylab, angle = 90, ...) +
      theme_void()
  }
  if (!is.null(Xlab)) {
    xlabgrob <- ggplot() +
      geom_text(aes(x = .5, y = .5), label = Xlab, ...) +
      theme_void()
  }
  if (!is.null(Ylab) & is.null(Xlab)) {
    return((ylabgrob + patchworkGrob(pwobj)) + 
             patchwork::plot_layout(widths = 100 * c(Ygap, 1 - Ygap)))
  }
  if (is.null(Ylab) & !is.null(Xlab)) {
    return((ylabgrob + pwobj) + 
             (xlabgrob) +
             patchwork::plot_layout(heights = 100 * c(1 - Xgap, Xgap),
                                    widths = c(0, 100),
                                    design = "
                                    AB
                                    CC
                                    "
             ))
  }
  if (!is.null(Ylab) & !is.null(Xlab)) {
    return((ylabgrob + pwobj) + 
             (xlabgrob) +
             patchwork::plot_layout(heights = 100 * c(1 - Xgap, Xgap),
                                    widths = 100 * c(Ygap, 1 - Ygap),
                                    design = "
                                    AB
                                    CC
                                    "
             ))
  }
  return(pwobj)
}


# x-axis scale
scaleFUN <- function(x) sprintf("%1.3f%%", 100*as.numeric(x))

# core microbiome calculation and heatmap
# src library(microbiome)
docoresV2 <- function(pseq, plotTitle, detection, prevalence){
  # pseq <- microbiome::aggregate_taxa(pseq, "Genus")
  pseq.rel2 <- microbiome::transform(pseq, "compositional")
  
  # core taxa: taxa that exceed the given prevalence and detection thresholds
  # core.taxa.standard2 <- core_members(pseq.rel2, detection = 0.001, prevalence = 80/100)
  
  # A full phyloseq object of the core microbiota is obtained as follows:
  # pseq.core2 <- core(pseq.rel2, detection = 0.001, prevalence = 80/100)
  pseq.core2 <- core(pseq.rel2, detection = detection, prevalence = prevalence)
  
  # get the taxonomy data
  tax.mat2 <- tax_table(pseq.core2)
  tax.df2 <- as.data.frame(tax.mat2[,c(1:6)])
  
  # add the OTus to last column
  tax.df2$OTU <- rownames(tax.df2)
  
  # plotting variables
  prevalences <- seq(.05, 1, .05)
  detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)
  gray <- gray(seq(0,1,length=5))
  
  # core heatmap
  myplot <- plot_core(pseq.rel2, 
                      plot.type = "heatmap", 
                      colours = gray,
                      prevalences = prevalences, 
                      detections = detections, 
                      min.prevalence = prevalence) 
  
  myplot <- myplot  + scale_fill_viridis(option = "E") +
    scale_x_discrete(labels=scaleFUN) +
    theme(axis.text.x = element_text(angle = 45, size = 8)) +
    labs(title = plotTitle, x = "")
  
  # Retrieving the associated taxa names from the phyloseq object:
  core.taxa2 <- taxa(pseq.core2)
  
  cat("Number of taxa in the whole dataset:\n")
  print(length(taxa(pseq)))
  cat("Number of taxa in core:\n")
  print(length(core.taxa2))
  # cat("Taxa:\n")
  # print(core.taxa2)
  
  # show number of all taxa
  tempdata <- taxa(pseq)
  length(tempdata)
  
  result <- list(core.taxa2, myplot)
  return(result)
  #return(myplot)
}


# standard shannon
doShannon <- function(temp.ps){
  # tempTitle <- "Richness (Shannon alpha diversity)"
  Richness_plot_group <- plot_richness(temp.ps, 
                                       x = "group", 
                                       color = "group", 
                                       measures = c("Shannon")) 
  alpha.diversity <- estimate_richness(temp.ps, 
                                       measures = c("Shannon"))
  tempdata <- cbind(sample_data(temp.ps), 
                    alpha.diversity) 
  ps.kruskal_group <- kruskal.test(Shannon ~ group, 
                                   data=tempdata) # Kruskal-Wallis test
  tempSubtitle <- paste0("Kruskal-Wallis p-value < ", 
                         round(ps.kruskal_group$p.value, 
                               digits = 4))
  Richness_plot_group <- Richness_plot_group + 
    geom_boxplot(data = Richness_plot_group$data, 
                 aes(x = group, y = value, color = group), alpha = 0.1) +
    labs(subtitle = tempSubtitle, y = "Diversity Measure") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    ylim(1, 5)
  return(Richness_plot_group)
}


