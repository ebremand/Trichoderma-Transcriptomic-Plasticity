###############################################################################
# Script: RNAseq_DE_AskoR.R
# Description:
#   Perform differential expression analysis using AskoR on RNA-seq count data.
#   Generates DE results, Upset plots, GO enrichment, and co-expression analyses.
###############################################################################

library(devtools)
# devtools::install_github("askomics/askoR", force = TRUE)
library(askoR)
library(dplyr)

# Clear workspace
rm(list = ls())

# Source AskoR functions
source("/path/to/AskoR.R")

# Set working directory
setwd("/path/to/output_directory")

# -------------------------
# === PARAMETERS ===
# -------------------------
parameters <- Asko_start()

parameters$analysis_name   <- "AskoR_results"
parameters$fileofcount     <- "counts_matrix.txt"
parameters$sep             <- "\t"
parameters$annotation      <- "gene_annotations.txt"
parameters$geneID2GO_file  <- "GO_annotations.txt"
parameters$contrast_file   <- "contrast_matrix.txt"
parameters$sample_file     <- "samples_description.txt"
parameters$rm_sample       <- c()

parameters$threshold_cpm   <- 0.5
parameters$replicate_cpm   <- 3
parameters$threshold_FDR   <- 0.05
parameters$threshold_logFC <- 1
parameters$normal_method   <- "TMM"
parameters$p_adj_method    <- "fdr"
parameters$glm             <- "qlf"

parameters$CompleteHeatmap <- FALSE
parameters$logFC           <- FALSE
parameters$FC              <- TRUE
parameters$logCPM          <- FALSE
parameters$FDR             <- TRUE
parameters$LR              <- FALSE
parameters$Sign            <- TRUE
parameters$Expression      <- TRUE
parameters$mean_counts     <- TRUE
parameters$norm_counts     <- TRUE

parameters$densbotmar      <- 20
parameters$densinset       <- 0.45
parameters$legendcol       <- 6

parameters$plotMD          <- FALSE
parameters$plotVO          <- FALSE
parameters$glimMD          <- FALSE
parameters$glimVO          <- FALSE

# -------------------------
# === LOAD DATA ===
# -------------------------
cat("\nLoading data...\n")
data <- loadData(parameters)

cat("\nSamples, contrast, design:\n")
data$samples
data$contrast
data$design
cat("Total genes:", dim(data$dge$counts)[1], "\n")
cat("Total samples:", dim(data$dge$counts)[2], "\n")
cat("CPM summary:\n")
summary(edgeR::cpm(data$dge))

# Prepare AskoR objects
asko_data <- asko3c(data, parameters)

# Filter genes
cat("\nFiltering genes...\n")
asko_filt <- GEfilt(data, parameters)
cat("Filtered genes:", dim(asko_filt$counts)[1], "\n")

# Normalize
asko_norm <- GEnorm(asko_filt, asko_data, data, parameters)

# Correlation
GEcorr(asko_norm, parameters)

# Differential expression
cat("\nDifferential expression analysis...\n")
resDEG <- DEanalysis(asko_norm, data, asko_data, parameters)

# -------------------------
# === Upset Plots ===
# -------------------------
parameters$upset_basic <- "all"
parameters$upset_type  <- "all"

# Examples of contrasts
parameters$upset_list <- c("LP_b_AbvsHP_b_Ab","LP_b_RsvsHP_b_Rs","LP_b_PuvsHP_b_Pu","LP_b_TavsHP_b_Ta")
UpSetGraph(resDEG, data, parameters)

parameters$upset_list <- c("LP_a_AbvsHP_a_Ab","LP_a_RsvsHP_a_Rs","LP_a_PuvsHP_a_Pu","LP_a_TavsHP_a_Ta")
UpSetGraph(resDEG, data, parameters)

# -------------------------
# === GO enrichment ===
# -------------------------
parameters$GO_threshold      <- 0.05
parameters$GO_min_num_genes  <- 10
parameters$GO                <- "both"
parameters$GO_algo           <- "weight01"
parameters$GO_stats          <- "fisher"

parameters$Ratio_threshold   <- 2
parameters$GO_max_top_terms  <- 10
parameters$GO_min_sig_genes  <- 2

GOenrichment(resDEG, data, parameters)

# -------------------------
# === Co-expression analysis ===
# -------------------------
parameters$coseq_data               <- "ExpressionProfiles"
parameters$coseq_ClustersNb         <- 15
parameters$coseq_HeatmapOrderSample <- TRUE

clust <- ClustAndGO(asko_norm, resDEG, parameters, data)
IncludeNonDEgenes_InClustering(data, asko_norm, resDEG, parameters, clust)

# -------------------------
# === Export DE results ===
# -------------------------
path <- "/path/to/DE_results_folder"
files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

read_results <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  comp_name <- tools::file_path_sans_ext(basename(file))
  df %>%
    select(is.gene, logFC, PValue, FDR) %>%
    rename(GeneID = is.gene,
           !!paste0("logFC_", comp_name) := logFC,
           !!paste0("Pval_", comp_name) := PValue,
           !!paste0("FDR_", comp_name) := FDR)
}

all_res <- lapply(files, read_results)
merged_res <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), all_res)

logFC_table <- merged_res %>% select(GeneID, starts_with("logFC_"))
pval_table  <- merged_res %>% select(GeneID, starts_with("Pval_"))
fdr_table   <- merged_res %>% select(GeneID, starts_with("FDR_"))

write.table(logFC_table, file = "merged_logFC.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pval_table,  file = "merged_PValue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fdr_table,   file = "merged_FDR.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
