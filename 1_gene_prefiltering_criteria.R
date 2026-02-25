# =============================================================================
# STEP 1: Load RNA-seq Data and Pre-filter Genes
# =============================================================================
# Purpose : Remove lowly expressed and ubiquitously high-expressed genes
#           before NMF analysis.
# Input   : TARGET-OS.htseq_fpkm.tsv
# Output  : OS_TARGET_gene_filtered_1.xlsx
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readr)
library(readxl)
library(dplyr)
library(openxlsx)
library(NMF)

# --- 1. Load Expression Data --------------------------------------------------
TARGET_OS_rna_seq <- read_tsv("TARGET-OS.htseq_fpkm.tsv")

# Trim Ensembl IDs to 15 characters (remove version suffix, e.g. ".1")
TARGET_OS_rna_seq$ENSGID <- substr(TARGET_OS_rna_seq$ENSGID, 1, 15)
# NOTE: If needed, match ENSGID to gene symbols and remove duplicate genes here.

# --- 2. Convert to Data Frame and Set Row Names --------------------------------
T_OS_df <- data.frame(TARGET_OS_rna_seq)  # 60,483 genes
rownames(T_OS_df) <- T_OS_df$ENSGID
T_OS_df <- T_OS_df[, -1]                  # Remove the ENSGID column

# --- 3. Filter 1: Remove Lowly Expressed Genes --------------------------------
# Keep genes with FPKM > 1 in at least 22 samples (~25% of samples)
T_OS_df_filtered <- T_OS_df[rowSums(T_OS_df > 1) >= 22, ]

# --- 4. Filter 2: Remove Ubiquitously High-Expressed Genes -------------------
# Exclude genes with total rowSum >= 280,000 (e.g. mitochondrial housekeeping genes)
T_OS_df_filtered <- T_OS_df_filtered[rowSums(T_OS_df_filtered) < 280000, ]

# --- 5. Save Filtered Data ----------------------------------------------------
write.xlsx(T_OS_df_filtered, "OS_TARGET_gene_filtered_1.xlsx", rowNames = TRUE)

message("Step 1 complete. Filtered gene matrix saved to OS_TARGET_gene_filtered_1.xlsx")
