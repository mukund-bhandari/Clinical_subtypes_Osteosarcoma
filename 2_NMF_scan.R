# =============================================================================
# STEP 2: NMF Rank Survey (K Scan)
# =============================================================================
# Purpose : Run NMF across a range of ranks (K = 2–6) to identify the optimal
#           number of metagenes/clusters using cophenetic correlation and
#           dispersion metrics.
# Input   : OS_TARGET_gene_filtered_1.xlsx
# Output  : OS_TARGET_gene_filtered_Kscan_1.rds
#           OS_TARGET_gene_filtered_Kscan_1.pdf
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(NMF)
library(dplyr)
library(openxlsx)

# --- Configuration ------------------------------------------------------------
WORKING_DIR  <- "/home/CBBI/bhandarim/NNMF"   # Update path as needed
INPUT_FILE   <- "OS_TARGET_gene_filtered_1.xlsx"
OUTPUT_RDS   <- "OS_TARGET_gene_filtered_Kscan_1.rds"
OUTPUT_PDF   <- "OS_TARGET_gene_filtered_Kscan_1.pdf"
NMF_SEED     <- 123456
NMF_NRUN     <- 100
K_RANGE      <- 2:6

# --- Set Working Directory ----------------------------------------------------
setwd(WORKING_DIR)

# --- 1. Load Filtered Expression Matrix ---------------------------------------
data_mat <- read.xlsx(INPUT_FILE, rowNames = TRUE)

# --- 2. Run NMF Rank Survey ---------------------------------------------------
# .opt = "vp25": verbose output, 25 parallel cores
nmf_scan <- nmf(data_mat, K_RANGE, seed = NMF_SEED, nrun = NMF_NRUN, .opt = "vp25")

# --- 3. Save NMF Result Object ------------------------------------------------
saveRDS(nmf_scan, OUTPUT_RDS)

# --- 4. Export Diagnostic Plots -----------------------------------------------
pdf(OUTPUT_PDF)
  plot(nmf_scan)          # Rank survey metrics (cophenetic, RSS, etc.)
  consensusmap(nmf_scan)  # Consensus matrices for each K
dev.off()

message("Step 2 complete. K-scan results saved to ", OUTPUT_RDS, " and ", OUTPUT_PDF)
