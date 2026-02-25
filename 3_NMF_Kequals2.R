# =============================================================================
# STEP 3: NMF at Selected Rank (K = 2)
# =============================================================================
# Purpose : Run NMF at the optimal rank chosen from the K scan (K = 2) to
#           derive final sample cluster assignments and metagene signatures.
# Input   : OS_TARGET_gene_filtered_1.xlsx
# Output  : OS_TARGET_gene_filtered_K2_1.rds
#           OS_TARGET_gene_filtered_K2_1.pdf
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(NMF)
library(dplyr)
library(openxlsx)

# --- Configuration ------------------------------------------------------------
WORKING_DIR  <- "/home/CBBI/bhandarim/NNMF"   # Update path as needed
INPUT_FILE   <- "OS_TARGET_gene_filtered_1.xlsx"
OUTPUT_RDS   <- "OS_TARGET_gene_filtered_K2_1.rds"
OUTPUT_PDF   <- "OS_TARGET_gene_filtered_K2_1.pdf"
NMF_SEED     <- 123456
NMF_NRUN     <- 100
K_SELECTED   <- 2          # Change this value if a different K is chosen

# --- Set Working Directory ----------------------------------------------------
setwd(WORKING_DIR)

# --- 1. Load Filtered Expression Matrix ---------------------------------------
data_mat <- read.xlsx(INPUT_FILE, rowNames = TRUE)

# --- 2. Run NMF at Selected Rank ----------------------------------------------
# .opt = "vp25": verbose output, 25 parallel cores
nmf_result <- nmf(data_mat, K_SELECTED, seed = NMF_SEED, nrun = NMF_NRUN, .opt = "vp25")

# --- 3. Save NMF Result Object ------------------------------------------------
saveRDS(nmf_result, OUTPUT_RDS)

# --- 4. Export Diagnostic Plots -----------------------------------------------
pdf(OUTPUT_PDF)
  plot(nmf_result)          # Model fit metrics
  consensusmap(nmf_result)  # Consensus map for K = 2
dev.off()

message("Step 3 complete. NMF K=", K_SELECTED, " result saved to ", OUTPUT_RDS, " and ", OUTPUT_PDF)
