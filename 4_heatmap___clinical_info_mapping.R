# =============================================================================
# STEP 4: Signature Gene Identification and Annotated Heatmap
# =============================================================================
# Purpose : Identify differentially expressed signature genes for each NMF
#           cluster (Group A vs B), then generate an annotated heatmap with
#           clinical metadata.
# Input   : OS_TARGET_gene_filtered_1.xlsx
#           OS_gene_list.xlsx  (sheets 3, 10, 13)
#           TargetOS_Genes_withEnsembelIDS.txt
# Output  : Sig_genes_OS_2grp.xlsx
#           Sig_gene_names_OS_2grp.xlsx
#           Heatmap plots (drawn to current device)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readr)
library(readxl)
library(dplyr)
library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(grid)   # for gpar() / unit()

# =============================================================================
# PART A: Identify Signature Genes
# =============================================================================

# --- 1. Load Filtered Expression Data ----------------------------------------
exp_df <- as.data.frame(read.xlsx("OS_TARGET_gene_filtered_1.xlsx", rowNames = TRUE))

# --- 2. Load Group Assignments and Order Samples Accordingly -----------------
# Sheet 10: TARGET_ID in column 1, Group label (A/B) in subsequent columns
group_data <- as.data.frame(read_excel("OS_gene_list.xlsx", sheet = 10))

exp_df_ordered <- exp_df[, order(match(colnames(exp_df), group_data[, 1]))]

# --- 3. Define Group Column Indices ------------------------------------------
# Group B: columns 1–46  |  Group A: columns 47–88
# Update these ranges if sample sizes change.
idx_B <- 1:46
idx_A <- 47:88

# --- Helper: Differential Expression Analysis --------------------------------
# Returns a filtered data frame of genes with |log2FC| >= 1 and FDR < 0.05
compute_DE_genes <- function(expr_mat, idx_num, idx_denom, fc_col, fdr_col) {
  mean_num   <- apply(expr_mat[, idx_num],   1, mean)
  mean_denom <- apply(expr_mat[, idx_denom], 1, mean)

  log2_fc  <- log2(mean_num / mean_denom)
  p_values <- apply(expr_mat, 1, function(x) t.test(x[idx_num], x[idx_denom])$p.value)
  p_adj    <- p.adjust(p_values, method = "hochberg")

  result_df <- cbind(
    expr_mat,
    setNames(data.frame(p_values), paste0("pval_", fc_col)),
    setNames(data.frame(p_adj),    fdr_col),
    setNames(data.frame(log2_fc),  fc_col)
  )

  result_df %>%
    dplyr::filter(.data[[fc_col]] >= 1,
                  .data[[fdr_col]] < 0.05)
}

# --- 4. Find Group B Signature Genes (B vs A) --------------------------------
gene_matrix_B <- compute_DE_genes(exp_df_ordered, idx_B, idx_A,
                                   fc_col  = "log2_BvsA",
                                   fdr_col = "FDR_BvsA")

heatmap(data.matrix(gene_matrix_B[, seq_len(length(idx_B) + length(idx_A))]),
        main = "Group B Signature Genes")

genes_B <- data.frame(genes = rownames(gene_matrix_B))

# --- 5. Find Group A Signature Genes (A vs B) --------------------------------
gene_matrix_A <- compute_DE_genes(exp_df_ordered, idx_A, idx_B,
                                   fc_col  = "log2_AvsB",
                                   fdr_col = "FDR_AvsB")

heatmap(data.matrix(gene_matrix_A[, seq_len(length(idx_B) + length(idx_A))]),
        main = "Group A Signature Genes")

genes_A <- data.frame(genes = rownames(gene_matrix_A))

# --- 6. Combine and Save Signature Genes -------------------------------------
sig_genes <- rbind(genes_A, genes_B)
write.xlsx(sig_genes, "Sig_genes_OS_2grp.xlsx", rowNames = FALSE)

message("Part A complete. Signature genes saved to Sig_genes_OS_2grp.xlsx")

# =============================================================================
# PART B: Annotated Heatmap with Clinical Metadata
# =============================================================================

# --- 7. Subset Expression Data to Signature Genes ----------------------------
selected_exp_df <- exp_df_ordered[rownames(exp_df_ordered) %in% sig_genes$genes, ]

# --- 8. Transpose and Prepare Sample-Level Data Frame ------------------------
selected_exp_t <- as.data.frame(t(selected_exp_df))
rownames(selected_exp_t) <- substr(rownames(selected_exp_t), 1, 16)
selected_exp_t$TARGET_ID  <- rownames(selected_exp_t)

# --- 9. Load and Harmonise Clinical Metadata ---------------------------------
# Replace "-" with "." in TARGET_IDs to match expression matrix naming

load_clinical_sheet <- function(file, sheet, cols = NULL) {
  df <- as.data.frame(read_excel(file, sheet = sheet))
  if (!is.null(cols)) df <- df[, cols]
  df$TARGET_ID <- gsub("\\-", ".", df$TARGET_ID)
  df
}

status_data       <- load_clinical_sheet("OS_gene_list.xlsx", sheet = 3)
grp_data          <- load_clinical_sheet("OS_gene_list.xlsx", sheet = 10)
grp_data$TARGET_ID <- substr(rownames(selected_exp_t), 1, 16)   # align IDs

gender_data       <- load_clinical_sheet("OS_gene_list.xlsx", sheet = 13, cols = c(2, 4))
race_ethn_data    <- load_clinical_sheet("OS_gene_list.xlsx", sheet = 13, cols = c(2, 5, 6))
meta_site_data    <- load_clinical_sheet("OS_gene_list.xlsx", sheet = 13, cols = c(2, 10, 16, 17, 18))

# --- 10. Join All Metadata to Transposed Expression Data ---------------------
selected_exp_t <- selected_exp_t %>%
  inner_join(status_data,    by = "TARGET_ID") %>%
  inner_join(grp_data,       by = "TARGET_ID") %>%
  inner_join(gender_data,    by = "TARGET_ID") %>%
  inner_join(race_ethn_data, by = "TARGET_ID") %>%
  inner_join(meta_site_data, by = "TARGET_ID")

rownames(selected_exp_t) <- selected_exp_t$TARGET_ID

# --- 11. Separate Expression Columns from Metadata ---------------------------
# Metadata occupies the last 10 columns (columns 62–71); adjust if needed.
n_meta_cols <- 10
expr_cols   <- seq_len(ncol(selected_exp_t) - n_meta_cols)
tidy_expr   <- as.data.frame(selected_exp_t[, expr_cols])

# --- 12. Scale and Transpose for Heatmap -------------------------------------
tidy_scaled    <- scale(tidy_expr, center = TRUE, scale = TRUE)
factors_matrix <- t(tidy_scaled)

# --- 13. Reorder Rows to Match Signature Gene Order --------------------------
factors_ordered <- factors_matrix[
  order(match(rownames(factors_matrix), sig_genes[, 1])), 
]

# --- 14. Map Ensembl IDs to Gene Symbols -------------------------------------
sig_gene_mat  <- read.xlsx("Sig_genes_OS_2grp.xlsx")
colnames(sig_gene_mat) <- "ENSG"

gene_list <- read.delim2("TargetOS_Genes_withEnsembelIDS.txt")
colnames(gene_list) <- c("HGNC", "ENSG")

sig_gene_annotated <- left_join(sig_gene_mat, gene_list, by = "ENSG")
write.xlsx(sig_gene_annotated, "Sig_gene_names_OS_2grp.xlsx", rowNames = FALSE)

# --- 15. Extract Annotation Vectors ------------------------------------------
Status            <- selected_exp_t$STATUS
Group             <- selected_exp_t$Group
Gender            <- selected_exp_t$Gender
Race              <- selected_exp_t$Race
Ethnicity         <- selected_exp_t$Ethnicity
Metastasis_site   <- selected_exp_t$`Metastasis site`
Primary_tumor_site  <- selected_exp_t$`Primary tumor site`
Specific_tumor_site <- selected_exp_t$`Specific tumor site`
Vital_Status      <- selected_exp_t$`Vital Status`

# --- 16. Define Colour Scale -------------------------------------------------
col_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#FFFFFF", "#CD3333"))

# --- 17. Build Top Annotation (Clinical / Disease) ---------------------------
annotation_legend_params <- list(
  ncol          = 1,
  labels_gp     = gpar(fontsize = 5),
  grid_width    = unit(0.2, "cm"),
  title_gp      = gpar(fontsize = 4),
  legend_height = unit(0.2, "cm"),
  legend_width  = unit(0.3, "cm"),
  direction     = "horizontal",
  border        = "black"
)

ha_top <- HeatmapAnnotation(
  Group               = Group,
  Status              = Status,
  Vital_Status        = Vital_Status,
  Metastasis_site     = Metastasis_site,
  Primary_tumor_site  = Primary_tumor_site,
  Specific_tumor_site = Specific_tumor_site,
  col = list(
    Group   = c("A" = "darkred",    "B" = "darkblue"),
    Status  = c("Metastatic" = "firebrick2", "Non-metastatic" = "tan1"),
    Metastasis_site = c(
      "Bone and lung" = "tan",
      "Bone only"     = "white",
      "Lung only"     = "pink"
    ),
    Primary_tumor_site = c(
      "Arm/hand" = "magenta",
      "Leg/Foot" = "tan",
      "Pelvis"   = "pink"
    ),
    Vital_Status = c("Alive" = "forestgreen", "Dead" = "firebrick2")
  ),
  annotation_name_gp       = gpar(fontsize = 6),
  annotation_legend_param  = annotation_legend_params
)

# --- 18. Build Bottom Annotation (Demographics) ------------------------------
ha_bottom <- HeatmapAnnotation(
  Group     = Group,
  Gender    = Gender,
  Race      = Race,
  Ethnicity = Ethnicity,
  col = list(
    Group     = c("A" = "darkred", "B" = "darkblue"),
    Gender    = c("Male" = "black",  "Female" = "pink"),
    Race      = c(
      "Asian"                      = "tan",
      "Black or African American"  = "black",
      "Unknown"                    = "cornsilk",
      "White"                      = "white"
    ),
    Ethnicity = c(
      "Hispanic or Latino"     = "bisque4",
      "Not Hispanic or Latino" = "azure",
      "Unknown"                = "cornsilk"
    )
  ),
  annotation_name_gp      = gpar(fontsize = 6),
  annotation_legend_param = annotation_legend_params
)

# --- 19. Shared Heatmap Parameters -------------------------------------------
heatmap_params <- list(
  cluster_rows          = FALSE,
  cluster_row_slices    = TRUE,
  cluster_columns       = FALSE,
  clustering_distance_rows = "pearson",
  clustering_method_rows   = "centroid",
  show_row_names        = TRUE,
  show_column_names     = TRUE,
  show_row_dend         = TRUE,
  show_column_dend      = FALSE,
  name                  = "Scaled_Exp(FPKM)",
  col                   = col_fun,
  width                 = unit(15, "cm"),
  height                = unit(13, "cm"),
  row_title_gp          = gpar(fontsize = 14),
  column_names_gp       = gpar(fontsize = 5, fontface = "bold"),
  row_names_gp          = gpar(fontsize = 4, fontface = "bold"),
  top_annotation        = ha_top,
  bottom_annotation     = ha_bottom,
  heatmap_legend_param  = list(
    labels_gp     = gpar(fontsize = 7),
    grid_width    = unit(2, "cm"),
    title_gp      = gpar(fontsize = 10),
    legend_height = unit(1, "cm"),
    legend_width  = unit(2, "cm"),
    direction     = "horizontal",
    border        = "black"
  )
)

# --- 20. Draw Full Heatmap (All Signature Genes) -----------------------------
H_full <- do.call(Heatmap, c(
  list(
    as.matrix(factors_ordered),
    row_labels = sig_gene_annotated[, 2]
  ),
  heatmap_params
))

draw(H_full, heatmap_legend_side = "top", annotation_legend_side = "right")

# --- 21. Draw Refined Heatmap (Remove Rows 24 & 25) -------------------------
# These two rows were identified as outliers / artefacts — remove and redraw.
rows_to_remove      <- c(24, 25)
sig_gene_annotated2 <- sig_gene_annotated[-rows_to_remove, ]

H_refined <- do.call(Heatmap, c(
  list(
    as.matrix(factors_ordered[-rows_to_remove, ]),
    row_labels = sig_gene_annotated2[, 2]
  ),
  heatmap_params
))

draw(H_refined, heatmap_legend_side = "top", annotation_legend_side = "right")

message("Step 4 complete. Heatmaps drawn and gene name table saved to Sig_gene_names_OS_2grp.xlsx")
