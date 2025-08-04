# ===============================================================================
# COMPLETE LDA ANALYSIS: ALL GENES vs BURSTEIN GENES
# ===============================================================================
# 
# Starting point: burstein_dataset object loaded
# Performs LDA classification using:
# 1. All common genes across datasets
# 2. 80 Burstein subtype marker genes
# Uses raw log-transformed data (no z-scoring)
# ===============================================================================

library(MASS)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

cat("=== COMPLETE LDA ANALYSIS ===\n\n")

# ===============================================================================
# STEP 1: PREPARE SAMPLE IDENTIFIERS
# ===============================================================================

cat("STEP 1: Preparing sample identifiers...\n")

# Get sample IDs with clinical data
tcga_samples <- intersect(colnames(tcga_common), 
                         combined_clinical$expression_id[combined_clinical$dataset == "TCGA"])
jiang_samples <- intersect(colnames(jiang_common), 
                          combined_clinical$expression_id[combined_clinical$dataset == "FUSCC"])

cat("Available samples:\n")
cat("  TCGA:", length(tcga_samples), "\n")
cat("  FUSCC:", length(jiang_samples), "\n")
cat("  Burstein:", ncol(burstein_dataset$expression_matrix), "\n\n")

# ===============================================================================
# STEP 2: PREPARE PATIENT DATA
# ===============================================================================

cat("STEP 2: Preparing patient data...\n")

# Log transform patient data
patient_log <- log2(expression_data_corrected$fpkm + 1)
names(patient_log) <- expression_data_corrected$gene_symbol

cat("Patient genes:", length(patient_log), "\n\n")

# ===============================================================================
# STEP 3: FIND COMMON GENES AND COMBINE DATASETS
# ===============================================================================

cat("STEP 3: Finding common genes and combining datasets...\n")

# Extract Burstein expression matrix
burstein_expression <- burstein_dataset$expression_matrix

# Find common genes across all datasets
all_common_genes <- intersect(intersect(intersect(rownames(tcga_common), 
                                                 rownames(jiang_common)), 
                                       rownames(burstein_expression)), 
                             names(patient_log))

cat("Gene availability:\n")
cat("  TCGA genes:", nrow(tcga_common), "\n")
cat("  FUSCC genes:", nrow(jiang_common), "\n")
cat("  Burstein genes:", nrow(burstein_expression), "\n")
cat("  Patient genes:", length(patient_log), "\n")
cat("  Common to all:", length(all_common_genes), "\n\n")

# Extract data for common genes (raw log data - NO Z-SCORING)
tcga_log <- tcga_common[all_common_genes, tcga_samples]
jiang_log <- jiang_common[all_common_genes, jiang_samples] 
burstein_log <- burstein_expression[all_common_genes, ]
patient_log_subset <- patient_log[all_common_genes]

# Combine all datasets
combined_log_raw <- cbind(tcga_log, jiang_log, burstein_log, patient_log_subset)
colnames(combined_log_raw)[ncol(combined_log_raw)] <- "Patient"

cat("Combined dataset:", nrow(combined_log_raw), "genes x", ncol(combined_log_raw), "samples\n\n")

# ===============================================================================
# STEP 4: PREPARE SAMPLE LABELS
# ===============================================================================

cat("STEP 4: Preparing sample labels...\n")

# Get Burstein metadata and clean subtype labels
metaburs <- burstein_dataset$sample_metadata
burstein_subtypes <- metaburs$`tnbc subtype:ch1`

# Extract abbreviations and map to match TCGA/FUSCC labels
burstein_subtypes_clean <- gsub(".*\\((.*)\\)", "\\1", burstein_subtypes)
burstein_subtypes_mapped <- burstein_subtypes_clean
burstein_subtypes_mapped[burstein_subtypes_mapped == "BLIA"] <- "IM"

cat("Subtype distributions:\n")
cat("TCGA/FUSCC:\n")
print(table(combined_clinical$mRNA_subtype))
cat("Burstein (mapped):\n")
print(table(burstein_subtypes_mapped))

# Create sample labels for combined dataset
all_sample_labels <- c(
  combined_clinical$mRNA_subtype[match(tcga_samples, combined_clinical$expression_id)],
  combined_clinical$mRNA_subtype[match(jiang_samples, combined_clinical$expression_id)],
  burstein_subtypes_mapped,
  NA  # Patient
)

# Create platform labels
all_platform_labels <- c(
  rep("TCGA", length(tcga_samples)),
  rep("FUSCC", length(jiang_samples)), 
  rep("Burstein", length(burstein_subtypes_mapped)),
  "Patient"
)

cat("\nCombined dataset labels:\n")
print(table(all_sample_labels, useNA = "ifany"))
cat("\n")

# ===============================================================================
# STEP 5: LDA WITH ALL COMMON GENES
# ===============================================================================

cat("STEP 5: LDA with ALL common genes...\n")

# Prepare data for LDA
training_mask_all <- !is.na(all_sample_labels)
lda_input_all <- t(combined_log_raw)  # samples in rows, genes in columns
training_data_all <- lda_input_all[training_mask_all, ]
training_subtypes_all <- all_sample_labels[training_mask_all]

cat("Training samples:", nrow(training_data_all), "\n")
cat("Training features (genes):", ncol(training_data_all), "\n")

# Train LDA
lda_all_genes <- lda(training_data_all, training_subtypes_all)
lda_all_pred <- predict(lda_all_genes, lda_input_all)

# Extract discriminant scores
discriminant_all <- data.frame(
  sample_id = rownames(lda_input_all),
  LD1 = lda_all_pred$x[, 1],
  LD2 = lda_all_pred$x[, 2],
  LD3 = lda_all_pred$x[, 3],
  platform = all_platform_labels,
  subtype = c(all_sample_labels[!is.na(all_sample_labels)], "Patient"),
  stringsAsFactors = FALSE
)

cat("Patient prediction (all genes):", as.character(lda_all_pred$class[nrow(lda_input_all)]), "\n")

# Patient scores
patient_all <- discriminant_all[discriminant_all$platform == "Patient", ]
cat("Patient discriminant scores (all genes):\n")
cat("  LD1:", round(patient_all$LD1, 3), "\n")
cat("  LD2:", round(patient_all$LD2, 3), "\n")
cat("  LD3:", round(patient_all$LD3, 3), "\n\n")

# ===============================================================================
# STEP 6: LDA WITH BURSTEIN GENES ONLY
# ===============================================================================

cat("STEP 6: LDA with Burstein marker genes only...\n")

# Define 80 Burstein marker genes
burstein_genes <- c("DHRS2", "GABRP", "AGR2", "PIP", "FOXA1", "PROM1", "TFF1", "NAT1", 
                   "BCL11A", "ESR1", "FOXC1", "CA12", "TFF3", "SCUBE2", "SFRP1", "ERBB4", 
                   "SIDT1", "PSAT1", "CHI3L1", "AR", "CD36", "OGN", "ABCA8", "CFD", 
                   "IGF1", "HBB", "CDH1", "MEOX2", "GPX3", "SCARA5", "PDK4", "ENPP2", 
                   "AGTR1", "LEPL", "PLDT", "TIMP4", "FHL1", "SRPX", "EDNRB", "SERPINB5", 
                   "SOX10", "IRX1", "MIA", "DSC2", "TTYH1", "COL9A3", "FGL2", "RARRES3", 
                   "PDE9A", "BST2", "PTGER4", "KCNK5", "PSMB9", "HLA-DMA", "EPHB3", 
                   "IGSF6", "ST3GAL6", "RHO", "CXCL9", "CXCL11", "GBP5", "GZMB", 
                   "LAMP3", "GBP1", "ADAMDEC1", "CCL5", "SPON1", "PBK", "STAT1", 
                   "EZH2", "PLAT", "TAP2", "SLAMF7", "HERC5", "SPOCK1", "TAP1", "CD2", "AIM2")

# Check availability and filter
available_burstein_genes <- intersect(burstein_genes, rownames(combined_log_raw))
cat("Burstein genes available:", length(available_burstein_genes), "out of", length(burstein_genes), "\n")

# Filter combined data to Burstein genes only
combined_log_burstein <- combined_log_raw[available_burstein_genes, ]

# Prepare data for Burstein LDA
lda_input_burstein <- t(combined_log_burstein)
training_data_burstein <- lda_input_burstein[training_mask_all, ]

cat("Training samples:", nrow(training_data_burstein), "\n")
cat("Training features (Burstein genes):", ncol(training_data_burstein), "\n")

# Train LDA on Burstein genes
lda_burstein_genes <- lda(training_data_burstein, training_subtypes_all)
lda_burstein_pred <- predict(lda_burstein_genes, lda_input_burstein)

# Extract discriminant scores
discriminant_burstein <- data.frame(
  sample_id = rownames(lda_input_burstein),
  LD1 = lda_burstein_pred$x[, 1],
  LD2 = lda_burstein_pred$x[, 2],
  LD3 = lda_burstein_pred$x[, 3],
  platform = all_platform_labels,
  subtype = c(all_sample_labels[!is.na(all_sample_labels)], "Patient"),
  stringsAsFactors = FALSE
)

cat("Patient prediction (Burstein genes):", as.character(lda_burstein_pred$class[nrow(lda_input_burstein)]), "\n")

# Patient scores
patient_burstein <- discriminant_burstein[discriminant_burstein$platform == "Patient", ]
cat("Patient discriminant scores (Burstein genes):\n")
cat("  LD1:", round(patient_burstein$LD1, 3), "\n")
cat("  LD2:", round(patient_burstein$LD2, 3), "\n")
cat("  LD3:", round(patient_burstein$LD3, 3), "\n\n")

# ===============================================================================
# STEP 7: CREATE VISUALIZATIONS
# ===============================================================================

cat("STEP 7: Creating visualizations...\n")

# Set up colors
unique_subtypes <- unique(discriminant_all$subtype[discriminant_all$subtype != "Patient"])
n_subtypes <- length(unique_subtypes)

if (n_subtypes <= 8) {
  subtype_colors <- brewer.pal(max(3, n_subtypes), "Set2")
} else {
  subtype_colors <- rainbow(n_subtypes)
}
names(subtype_colors) <- unique_subtypes
subtype_colors["Patient"] <- "red"

# Set up shapes
platform_shapes <- c("TCGA" = 16, "FUSCC" = 17, "Burstein" = 15, "Patient" = 8)

# Create output folder
output_folder <- "complete_lda_analysis"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

pdf(file.path(output_folder, "lda_comparison_analysis.pdf"), width = 16, height = 12)

# Plot 1: All genes LDA
p1 <- ggplot(discriminant_all, aes(x = LD1, y = LD2, color = subtype, shape = platform)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = subtype_colors) +
  scale_shape_manual(values = platform_shapes) +
  labs(title = "LDA: All Common Genes",
       subtitle = paste("Patient predicted as:", lda_all_pred$class[nrow(lda_input_all)]),
       x = "LD1", y = "LD2") +
  theme_bw() +
  theme(legend.position = "bottom")

# Plot 2: Burstein genes LDA
p2 <- ggplot(discriminant_burstein, aes(x = LD1, y = LD2, color = subtype, shape = platform)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = subtype_colors) +
  scale_shape_manual(values = platform_shapes) +
  labs(title = "LDA: Burstein Marker Genes",
       subtitle = paste("Patient predicted as:", lda_burstein_pred$class[nrow(lda_input_burstein)]),
       x = "LD1", y = "LD2") +
  theme_bw() +
  theme(legend.position = "bottom")

# Combined plot
library(gridExtra)
grid.arrange(p1, p2, ncol = 2, top = "LDA Classification Comparison")

# Individual plots
print(p1)
print(p2)

dev.off()

# ===============================================================================
# STEP 8: SAVE RESULTS AND SUMMARY
# ===============================================================================

cat("STEP 8: Saving results...\n")

# Save comprehensive results
save(lda_all_genes, lda_burstein_genes, 
     discriminant_all, discriminant_burstein,
     lda_all_pred, lda_burstein_pred,
     combined_log_raw, combined_log_burstein,
     all_sample_labels, all_platform_labels,
     available_burstein_genes, all_common_genes,
     file = file.path(output_folder, "complete_lda_results.RData"))

# Save coordinate files
write.csv(discriminant_all, file.path(output_folder, "lda_all_genes_coordinates.csv"), row.names = FALSE)
write.csv(discriminant_burstein, file.path(output_folder, "lda_burstein_genes_coordinates.csv"), row.names = FALSE)

# Create comparison summary
sink(file.path(output_folder, "lda_comparison_report.txt"))

cat("=================================================================\n")
cat("COMPLETE LDA ANALYSIS COMPARISON REPORT\n")
cat("=================================================================\n\n")

cat("Analysis Date:", as.character(Sys.Date()), "\n\n")

cat("DATASETS INCLUDED:\n")
cat("  TCGA samples:", length(tcga_samples), "\n")
cat("  FUSCC samples:", length(jiang_samples), "\n")
cat("  Burstein samples:", ncol(burstein_expression), "\n")
cat("  Patient: 1 sample\n")
cat("  Total samples:", ncol(combined_log_raw), "\n\n")

cat("GENE SETS:\n")
cat("  All common genes:", length(all_common_genes), "\n")
cat("  Burstein marker genes available:", length(available_burstein_genes), "out of 80\n\n")

cat("CLASSIFICATION RESULTS:\n")
cat("All genes LDA:\n")
cat("  Patient prediction:", as.character(lda_all_pred$class[nrow(lda_input_all)]), "\n")
cat("  LD1 score:", round(patient_all$LD1, 3), "\n")
cat("  LD2 score:", round(patient_all$LD2, 3), "\n")
cat("  LD3 score:", round(patient_all$LD3, 3), "\n\n")

cat("Burstein genes LDA:\n")
cat("  Patient prediction:", as.character(lda_burstein_pred$class[nrow(lda_input_burstein)]), "\n")
cat("  LD1 score:", round(patient_burstein$LD1, 3), "\n")
cat("  LD2 score:", round(patient_burstein$LD2, 3), "\n")
cat("  LD3 score:", round(patient_burstein$LD3, 3), "\n\n")

cat("CONSISTENCY CHECK:\n")
all_prediction <- as.character(lda_all_pred$class[nrow(lda_input_all)])
burstein_prediction <- as.character(lda_burstein_pred$class[nrow(lda_input_burstein)])

if (all_prediction == burstein_prediction) {
  cat("  CONSISTENT: Both methods predict", all_prediction, "\n")
} else {
  cat("  INCONSISTENT: All genes =", all_prediction, ", Burstein genes =", burstein_prediction, "\n")
}

cat("\nMETHOD NOTES:\n")
cat("  - Raw log2-transformed data (no z-scoring)\n")
cat("  - Combined TCGA, FUSCC, and Burstein datasets\n")
cat("  - 4-group classification: BLIS, IM, LAR, MES\n")
cat("  - Patient compared to population-based subtypes\n")

sink()

cat("Complete LDA analysis finished!\n")
cat("Files created in folder:", output_folder, "\n")
cat("  - lda_comparison_analysis.pdf\n")
cat("  - lda_all_genes_coordinates.csv\n")
cat("  - lda_burstein_genes_coordinates.csv\n")
cat("  - complete_lda_results.RData\n")
cat("  - lda_comparison_report.txt\n")

# ===============================================================================
# FINAL SUMMARY
# ===============================================================================

cat("\n=== FINAL SUMMARY ===\n")
cat("Patient Classification Results:\n")
cat("  All genes:", as.character(lda_all_pred$class[nrow(lda_input_all)]), "\n")
cat("  Burstein genes:", as.character(lda_burstein_pred$class[nrow(lda_input_burstein)]), "\n")

if (all_prediction == burstein_prediction) {
  cat("  CONSENSUS PREDICTION:", all_prediction, "\n")
} else {
  cat("  NO CONSENSUS - further investigation needed\n")
}

cat("=== ANALYSIS COMPLETE ===\n")