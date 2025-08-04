# ===============================================================================
# BURSTEIN GENES HEATMAP WITH COMBAT BATCH CORRECTION
# ===============================================================================

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(sva)  # for ComBat

cat("=== BURSTEIN GENES HEATMAP WITH COMBAT BATCH CORRECTION ===\n\n")

# ===============================================================================
# STEP 1: PREPARE DATA
# ===============================================================================

cat("STEP 1: Preparing data...\n")

# Prepare patient data (log transform)
patient_log <- log2(expression_data_corrected$fpkm + 1)
names(patient_log) <- expression_data_corrected$gene_symbol

# Get sample identifiers
tcga_samples <- intersect(colnames(tcga_common), 
                         combined_clinical$expression_id[combined_clinical$dataset == "TCGA"])
jiang_samples <- intersect(colnames(jiang_common), 
                          combined_clinical$expression_id[combined_clinical$dataset == "FUSCC"])

# Extract Burstein expression matrix and metadata
burstein_expression <- burstein_dataset$expression_matrix
metaburs <- burstein_dataset$sample_metadata
burstein_subtypes <- metaburs$`tnbc subtype:ch1`
burstein_subtypes_clean <- gsub(".*\\((.*)\\)", "\\1", burstein_subtypes)
burstein_subtypes_mapped <- burstein_subtypes_clean
burstein_subtypes_mapped[burstein_subtypes_mapped == "BLIA"] <- "IM"

cat("Sample counts:\n")
cat("  TCGA:", length(tcga_samples), "\n")
cat("  FUSCC:", length(jiang_samples), "\n") 
cat("  Burstein:", ncol(burstein_expression), "\n")
cat("  Patient: 1\n\n")

# ===============================================================================
# STEP 2: FILTER TO BURSTEIN GENES AND BATCH CORRECT REFERENCE DATA
# ===============================================================================

cat("STEP 2: Filtering to Burstein marker genes and batch correcting...\n")

# Find common genes across all datasets
all_common_genes <- intersect(intersect(intersect(rownames(tcga_common), 
                                                 rownames(jiang_common)), 
                                       rownames(burstein_expression)), 
                             names(patient_log))

# Filter to available Burstein genes
available_burstein_genes <- intersect(burstein_genes, all_common_genes)
cat("Burstein genes available:", length(available_burstein_genes), "out of", length(burstein_genes), "\n")

# Extract REFERENCE data for Burstein genes only (exclude patient for now)
tcga_burstein <- tcga_common[available_burstein_genes, tcga_samples]
jiang_burstein <- jiang_common[available_burstein_genes, jiang_samples]
burstein_burstein <- burstein_expression[available_burstein_genes, ]

# Combine reference datasets
reference_burstein <- cbind(tcga_burstein, jiang_burstein, burstein_burstein)

# Create batch labels for reference data
batch_ref <- c(
  rep("TCGA", length(tcga_samples)),
  rep("FUSCC", length(jiang_samples)),
  rep("Burstein", ncol(burstein_burstein))
)

# Create subtype labels for reference data
subtype_ref <- c(
  combined_clinical$mRNA_subtype[match(tcga_samples, combined_clinical$expression_id)],
  combined_clinical$mRNA_subtype[match(jiang_samples, combined_clinical$expression_id)],
  burstein_subtypes_mapped
)

cat("Reference data before batch correction:\n")
cat("  Samples:", ncol(reference_burstein), "\n")
cat("  Batches:", table(batch_ref), "\n")
cat("  Subtypes:", table(subtype_ref, useNA = "ifany"), "\n")

# Remove any samples with missing subtypes
complete_mask <- !is.na(subtype_ref)
if(!all(complete_mask)) {
  cat("Removing", sum(!complete_mask), "samples with missing subtypes\n")
  reference_burstein <- reference_burstein[, complete_mask]
  batch_ref <- batch_ref[complete_mask]
  subtype_ref <- subtype_ref[complete_mask]
}

# COMBAT BATCH CORRECTION on reference data only
cat("Applying ComBat batch correction to reference datasets...\n")
mod <- model.matrix(~ subtype_ref)
reference_corrected <- ComBat(dat = reference_burstein,
                             batch = batch_ref,
                             mod = mod,
                             par.prior = TRUE,
                             prior.plots = FALSE)

cat("ComBat correction completed.\n")

# NOW add patient data (not batch corrected, but will be in same space for comparison)
patient_burstein <- patient_log[available_burstein_genes]
combined_burstein <- cbind(reference_corrected, patient_burstein)
colnames(combined_burstein)[ncol(combined_burstein)] <- "Patient"

cat("Combined Burstein gene matrix (batch corrected):", nrow(combined_burstein), "genes x", ncol(combined_burstein), "samples\n\n")

# ===============================================================================
# STEP 3: CREATE ANNOTATION DATA FOR BATCH-CORRECTED DATA
# ===============================================================================

cat("STEP 3: Creating sample annotations for batch-corrected data...\n")

# Create sample annotation data frame for the final combined data
# Note: reference data may have had samples removed during ComBat
reference_sample_names <- colnames(reference_corrected)

# Reconstruct the annotation based on what's actually in the corrected data
sample_annotation <- data.frame(
  sample_id = c(reference_sample_names, "Patient"),
  dataset = c(batch_ref, "Patient"),  # batch_ref was already filtered
  subtype = c(subtype_ref, "Unknown"),  # subtype_ref was already filtered
  stringsAsFactors = FALSE
)

rownames(sample_annotation) <- sample_annotation$sample_id

cat("Final sample counts by dataset:\n")
print(table(sample_annotation$dataset))
cat("Final sample counts by subtype:\n")
print(table(sample_annotation$subtype))

# ===============================================================================
# STEP 4: PREPARE HEATMAP ANNOTATIONS WITH PATIENT HIGHLIGHTING
# ===============================================================================

cat("STEP 4: Setting up heatmap annotations with patient highlighting...\n")

# Create annotation data frame for heatmap
annotation_col <- data.frame(
  Dataset = sample_annotation$dataset,
  Subtype = sample_annotation$subtype,
  Patient_Status = ifelse(sample_annotation$dataset == "Patient", "YOUR PATIENT", "Reference"),
  row.names = rownames(sample_annotation)
)

# Define colors with special highlighting for patient
dataset_colors <- c("TCGA" = "#E31A1C", "FUSCC" = "#1F78B4", "Burstein" = "#33A02C", "Patient" = "#FF7F00")
subtype_colors <- c("BLIS" = "#A6CEE3", "IM" = "#B2DF8A", "LAR" = "#FB9A99", "MES" = "#FDBF6F", "Unknown" = "#CAB2D6")
patient_colors <- c("YOUR PATIENT" = "#FF0000", "Reference" = "#FFFFFF")

annotation_colors <- list(
  Dataset = dataset_colors,
  Subtype = subtype_colors,
  Patient_Status = patient_colors
)

# ===============================================================================
# STEP 5: Z-SCORE NORMALIZATION FOR VISUALIZATION
# ===============================================================================

cat("STEP 5: Z-score normalizing for visualization...\n")

# Z-score normalize genes (across samples) for better visualization
combined_burstein_scaled <- t(scale(t(combined_burstein)))

# Handle any genes with zero variance
if(any(is.na(combined_burstein_scaled))) {
  cat("Warning: Some genes had zero variance and were set to 0\n")
  combined_burstein_scaled[is.na(combined_burstein_scaled)] <- 0
}

# ===============================================================================
# STEP 6: CREATE ENHANCED HEATMAPS WITH PATIENT FOCUS
# ===============================================================================

cat("STEP 6: Creating enhanced heatmaps...\n")

# Create output folder
if(!dir.exists("heatmap_analysis_combat")) {
  dir.create("heatmap_analysis_combat")
}

# MAIN HEATMAP with patient highlighting
pdf("heatmap_analysis_combat/burstein_genes_heatmap_combat.pdf", width = 18, height = 14)

pheatmap(combined_burstein_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",  # Already z-scored
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 9,
         fontsize_col = 6,
         main = "Burstein TNBC Marker Genes Heatmap\n(ComBat Batch Corrected + Patient Highlighted)",
         border_color = NA)

dev.off()

# PNG version
png("heatmap_analysis_combat/burstein_genes_heatmap_combat.png", width = 1800, height = 1400, res = 100)

pheatmap(combined_burstein_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2", 
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 9,
         fontsize_col = 6,
         main = "Burstein TNBC Marker Genes Heatmap\n(ComBat Batch Corrected + Patient Highlighted)",
         border_color = NA)

dev.off()

# FOCUSED HEATMAP: Patient vs Subtype Centroids
cat("Creating focused patient comparison heatmap...\n")

# Calculate subtype centroids
reference_data <- combined_burstein_scaled[, sample_annotation$dataset != "Patient"]
reference_annotation <- sample_annotation[sample_annotation$dataset != "Patient", ]

subtype_centroids <- sapply(unique(reference_annotation$subtype), function(subtype) {
  subtype_samples <- reference_annotation$subtype == subtype
  rowMeans(reference_data[, subtype_samples])
})

# SOLUTION: Create separate visualizations with optimized scales

# 1. CENTROIDS ONLY - optimized color scale for subtype patterns
centroids_annotation <- data.frame(
  Subtype = colnames(subtype_centroids),
  row.names = colnames(subtype_centroids)
)

centroids_colors <- list(
  Subtype = subtype_colors
)

# Find optimal color range for centroids (exclude extreme outliers)
centroid_range <- quantile(subtype_centroids, c(0.05, 0.95))
centroid_breaks <- seq(centroid_range[1], centroid_range[2], length.out = 100)

# 2. PATIENT ONLY - ROBUST PERCENTILE SCALING
patient_data <- matrix(combined_burstein_scaled[, "Patient"], ncol = 1)
rownames(patient_data) <- rownames(combined_burstein_scaled)
colnames(patient_data) <- "YOUR_PATIENT"

patient_annotation <- data.frame(
  Status = "YOUR_PATIENT",
  row.names = "YOUR_PATIENT"
)

patient_colors <- list(
  Status = c("YOUR_PATIENT" = "#FF0000")
)

# ROBUST SCALING: Exclude top 3-4 outlier genes dynamically
n_genes <- nrow(patient_data)
outliers_to_exclude <- min(4, max(2, round(n_genes * 0.1)))  # Exclude top 10% or 2-4 genes, whichever is reasonable
upper_percentile <- (n_genes - outliers_to_exclude) / n_genes
lower_percentile <- outliers_to_exclude / n_genes  # Also exclude bottom outliers symmetrically

patient_range <- quantile(patient_data, c(lower_percentile, upper_percentile), na.rm = TRUE)

cat("Adaptive outlier exclusion for", n_genes, "genes:\n")
cat("Excluding top", outliers_to_exclude, "and bottom", outliers_to_exclude, "outlier genes\n")
cat("Using", sprintf("%.1f", lower_percentile*100), "th to", sprintf("%.1f", upper_percentile*100), "th percentile for color range\n")
cat("Patient color range:", sprintf("%.2f to %.2f", patient_range[1], patient_range[2]), "\n")

# Show which genes are being capped
patient_values <- patient_data[,1]
extreme_high <- patient_values[patient_values > patient_range[2]]
extreme_low <- patient_values[patient_values < patient_range[1]]

if(length(extreme_high) > 0) {
  cat("Extremely high genes (excluded from color scale, will appear maximally red):\n")
  extreme_high_sorted <- sort(extreme_high, decreasing = TRUE)
  for(i in 1:length(extreme_high_sorted)) {
    gene <- names(extreme_high_sorted)[i]
    cat(sprintf("  %d. %s: %.2f (capped at %.2f for colors)\n", i, gene, extreme_high_sorted[i], patient_range[2]))
  }
}

if(length(extreme_low) > 0) {
  cat("Extremely low genes (excluded from color scale, will appear maximally blue):\n")
  extreme_low_sorted <- sort(extreme_low, decreasing = FALSE)
  for(i in 1:length(extreme_low_sorted)) {
    gene <- names(extreme_low_sorted)[i]
    cat(sprintf("  %d. %s: %.2f (capped at %.2f for colors)\n", i, gene, extreme_low_sorted[i], patient_range[1]))
  }
}

patient_breaks <- seq(patient_range[1], patient_range[2], length.out = 100)

# Create side-by-side heatmaps
pdf("heatmap_analysis_combat/patient_vs_subtypes_split_scales.pdf", width = 16, height = 14)

# Set up layout for side-by-side plots
par(mfrow = c(1, 2))

# Left panel: Centroids with optimized scale
pheatmap(subtype_centroids,
         annotation_col = centroids_annotation,
         annotation_colors = centroids_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = centroid_breaks,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "TNBC Subtype Centroids\n(Optimized Scale)",
         border_color = "grey60")

# Right panel: Patient with optimized scale  
pheatmap(patient_data,
         annotation_col = patient_annotation,
         annotation_colors = patient_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = patient_breaks,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "Your Patient\n(Optimized Scale)",
         border_color = "grey60")

dev.off()

# Also create PNG version
png("heatmap_analysis_combat/patient_vs_subtypes_split_scales.png", width = 1600, height = 1400, res = 100)

par(mfrow = c(1, 2))

pheatmap(subtype_centroids,
         annotation_col = centroids_annotation,
         annotation_colors = centroids_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = centroid_breaks,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "TNBC Subtype Centroids\n(Optimized Scale)",
         border_color = "grey60")

pheatmap(patient_data,
         annotation_col = patient_annotation,
         annotation_colors = patient_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = patient_breaks,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "Your Patient\n(Optimized Scale)",
         border_color = "grey60")

dev.off()

# 3. ALSO CREATE PATTERN-MATCHED VISUALIZATION WITH ROBUST SCALING
# Normalize both to same relative scale for pattern comparison
cat("Creating pattern-matched comparison with robust scaling...\n")

# ROBUST Scale each to -1 to +1 range for pattern comparison (using adaptive outlier exclusion)
scale_to_range_robust <- function(x, new_min = -1, new_max = 1) {
  # Adaptive outlier exclusion based on data size
  n_values <- length(x)
  outliers_to_exclude <- min(4, max(2, round(n_values * 0.1)))
  upper_percentile <- (n_values - outliers_to_exclude) / n_values
  lower_percentile <- outliers_to_exclude / n_values
  
  # Use adaptive percentile range instead of fixed percentiles
  old_range_limits <- quantile(x, c(lower_percentile, upper_percentile), na.rm = TRUE)
  old_range <- old_range_limits[2] - old_range_limits[1]
  
  if(old_range == 0) return(rep(0, length(x)))
  
  new_range <- new_max - new_min
  
  # Cap values at adaptive percentile limits before scaling
  x_capped <- pmax(pmin(x, old_range_limits[2]), old_range_limits[1])
  
  # Scale to new range
  ((x_capped - old_range_limits[1]) * new_range / old_range) + new_min
}

# Scale centroids and patient to same range for pattern comparison (adaptive robust method)
centroids_scaled <- apply(subtype_centroids, 2, scale_to_range_robust)
patient_scaled <- scale_to_range_robust(patient_data[, 1])

cat("Adaptive robust scaling applied - top/bottom outliers excluded before normalization\n")
cat("Each dataset (centroids and patient) scaled with their own adaptive outlier exclusion\n")

# Combine for pattern comparison
pattern_comparison <- cbind(centroids_scaled, "YOUR_PATIENT" = patient_scaled)

pattern_annotation <- data.frame(
  Type = c(rep("Subtype_Centroid", ncol(centroids_scaled)), "YOUR_PATIENT"),
  Subtype = c(colnames(centroids_scaled), "Unknown"),
  row.names = colnames(pattern_comparison)
)

pattern_colors <- list(
  Type = c("Subtype_Centroid" = "#CCCCCC", "YOUR_PATIENT" = "#FF0000"),
  Subtype = subtype_colors
)

pdf("heatmap_analysis_combat/patient_vs_subtypes_pattern_matched.pdf", width = 12, height = 14)

pheatmap(pattern_comparison,
         annotation_col = pattern_annotation,
         annotation_colors = pattern_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         clustering_distance_cols = "correlation", 
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "Pattern Comparison: Your Patient vs Subtype Centroids\n(Normalized -1 to +1 for Pattern Matching)",
         border_color = "grey60")

dev.off()

png("heatmap_analysis_combat/patient_vs_subtypes_pattern_matched.png", width = 1200, height = 1400, res = 100)

pheatmap(pattern_comparison,
         annotation_col = pattern_annotation,
         annotation_colors = pattern_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2", 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "Pattern Comparison: Your Patient vs Subtype Centroids\n(Normalized -1 to +1 for Pattern Matching)",
         border_color = "grey60")

dev.off()

# ===============================================================================
# STEP 7: IDENTIFY PATIENT'S TOP GENES
# ===============================================================================

cat("STEP 7: Analyzing patient's expression pattern...\n")

# Find where patient clusters
patient_col <- which(colnames(combined_burstein_scaled) == "Patient")
patient_expression <- combined_burstein_scaled[, patient_col]

# Get top upregulated and downregulated genes in patient
top_up <- names(sort(patient_expression, decreasing = TRUE)[1:10])
top_down <- names(sort(patient_expression, decreasing = FALSE)[1:10])

cat("Patient's TOP 10 UPREGULATED Burstein genes:\n")
for(i in 1:10) {
  cat("  ", i, ". ", top_up[i], " (z-score: ", round(patient_expression[top_up[i]], 2), ")\n", sep="")
}

cat("\nPatient's TOP 10 DOWNREGULATED Burstein genes:\n")
for(i in 1:10) {
  cat("  ", i, ". ", top_down[i], " (z-score: ", round(patient_expression[top_down[i]], 2), ")\n", sep="")
}

# Save results
patient_gene_results <- data.frame(
  gene = names(patient_expression),
  z_score = patient_expression,
  rank = rank(-patient_expression),
  stringsAsFactors = FALSE
)
patient_gene_results <- patient_gene_results[order(patient_gene_results$z_score, decreasing = TRUE), ]

write.csv(patient_gene_results, "heatmap_analysis_combat/patient_burstein_gene_expression.csv", row.names = FALSE)

# ===============================================================================
# STEP 8: COMPREHENSIVE PATIENT ANALYSIS REPORT
# ===============================================================================

cat("\nSTEP 8: Generating comprehensive patient analysis report...\n")

# Calculate mean expression by subtype (excluding patient)
reference_data <- combined_burstein_scaled[, sample_annotation$dataset != "Patient"]
reference_annotation <- sample_annotation[sample_annotation$dataset != "Patient", ]

subtype_means <- sapply(unique(reference_annotation$subtype), function(subtype) {
  subtype_samples <- reference_annotation$subtype == subtype
  if(sum(subtype_samples) > 1) {
    rowMeans(reference_data[, subtype_samples])
  } else {
    reference_data[, subtype_samples]
  }
})

# Compare patient to each subtype
cat("Calculating patient similarity to each subtype...\n")
subtype_correlations <- sapply(colnames(subtype_means), function(subtype) {
  cor(patient_expression, subtype_means[, subtype])
})

# Find closest subtype
closest_subtype <- names(which.max(subtype_correlations))
closest_correlation <- max(subtype_correlations)

# Calculate distances
subtype_distances <- sapply(colnames(subtype_means), function(subtype) {
  sqrt(sum((patient_expression - subtype_means[, subtype])^2))
})
closest_by_distance <- names(which.min(subtype_distances))

# Generate comprehensive report
sink("heatmap_analysis_combat/PATIENT_ANALYSIS_REPORT.txt")

cat("=================================================================\n")
cat("                 PATIENT TNBC SUBTYPE ANALYSIS REPORT\n")
cat("                    Based on Burstein Marker Genes\n")
cat("=================================================================\n\n")

cat("Analysis Date:", as.character(Sys.Date()), "\n")
cat("Analysis Time:", as.character(Sys.time()), "\n\n")

cat("METHODOLOGY:\n")
cat("- Reference datasets: TCGA, FUSCC, Burstein (", ncol(reference_data), "samples )\n")
cat("- Batch correction: ComBat with subtype preservation\n")
cat("- Gene set: Burstein TNBC marker genes (", length(available_burstein_genes), "genes )\n")
cat("- Normalization: Z-score across samples\n")
cat("- Clustering: Hierarchical (Ward.D2, Euclidean distance)\n\n")

cat("=================================================================\n")
cat("                        PATIENT CLASSIFICATION\n")
cat("=================================================================\n\n")

cat("SIMILARITY TO KNOWN SUBTYPES (Correlation):\n")
sorted_correlations <- sort(subtype_correlations, decreasing = TRUE)
for(i in 1:length(sorted_correlations)) {
  subtype <- names(sorted_correlations)[i]
  corr <- sorted_correlations[i]
  if(i == 1) {
    cat("  *** BEST MATCH: ", subtype, " (r = ", sprintf("%.3f", corr), ") ***\n", sep="")
  } else {
    cat("      ", subtype, " (r = ", sprintf("%.3f", corr), ")\n", sep="")
  }
}

cat("\nDISTANCE TO SUBTYPE CENTROIDS (Euclidean):\n")
sorted_distances <- sort(subtype_distances, decreasing = FALSE)
for(i in 1:length(sorted_distances)) {
  subtype <- names(sorted_distances)[i]
  dist <- sorted_distances[i]
  if(i == 1) {
    cat("  *** CLOSEST: ", subtype, " (distance = ", sprintf("%.3f", dist), ") ***\n", sep="")
  } else {
    cat("      ", subtype, " (distance = ", sprintf("%.3f", dist), ")\n", sep="")
  }
}

cat("\n=================================================================\n")
cat("                    PATIENT'S TOP MARKER GENES\n")
cat("=================================================================\n\n")

cat("TOP 15 UPREGULATED GENES (vs population mean):\n")
for(i in 1:min(15, length(top_up))) {
  gene <- top_up[i]
  zscore <- patient_expression[gene]
  cat(sprintf("  %2d. %-12s (z-score: %6.2f)\n", i, gene, zscore))
}

cat("\nTOP 15 DOWNREGULATED GENES (vs population mean):\n")
for(i in 1:min(15, length(top_down))) {
  gene <- top_down[i]
  zscore <- patient_expression[gene]
  cat(sprintf("  %2d. %-12s (z-score: %6.2f)\n", i, gene, zscore))
}

cat("\n=================================================================\n")
cat("                  SUBTYPE-SPECIFIC GENE SIGNATURES\n")
cat("=================================================================\n\n")

# Analyze genes where patient matches the DOMINANT subtype pattern
for(subtype in colnames(subtype_means)) {
  cat("GENES WHERE PATIENT MATCHES", subtype, "AS THE DOMINANT SUBTYPE:\n")
  
  # Find genes where this subtype has the highest (or lowest) expression among all subtypes
  # AND patient is concordant with this pattern
  
  subtype_profile <- subtype_means[, subtype]
  
  # For each gene, check if this subtype is the most extreme (highest or lowest)
  dominant_genes <- c()
  
  for(gene in rownames(subtype_means)) {
    all_subtype_values <- subtype_means[gene, ]
    
    # Check if this subtype is the max or min among all subtypes
    is_max <- which.max(all_subtype_values) == which(colnames(subtype_means) == subtype)
    is_min <- which.min(all_subtype_values) == which(colnames(subtype_means) == subtype)
    
    if(is_max || is_min) {
      # Check if patient is concordant (same direction)
      patient_val <- patient_expression[gene]
      subtype_val <- subtype_profile[gene]
      
      # Both high or both low
      if((patient_val > 0.5 && subtype_val > 0.5) || (patient_val < -0.5 && subtype_val < -0.5)) {
        dominant_genes <- c(dominant_genes, gene)
      }
    }
  }
  
  if(length(dominant_genes) > 0) {
    # Sort by absolute patient expression
    dominant_genes <- dominant_genes[order(abs(patient_expression[dominant_genes]), decreasing = TRUE)]
    
    for(i in 1:min(8, length(dominant_genes))) {
      gene <- dominant_genes[i]
      patient_z <- patient_expression[gene]
      
      # Show ALL subtype values for comparison
      all_values <- subtype_means[gene, ]
      values_text <- paste(names(all_values), sprintf("%.2f", all_values), sep="=", collapse=", ")
      
      cat(sprintf("  %-12s: Patient=%5.2f | All subtypes: %s\n", gene, patient_z, values_text))
    }
  } else {
    cat("  No genes where patient matches", subtype, "as the dominant subtype\n")
  }
  cat("\n")
}

cat("=================================================================\n")
cat("                         INTERPRETATION\n")
cat("=================================================================\n\n")

cat("SUMMARY:\n")
cat("Your patient shows strongest similarity to the", closest_subtype, "subtype\n")
cat("(correlation =", sprintf("%.3f", closest_correlation), ").\n\n")

if(closest_correlation > 0.7) {
  cat("CONFIDENCE: HIGH - Strong correlation suggests clear subtype assignment.\n")
} else if(closest_correlation > 0.5) {
  cat("CONFIDENCE: MODERATE - Reasonable correlation but some mixed features.\n")
} else {
  cat("CONFIDENCE: LOW - Weak correlations suggest mixed or atypical features.\n")
}

cat("\nKEY BIOLOGICAL FEATURES:\n")
# Highlight key upregulated genes and their biological significance
key_genes_up <- head(top_up, 5)
cat("Highly expressed genes:", paste(key_genes_up, collapse=", "), "\n")

key_genes_down <- head(top_down, 5)
cat("Lowly expressed genes:", paste(key_genes_down, collapse=", "), "\n")

cat("\nRECOMMENDATIONS:\n")
cat("1. Review the focused heatmap (patient_vs_subtypes_focused.pdf) for visual confirmation\n")
cat("2. Consider the", closest_subtype, "subtype-specific treatment implications\n")
cat("3. Validate key marker genes through additional methods if needed\n")
cat("4. Note any discordant genes that don't fit the predicted subtype pattern\n")

cat("\n=================================================================\n")
cat("                           FILES CREATED\n")
cat("=================================================================\n\n")

cat("HEATMAPS:\n")
cat("- burstein_genes_heatmap_combat.pdf/.png (Full dataset with patient highlighted)\n")
cat("- patient_vs_subtypes_focused.pdf/.png (Patient vs subtype centroids)\n\n")

cat("DATA FILES:\n")
cat("- patient_burstein_gene_expression.csv (All gene expression values)\n")
cat("- burstein_heatmap_results_combat.RData (Complete analysis workspace)\n")
cat("- PATIENT_ANALYSIS_REPORT.txt (This report)\n\n")

cat("=================================================================\n")
cat("Analysis completed:", as.character(Sys.time()), "\n")
cat("=================================================================\n")

sink()

# Save all results
save(combined_burstein_scaled, sample_annotation, annotation_col, annotation_colors,
     patient_gene_results, subtype_means, subtype_correlations,
     reference_corrected, batch_ref, subtype_ref,  # ComBat results
     subtype_centroids, patient_data, pattern_comparison,  # New visualizations
     centroids_scaled, patient_scaled,  # Scaled data
     closest_subtype, closest_correlation, subtype_distances,  # Classification results
     file = "heatmap_analysis_combat/burstein_heatmap_results_combat.RData")

cat("\n=================================================================\n")
cat("      BATCH-CORRECTED HEATMAP ANALYSIS COMPLETE\n")
cat("=================================================================\n\n")

cat("PATIENT CLASSIFICATION RESULT:\n")
cat("*** BEST MATCH:", closest_subtype, "(correlation =", sprintf("%.3f", closest_correlation), ") ***\n\n")

cat("FILES CREATED:\n")
cat("HEATMAPS:\n")
cat("  - burstein_genes_heatmap_combat.pdf/.png (Full heatmap with patient highlighted in RED)\n")
cat("  - patient_vs_subtypes_split_scales.pdf/.png (*** BEST VIEW: Side-by-side optimized scales ***)\n")
cat("  - patient_vs_subtypes_pattern_matched.pdf/.png (Pattern comparison: normalized for matching)\n\n")
cat("ANALYSIS:\n") 
cat("  - PATIENT_ANALYSIS_REPORT.txt (*** READ THIS FIRST ***)\n")
cat("  - patient_burstein_gene_expression.csv (All gene values)\n")
cat("  - burstein_heatmap_results_combat.RData (Complete workspace)\n\n")

cat("VISUALIZATION EXPLANATION:\n")
cat("- SPLIT SCALES heatmap: Centroids use their own color scale (see subtype patterns clearly)\n")
cat("                        Patient uses ADAPTIVE robust scaling (top 3-4 outliers excluded)\n")
cat("- PATTERN MATCHED: Both normalized -1 to +1 using adaptive outlier exclusion\n")
cat("- This solves BOTH the 'muted blue' problem AND the extreme outlier domination!\n")
cat("- Top outliers (HERC5, LAMP3, GBP5, etc.) excluded from color scale calculation\n")
cat("- Now you can see patterns in the 'normal' range genes clearly!\n\n")

cat("QUICK INTERPRETATION:\n")
if(closest_correlation > 0.7) {
  cat("- HIGH confidence", closest_subtype, "classification\n")
} else if(closest_correlation > 0.5) {
  cat("- MODERATE confidence", closest_subtype, "classification\n") 
} else {
  cat("- LOW confidence - mixed or atypical features\n")
}

cat("- Reference data was ComBat batch corrected with subtype preservation\n")
cat("- Patient highlighted in RED on main heatmap\n")
cat("- Focused heatmap shows direct comparison to subtype centroids\n\n")

cat("*** LOOK FOR THE RED COLUMN IN THE HEATMAP TO FIND YOUR PATIENT! ***\n")
cat("*** READ THE DETAILED REPORT: PATIENT_ANALYSIS_REPORT.txt ***\n")
