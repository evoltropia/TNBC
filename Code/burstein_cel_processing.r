# ===============================================================================
# PROCESSING BURSTEIN TNBC CEL FILES
# ===============================================================================

# ===============================================================================
# PROCESSING BURSTEIN TNBC CEL FILES
# ===============================================================================

# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("affy", "GEOquery", "hgu133plus2.db", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("affy", "GEOquery", "hgu133plus2.db")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load required libraries
library(affy)          # For reading CEL files and RMA normalization
library(GEOquery)      # For downloading sample metadata
library(hgu133plus2.db) # Annotation for Affymetrix U133 Plus 2.0 (likely platform)
library(dplyr)

# ===============================================================================
# STEP 1: READ ALL CEL FILES
# ===============================================================================

# Set the directory containing your CEL files
cel_dir <- "GSE76124_RAW"  # Your downloaded folder

# Look for uncompressed CEL files only (skip automatic decompression)
cel_files <- list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE)

# Check if we found any CEL files
if (length(cel_files) == 0) {
  cat("No uncompressed .CEL files found.\n")
  cat("Please manually extract all .CEL.gz files in the", cel_dir, "folder.\n")
  cat("Then re-run this script.\n")
  stop("No CEL files ready for processing.")
}
cat("Found", length(cel_files), "CEL files\n")

# Read all CEL files into an AffyBatch object
cat("Reading CEL files... (this may take several minutes)\n")
raw_data <- ReadAffy(filenames = cel_files)

# Check the data
cat("Array dimensions:", dim(raw_data), "\n")
cat("Sample names (first 5):", sampleNames(raw_data)[1:5], "\n")

# ===============================================================================
# STEP 2: RMA NORMALIZATION
# ===============================================================================

cat("Performing RMA normalization... (this may take 10-15 minutes)\n")
# RMA performs background correction, quantile normalization, and summarization
normalized_data <- rma(raw_data)

# Extract the expression matrix
expr_matrix <- exprs(normalized_data)
cat("Normalized expression matrix dimensions:", dim(expr_matrix), "\n")
cat("Expression values are log2-transformed\n")

# ===============================================================================
# STEP 3: GET SAMPLE METADATA FROM GEO
# ===============================================================================

cat("Downloading sample metadata from GEO...\n")

# The GSE series for the Burstein study
gse_id <- "GSE76124"  # Burstein TNBC dataset

# Download the series matrix
gse_data <- getGEO(gse_id, GSEMatrix = TRUE)
if (length(gse_data) > 1) {
  # If multiple platforms, select the appropriate one
  gse_data <- gse_data[[1]]  # Usually the first one
} else {
  gse_data <- gse_data[[1]]
}

# Extract sample metadata
sample_metadata <- pData(gse_data)

# Clean up sample names to match CEL file names
# Extract GSM IDs from CEL filenames
cel_gsm_ids <- gsub(".*/(GSM\\d+)_.*\\.CEL$", "\\1", cel_files)
sample_metadata$gsm_id <- rownames(sample_metadata)

# Match metadata to your samples
sample_metadata_matched <- sample_metadata[sample_metadata$gsm_id %in% cel_gsm_ids, ]

cat("Matched", nrow(sample_metadata_matched), "samples with metadata\n")

# ===============================================================================
# STEP 4: GENE ANNOTATION
# ===============================================================================

cat("Adding gene annotations...\n")

# Get probe IDs
probe_ids <- rownames(expr_matrix)

# Map probes to gene symbols
gene_symbols <- select(hgu133plus2.db, 
                      keys = probe_ids, 
                      columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                      keytype = "PROBEID")

# Remove probes with no gene symbol
gene_symbols_clean <- gene_symbols[!is.na(gene_symbols$SYMBOL), ]

# Filter expression matrix to only annotated probes
expr_matrix_annotated <- expr_matrix[gene_symbols_clean$PROBEID, ]

# ===============================================================================
# STEP 5: HANDLE MULTIPLE PROBES PER GENE
# ===============================================================================

cat("Collapsing multiple probes per gene...\n")

# Create a data frame for easier manipulation
expr_df <- data.frame(
  PROBEID = gene_symbols_clean$PROBEID,
  SYMBOL = gene_symbols_clean$SYMBOL,
  expr_matrix_annotated,
  stringsAsFactors = FALSE
)

# For genes with multiple probes, keep the one with highest mean expression
expr_collapsed <- expr_df %>%
  group_by(SYMBOL) %>%
  # Calculate mean expression across samples for each probe
  mutate(mean_expr = rowMeans(select(., starts_with("GSM")), na.rm = TRUE)) %>%
  # Keep probe with highest mean expression
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-mean_expr, -PROBEID)

# Convert back to matrix with gene symbols as rownames
final_expr_matrix <- as.matrix(expr_collapsed[, -1])
rownames(final_expr_matrix) <- expr_collapsed$SYMBOL

cat("Final expression matrix:", dim(final_expr_matrix), "\n")

# ===============================================================================
# STEP 6: EXTRACT BURSTEIN 80-GENE SIGNATURE
# ===============================================================================

# Load your Burstein signature genes (from your existing file)
burstein_genes <- read.csv("burstein_signatures.csv", stringsAsFactors = FALSE)
signature_genes <- unique(burstein_genes$Gene)

# Extract expression for signature genes
signature_expr <- final_expr_matrix[rownames(final_expr_matrix) %in% signature_genes, ]

cat("Extracted", nrow(signature_expr), "of", length(signature_genes), "signature genes\n")

# ===============================================================================
# STEP 7: CREATE FINAL DATASET
# ===============================================================================

# Combine everything into a list for easy access
burstein_dataset <- list(
  expression_matrix = final_expr_matrix,
  signature_expression = signature_expr,
  sample_metadata = sample_metadata_matched,
  gene_annotation = gene_symbols_clean
)

# Save the dataset
save(burstein_dataset, file = "burstein_processed_data.RData")
cat("Dataset saved as 'burstein_processed_data.RData'\n")

# ===============================================================================
# STEP 8: QUICK QUALITY CHECK
# ===============================================================================

cat("\n=== QUALITY CHECK ===\n")

# Check if we have the expected subtypes
if ("tnbc.subtype" %in% colnames(sample_metadata_matched)) {
  subtype_counts <- table(sample_metadata_matched$tnbc.subtype)
  cat("Subtype distribution:\n")
  print(subtype_counts)
}

# Check expression levels
cat("Expression statistics (log2 scale):\n")
cat("Min:", min(final_expr_matrix), "\n")
cat("Max:", max(final_expr_matrix), "\n")
cat("Median:", median(final_expr_matrix), "\n")

# Simple correlation check between samples
sample_correlations <- cor(final_expr_matrix[, 1:min(10, ncol(final_expr_matrix))])
cat("Sample correlation range:", range(sample_correlations[lower.tri(sample_correlations)]), "\n")

cat("\n=== PROCESSING COMPLETE ===\n")
cat("You can now load the data with: load('burstein_processed_data.RData')\n")
cat("Access components with:\n")
cat("  - burstein_dataset$expression_matrix (all genes)\n")
cat("  - burstein_dataset$signature_expression (80-gene signature)\n")
cat("  - burstein_dataset$sample_metadata (clinical data)\n")
