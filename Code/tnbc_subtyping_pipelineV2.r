# ===============================================================================
# TNBC Subtyping Pipeline: From BAM File to Molecular Subtype Classification
# ===============================================================================
# 
# This script takes a single RNA-seq BAM file from Tempus xR and performs
# TNBC (Triple Negative Breast Cancer) molecular subtyping using established
# Lehmann signatures.
#
# Input: BAM file from RNA-seq
# Output: TNBC subtype classification and associated analyses
#
# 
# ===============================================================================

# ===============================================================================
# SECTION 1: PACKAGE INSTALLATION AND LOADING
# ===============================================================================

# Check if required packages are installed, install if needed
required_packages <- c(
  "Rsubread",      # For reading BAM files and counting features
  "GenomicFeatures", # For handling genomic annotations
  "DESeq2",        # For normalization and data handling
  "edgeR",         # For TPM calculation
  "TNBCtype",      # For TNBC subtyping (if available)
  "biomaRt",       # For gene annotation
  "ggplot2",       # For plotting
  "pheatmap",      # For heatmaps
  "corrplot",      # For correlation plots
  "dplyr",         # For data manipulation
  "readr",         # For reading files
  "tibble"         # For data structures
)

# Install Bioconductor packages if not present
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("Rsubread", "GenomicFeatures", "DESeq2", "edgeR", "biomaRt")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load all required libraries
library(Rsubread)
library(GenomicFeatures)
library(DESeq2)
library(edgeR)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(dplyr)
library(readr)
library(tibble)

# ===============================================================================
# SECTION 2: FILE PATHS AND PARAMETERS SETUP
# ===============================================================================

# MODIFY THESE PATHS TO MATCH YOUR FILES
bam_file <- "jennybams_RNAseq_TL-25-6Z9K9WW8PF_T_sorted.bam"          # Your BAM file path
gtf_file <- "gencode.v19.annotation.gtf_withproteinids"           # GTF annotation file
sample_name <- "JennyTNBC"                # Name for your sample
output_dir <- "TNBC_Analysis_Results"          # Output directory

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set parameters for analysis
min_count_threshold <- 10    # Minimum read count to consider gene as expressed
tpm_threshold <- 0.01          # TPM threshold for expressed genes

cat("=== TNBC Subtyping Pipeline Started ===\n")
cat("Sample:", sample_name, "\n")
cat("BAM file:", bam_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# ===============================================================================
# SECTION 3: GENE EXPRESSION QUANTIFICATION FROM BAM FILE
# ===============================================================================

cat("STEP 1: Quantifying gene expression from BAM file...\n")

# Count reads mapping to genes using featureCounts
# This is the core step that converts aligned reads to gene-level counts
feature_counts <- featureCounts(
  files = bam_file,                    # Input BAM file
  annot.ext = gtf_file,               # GTF annotation file
  isGTFAnnotationFile = TRUE,         # Specify that annotation is GTF format
  GTF.featureType = "exon",           # Count reads mapping to exons
  GTF.attrType = "gene_id",           # Use gene_id attribute for grouping
  useMetaFeatures = TRUE,             # Sum counts across all exons of a gene
  allowMultiOverlap = FALSE,          # Don't count reads mapping to multiple features
  minOverlap = 1,                     # Minimum overlap between read and feature
  countMultiMappingReads = FALSE,     # Don't count multi-mapping reads
  strandSpecific = 0,                 # Assume unstranded library (0=unstranded, 1=stranded, 2=reverse)
  isPairedEnd = FALSE,                 # Assume paired-end sequencing
  nthreads = 4                        # Use 4 threads for faster processing
)

# Extract the count matrix and gene information
raw_counts <- feature_counts$counts
gene_info <- feature_counts$annotation

# Check the results
cat("Number of genes quantified:", nrow(raw_counts), "\n")
cat("Total mapped reads:", sum(raw_counts), "\n")
cat("Genes with >", min_count_threshold, "counts:", sum(raw_counts > min_count_threshold), "\n\n")

# ===============================================================================
# SECTION 4: GENE ANNOTATION AND CLEANUP
# ===============================================================================

cat("STEP 2: Cleaning gene names and adding annotations...\n")

# Clean up gene IDs (remove version numbers if present)
# Ensembl IDs often have version numbers (e.g., ENSG00000000003.14)
clean_gene_ids <- sub("\\..*", "", rownames(raw_counts))
rownames(raw_counts) <- clean_gene_ids

# Load gene symbol mapping from local file (instead of biomaRt)
cat("Loading gene mapping from local file...\n")
gene_mapping <- read.delim("ensembl_gene_mapping_grch37.txt", stringsAsFactors = FALSE)

# Ensure column names match what the rest of the script expects
names(gene_mapping) <- c("ensembl_gene_id", "external_gene_name", "gene_biotype")

# Filter mapping to only genes in our data
gene_mapping <- gene_mapping[gene_mapping$ensembl_gene_id %in% clean_gene_ids, ]

cat("Successfully loaded", nrow(gene_mapping), "gene mappings\n")

# Merge count data with gene symbols
gene_info_extended <- merge(
  data.frame(ensembl_gene_id = clean_gene_ids, stringsAsFactors = FALSE),
  gene_mapping,
  by = "ensembl_gene_id",
  all.x = TRUE
)

# Filter for protein-coding genes only (removes pseudogenes, etc.)
protein_coding_genes <- gene_info_extended$gene_biotype == "protein_coding" & 
                       !is.na(gene_info_extended$gene_biotype)

raw_counts_pc <- raw_counts[protein_coding_genes, , drop = FALSE]
gene_info_pc <- gene_info_extended[protein_coding_genes, ]

cat("Protein-coding genes retained:", nrow(raw_counts_pc), "\n\n")

# ===============================================================================
# SECTION 5: NORMALIZATION AND TPM CALCULATION
# ===============================================================================

cat("STEP 3: Calculating normalized expression values (TPM)...\n")

# Calculate gene lengths from the GTF annotation
# Gene length is needed for TPM calculation
gene_lengths <- gene_info$Length[protein_coding_genes]
names(gene_lengths) <- gene_info_pc$ensembl_gene_id

# Calculate TPM (Transcripts Per Million)
# TPM accounts for both library size and gene length
# Formula: TPM = (raw_count / gene_length) / sum(raw_count / gene_length) * 1e6

# Calculate FPKM instead of TPM
total_mapped_reads <- sum(raw_counts_pc[, 1])
fpkm_values <- (raw_counts_pc[, 1] * 1e9) / (gene_lengths * total_mapped_reads)

# Create a data frame with gene information and expression
expression_data <- data.frame(
  ensembl_gene_id = gene_info_pc$ensembl_gene_id,
  gene_symbol = gene_info_pc$external_gene_name,
  raw_counts = raw_counts_pc[, 1],
  fpkm = fpkm_values,
  stringsAsFactors = FALSE
)

# Remove genes with missing symbols and low expression
expression_data <- expression_data[
  !is.na(expression_data$gene_symbol) & 
  expression_data$gene_symbol != "" & 
  expression_data$fpkm >= tpm_threshold,
]

# Handle duplicate gene symbols by keeping the highest expressed
expression_data <- expression_data %>%
  group_by(gene_symbol) %>%
  slice_max(tpm, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Final gene count for analysis:", nrow(expression_data), "\n")
cat("Median TPM:", median(expression_data$tpm), "\n")
cat("Mean TPM:", mean(expression_data$tpm), "\n\n")

# ===============================================================================
# SECTION 6: QUALITY CONTROL CHECKS
# ===============================================================================

cat("STEP 4: Performing quality control checks...\n")

# Check expression of key breast cancer markers
# These should be low/absent in TNBC samples
key_markers <- c("ESR1", "PGR", "ERBB2")  # ER, PR, HER2
marker_expression <- expression_data[expression_data$gene_symbol %in% key_markers, ]

cat("Key marker expression levels:\n")
if (nrow(marker_expression) > 0) {
  for (i in 1:nrow(marker_expression)) {
    cat(sprintf("  %s: %.2f TPM\n", 
                marker_expression$gene_symbol[i], 
                marker_expression$tpm[i]))
  }
} else {
  cat("  Warning: Key markers not found in data\n")
}

# Check housekeeping gene expression (should be moderate to high)
housekeeping <- c("ACTB", "GAPDH", "B2M", "HPRT1")
hk_expression <- expression_data[expression_data$gene_symbol %in% housekeeping, ]

cat("\nHousekeeping gene expression:\n")
if (nrow(hk_expression) > 0) {
  for (i in 1:nrow(hk_expression)) {
    cat(sprintf("  %s: %.2f TPM\n", 
                hk_expression$gene_symbol[i], 
                hk_expression$tpm[i]))
  }
} else {
  cat("  Warning: Housekeeping genes not found\n")
}

# Create expression vector for downstream analysis
tpm_vector <- setNames(expression_data$tpm, expression_data$gene_symbol)

cat("\nQuality control completed.\n\n")

# ===============================================================================
# SECTION 7: TNBC SUBTYPE GENE SIGNATURES
# ===============================================================================

cat("STEP 5: Defining TNBC subtype gene signatures...\n")

# Lehmann et al. TNBC subtype signatures
# These are the established gene signatures for TNBC molecular subtypes
tnbc_signatures <- list(
  
  # BL1 (Basal-like 1): DNA damage response and cell cycle
  BL1 = c("CCNE1", "CDC25A", "CDC6", "CHEK1", "RRM2", "E2F1", "CCNE2", 
          "PCNA", "CDC25B", "CDK1", "TOP2A", "CCNB1", "BUB1", "PLK1",
          "AURKA", "AURKB", "CENPE", "CENPF"),
  
  # BL2 (Basal-like 2): Growth factor signaling
  BL2 = c("EGFR", "MET", "NRAS", "KRAS", "IGF1R", "EPHA2", "TNF", "IL6",
          "VEGFA", "PDGFRA", "FGFR1", "FGFR2", "ERBB3", "PIK3CA"),
  
  # M (Mesenchymal): EMT and cell motility
  M = c("VIM", "SNAI1", "TWIST1", "ZEB1", "CDH2", "FN1", "FOXC1", "SPARC",
        "COL1A1", "COL3A1", "COL5A1", "ACTA2", "TAGLN", "MMP2", "MMP9"),
  
  # MSL (Mesenchymal Stem-like): Stem cell characteristics
  MSL = c("ALDH1A1", "CD44", "SOX2", "NANOG", "POU5F1", "KLF4", "MYC",
          "NOTCH1", "NOTCH3", "WNT3A", "WNT5A", "DVL1", "CTNNB1"),
  
  # IM (Immunomodulatory): Immune signaling
  IM = c("STAT1", "IRF1", "TAP1", "CXCL9", "CCL5", "IDO1", "PDCD1LG2",
         "CD274", "IFNG", "TNF", "IL2RG", "CXCR3", "CCR5", "CD8A", "CD4"),
  
  # LAR (Luminal Androgen Receptor): Androgen receptor pathway
  LAR = c("AR", "KRT18", "XBP1", "GATA3", "FOXA1", "TFF3", "AGR2", "SPDEF",
          "SCUBE2", "ANKRD30A", "ERBB4", "TFF1", "PRLR", "TFAP2C")
)

# Count how many signature genes are present in our data
for (subtype in names(tnbc_signatures)) {
  genes_present <- sum(tnbc_signatures[[subtype]] %in% names(tpm_vector))
  total_genes <- length(tnbc_signatures[[subtype]])
  cat(sprintf("%s signature: %d/%d genes present (%.1f%%)\n", 
              subtype, genes_present, total_genes, 
              (genes_present/total_genes)*100))
}

cat("\n")

# ===============================================================================
# SECTION 8: SIGNATURE SCORE CALCULATION
# ===============================================================================

cat("STEP 6: Calculating signature scores for TNBC subtypes...\n")

# Calculate signature scores using mean expression of available genes
signature_scores <- sapply(tnbc_signatures, function(signature_genes) {
  # Find genes in signature that are present in our data
  available_genes <- signature_genes[signature_genes %in% names(tpm_vector)]
  
  if (length(available_genes) == 0) {
    return(NA)  # No genes available for this signature
  }
  
  # Calculate mean log2(TPM+1) for available genes
  mean(log2(tpm_vector[available_genes] + 1), na.rm = TRUE)
})

# Remove any signatures with missing scores
signature_scores <- signature_scores[!is.na(signature_scores)]

# Normalize scores to z-scores (subtract mean, divide by SD)
signature_scores_norm <- scale(signature_scores)[, 1]

cat("Signature scores (normalized):\n")
for (i in 1:length(signature_scores_norm)) {
  cat(sprintf("  %s: %.3f\n", names(signature_scores_norm)[i], signature_scores_norm[i]))
}

# Determine predicted subtype (highest scoring signature)
predicted_subtype <- names(signature_scores_norm)[which.max(signature_scores_norm)]
max_score <- max(signature_scores_norm)

cat(sprintf("\nPredicted TNBC Subtype: %s (score: %.3f)\n", predicted_subtype, max_score))

# ===============================================================================
# SECTION 9: CONFIDENCE ASSESSMENT
# ===============================================================================

cat("\nSTEP 7: Assessing prediction confidence...\n")

# Calculate confidence metrics
scores_sorted <- sort(signature_scores_norm, decreasing = TRUE)
top_score <- scores_sorted[1]
second_score <- scores_sorted[2]
score_difference <- top_score - second_score

# Assess confidence level
if (score_difference > 1.0) {
  confidence <- "High"
} else if (score_difference > 0.5) {
  confidence <- "Medium"
} else {
  confidence <- "Low"
}

cat(sprintf("Confidence level: %s\n", confidence))
cat(sprintf("Score difference between top subtypes: %.3f\n", score_difference))
cat(sprintf("Top 2 subtypes: %s (%.3f), %s (%.3f)\n", 
            names(scores_sorted)[1], scores_sorted[1],
            names(scores_sorted)[2], scores_sorted[2]))

# ===============================================================================
# SECTION 10: VISUALIZATION AND OUTPUT
# ===============================================================================

cat("\nSTEP 8: Creating visualizations and saving results...\n")

# 1. Signature score barplot
pdf(file.path(output_dir, "TNBC_signature_scores.pdf"), width = 10, height = 6)
barplot(signature_scores_norm, 
        main = paste("TNBC Subtype Signature Scores -", sample_name),
        ylab = "Normalized Signature Score",
        xlab = "TNBC Subtype",
        col = rainbow(length(signature_scores_norm)),
        las = 2)
abline(h = 0, lty = 2)
text(x = which.max(signature_scores_norm), 
     y = max(signature_scores_norm) + 0.1, 
     labels = "PREDICTED", 
     col = "red", font = 2)
dev.off()

# 2. Expression heatmap of signature genes
signature_genes_all <- unique(unlist(tnbc_signatures))
available_sig_genes <- signature_genes_all[signature_genes_all %in% names(tpm_vector)]

if (length(available_sig_genes) > 10) {  # Only create heatmap if enough genes
  
  # Create expression matrix for heatmap
  heatmap_data <- log2(tpm_vector[available_sig_genes] + 1)
  heatmap_matrix <- matrix(heatmap_data, ncol = 1)
  rownames(heatmap_matrix) <- available_sig_genes
  colnames(heatmap_matrix) <- sample_name
  
  # Create annotation for genes by subtype
  gene_annotations <- data.frame(
    Subtype = character(length(available_sig_genes)),
    stringsAsFactors = FALSE
  )
  
  for (subtype in names(tnbc_signatures)) {
    gene_annotations$Subtype[available_sig_genes %in% tnbc_signatures[[subtype]]] <- subtype
  }
  
  rownames(gene_annotations) <- available_sig_genes
  
  # Create heatmap
  pdf(file.path(output_dir, "TNBC_signature_genes_heatmap.pdf"), width = 8, height = 12)
  pheatmap(heatmap_matrix,
           main = paste("TNBC Signature Gene Expression -", sample_name),
           show_rownames = TRUE,
           show_colnames = TRUE,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_row = gene_annotations,
           color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
}

# ===============================================================================
# SECTION 11: DETAILED RESULTS SUMMARY
# ===============================================================================

# Create comprehensive results summary
results_summary <- list(
  sample_name = sample_name,
  analysis_date = Sys.Date(),
  
  # Expression summary
  total_genes_quantified = nrow(raw_counts),
  protein_coding_genes = nrow(raw_counts_pc),
  expressed_genes = nrow(expression_data),
  median_tpm = median(expression_data$tpm),
  
  # Key marker expression
  key_markers = marker_expression,
  housekeeping_genes = hk_expression,
  
  # Subtype prediction
  predicted_subtype = predicted_subtype,
  confidence_level = confidence,
  signature_scores = signature_scores_norm,
  score_difference = score_difference,
  
  # Signature gene coverage
  signature_coverage = sapply(tnbc_signatures, function(sig) {
    sum(sig %in% names(tpm_vector)) / length(sig)
  })
)

# Save results summary
saveRDS(results_summary, file.path(output_dir, "TNBC_analysis_results.rds"))

# Create human-readable summary report
sink(file.path(output_dir, "TNBC_analysis_report.txt"))

cat("=================================================================\n")
cat("TNBC MOLECULAR SUBTYPING ANALYSIS REPORT\n")
cat("=================================================================\n\n")

cat("Sample Information:\n")
cat("  Sample Name:", sample_name, "\n")
cat("  Analysis Date:", as.character(Sys.Date()), "\n")
cat("  BAM File:", bam_file, "\n\n")

cat("Expression Data Summary:\n")
cat("  Total genes quantified:", nrow(raw_counts), "\n")
cat("  Protein-coding genes:", nrow(raw_counts_pc), "\n")
cat("  Genes passing filters:", nrow(expression_data), "\n")
cat("  Median TPM:", round(median(expression_data$tpm), 2), "\n\n")

cat("Key Marker Expression (should be low in TNBC):\n")
if (nrow(marker_expression) > 0) {
  for (i in 1:nrow(marker_expression)) {
    cat(sprintf("  %s: %.2f TPM\n", 
                marker_expression$gene_symbol[i], 
                marker_expression$tpm[i]))
  }
} else {
  cat("  No key markers detected\n")
}

cat("\nTNBC SUBTYPE PREDICTION:\n")
cat("  Predicted Subtype:", predicted_subtype, "\n")
cat("  Confidence Level:", confidence, "\n")
cat("  Top Score:", round(max(signature_scores_norm), 3), "\n")
cat("  Score Difference:", round(score_difference, 3), "\n\n")

cat("All Signature Scores:\n")
for (i in 1:length(signature_scores_norm)) {
  cat(sprintf("  %s: %.3f\n", names(signature_scores_norm)[i], signature_scores_norm[i]))
}

cat("\nSignature Gene Coverage:\n")
for (subtype in names(tnbc_signatures)) {
  coverage <- sum(tnbc_signatures[[subtype]] %in% names(tpm_vector)) / length(tnbc_signatures[[subtype]])
  cat(sprintf("  %s: %.1f%% (%d/%d genes)\n", 
              subtype, coverage*100,
              sum(tnbc_signatures[[subtype]] %in% names(tpm_vector)),
              length(tnbc_signatures[[subtype]])))
}

cat("\n=================================================================\n")
cat("INTERPRETATION NOTES:\n")
cat("=================================================================\n\n")

if (predicted_subtype == "BL1") {
  cat("BL1 (Basal-like 1) subtype is characterized by:\n")
  cat("- High DNA damage response signaling\n")
  cat("- Cell cycle dysregulation\n")
  cat("- May respond well to DNA-damaging agents\n")
} else if (predicted_subtype == "BL2") {
  cat("BL2 (Basal-like 2) subtype is characterized by:\n")
  cat("- Growth factor receptor signaling\n")
  cat("- May respond to targeted therapies (EGFR, etc.)\n")
} else if (predicted_subtype == "M") {
  cat("M (Mesenchymal) subtype is characterized by:\n")
  cat("- Epithelial-mesenchymal transition\n")
  cat("- Cell motility and invasion\n")
  cat("- May be more aggressive\n")
} else if (predicted_subtype == "MSL") {
  cat("MSL (Mesenchymal Stem-like) subtype is characterized by:\n")
  cat("- Stem cell-like properties\n")
  cat("- May be resistant to conventional therapy\n")
} else if (predicted_subtype == "IM") {
  cat("IM (Immunomodulatory) subtype is characterized by:\n")
  cat("- High immune cell infiltration\n")
  cat("- May respond well to immunotherapy\n")
} else if (predicted_subtype == "LAR") {
  cat("LAR (Luminal Androgen Receptor) subtype is characterized by:\n")
  cat("- Androgen receptor pathway activation\n")
  cat("- May respond to anti-androgen therapy\n")
}

cat("\nConfidence Interpretation:\n")
if (confidence == "High") {
  cat("- High confidence: Clear subtype assignment\n")
} else if (confidence == "Medium") {
  cat("- Medium confidence: Likely subtype but consider alternatives\n")
} else {
  cat("- Low confidence: Uncertain classification, may be mixed subtype\n")
}

sink()

# Save expression data
write.csv(expression_data, file.path(output_dir, "gene_expression_data.csv"), row.names = FALSE)

cat("Analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
cat("Files created:\n")
cat("  - TNBC_analysis_results.rds (R data object)\n")
cat("  - TNBC_analysis_report.txt (human-readable report)\n")
cat("  - TNBC_signature_scores.pdf (barplot visualization)\n")
cat("  - gene_expression_data.csv (expression data)\n")
if (length(available_sig_genes) > 10) {
  cat("  - TNBC_signature_genes_heatmap.pdf (heatmap visualization)\n")
}

cat("\n=== PIPELINE COMPLETED SUCCESSFULLY ===\n")

# ===============================================================================
# END OF SCRIPT
# ===============================================================================
