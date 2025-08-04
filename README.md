# TNBC Single-Patient Molecular Subtyping Analysis

Analysis pipeline for molecular subtyping of triple-negative breast cancer using patient RNA-seq data integrated with reference cohorts.

## Overview

This repository contains the analysis for a single TNBC patient case where molecular subtyping was performed using RNA-seq data originally obtained from Tempus for clinical testing. The raw sequencing data, which was only used to detect fusions in the original clinical report, was requested by the patient under HIPAA rights and subjected to comprehensive molecular subtyping analysis.

## Why TNBC Subtyping Matters

Triple-negative breast cancer represents 15-20% of all breast cancers but accounts for disproportionate mortality. Recent research has identified distinct molecular subtypes (BLIS, IM, LAR, MES) with significantly different prognoses, treatment responses, and clinical trial eligibility. Immunotherapy efficacy varies 3-fold between subtypes, and many new clinical trials now stratify enrollment by molecular subtype. However, most patients never receive this analysis because it's not standard of care and existing genomic data from clinical tests often goes unused for secondary analyses.

## Project Context

This analysis demonstrates the principle of patients owning their raw and comprehensive molecular data for ongoing secondary analyses as methods and knowledge improve. Rather than letting valuable RNA-seq data sit unused after the initial clinical report, this patient requested their complete dataset and collaborated on this analysis. The patient is a co-author on the resulting manuscript, showcasing patient-partnered research principles.

## Repository Structure

```
TNBC_single_patient_subtyping/
├── Code/
│   ├── patient_bam_processing.R         # Processes patient RNA-seq BAM file from Tempus
│   ├── burstein_cel_processing.R        # Processes Burstein reference dataset CEL files  
│   ├── complete_lda_analysis.R         # LDA classification with ComBat batch correction
│   └── burstein_heatmap_combat4.R               # Generates visualizations and result summaries
├── Manuscript/
│   └── [manuscript draft]               # First draft of publication
└── Results/
    ├── [folder from complete_lda_analysis]     # LDA classification results
    └── [folder from heatmap_analysis_combat]           # Visualization outputs
```

## Files

### Data Processing Scripts (Code/)
- `patient_bam_processing.R` - Processes patient RNA-seq BAM file from Tempus to extract gene expression
- `burstein_cel_processing.R` - Processes Burstein reference dataset CEL files and performs normalization

### Analysis Scripts (Code/)
- `lda_combat_analysis_v4.R` - Performs LDA classification with ComBat batch correction using both all genes and Burstein marker genes on Burstein HUSCC and TCGA Datasets
- `burstein_heatmap.R` - Generates visualizations and result summaries

## Results

The analysis scripts generate output folders in Results/:
- Results from `lda_combat_analysis_v4.R` - Classification outcomes, distance metrics, LDA plots
- Results from `burstein_heatmap.R` - Heatmaps, correlation analyses, expression profiles

## Methods Summary

1. **Data Integration**: Patient RNA-seq integrated with 716 reference TNBC samples from TCGA, FUSCC, and Burstein cohorts
2. **Batch Correction**: ComBat normalization applied to harmonize data across platforms
3. **Classification**: Linear discriminant analysis using both genome-wide (4,887 genes) and targeted (26 Burstein markers) approaches
4. **Validation**: Correlation analysis and hierarchical clustering for classification confidence assessment

## Key Results

- Genome-wide LDA classified patient as IM subtype
- Targeted LDA classified patient as BLIS subtype  
- Correlation analysis showed strongest similarity to IM centroid (r=0.392)
- Results indicate mixed molecular features with classification uncertainty

## Requirements

```
R packages:
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
- sva (ComBat batch correction)
- pheatmap (visualization)
- RColorBrewer
- dplyr
- MASS (LDA)
- ggplot2
- RColorBrewer
- dplyr
```

## Data Sources

- Patient RNA-seq: Tempus xR clinical sequencing
- Reference cohorts: TCGA (n=158), FUSCC (n=360), Burstein (n=198)
- Gene signatures: 78 Burstein TNBC subtype marker genes
  Burstein et al. 2015:


"Comprehensive Genomic Analysis Identifies Novel Subtypes and Targets of Triple-Negative Breast Cancer"

Published in: Clinical Cancer Research
Link: https://aacrjournals.org/clincancerres/article/21/7/1688/248479/Comprehensive-Genomic-Analysis-Identifies-Novel

Jiang et al. 2019:
"Genomic and Transcriptomic Landscape of Triple-Negative Breast Cancers: Subtypes and Treatment Strategies"

Published in: Cancer Cell
Link: https://www.cell.com/cancer-cell/fulltext/S1535-6108(19)30096-0

## Publication

Manuscript in preparation with patient as co-author.

## Contact

CEtHI (Community Empowerment through Health Information) https://cethi.org - Providing infrastructure, culture and path towards patient-owned healthcare.
