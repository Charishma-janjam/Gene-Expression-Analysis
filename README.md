**From Developmental Transition to Functional Divergence: Transcriptomic Profiling of Placenta and Bone Marrow Identifies Bone Marrow–Specific Links to Hematologic and Immune Disorders**

**Author**: Charishma Sruthi Janjam

## Overview
This repository presents a **differential gene expression (DGE) analysis comparing human bone marrow and placenta** using RNA-seq data. The primary goal of this project is to uncover how **gene dysregulation in these tissues may contribute to hematologic and immune disorders.**
Using pathway enrichment and STRING-based network analysis, **key bone marrow genes like MPO, BPI, and ELANE were linked to diseases such as leukemia, SCID, juvenile arthritis, thalassemia, lymphoma, multiple myeloma, Diamond-Blackfan anemia and SLE.** Placenta-upregulated genes were primarily involved in structural development, reflecting its limited hematopoietic role in late gestation.

## Repository Structure
├── gene_counts/               
├── report/                   
├── results/                   
│   ├── plots/                 
│   ├── gene_lists/            
│   └── dge_results_all_genes.csv  
├── scripts/                   
└── README.md                   

## Data & Methods
**Samples**: 3 placenta + 3 bone marrow replicates (sourced from Human Protein Atlas)

**Tools**: FastQC, Trim Galore, HISAT2, HTSeq, edgeR, DESeq2, pheatmap, Enrichr, DAVID, STRINGdb, Cytoscape

**Analysis pipeline**:
Quality control & trimming,
read mapping to GRCh38,
gene quantification,
differential expression using edgeR,
visualization (MDS, PCA, smear plot, BCV Plot, pheatmap, heatmap),
functional enrichment (DAVID, Enrichr),
PPI & disease network analysis,

## Key Results
**Total Differentially Expressed Genes (DEGs)**:

Upregulated in Placenta: 6,186

Upregulated in Bone Marrow: 4,496

Not Significant: 66,629

**Bone Marrow Upregulated Genes (e.g., MPO, AZU1, CTSG)**

→ Enriched for immune processes, T-cell signaling, myeloid differentiation

→ Linked to diseases: leukemia, SCID, juvenile arthritis, SLE

**Placenta Upregulated Genes (e.g., DLK1, EGFL6, COL4A1)**

→ Involved in extracellular matrix organization, development, angiogenesis

**Insights**

The placenta’s role in hematopoiesis appears limited to early gestation; late-gestation samples lacked immune/hematopoietic gene activity.
Bone marrow is a transcriptional hub for genes vital to lifelong immune function.
Disease linkage reveals that dysregulation in bone marrow genes can contribute to hematologic and autoimmune conditions.

## Visualizations

MDS_plot.png: Tissue-level clustering of samples

PCA_plot.png: Variance explained across samples

heatmap.png: Top 50 differentially expressed genes

smear_plot.png: Highlights DEGs (logFC vs. average expression)

bcv_plot.png: Gene expression variability across samples

pheatmap.png: Sample-to-sample correlation matrix

## How to Reproduce

Scripts are in the scripts/ folder. The analysis was conducted in R using the following packages:

library(edgeR)

library(DESeq2)

library(pheatmap)

library(RColorBrewer)

library(genefilter)

**Steps**:

Reading count files

Normalizing with calcNormFactors()

Running exactTest() and filtering with FDR

Visualizing with MDS, PCA, and heatmaps

## Report

The full PDF report is available in the /report folder: DGE Analysis of Bone Marrow and Placenta.pdf
