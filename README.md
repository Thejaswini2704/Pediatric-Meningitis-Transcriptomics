Pediatric Meningitis Transcriptomics

This project investigates transcript-level molecular signatures associated with pediatric meningitis, with the primary objective of identifying transcript/isoform candidates that can distinguish **bacterial meningitis from viral meningitis**.

The analysis integrates RNA-seq transcript quantification, transcript-level differential expression, and repeated machine-learning-based stability selection to prioritize robust transcript biomarker candidates.

## Study Design

The dataset contains **15 samples** across three study groups:

* 5 Bacterial meningitis
* 5 Viral meningitis
* 5 Healthy controls

The primary biomarker discovery analysis focuses on the **Bacterial vs Viral** comparison.

## Objectives

1. Perform transcriptome profiling of pediatric meningitis samples.
2. Identify transcript-level expression differences between bacterial and viral meningitis.
3. Investigate infection-associated transcript/isoform patterns.
4. Identify stable transcript candidates using repeated machine-learning-based stability selection.
5. Annotate prioritized transcripts with corresponding gene information.
6. Explore expression patterns and relationships among the prioritized transcript biomarkers.

## Analysis Workflow

```text
RNA-seq FASTQ
      ↓
Quality Control
      ↓
HISAT2 Alignment
      ↓
StringTie Transcript Quantification
      ↓
Transcript Count Matrix
      ↓
Transcript-level Differential Expression
      ↓
Bacterial vs Viral Candidate Selection
      ↓
Repeated Stability Selection
      ↓
Top 20 Stable Transcript Candidates
      ↓
Transcript Annotation
      ↓
Expression / Correlation Analysis
      ↓
Heatmap and PCA Visualization
```

## Key Result

Repeated stability selection identified **20 stable transcript candidates** from the Bacterial vs Viral comparison.

The final prioritized candidates include transcripts associated with genes such as:

* **ST3GAL3**
* **PKP2**
* **GOT1-DT**
* **MAP2K1**
* **HEXA**
* **ANKRD24**
* **RANGAP1**
* **LINC00299**
* **HOMER3**
* **CCRL2**
* **PRMT1**
* **OAT**
* **DAZAP1**

Some transcript identifiers are assembled StringTie transcripts (`MSTRG.*`) and therefore do not currently have a standard gene symbol annotation.

## Stability Selection

The candidate transcript set was evaluated using repeated machine-learning-based stability selection.

The analysis produced **53 successful model iterations**. The top 20 transcripts were prioritized according to their selection frequency across successful iterations.

The highest-ranked 15 transcripts were selected in **100% of successful iterations**, while the remaining candidates showed selection frequencies ranging from approximately **66% to 92%**.

## Repository Structure

```text
Pediatric-Meningitis-Transcriptomics/
│
├── metadata/
│   └── samples.csv
│
├── scripts/
│   ├── 02_transcript_biomarker_analysis.R
│   └── 03_ml_stability_selection.R
│
├── results/
│   └── biomarkers/
│       ├── Final_20_Bacterial_vs_Viral_Biomarkers.csv
│       ├── Bacterial_vs_Viral_Stable_Top20.csv
│       ├── Stable_Top20_Transcript_Annotation.csv
│       ├── Stable_Top20_Expression_Matrix.csv
│       ├── Stable_Top20_Correlation_Matrix.csv
│       └── README.md
│
└── figures/
    ├── Stable_Top20_Biomarker_Heatmap.png
    └── Stable_Top20_PCA.png
```

## Reproducibility

Large sequencing files, genome annotations, alignment files, intermediate StringTie outputs, and other computationally intensive intermediate results are intentionally excluded from the GitHub repository.

The repository instead contains the analysis scripts, sample metadata, key biomarker results, annotations, and final visualizations required to document the biomarker-discovery workflow.

## Tools

* FastQC
* MultiQC
* HISAT2
* StringTie
* R
* edgeR
* limma
* glmnet
* caret
* ggplot2
* pheatmap
* rtracklayer

## Current Status

The transcript-level biomarker discovery workflow has been completed for the current 15-sample dataset.

The current result is a **20-transcript candidate biomarker panel for distinguishing bacterial from viral meningitis**. These candidates should be considered research candidates requiring validation in an independent cohort before clinical interpretation or diagnostic use.
****
