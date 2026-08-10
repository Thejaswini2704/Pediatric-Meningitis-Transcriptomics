# Pediatric Meningitis Transcriptomics

## Overview

This project investigates transcript-level molecular signatures associated with pediatric meningitis, with the primary objective of identifying transcript/isoform candidates that can distinguish **bacterial meningitis from viral meningitis**.

The analysis integrates RNA-seq transcript quantification, transcript-level differential expression, and repeated stability selection using machine learning to prioritize robust transcript biomarkers.

## Study Design

The dataset contains **15 samples** across three study groups:

- 5 Bacterial meningitis
- 5 Viral meningitis
- 5 Healthy controls

The primary biomarker analysis focuses on the **Bacterial vs Viral** comparison.

## Objectives

1. Perform transcriptome profiling of pediatric meningitis samples.
2. Identify transcript-level expression differences between bacterial and viral meningitis.
3. Investigate infection-associated isoform/transcript patterns.
4. Identify stable transcript candidates using repeated machine-learning-based stability selection.
5. Annotate prioritized transcripts with their corresponding gene information.
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
