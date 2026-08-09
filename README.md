# Pediatric Meningitis Transcriptomics


## Overview

This project aims to identify transcript-level molecular biomarkers that can distinguish bacterial meningitis, viral meningitis, and healthy controls using RNA-seq data.

## Study Groups

- 5 Viral meningitis samples
- 5 Bacterial meningitis samples
- 5 Healthy controls

## Planned Workflow
FASTQ
→ Quality Control
→ HISAT2 Alignment
→ StringTie Transcript Quantification
→ Differential Expression
→ Isoform Switching Analysis
→ Machine Learning
→ Biomarker Selection
→ Functional Annotation

## Tools

- FastQC
- MultiQC
- HISAT2
- StringTie
- R
- edgeR
- limma
- IsoformSwitchAnalyzeR
- glmnet
- caret
- pROC
- ggplot2
- pheatmap
- biomaRt

## Current Status

The project is being developed as an end-to-end RNA-seq analysis workflow, starting from publicly available GEO/SRA sequencing data.
