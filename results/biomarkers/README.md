# Bacterial vs Viral Transcript Biomarkers

This directory contains transcript-level candidate biomarkers identified
from the pediatric meningitis transcriptomic dataset.

## Study groups

- Bacterial meningitis: 5 samples
- Viral meningitis: 5 samples
- Healthy controls: 5 samples

## Biomarker selection strategy

Transcript-level differential expression was performed using limma/voom.

For bacterial-versus-viral discrimination, the top 100 candidate transcripts
were subjected to repeated stratified subsampling followed by multinomial
elastic-net regression.

The analysis produced 53 successful models from 100 bootstrap iterations.

The 20 transcripts with the highest selection stability were retained as
candidate bacterial-versus-viral biomarkers.

## Final biomarker file

`Final_20_Bacterial_vs_Viral_Biomarkers.csv`

## Supporting files

- `Bacterial_vs_Viral_Stable_Top20.csv`
- `Stable_Top20_Transcript_Annotation.csv`
- `Stable_Top20_Expression_Matrix.csv`
- `Stable_Top20_Correlation_Matrix.csv`
- `Stable_Top20_Biomarker_Heatmap.png`
- `Stable_Top20_PCA.png`

## Important interpretation

These 20 transcripts are candidate biomarkers identified through stability
selection. They should not be described as clinically validated biomarkers.
Independent validation is required before clinical diagnostic claims.
