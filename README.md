# Placental metabolomics and lipidomics reproducibility code

This repository contains the reproducible preprocessing and batch-correction code supporting the manuscript *Chemometrically optimized sequential dual extraction for integrated metabolomics and lipidomics of human placenta by ultra-high performance liquid chromatography high-resolution mass spectrometry*.

The workflow was applied to placental samples from 371 participants in the ENVIRONAGE birth cohort. It supports the targeted metabolomics and lipidomics preprocessing steps described in the manuscript; it is not a complete end-to-end reproduction of the study.

## Repository contents

- `GapFilling.ipynb` - Python notebook that imputes missing targeted values feature-wise. Each imputed value is sampled uniformly between 50% and 100% of the observed minimum for that analyte.
- `ComBat_normalization.R` - R script for ComBat batch correction using collection year as the batch factor. It protects SGA status in the SGA analysis or low 1-minute APGAR score (<7) in the APGAR analysis.
