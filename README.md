# Blood-Based Machine Learning Biomarkers for Parkinson's Disease

A leakage-free reanalysis of blood transcriptomic biomarkers for Parkinson's disease (PD), using nested cross-validation, WGCNA network analysis, and cross-tissue comparison.

## Overview

This repository contains the complete analysis code for the manuscript:

> **"Machine Learning-Derived Blood Transcriptomic Biomarkers for Parkinson's Disease Show Limited Cross-Cohort Generalizability"**

The study evaluates whether blood-based ML biomarkers for PD generalize across independent cohorts and tissues, using a leakage-free nested pipeline that prevents the information leakage common in transcriptomic biomarker studies.

## Datasets

| Dataset | Tissue | Role | Platform | Samples |
|---------|--------|------|----------|---------|
| [GSE99039](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99039) | Peripheral blood | Discovery | GPL570 | 438 (205 PD, 233 control) |
| [GSE6613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6613) | Peripheral blood | Independent validation | GPL96 | 72 (50 PD, 22 control) |
| [GSE7621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7621) | Substantia nigra | Cross-tissue comparison | GPL570 | 25 (16 PD, 9 control) |
| [GSE54536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54536) | Peripheral blood | Exploratory only | GPL10558 | 10 (5 PD, 5 control) |

## Key Results

| Analysis | AUC | 95% CI |
|----------|-----|--------|
| Repeated 5-fold CV (GSE99039, discovery) | 0.682 ± 0.049 | — |
| Independent validation (GSE6613, blood) | 0.618 | 0.470 – 0.767 |
| Cross-tissue comparison (GSE7621, brain) | 0.549 | 0.299 – 0.798 |
| WGCNA + LASSO validation | 0.543 | 0.394 – 0.692 |

**Cross-tissue direction concordance:** 51% (chance = 50%)

**Conclusion:** Neither single-gene ML nor WGCNA-based feature selection yields a blood biomarker panel for PD that generalizes to independent cohorts or tissues. Top enriched pathways are immune-related (NK cells, leukocyte activation), not neuron-specific.

## Repository Structure
