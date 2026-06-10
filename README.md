# TAOLS Multicointegration ECS Estimation — Replication Package

This repository contains all data, code, and figures required to fully reproduce the empirical results presented in the accompanying manuscript. It is provided in direct response to reviewer requests for full reproducibility, including all preprocessing steps, model specifications, and sensitivity analyses.

---
## 1. Overview

This project implements the **TAOLS (Transformed and Augmented Ordinary Least Squares)** estimator within a multicointegration framework to estimate equilibrium climate sensitivity (ECS) from historical radiative forcing and surface temperature data.

The replication package allows full reconstruction of:

- Data preprocessing and alignment
- Unit root and stationarity tests (ADF, KPSS)
- Cointegration and multicointegration diagnostics
- TAOLS estimation across K ∈ [10, 150]
- Confidence interval construction across all specifications
- Sensitivity analysis across datasets and model pairings
This repository contains two subfolders: one for MATLAB and one for Python. Each folder contains its own README file with software-specific instructions for reproducing the analyses.---
## 2. Data Sources

All datasets are publicly available and included in processed form for reproducibility:

- Radiative forcing series (total and partial efficacy)
- Surface temperature series (Berkeley Earth, HadCRUT)
- Ocean heat content series (for multicointegration structure)

No proprietary or restricted data sources are used.

---
