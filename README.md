# Planning CT Radiomics Model for Grade ≥ 2 Leukopenia After Pelvic Radiotherapy Without Concurrent Chemotherapy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![R 4.0+](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

> **Public release of analysis code for the radiomics prediction model described in our paper submitted to *Radiotherapy and Oncology* (IF ≈ 6).**
> This repository contains **de-identified summary data and analysis scripts only**. No individual patient data, raw imaging, or protected health information (PHI) are included.

---

## 📋 Table of Contents

1. [Study Summary](#-study-summary)
2. [Repository Structure](#-repository-structure)
3. [Reproducing the Analysis](#-reproducing-the-analysis)
4. [File Map](#-file-map)
5. [De-identification Policy](#-de-identification-policy)
6. [Citing this Work](#-citing-this-work)
7. [License](#-license)
8. [Contact](#-contact)

---

## 🔬 Study Summary

| Item | Value |
|---|---|
| **Population** | 60 patients with pelvic cancer (cervical, endometrial, rectal) treated with pelvic radiotherapy |
| **Key exclusion** | Patients receiving concurrent chemotherapy or with prior radiotherapy / baseline bone-marrow disease |
| **Prescription** | Uniform 45 Gy / 25 fractions, 5 weeks |
| **Imaging** | Planning CT (Siemens SOMATOM Definition AS, 120 kVp, 3 mm slice) |
| **Feature extraction** | AccuContour (Manteia) radiomics module, 1,409 features, bin width 25 HU |
| **Modeling** | LASSO logistic (λ.1se, 10-fold CV) → unpenalized IRLS refit |
| **Validation** | Full-pipeline bootstrap optimism correction (B = 500, including LASSO re-screening) |
| **Final model** | 2 features: MeshVolume (shape) + DependenceVariance (GLCM) |
| **Apparent AUC** | 0.793 (95% CI 0.69–0.90) |
| **Optimism-corrected AUC** | 0.695 |
| **Calibration** | Hosmer–Lemeshow P = 0.535 |

The full paper, when published, will appear in *Radiotherapy and Oncology*.

---

## 📁 Repository Structure

```
pelvic-radiomics-grade2-leukopenia/
├── README.md                 ← this file
├── LICENSE                   ← MIT License
├── CITATION.cff              ← academic citation format
├── .gitignore                ← Python + Node + R + PHI defense
├── requirements.txt          ← Python dependencies
├── requirements_R.txt        ← R dependencies
├── GITHUB_SETUP_GUIDE.md     ← (Chinese) how to register a GitHub account
├── PUSH_GUIDE.md             ← (Chinese) how to push this folder to GitHub
│
├── scripts/
│   ├── python/               ← analysis scripts in Python
│   │   ├── make_figures.py
│   │   ├── nomogram.py
│   │   ├── make_master_table_DEID.py      ← DE-IDENTIFIED
│   │   ├── make_master_table_v2_DEID.py   ← DE-IDENTIFIED
│   │   ├── make_master_table_v3_DEID.py   ← DE-IDENTIFIED
│   │   ├── compare_dose.py
│   │   └── recalc_table2.py
│   ├── R/
│   │   └── radiomics_bm_pipeline_validation.R
│   └── js/                   ← Node.js cross-validation kernels
│       ├── validation_pipeline.js          ← main LASSO + IRLS + bootstrap kernel
│       ├── roc_curves.js
│       ├── model_comparison.js
│       ├── model_comparison_v2_DEID.js     ← DE-IDENTIFIED
│       ├── nri_idi.js
│       └── nri_idi_36win.js
│
├── data/                     ← summary statistics only (no PHI)
│   ├── README.md
│   ├── nomogram_params.json
│   ├── model_comparison_results.json
│   ├── model_comparison_v2_results.json
│   ├── nri_idi_result.json
│   ├── nri_idi_36win.json
│   ├── nri_idi_48all.json
│   ├── roc_curves.json
│   ├── ambiguous_resolved.json
│   ├── dosi_group_comparison.csv
│   ├── dosi_group_comparison_v2.csv
│   └── 剂量学组间比较.xlsx
│
├── results/
│   ├── v2_final/
│   │   ├── results.json
│   │   ├── full_pipeline_validation_summary.csv
│   │   └── feature_selection_frequency.csv
│   ├── figures/              ← 8 paper figures (SVG, vector)
│   │   ├── fig1_dca.svg
│   │   ├── fig2_feature_freq.svg
│   │   ├── fig3_baseline_boxplot.svg
│   │   ├── fig4_dose_boxplot.svg
│   │   ├── fig5_auc_bar.svg
│   │   ├── fig6_roc.svg
│   │   ├── fig7_calibration.svg
│   │   └── fig8_nomogram.svg
│   └── tables/
│       └── (Table 1–4 in supplementary CSV form — see v2_final/ for summary stats)
│
└── paper/                    ← (empty at this time)
    └── README.md             ← RoO policy: full text added after acceptance
```

---

## 🔁 Reproducing the Analysis

The original analysis was performed in two parallel environments to cross-validate:

### Python + Node.js path (recommended for non-R users)

```bash
# 1. Install Python dependencies
pip install -r requirements.txt          # only openpyxl

# 2. Verify Node.js is installed (built-in libraries used; no npm packages)
node --version                           # should be v18+

# 3. Run the LASSO + IRLS + bootstrap kernel
node scripts/js/validation_pipeline.js

# 4. Run the ROC comparison
node scripts/js/roc_curves.js

# 5. Run the model comparison
node scripts/js/model_comparison.js
node scripts/js/model_comparison_v2_DEID.js

# 6. Run the NRI/IDI bootstrap
node scripts/js/nri_idi.js
node scripts/js/nri_idi_36win.js

# 7. Generate the figures
python scripts/python/make_figures.py

# 8. Build the nomogram
python scripts/python/nomogram.py
```

### R path (canonical statistical environment)

```r
# Install dependencies
install.packages(c("readxl", "glmnet", "pROC", "rms"))

# Run the canonical R pipeline
source("scripts/R/radiomics_bm_pipeline_validation.R")
```

The R script and the Node.js kernels produce **numerically identical** results on the same input data, providing a strong cross-validation.

---

## 📂 File Map

### Core analysis kernels

| File | What it does | Lines |
|---|---|---|
| `scripts/js/validation_pipeline.js` | LASSO + IRLS + full-pipeline bootstrap (B=500) | 509 |
| `scripts/R/radiomics_bm_pipeline_validation.R` | R canonical version of the same | 170 |
| `scripts/python/nomogram.py` | Harrell nomogram construction | 188 |
| `scripts/python/make_figures.py` | 8 paper figures (SVG) | 385 |
| `scripts/js/roc_curves.js` | 4-model ROC coordinates | 175 |
| `scripts/js/model_comparison.js` | 4-model AUC comparison | 166 |
| `scripts/js/nri_idi.js` | NRI/IDI bootstrap (48-case main) | 221 |
| `scripts/js/nri_idi_36win.js` | NRI/IDI bootstrap (36-case strict window) | 221 |

### Data inputs (summary statistics only)

| File | Contents |
|---|---|
| `data/nomogram_params.json` | Nomogram coefficients (β0, βDV, βMV) and z-score parameters |
| `data/model_comparison_results.json` | 4-model AUC and DeLong tests |
| `data/nri_idi_result.json` | 48-case NRI/IDI bootstrap results |
| `data/nri_idi_36win.json` | 36-case strict-window NRI/IDI bootstrap results |
| `data/roc_curves.json` | 4-model ROC coordinates and calibration data |
| `data/ambiguous_resolved.json` | Feature disambiguation map (Contrast→GLCM etc.) |
| `data/dosi_group_comparison.csv` | Dose-volume metrics by event status |

### Results

| File | Contents |
|---|---|
| `results/v2_final/results.json` | Full bootstrap validation summary |
| `results/v2_final/feature_selection_frequency.csv` | Per-feature selection frequency over B=500 resamples |
| `results/figures/fig*.svg` | 8 paper-quality figures (vector SVG) |

---

## 🔒 De-identification Policy

This repository follows a strict **de-identification standard**:

1. **No individual patient data.** All 60 patients are represented only by summary statistics (model coefficients, AUC, p-values) and aggregate outputs.
2. **Patient names → P001–P060 codes.** Where any clinical mapping is needed (e.g., in the master-table generator scripts), patients are referenced as `P001` … `P060`. The mapping from real names to P-codes is held by the corresponding author and is **not** included in this repository.
3. **2 excluded patients** (originally in the 62-patient candidate set but dropped due to dose heterogeneity) are labeled `P061` and `P062` in the de-identified scripts.
4. **No raw imaging data.** Planning CT, dose distributions, and structure sets are not included.
5. **No dates of treatment.** Time references are aggregated to "days before radiotherapy" or to calendar quarters, never to specific calendar dates.
6. **No geographic identifiers.** The hospital is identified only as "[Single-center, China]."

The corresponding author maintains the master key linking P-codes to original records, in compliance with the original IRB approval.

---

## 📑 Citing this Work

A `CITATION.cff` file in this repository provides machine-readable citation metadata. When the paper is published, the canonical citation will be:

> [Author et al.]. A Planning CT Radiomics Model for Predicting Grade ≥ 2 Leukopenia After Pelvic Radiotherapy Without Concurrent Chemotherapy. *Radiotherapy and Oncology*. 2026. (In press)

Until publication, please cite the preprint / GitHub repository as:

> [Author]. (2026). *pelvic-radiomics-grade2-leukopenia* (Version 1.0.0) [Computer software]. https://github.com/[username]/pelvic-radiomics-grade2-leukopenia

---

## 📜 License

This repository is released under the [MIT License](LICENSE). You are free to use, modify, and redistribute the code with attribution. The methodology and clinical findings described in the associated paper remain the intellectual property of the authors.

```
MIT License

Copyright (c) 2026 [Author]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND...
```

---

## ✉️ Contact

For questions about the analysis, the de-identification policy, or to request access to the master de-identification key for research collaboration, please contact the corresponding author via the email address in the published paper.

---

*Last updated: 2026-10-13 · Repository version 1.0.0 · No PHI included.*
