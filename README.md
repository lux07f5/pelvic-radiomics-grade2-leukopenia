# Analysis code for *Pelvic bone volume predicts grade 2 or worse leukopenia after pelvic radiotherapy alone, but radiomic texture analysis does not*

This directory reproduces every number, table and figure in the manuscript from
two de-identified CSV files. It is the code that produced the submitted results;
the numbers it prints at each step are the numbers that appear in the paper.

The pipeline is deliberately boring. It uses the Python standard library only —
no numpy, no scipy, no scikit-learn, no R — so a reader who has Python can run it
without installing anything, and the statistical estimator is fully visible in
`lib/stats_lib.py` instead of being hidden inside a compiled dependency.

---

## Requirements

| item | version |
|---|---|
| Python | 3.8 or later (developed and tested on 3.12.9) |
| third-party packages | none |
| disk | about 20 MB for inputs and outputs |

Check the interpreter:

```bash
python3 --version
```

There is nothing to install. `pip install -r requirements.txt` is a no-op kept
for people who expect it to exist.

---

## Quick start

```bash
# 1. put the two de-identified CSV files in data/ (see "Input data" below)
ls data/cohort.csv data/radiomics.csv

# 2. see the steps and their order
python3 run_all.py --list

# 3. run everything
python3 run_all.py
```

Useful variants:

```bash
python3 run_all.py --only 03          # a single step
python3 run_all.py --from 04          # resume after a failure
python3 run_all.py --skip-icc         # when the second-observer tables are absent
python3 run_all.py --repeats 20       # shorter stability analysis (default 100)
python3 run_all.py --jobs 4           # stability step on 4 cores, about 30 min
```

Every step writes its console output to `outputs/logs/`. A step that fails stops
the run and prints the tail of its log, so the cause is visible immediately.

Inputs and outputs can live elsewhere:

```bash
RAD_DATA=/path/to/data RAD_OUT=/path/to/out python3 run_all.py
```

No absolute path from the original workstation remains anywhere in this
directory. All paths are derived from the location of `config.py`.

---

## Input data

The patient-level data are **not** included. They remain with the authors under
the institutional review board approval covering the study and can be shared
with a named researcher on reasonable request. `data/schema/` holds templates
with the exact expected headers; `data/README.md` states the rules.

### `data/cohort.csv` — 60 rows, one per patient

| column | type | meaning |
|---|---|---|
| `patient_id` | text | pseudonym, the only key linking to the feature table |
| `group` | 0 / 1 | 1 = grade ≥2 leukopenia within one month of radiotherapy |
| `tumor` | text | 直肠癌 / 宫颈癌 / 子宫内膜癌 (rectal, cervical, endometrial) |
| `age` | number | years |
| `sex` | text | 男 / 女 (male / female) |
| `stage` | 1–4 | FIGO or TNM stage, collapsed to I–II versus III–IV for Table 1 |
| `base_wbc` | number, may be empty | baseline white cell count, ×10⁹/L |
| `base_hb` | number, may be empty | baseline haemoglobin, g/L |
| `base_lym` | number, may be empty | baseline lymphocyte count, ×10⁹/L |
| `bmi` | number, may be empty | kg/m² |
| `bleed_rt_gap_days` | integer, may be empty | days from baseline blood draw to first fraction |
| `pelvis_volume` | number | contoured bone volume, cm³ |
| `min_dose`, `max_dose`, `mean_dose` | number | cGy |
| `V5`, `V10`, `V20`, `V35`, `V40` | number | percentage of the bone volume receiving at least that dose |

Baseline laboratory values are missing for 12 of the 60 patients. They are not
imputed: analyses that use them run on the 48 complete cases, which is why the
denominators in Table 1 differ between rows.

### `data/radiomics.csv` — 60 rows, one column per feature

```
patient_id,original:shape:MeshVolume,original:shape:SurfaceVolumeRatio,...,wavelet-HHH:glrlm:RunVariance
SAMPLE001,1095.42,1.241,...,0.842
```

* One row per patient, `patient_id` first, then **1408 feature columns**
  (1409 as exported, minus one constant feature dropped in step 01).
* Column names are `<image type>:<feature class>:<feature name>`. The three-part
  name is mandatory: 1307 of the feature names repeat across image types, so
  `Elongation` alone is ambiguous.
* Values are extracted by PyRadiomics through MANTEIA AccuContour 4.0, fixed bin
  width 25 HU, Laplacian of Gaussian σ = 1.0 mm, on the original image and 14
  derived images.
* No missing values. A patient present in `cohort.csv` but absent here is a
  hard error, since the models are fitted on the full 60.

### `data/icc/` — only needed for steps 06 and 07

| file | content |
|---|---|
| `observer1.csv` | `patient_id` + feature columns, first observer |
| `observer2.csv` | same layout, second observer |
| `subset.csv` | `patient_id,set` with `set` = `preset20` or `extra6` |

If these three files are absent, `run_all.py` skips steps 06 and 07 and says so.
All other steps run unaffected.

---

## Running order

The numbering is the running order; it is not cosmetic. Step 03 must follow 01,
step 04 must follow 03, and the figures must follow 08 and 09.

| # | script | consumes | produces | typical runtime |
|---|---|---|---|---|
| 01 | `01_build_dataset.py` | `cohort.csv`, `radiomics.csv` | `merged_dataset.csv`, `radiomics_matrix.json` | < 1 s |
| 02 | `02_table1_univariate.py` | `merged_dataset.csv` | `phase1_results.json` | 2 s |
| 03 | `03_radiomics_lasso.py` | `radiomics_matrix.json` | `phase2_results.json` | 3–5 min |
| 04 | `04_model_comparison.py` | 02, 03 | `phase3_results.json` | 1–2 min |
| 05 | `05_calibration_dca.py` | 02, 03 | `phase4_results.json` | < 1 min |
| 06 | `06_icc.py` | `data/icc/*` | `icc_results.json` | < 1 min |
| 07 | `07_icc_sensitivity.py` | 03, 06 | `phase5_all26.json` | 1–2 min |
| 08 | `08_mainline.py` | 03, `merged_dataset.csv` | `phase6_mainline.json` | 1–2 min |
| 09 | `09_stability.py` | 03 | `stability_results.json` | 1.5–3 h on one core |
| 10 | `10_figures_main.py` | `merged_dataset.csv`, 08, 09 | `figures/Figure_1..3*.svg` | < 1 s |
| 11 | `11_figures_supp.py` | 08 | `figures/Supplementary_Figure_S1,S2*.svg` | < 1 s |

Step 09 is the slow one: 100 repeats × 10-fold cross-validated LASSO, each fold
fitted by coordinate descent in pure Python, about 100 s per repeat. Run it with
`--jobs 4` and the same 100 repeats finish in roughly half an hour; the workers
take disjoint seed ranges, so the result is identical to the serial run. Each
repeat is appended to a JSONL file as it completes, so an interrupted step
resumes rather than starting over. `--repeats` is for smoke tests; the
manuscript uses 100.

---

## Outputs, and where they land in the paper

| output | manuscript item |
|---|---|
| `phase1_results.json` | Table 1 (baseline), Table 2 (dosimetric univariable) |
| `phase3_results.json` | discriminative performance of the three models, Figure 2A |
| `phase4_results.json` | calibration, Figure 2B cut-off, Supplementary Figure S2 (decision curves) |
| `icc_results.json` | Supplementary Table S2, inter-observer agreement |
| `phase5_all26.json` | ICC-filtered sensitivity analysis (Discussion) |
| `phase6_mainline.json` | Table 3, Figures 1 and 2 |
| `stability_results.json` | Figure 3, selection frequency across 100 repeats |
| `figures/*.svg` | final line art at 190 mm width, 7 pt type, ready for 1000 dpi export |

---

## Settings that determine the numbers

Repeating the analysis with these values reproduces the published results
exactly; changing any of them changes the results, and several are choices
rather than facts.

| setting | value | where |
|---|---|---|
| univariable screen | Mann-Whitney U, P < 0.05 | step 03 |
| redundancy pruning | \|Spearman ρ\| > 0.9, keep the smaller P | step 03 |
| candidate cap | 30 features (twice the 15 events) | step 03 |
| LASSO | 10-fold CV, 100 λ values, λ.1se, seed 42 | step 03 |
| AUC confidence interval | bootstrap, B = 2000 | steps 02–08 |
| optimism correction | bootstrap, B = 1000, seed 42 | step 03 |
| stability analysis | 100 repeats, seeds 1000…1099 | step 09 |
| ICC | ICC(A,1), two-way ANOVA, absolute agreement | step 06 |
| ICC threshold | 0.75 | step 07 |

Standardisation is fitted on the full sample before the LASSO, as in the
original analysis; it is not refitted inside each cross-validation fold. This
is a known source of mild optimism and is one of the reasons the Rad-score is
presented as exploratory rather than as a validated model.

---

## Verification

Before release the pipeline was rerun end to end on the study cohort, starting
from a de-identified export rather than from the original spreadsheets, and
compared field by field with the results behind the manuscript:

| comparison | outcome |
|---|---|
| 60 × 1408 feature matrix, rebuilt from CSV | identical to the last bit |
| `phase1`, `phase3`, `phase4`, `phase6` output | identical to 1e-9 |
| inter-observer ICC (26 patients) | median 0.947, 1226 of 1408 features ≥ 0.75, as reported |
| Supplementary Figures S1 and S2 | byte-identical |
| Figures 2 and 3 | byte-identical |
| Figure 1 | identical except the scatter jitter, now seeded on `patient_id` instead of the patient name |

The only intentional differences from the original working scripts are the
removal of absolute paths, the replacement of names by pseudonyms, the change of
two Chinese log labels to English, and one hard-coded feature index in step 08
that is now looked up by name.


---

## De-identification

Before `data/` is copied anywhere, run:

```bash
python3 check_phi.py data/
```

It fails on a direct-identifier column name (`name`, `姓名`, `住院号`, `mrn`, …),
on DICOM keywords and tag numbers, on a `.dcm` file, and on any calendar date in
the text. It warns on runs of six or more digits and on Chinese name-like values.
Add `--strict` to turn warnings into failures.

The screen is not a proof of anonymity. Read the headers once by eye and have
the data manager sign off; the screen only catches what is easy to forget.

---

## What this code does not do

* **No external validation.** Every number is apparent and internally
  cross-validated performance in a single-centre cohort of 60 patients with 15
  events. It is not expected performance elsewhere.
* **No marrow imaging.** The segmented structure is bone, not haematopoietic
  marrow; no MRI, PET or fat-fraction imaging was available.
* **The inter-observer coefficient is not a reproducibility study of
  contouring.** Both observers started from the same automatic segmentation, so
  the ICC describes the editing step only. Steps 06 and 07 are written to make
  that limitation explicit rather than to hide it.
* **No modelling of the treatment plan.** Dose was summarised by mean, maximum,
  minimum and V5–V40; nothing about beam arrangement or fractionation enters.

---

## License

MIT. See `LICENSE`. Use it, change it, publish with it; there is no warranty
that it will reproduce your own data.

## Contact

Zhang HL, Department of Radiation Oncology, Suzhou Ninth People's Hospital,
Suzhou, Jiangsu, China. zhl520@163.com

## Citation

If this code contributes to published work, please cite the article it
accompanies rather than the repository.

If the article is not yet published, cite this repository:

> Lu X, Liu J, Du M, Li Z, Hong Y, Zhang HL. Analysis code for: Pelvic bone
> volume predicts grade 2 or worse leukopenia after pelvic radiotherapy alone,
> but radiomic texture analysis does not (version 2.0.0) [Computer software].
> https://github.com/lux07f5/pelvic-radiomics-grade2-leukopenia

Machine-readable metadata is in `CITATION.cff`.

## Version history

| tag | date | contents |
|---|---|---|
| v1.0.0 | 2026-10-14 | First release, from an earlier stage of the study when it was framed as a prediction model with two retained features. Its scripts read intermediate files that were never published, and its results do not correspond to the manuscript. Superseded. |
| v2.0.0 | 2026-11-19 | The analysis behind the submitted manuscript: dataset assembly, feature selection, model comparison, calibration, decision curves, inter-observer agreement, selection stability and figures, plus the de-identification report. |
