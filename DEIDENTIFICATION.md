# De-identification record

The repository is published; the patient data are not. This file records what
was checked before the code was released, file by file, so that the claim does
not rest on a memory of having looked.

## What was checked

1. **Direct identifiers.** Every file was searched for the 60 patient names of
   the study cohort (the full list held in the authors' local dataset), not just
   for a pattern. Match count: **0**.
2. **Identifier-shaped content.** `check_phi.py` was run over the whole
   directory: column names such as `name`, `姓名`, `住院号`, `mrn`, `accession`;
   DICOM keywords (`PatientName`, `PatientID`, `StudyInstanceUID`, …); DICOM tag
   numbers (`(0010,0010)`); calendar dates; runs of six or more digits; Chinese
   two- and three-character values that look like personal names. Result:
   **0 failures**, 1 warning, explained below.
3. **No DICOM objects.** The file-type screen reports any `.dcm`, `.ima` or
   `.dicom`, and none is present. No image header, no RTSTRUCT, no RTPLAN and no
   dose grid enters this repository.
4. **No hospital numbers.** The pipeline never reads a hospital number: the only
   key is `patient_id`, which the template defines as a locally assigned
   pseudonym (`SAMPLE001`, `SAMPLE002`, …). The scripts contain no lookup table,
   no cross-walk and no comment naming a patient.
5. **Figure code.** The scatter jitter in Figure 1 was seeded on the patient
   name in the original working script. It now seeds on `patient_id`, so the
   figure can be rebuilt without the name ever entering the process.

The single warning is `1000003` in `steps/10_figures_main.py`, a figure
coordinate, not an identifier.

Two files are exempt from the screen by name: `check_phi.py` and this file. Both
exist to write the identifiers down, so scanning them would only report their own
vocabulary. Everything else is scanned.

## File-by-file

Paths are relative to the repository root.

| file | content | identifiers |
|---|---|---|
| `README.md` | usage, input specification, running order | none |
| `LICENSE` | MIT text, author initials only | none |
| `requirements.txt` | comment stating that no package is needed | none |
| `.gitignore` | excludes `data/`, `outputs/`, `*.xlsx`, `*.dcm` | none |
| `config.py` | path resolution | none |
| `run_all.py` | step runner | none |
| `check_phi.py` | the screening tool itself | none (keyword list is inert text) |
| `lib/__init__.py` | empty package marker | none |
| `lib/stats_lib.py` | estimators: tests, logistic, LASSO, AUC, bootstrap | none |
| `steps/01_build_dataset.py` | builds the two analysis files from CSV | none |
| `steps/02_table1_univariate.py` | Table 1, dosimetric univariable analysis | none |
| `steps/03_radiomics_lasso.py` | filter, LASSO, Rad-score | none |
| `steps/04_model_comparison.py` | three-model comparison, ROC | none |
| `steps/05_calibration_dca.py` | calibration, decision curves, Firth refit | none |
| `steps/06_icc.py` | inter-observer ICC from two observer tables | none |
| `steps/07_icc_sensitivity.py` | ICC-filtered sensitivity analysis | none |
| `steps/08_mainline.py` | pelvis-volume-first analysis, DeLong, NRI/IDI | none |
| `steps/09_stability.py` | 100 repeats of cross-validated LASSO | none |
| `steps/09_stability_worker.py` | parallel worker, seed range as argument | none |
| `steps/09_stability_combine.py` | merges worker output | none |
| `steps/10_figures_main.py` | Figures 1–3 as SVG | none |
| `steps/11_figures_supp.py` | Supplementary Figures S1–S2 as SVG | none |
| `data/README.md` | what belongs in `data/` and what must not | none |
| `data/schema/cohort_template.csv` | header plus two invented rows | none (fictional `SAMPLE…`) |
| `data/schema/radiomics_template.csv` | header plus two invented rows | none (fictional `SAMPLE…`) |
| `data/schema/subset_template.csv` | header plus three invented rows | none (fictional `SAMPLE…`) |

## What is not in the repository

`data/cohort.csv`, `data/radiomics.csv` and `data/icc/*.csv` are absent by
design: they hold patient-level measurements and remain with the authors under
the ethics approval. `.gitignore` blocks them, so an accidental `git add .`
cannot commit them.

If the data are ever to be shared, de-identify first and then run:

```bash
python3 check_phi.py /path/to/shared/data --strict
```

A clean screen is a floor, not a ceiling. It cannot detect a re-identification
risk from a rare combination of tumour type, age and dosimetry in a small
cohort, which is the reason the data are shared only under a data-use agreement
rather than posted alongside the code.
