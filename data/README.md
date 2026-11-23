# Input data

The patient-level data are not distributed with this repository. They are
retained by the authors under the institutional review board approval that
covers the study (KYLW2026-13) and can be made available to a named researcher
on reasonable request, subject to a data-use agreement.

The pipeline reads two files from this directory, and optionally three more for
the inter-observer step. Templates with the exact header layout are in
`schema/`.

```
data/
├── cohort.csv            one row per patient - outcome, covariates, dosimetry
├── radiomics.csv         one row per patient, one column per radiomic feature
└── icc/                  only needed for steps 06 and 07
    ├── observer1.csv
    ├── observer2.csv
    └── subset.csv
```

## De-identification

Nothing in `cohort.csv`, `radiomics.csv` or the files in `icc/` may carry a
name, a hospital number, an accession number, a study date or any DICOM header
field. `patient_id` is the only link between the files and should be a locally
assigned pseudonym (the shipped templates use `SAMPLE001`, `SAMPLE002`, ...).

`check_phi.py` in the repository root scans a directory for the residues that
usually survive a manual export - Chinese personal names, hospital-number
patterns, dates and DICOM keywords - and should be run on `data/` before the
directory is copied anywhere.

## Feature column names

Radiomic features must be named `<image type>:<feature class>:<feature name>`,
for example `original:shape:MeshVolume` or `wavelet-HHH:glrrm:RunVariance`.
The three-part name is required because 1307 of the 1408 feature names repeat
across image types; a bare feature name is not a unique key.
