# -*- coding: utf-8 -*-
"""Step 01 - build the analysis dataset from two de-identified CSV files.

Inputs (see ../data/schema for templates)
    cohort.csv      one row per patient: identifiers, outcome and dosimetry
    radiomics.csv   one row per patient, one column per feature

Outputs
    merged_dataset.csv      cohort.csv with the patient key renamed to patient_id
    radiomics_matrix.json   {names, group, feat_class, feat_names, X}

The radiomics column names must be fully qualified as
    <image type>:<feature class>:<feature name>
e.g. original:shape:MeshVolume, wavelet-HHH:glrrm:RunVariance
which is the only way to keep the 1408 feature names unambiguous: 1307 of them
repeat across image types.

Zero-variance features are dropped here, which is why the matrix ends up with
1408 columns rather than the 1409 originally exported.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

COHORT_COLS = [
    'patient_id', 'group', 'tumor', 'age', 'sex', 'stage',
    'base_wbc', 'base_hb', 'base_lym', 'bmi', 'bleed_rt_gap_days',
    'pelvis_volume', 'min_dose', 'max_dose', 'mean_dose',
    'V5', 'V10', 'V20', 'V35', 'V40',
]


def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh))


def main():
    cohort_path = config.data('cohort.csv')
    rad_path = config.data('radiomics.csv')
    for p in (cohort_path, rad_path):
        if not os.path.exists(p):
            sys.exit('missing input: %s\nsee data/schema/ for the expected layout' % p)

    rows = read_csv(cohort_path)
    if not rows:
        sys.exit('cohort.csv is empty')

    missing = [c for c in COHORT_COLS if c not in rows[0]]
    if missing:
        sys.exit('cohort.csv is missing columns: %s' % ', '.join(missing))

    ids = [r['patient_id'].strip() for r in rows]
    if len(set(ids)) != len(ids):
        sys.exit('cohort.csv has duplicate patient_id values')
    if any(not i for i in ids):
        sys.exit('cohort.csv has an empty patient_id')

    groups = [r['group'].strip() for r in rows]
    bad = sorted({g for g in groups} - {'0', '1'})
    if bad:
        sys.exit('group must be 0 or 1, found: %s' % bad)

    # ---- cohort ----
    with open(config.out('merged_dataset.csv'), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.DictWriter(fh, fieldnames=COHORT_COLS, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: (r.get(k) or '').strip() for k in COHORT_COLS})

    # ---- radiomics ----
    rad = read_csv(rad_path)
    rad_ids = [r['patient_id'].strip() for r in rad]
    if len(set(rad_ids)) != len(rad_ids):
        sys.exit('radiomics.csv has duplicate patient_id values')

    missing_ids = [i for i in ids if i not in set(rad_ids)]
    if missing_ids:
        sys.exit('radiomics.csv has no row for: %s' % ', '.join(missing_ids[:10]))

    by_id = {r['patient_id'].strip(): r for r in rad}
    feat_cols = [c for c in rad[0].keys() if c != 'patient_id']

    feat_class, feat_names, matrix = [], [], []
    skipped_empty = 0
    for c in feat_cols:
        parts = c.split(':')
        if len(parts) < 3:
            sys.exit('feature column "%s" is not <image>:<class>:<feature>' % c)
        image, cls, name = parts[0], parts[1], ':'.join(parts[2:])
        vals = []
        ok = True
        for i in ids:
            v = by_id[i].get(c, '')
            if v is None or str(v).strip() == '':
                ok = False
                break
            try:
                vals.append(float(v))
            except ValueError:
                ok = False
                break
        if not ok:
            skipped_empty += 1
            continue
        feat_class.append('%s:%s' % (image, cls))
        feat_names.append(name)
        matrix.append(vals)

    # zero-variance filter, identical to the original pipeline
    kept = []
    for fi, vals in enumerate(matrix):
        m = sum(vals) / len(vals)
        var = sum((x - m) ** 2 for x in vals) / (len(vals) - 1)
        if var > 1e-12:
            kept.append(fi)

    data = {
        'names': ids,
        'group': [int(g) for g in groups],
        'feat_class': [feat_class[i] for i in kept],
        'feat_names': [feat_names[i] for i in kept],
        'X': [[matrix[i][j] for i in kept] for j in range(len(ids))],
    }
    with open(config.out('radiomics_matrix.json'), 'w') as fh:
        json.dump(data, fh)

    summary = (
        'cohort rows        : %d\n'
        'events (group=1)   : %d\n'
        'feature columns in : %d\n'
        'dropped (empty)    : %d\n'
        'dropped (zero var) : %d\n'
        'features kept      : %d\n'
        % (len(ids), sum(data['group']), len(feat_cols), skipped_empty,
           len(matrix) - len(kept), len(kept))
    )
    with open(config.out('dataset_summary.txt'), 'w') as fh:
        fh.write(summary)
    print(summary)


if __name__ == '__main__':
    main()
