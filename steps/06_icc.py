# -*- coding: utf-8 -*-
"""Step 05 - inter-observer reproducibility of the contour-editing step.

Inputs
    icc/observer1.csv   patient_id + one column per feature, first observer
    icc/observer2.csv   same layout, second observer
    icc/subset.csv      patient_id,set   (set is preset20 or extra6)

Output
    icc_results.json    {'preset20': ..., 'all26': ..., 'kept_keys_20': ..., 'kept_keys_26': ...}

Agreement is ICC(A,1) from a two-way analysis of variance - absolute agreement,
single measurement - computed with the closed-form estimators below rather than
a package, because the analysis has to run on the same bare Python 3 the rest of
the pipeline uses.

Both observers re-contour from the same automatic segmentation, so the number
quantifies the editing step and not fully independent delineation. That caveat
is stated in the manuscript as well.
"""
import csv
import json
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

ICC_MIN = 0.75

# features quoted in the manuscript; printed for the reader of the log
TARGETS = [
    ('original', 'shape', 'MeshVolume'),
    ('original', 'shape', 'SurfaceVolumeRatio'),
    ('original', 'shape', 'Maximum2DDiameterRow'),
    ('wavelet-HHH', 'glrlm', 'RunVariance'),
    ('wavelet-HHL', 'ngtdm', 'Busyness'),
    ('log-sigma-1-mm-3D', 'ngtdm', 'Busyness'),
]


def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh))


def load_observer(path):
    """Return (patient_id -> {feature column -> float})."""
    rows = read_csv(path)
    if not rows:
        sys.exit('%s is empty' % path)
    if 'patient_id' not in rows[0]:
        sys.exit('%s has no patient_id column' % path)
    out = {}
    for r in rows:
        pid = (r.get('patient_id') or '').strip()
        if not pid:
            continue
        vals = {}
        for k, v in r.items():
            if k == 'patient_id' or v is None or str(v).strip() == '':
                continue
            try:
                vals[k] = float(v)
            except ValueError:
                continue
        out[pid] = vals
    return out


def icc_a1(x, y):
    """ICC(A,1): two-way random effects, absolute agreement, single rater."""
    n, k = len(x), 2
    if n < 3:
        return float('nan')
    grand = sum(x + y) / (n * k)
    ssr = k * sum(((x[i] + y[i]) / 2 - grand) ** 2 for i in range(n))
    ssc = n * (((sum(x) / n) - grand) ** 2 + ((sum(y) / n) - grand) ** 2)
    sst = sum((v - grand) ** 2 for v in x + y)
    sse = sst - ssr - ssc
    msr = ssr / (n - 1)
    msc = ssc / (k - 1)
    mse = sse / ((n - 1) * (k - 1))
    den = msr + (k - 1) * mse + k * (msc - mse) / n
    return (msr - mse) / den if den > 0 else float('nan')


def load_subset():
    path = config.data('icc', 'subset.csv')
    if not os.path.exists(path):
        sys.exit('missing %s' % path)
    groups = {'preset20': [], 'extra6': []}
    for r in read_csv(path):
        pid = (r.get('patient_id') or '').strip()
        s = (r.get('set') or '').strip().lower()
        if pid and s in groups:
            groups[s].append(pid)
    return groups


def run(obs1, obs2, ids, tag):
    ids = [i for i in ids if i in obs1 and i in obs2]
    if not ids:
        sys.exit('no patient with both observer tables in %s' % tag)

    common = set(obs1[ids[0]]) & set(obs2[ids[0]])
    for i in ids[1:]:
        common &= set(obs1[i]) & set(obs2[i])
    cols = [c for c in obs1[ids[0]] if c in common]

    res = []
    for c in cols:
        a = [obs1[i][c] for i in ids]
        b = [obs2[i][c] for i in ids]
        v = icc_a1(a, b)
        if v == v:  # not NaN
            parts = c.split(':')
            if len(parts) < 3:
                continue
            res.append({'image': parts[0], 'class': parts[1],
                        'feature': ':'.join(parts[2:]), 'icc': v})

    res.sort(key=lambda d: -d['icc'])
    iv = [r['icc'] for r in res]
    ge75 = sum(1 for v in iv if v >= ICC_MIN)
    ge90 = sum(1 for v in iv if v >= 0.90)

    print('\n===== %s: n=%d patients, %d features =====' % (tag, len(ids), len(res)))
    print('median ICC %.3f (range %.3f-%.3f)' % (st.median(iv), min(iv), max(iv)))
    print('>=0.90: %d (%.1f%%); >=0.75: %d (%.1f%%)'
          % (ge90, ge90 / len(iv) * 100, ge75, ge75 / len(iv) * 100))
    by_class = {}
    for r in res:
        by_class.setdefault(r['class'], []).append(r['icc'])
    print('  ' + ' | '.join('%s %.3f' % (c, st.median(v))
                            for c, v in sorted(by_class.items(), key=lambda t: -len(t[1]))))

    key_icc = {}
    for t in TARGETS:
        hit = [r for r in res if (r['image'], r['class'], r['feature']) == t]
        if hit:
            key_icc['%s:%s_%s' % t] = hit[0]['icc']
            print('    %s:%s_%s  ICC %.3f' % (t[0], t[1], t[2], hit[0]['icc']))

    kept = [r for r in res if r['icc'] >= ICC_MIN]
    return {'tag': tag, 'n_patients': len(ids), 'n_features': len(res),
            'median': st.median(iv), 'ge75': ge75, 'ge90': ge90,
            'key_features': key_icc,
            'kept_keys': [(r['image'], r['class'], r['feature']) for r in kept],
            'features': res}


def main():
    obs1_path = config.data('icc', 'observer1.csv')
    obs2_path = config.data('icc', 'observer2.csv')
    for p in (obs1_path, obs2_path):
        if not os.path.exists(p):
            sys.exit('missing input: %s' % p)

    obs1 = load_observer(obs1_path)
    obs2 = load_observer(obs2_path)
    groups = load_subset()
    preset, extra = groups['preset20'], groups['extra6']
    print('observer 1: %d patients | observer 2: %d patients | subset %d + %d'
          % (len(obs1), len(obs2), len(preset), len(extra)))

    a = run(obs1, obs2, preset, 'preset20')
    b = run(obs1, obs2, preset + extra, 'all26')

    with open(config.out('icc_results.json'), 'w', encoding='utf-8') as fh:
        json.dump({'preset20': {k: v for k, v in a.items() if k != 'kept_keys'},
                   'all26': {k: v for k, v in b.items() if k != 'kept_keys'},
                   'kept_keys_20': a['kept_keys'], 'kept_keys_26': b['kept_keys']},
                  fh, ensure_ascii=False)
    print('\nwritten icc_results.json')


if __name__ == '__main__':
    main()
