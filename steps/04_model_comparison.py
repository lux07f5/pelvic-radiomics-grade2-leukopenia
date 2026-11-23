# -*- coding: utf-8 -*-
"""Phase 3: model comparison (Rad-score / pelvis volume / combined), ROC, optimism correction.
Also traces selected radiomics features back to their row index in the source table."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, csv, math, random, sys
from lib.stats_lib import *

OUT = config.OUT
p2 = json.load(open(OUT + 'phase2_results.json'))
D = json.load(open(OUT + 'radiomics_matrix.json'))
rows = list(csv.DictReader(open(OUT + 'merged_dataset.csv', encoding='utf-8-sig')))
y = p2['y']
n = len(y)
names = D['names']
rowmap = {r['patient_id']: r for r in rows}

# ---- label rad features with source row (matrix_row saved by phase2) ----
for rf in p2['rad_features']:
    best = rf['matrix_row']
    rf['source_row'] = best + 2  # row in the source feature table (1-based, header included)
    rf['label'] = f"{rf['class']}_{rf['name']} [table row {best+2}]"

# recompute rad linear predictor from raw values + coefficients
rad_lin = []
for i in range(n):
    e = p2['rad_intercept'] + sum(rf['coef'] * D['X'][i][rf['matrix_row']] for rf in p2['rad_features'])
    rad_lin.append(e)
rad_prob = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in rad_lin]
A_rad = auc(rad_prob, y)
assert abs(A_rad - p2['rad_auc']) < 1e-6

v5 = [float(rowmap[nm]['V5']) for nm in names]
pv = [float(rowmap[nm]['pelvis_volume']) for nm in names]

def fit_auc(X, yv):
    beta, se, conv = logistic_fit(X, yv)
    pr = logit_predict(beta, X)
    A = auc(pr, yv)
    lo, hi = bootstrap_auc_ci(pr, yv, B=2000)
    return beta, pr, A, lo, hi

def optimism_correct(X, yv, B=500, seed=7):
    rnd = random.Random(seed)
    opts = []
    for _ in range(B):
        idx = [rnd.randrange(len(yv)) for _ in range(len(yv))]
        yy = [yv[i] for i in idx]
        if sum(yy) in (0, len(yy)): continue
        Xb = [X[i] for i in idx]
        bb, _, cv = logistic_fit(Xb, yy, max_iter=60)
        if not cv: continue
        opts.append(auc(logit_predict(bb, Xb), yy) - auc(logit_predict(bb, X), yv))
    return mean(opts)

models = []
# Model B: pelvis volume
betaB, prB, AB, loB, hiB = fit_auc([[v] for v in pv], y)
optB = optimism_correct([[v] for v in pv], y)
models.append({'label': '骨盆体积 (TPS)', 'nvar': 1, 'auc': AB, 'lo': loB, 'hi': hiB,
               'auc_adj': AB - optB, 'roc': roc_curve(prB, y)})
# Model A: rad-score
optA = p2['optimism']
_nrad = len(p2['rad_features'])
_lam_note = 'λ.1se' if len(p2['nz_1se']) >= 1 else 'λ.min(λ.1se为空模型)'
models.append({'label': f'Rad-score ({_nrad} 特征, {_lam_note})', 'nvar': _nrad, 'auc': A_rad, 'lo': p2['rad_auc_lo'], 'hi': p2['rad_auc_hi'],
               'auc_adj': p2['rad_auc_adj'], 'roc': roc_curve(rad_prob, y)})
# Model C: rad linear predictor + V5
v5z = [(v - mean(v5)) / sd(v5) for v in v5]
XC = [[rad_lin[i], v5z[i]] for i in range(n)]
betaC, prC, AC, loC, hiC = fit_auc(XC, y)
optC = optimism_correct(XC, y)
models.append({'label': '联合模型 (Rad-score + V5)', 'nvar': 2, 'auc': AC, 'lo': loC, 'hi': hiC,
               'auc_adj': AC - optC, 'roc': roc_curve(prC, y)})
# reference: V5 alone
betaV, prV, AV, loV, hiV = fit_auc([[v] for v in v5], y)
models.append({'label': 'V5 (参考)', 'nvar': 1, 'auc': AV, 'lo': loV, 'hi': hiV,
               'auc_adj': '-', 'roc': roc_curve(prV, y)})

print(f"{'模型':<26s} AUC (95%CI)            校正AUC")
for m in models:
    adj = f"{m['auc_adj']:.3f}" if isinstance(m['auc_adj'], float) else m['auc_adj']
    print(f"{m['label']:<26s} {m['auc']:.3f} ({m['lo']:.3f}-{m['hi']:.3f})   {adj}")
print('\n联合模型系数: rad_lin beta=%.4f, V5z beta=%.4f, intercept=%.4f' % (betaC[1], betaC[2], betaC[0]))
print('\nRad-score 特征溯源:')
for rf in p2['rad_features']:
    print(f"  {rf['label']}  coef={rf['coef']:+.6f}")

json.dump({'models': models, 'rad_features': p2['rad_features'],
           'combined': {'beta_rad': betaC[1], 'beta_v5z': betaC[2], 'intercept': betaC[0]}},
          open(OUT + 'phase3_results.json', 'w'), ensure_ascii=False, indent=1)
# also patch phase2 rad_features with labels
json.dump(p2, open(OUT + 'phase2_results.json', 'w'), ensure_ascii=False, indent=1)
print('\n保存 phase3_results.json (并回填 phase2 特征标签)')
