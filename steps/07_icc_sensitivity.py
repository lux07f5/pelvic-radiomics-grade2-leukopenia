# -*- coding: utf-8 -*-
"""Phase 5: 敏感性分析 —— 仅保留 ICC>=0.75 的特征，重跑 单因素→去冗余→LASSO→Rad-score
用法: phase5_sensitivity.py {preset20|all26}
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, sys, math, random
from lib.stats_lib import *

OUT = config.OUT
TAG = sys.argv[1] if len(sys.argv) > 1 else 'all26'
ICC_MIN = 0.75

D = json.load(open(OUT + 'radiomics_matrix.json'))
ICC = json.load(open(OUT + 'icc_results.json'))
names, y = D['names'], D['group']
feat_names, feat_class, X = D['feat_names'], D['feat_class'], D['X']
n, p = len(y), len(feat_names)

kept = {(k[0], k[1], k[2]) for k in ICC['kept_keys_20' if TAG == 'preset20' else 'kept_keys_26']}
keep_idx = []
for j in range(p):
    img, cls = feat_class[j].split(':', 1)
    if (img, cls, feat_names[j]) in kept:
        keep_idx.append(j)
print(f'[{TAG}] ICC>={ICC_MIN} 保留特征 {len(keep_idx)}/{p}')


def spearman(xa, xb):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    ra, rb = rank(xa), rank(xb)
    ma, mb = mean(ra), mean(rb)
    cov = sum((ra[i] - ma) * (rb[i] - mb) for i in range(len(ra)))
    va = sum((v - ma) ** 2 for v in ra)
    vb = sum((v - mb) ** 2 for v in rb)
    return cov / math.sqrt(va * vb) if va > 0 and vb > 0 else 0.0


cand = []
for j in keep_idx:
    a = [X[i][j] for i in range(n) if y[i] == 1]
    b = [X[i][j] for i in range(n) if y[i] == 0]
    _, pv = mann_whitney(a, b)
    if pv == pv and pv < 0.05:
        cand.append({'idx': j, 'class': feat_class[j], 'name': feat_names[j], 'p': pv,
                     'auc': auc([X[i][j] for i in range(n)], y)})
cand.sort(key=lambda d: d['p'])
print(f'单因素 P<0.05: {len(cand)}')

sel = []
for c in cand:
    if all(abs(spearman([X[i][c['idx']] for i in range(n)], [X[i][s['idx']] for i in range(n)])) <= 0.9
           for s in sel):
        sel.append(c)
print(f'去冗余后: {len(sel)}')
sel = sel[:30]
print(f'进入 LASSO: {len(sel)}')

Xs = [[X[i][c['idx']] for c in sel] for i in range(n)]
mu = [mean([r[j] for r in Xs]) for j in range(len(sel))]
sdd = [sd([r[j] for r in Xs]) or 1.0 for j in range(len(sel))]
Xs_std = [[(r[j] - mu[j]) / sdd[j] for j in range(len(sel))] for r in Xs]

cv = lasso_cv(Xs_std, y, nfolds=10, nlam=100, lam_min_ratio=0.01, seed=42)
b0_min, b_min = lasso_logistic_fit(Xs_std, y, cv['lam_min'], max_iter=500, tol=1e-7)
b0_1se, b_1se = lasso_logistic_fit(Xs_std, y, cv['lam_1se'], max_iter=500, tol=1e-7)
nz_min = [(sel[j], b_min[j]) for j in range(len(sel)) if abs(b_min[j]) > 1e-8]
nz_1se = [(sel[j], b_1se[j]) for j in range(len(sel)) if abs(b_1se[j]) > 1e-8]
print(f'lambda.min 非零 {len(nz_min)}; lambda.1se 非零 {len(nz_1se)}')

use = nz_1se if nz_1se else nz_min
use_idx = [sel.index(c) for c, b in use]
Xr = [[Xs[i][j] for j in use_idx] for i in range(n)]
beta_r, se_r, conv = logistic_fit(Xr, y)
rad = logit_predict(beta_r, Xr)
A = auc(rad, y)
ci = bootstrap_auc_ci(rad, y, B=2000)
rnd = random.Random(42)
opt = []
for _ in range(1000):
    idx = [rnd.randrange(n) for _ in range(n)]
    if sum(y[i] for i in idx) in (0, n):
        continue
    bb, _, cvv = logistic_fit([[Xr[i][j] for j in range(len(use_idx))] for i in idx],
                              [y[i] for i in idx], max_iter=60)
    if not cvv:
        continue
    opt.append(auc(logit_predict(bb, [[Xr[i][j] for j in range(len(use_idx))] for i in idx]),
                   [y[i] for i in idx]) - auc(logit_predict(bb, Xr), y))
A_adj = A - mean(opt)
print(f'\n[{TAG}] Rad-score AUC={A:.3f} ({ci[0]:.3f}-{ci[1]:.3f}), 校正 {A_adj:.3f}')
print('入选特征:')
for c, b in use:
    print(f"  {c['class']:<22s} {c['name']:<24s} beta_std={b:+.4f} P={c['p']:.5f}")

json.dump({'tag': TAG, 'n_kept': len(keep_idx), 'n_univariate': len(cand), 'n_after_dedup': len(sel),
           'lambda_min': cv['lam_min'], 'lambda_1se': cv['lam_1se'],
           'n_nz_min': len(nz_min), 'n_nz_1se': len(nz_1se),
           'features': [{'class': c['class'], 'name': c['name'], 'p': c['p'],
                         'beta_std': b, 'coef': beta_r[k + 1]} for k, (c, b) in enumerate(use)],
           'intercept': beta_r[0], 'auc': A, 'lo': ci[0], 'hi': ci[1], 'auc_adj': A_adj,
           'rad_scores': rad, 'y': y},
          open(OUT + f'phase5_{TAG}.json', 'w'), ensure_ascii=False)
print(f'保存 phase5_{TAG}.json')
