# -*- coding: utf-8 -*-
"""Phase 2: radiomics univariate filter -> LASSO (10-fold CV, lambda.1se) -> Rad-score."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, sys, math, random
from lib.stats_lib import *

OUT = config.OUT
with open(OUT + 'radiomics_matrix.json') as f:
    D = json.load(f)
names, y = D['names'], D['group']
feat_names, feat_class, X = D['feat_names'], D['feat_class'], D['X']
n, p = len(y), len(feat_names)
print(f'队列 n={n} (事件 {sum(y)}), 特征 p={p}')

# ---------- Step 1: univariate Mann-Whitney filter P<0.05 ----------
cand = []
for j in range(p):
    a = [X[i][j] for i in range(n) if y[i] == 1]
    b = [X[i][j] for i in range(n) if y[i] == 0]
    U, pv = mann_whitney(a, b)
    if pv == pv and pv < 0.05:  # not NaN
        A = auc([X[i][j] for i in range(n)], y)
        cand.append({'idx': j, 'class': feat_class[j], 'name': feat_names[j], 'p': pv, 'auc': A})
print(f'单因素 P<0.05: {len(cand)} 个特征')
cand.sort(key=lambda d: d['p'])
for c in cand[:15]:
    print(f"  {c['class']:<12s} {c['name']:<28s} P={c['p']:.5f} AUC={c['auc']:.3f}")

# ---------- Step 2: remove highly redundant features (|Spearman|>0.9 keep smaller p) ----------
def spearman(xa, xb):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
            for k in range(i, j + 1): r[order[k]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    ra, rb = rank(xa), rank(xb)
    ma, mb = mean(ra), mean(rb)
    cov = sum((ra[i] - ma) * (rb[i] - mb) for i in range(len(ra)))
    va = sum((v - ma) ** 2 for v in ra); vb = sum((v - mb) ** 2 for v in rb)
    return cov / math.sqrt(va * vb) if va > 0 and vb > 0 else 0.0

sel = []
for c in cand:  # cand sorted by p asc; keep first of each correlated cluster
    redundant = False
    for s in sel:
        r = spearman([X[i][c['idx']] for i in range(n)], [X[i][s['idx']] for i in range(n)])
        if abs(r) > 0.9:
            redundant = True
            break
    if not redundant:
        sel.append(c)
print(f'去冗余(|rho|>0.9)后: {len(sel)} 个特征')
for c in sel:
    print(f"  {c['class']:<12s} {c['name']:<28s} P={c['p']:.5f}")

# cap at 30 features for LASSO stability if too many
MAXF = 30
if len(sel) > MAXF:
    sel = sel[:MAXF]
    print(f'(按 P 值截断至前 {MAXF} 个进入 LASSO)')

Xs = [[X[i][c['idx']] for c in sel] for i in range(n)]
# standardize full-sample for final fit
mu = [mean([r[j] for r in Xs]) for j in range(len(sel))]
sdd = [sd([r[j] for r in Xs]) or 1.0 for j in range(len(sel))]
Xs_std = [[(r[j] - mu[j]) / sdd[j] for j in range(len(sel))] for r in Xs]

# ---------- Step 3: LASSO 10-fold CV ----------
cv = lasso_cv(Xs_std, y, nfolds=10, nlam=100, lam_min_ratio=0.01, seed=42)
print(f"\nLASSO CV: lambda.min={cv['lam_min']:.5f} (idx {cv['i_min']}), lambda.1se={cv['lam_1se']:.5f} (idx {cv['i_1se']})")

# fit at both lambdas on full data
b0_min, b_min = lasso_logistic_fit(Xs_std, y, cv['lam_min'], max_iter=500, tol=1e-7)
b0_1se, b_1se = lasso_logistic_fit(Xs_std, y, cv['lam_1se'], max_iter=500, tol=1e-7)
nz_min = [(sel[j], b_min[j]) for j in range(len(sel)) if abs(b_min[j]) > 1e-8]
nz_1se = [(sel[j], b_1se[j]) for j in range(len(sel)) if abs(b_1se[j]) > 1e-8]
print(f'lambda.min 非零特征 {len(nz_min)} 个; lambda.1se 非零特征 {len(nz_1se)} 个')
print('--- lambda.1se 选中特征 (标准化系数) ---')
for c, b in nz_1se:
    print(f"  {c['class']:<12s} {c['name']:<28s} beta={b:+.4f}  (单因素P={c['p']:.5f})")

# ---------- Step 4: Rad-score via plain logistic on selected features ----------
sel_use = nz_1se if len(nz_1se) >= 1 else nz_min
use_idx = [sel.index(c) for c, b in sel_use]
Xr = [[Xs[i][j] for j in use_idx] for i in range(n)]
beta_r, se_r, conv = logistic_fit(Xr, y)
print(f'\nRad-score logistic (原始量纲): converged={conv}')
rad_formula = []
for k, j in enumerate(use_idx):
    rad_formula.append({'class': sel[j]['class'], 'name': sel[j]['name'],
                        'coef': beta_r[k + 1], 'se': se_r[k + 1], 'matrix_row': sel[j]['idx']})
    print(f"  {sel[j]['name']:<28s} coef={beta_r[k+1]:+.6f} (SE {se_r[k+1]:.4f})")
print(f'  intercept = {beta_r[0]:+.6f}')
rad_scores = logit_predict(beta_r, Xr)
A_rad = auc(rad_scores, y)
ci_rad = bootstrap_auc_ci(rad_scores, y, B=2000)
print(f'Rad-score AUC = {A_rad:.3f} (95%CI {ci_rad[0]:.3f}-{ci_rad[1]:.3f})')

# ---------- Step 5: bootstrap optimism correction ----------
rnd = random.Random(42)
B = 1000
optimisms = []
nfail = 0
for bi in range(B):
    idx = [rnd.randrange(n) for _ in range(n)]
    if sum(y[i] for i in idx) in (0, len(idx)):
        nfail += 1
        continue
    Xb = [[Xr[i][j] for j in range(len(use_idx))] for i in idx]
    yb = [y[i] for i in idx]
    bb, _, cvv = logistic_fit(Xb, yb, max_iter=60)
    if not cvv:
        nfail += 1
        continue
    auc_boot = auc(logit_predict(bb, Xb), yb)
    auc_orig = auc(logit_predict(bb, Xr), y)
    optimisms.append(auc_boot - auc_orig)
opt = mean(optimisms)
A_adj = A_rad - opt
print(f'Bootstrap 乐观度 = {opt:.4f} ({len(optimisms)} 次成功, {nfail} 次失败), 校正后 AUC = {A_adj:.3f}')

# ---------- save ----------
out = {'univariate_hits': [{'class': c['class'], 'name': c['name'], 'p': c['p'], 'auc': c['auc'], 'matrix_row': c['idx']} for c in cand],
       'after_dedup': [{'class': c['class'], 'name': c['name'], 'p': c['p'], 'matrix_row': c['idx']} for c in sel],
       'lambda_min': cv['lam_min'], 'lambda_1se': cv['lam_1se'],
       'cv_lams': cv['lams'], 'cv_mean': cv['cv_mean'], 'i_min': cv['i_min'], 'i_1se': cv['i_1se'],
       'nz_min': [{'class': c['class'], 'name': c['name'], 'beta_std': b, 'matrix_row': c['idx']} for c, b in nz_min],
       'nz_1se': [{'class': c['class'], 'name': c['name'], 'beta_std': b, 'matrix_row': c['idx']} for c, b in nz_1se],
       'rad_features': rad_formula, 'rad_intercept': beta_r[0],
       'rad_auc': A_rad, 'rad_auc_lo': ci_rad[0], 'rad_auc_hi': ci_rad[1],
       'optimism': opt, 'rad_auc_adj': A_adj,
       'rad_scores': rad_scores, 'y': y, 'names': names,
       'selected_std': {'use_idx': use_idx, 'mu': [mu[j] for j in use_idx], 'sd': [sdd[j] for j in use_idx]}}
with open(OUT + 'phase2_results.json', 'w') as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print('\n保存 phase2_results.json')
