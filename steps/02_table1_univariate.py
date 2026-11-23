# -*- coding: utf-8 -*-
"""Phase 1: Table 1 (baseline) + dosimetric univariate analysis + ROC."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import csv, json, sys
from lib.stats_lib import *

OUT = config.OUT

rows = []
with open(OUT + 'merged_dataset.csv', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        rows.append(r)

def getf(r, k):
    v = r.get(k, '')
    if v in ('', None): return None
    return float(v)

ev = [r for r in rows if r['group'] == '1']
ne = [r for r in rows if r['group'] == '0']
print(f'事件组 n={len(ev)}, 对照组 n={len(ne)}')

results = {'table1': [], 'dose_univar': []}

# ---------- Table 1: continuous ----------
cont_vars = [('age', '年龄 (岁)'), ('bmi', 'BMI (kg/m²)'), ('base_wbc', '基线 WBC (×10⁹/L)'),
             ('base_hb', '基线 Hb (g/L)'), ('base_lym', '基线淋巴细胞 (×10⁹/L)'),
             ('bleed_rt_gap_days', '基线采血距放疗 (天)')]
for key, label in cont_vars:
    a = [getf(r, key) for r in ev]; a = [x for x in a if x is not None]
    b = [getf(r, key) for r in ne]; b = [x for x in b if x is not None]
    if len(a) < 3 or len(b) < 3:
        continue
    t, p_t = ttest_welch(a, b)
    U, p_mw = mann_whitney(a, b)
    q1a, q3a = iqr(a); q1b, q3b = iqr(b)
    results['table1'].append({
        'var': label, 'n_ev': len(a), 'n_ne': len(b),
        'ev': f'{mean(a):.2f}±{sd(a):.2f}', 'ne': f'{mean(b):.2f}±{sd(b):.2f}',
        'ev_med': f'{median(a):.2f} ({q1a:.2f}-{q3a:.2f})', 'ne_med': f'{median(b):.2f} ({q1b:.2f}-{q3b:.2f})',
        'p_t': p_t, 'p_mw': p_mw})

# ---------- Table 1: categorical ----------
def cat2x2(key, levels=None):
    from collections import Counter
    ce = Counter(r[key] for r in ev); cn = Counter(r[key] for r in ne)
    lv = levels or sorted(set(ce) | set(cn))
    return lv, ce, cn

# sex
lv, ce, cn = cat2x2('sex')
tab = [[ce.get('男', 0), ce.get('女', 0)], [cn.get('男', 0), cn.get('女', 0)]]
p, method = chi2_or_fisher(tab)
results['table1'].append({'var': '性别(男/女)', 'ev': f"{tab[0][0]}/{tab[0][1]}", 'ne': f"{tab[1][0]}/{tab[1][1]}",
                          'p_cat': p, 'method': method})
# tumor (3 levels -> chi2 on 2x3, or Fisher-Freeman-Halton too heavy; use chi2 if expected ok, else note)
lv, ce, cn = cat2x2('tumor')
obs = [[ce.get(l, 0) for l in lv], [cn.get(l, 0) for l in lv]]
N = sum(sum(r) for r in obs)
rowsum = [sum(r) for r in obs]; colsum = [sum(obs[i][j] for i in range(2)) for j in range(len(lv))]
exp = [[rowsum[i] * colsum[j] / N for j in range(len(lv))] for i in range(2)]
chi = sum((obs[i][j] - exp[i][j]) ** 2 / exp[i][j] for i in range(2) for j in range(len(lv)) if exp[i][j] > 0)
p_tumor = chi2_sf(chi, len(lv) - 1)
results['table1'].append({'var': '瘤种 (' + '/'.join(lv) + ')',
                          'ev': '/'.join(str(obs[0][j]) for j in range(len(lv))),
                          'ne': '/'.join(str(obs[1][j]) for j in range(len(lv))),
                          'p_cat': p_tumor, 'method': 'chi2' if min(min(e) for e in exp) >= 5 else 'chi2(期望频数偏小)'})
# stage: collapse I-II vs III-IV
def stage_bin(s):
    s = str(s).replace('Ⅰ', '1').replace('Ⅱ', '2').replace('Ⅲ', '3').replace('Ⅳ', '4')
    try: return 'III-IV' if int(float(s)) >= 3 else 'I-II'
    except ValueError: return None
se_ = [stage_bin(r['stage']) for r in ev]; sn = [stage_bin(r['stage']) for r in ne]
tab = [[sum(1 for x in se_ if x == 'III-IV'), sum(1 for x in se_ if x == 'I-II')],
       [sum(1 for x in sn if x == 'III-IV'), sum(1 for x in sn if x == 'I-II')]]
p, method = chi2_or_fisher(tab)
results['table1'].append({'var': '分期(III-IV/I-II)', 'ev': f"{tab[0][0]}/{tab[0][1]}", 'ne': f"{tab[1][0]}/{tab[1][1]}",
                          'p_cat': p, 'method': method})

# ---------- Dosimetric univariate + ROC ----------
dose_vars = [('pelvis_volume', '骨盆体积 (cm³)'), ('min_dose', '最小剂量 (cGy)'), ('max_dose', '最大剂量 (cGy)'),
             ('mean_dose', '平均剂量 (cGy)'), ('V5', 'V5 (%)'), ('V10', 'V10 (%)'),
             ('V20', 'V20 (%)'), ('V35', 'V35 (%)'), ('V40', 'V40 (%)')]
y_all = [int(r['group']) for r in rows]
for key, label in dose_vars:
    a = [getf(r, key) for r in ev]; b = [getf(r, key) for r in ne]
    t, p_t = ttest_welch(a, b)
    U, p_mw = mann_whitney(a, b)
    # univariate logistic
    X = [[getf(r, key)] for r in rows]
    beta, se_b, conv = logistic_fit(X, y_all)
    OR = math.exp(beta[1]); lo = math.exp(beta[1] - 1.96 * se_b[1]); hi = math.exp(beta[1] + 1.96 * se_b[1])
    z = beta[1] / se_b[1]; p_or = 2 * (1 - norm_cdf(abs(z)))
    scores = logit_predict(beta, X)
    A = auc(scores, y_all)
    ci = bootstrap_auc_ci(scores, y_all, B=2000)
    results['dose_univar'].append({'var': label, 'key': key,
        'ev': f'{mean(a):.2f}±{sd(a):.2f}', 'ne': f'{mean(b):.2f}±{sd(b):.2f}',
        'p_mw': p_mw, 'OR': OR, 'OR_lo': lo, 'OR_hi': hi, 'p_or': p_or,
        'auc': A, 'auc_lo': ci[0], 'auc_hi': ci[1]})

# ---------- baseline labs univariate logistic (complete cases) ----------
lab_vars = [('age', '年龄'), ('bmi', 'BMI'), ('base_wbc', '基线 WBC'), ('base_hb', '基线 Hb'), ('base_lym', '基线淋巴细胞')]
results['clin_univar'] = []
for key, label in lab_vars:
    sub = [(getf(r, key), int(r['group'])) for r in rows]
    sub = [(x, y) for x, y in sub if x is not None]
    X = [[x] for x, y in sub]; yy = [y for x, y in sub]
    beta, se_b, conv = logistic_fit(X, yy)
    OR = math.exp(beta[1]); lo = math.exp(beta[1] - 1.96 * se_b[1]); hi = math.exp(beta[1] + 1.96 * se_b[1])
    z = beta[1] / se_b[1]; p_or = 2 * (1 - norm_cdf(abs(z)))
    A = auc(logit_predict(beta, X), yy)
    results['clin_univar'].append({'var': label, 'key': key, 'n': len(sub),
        'OR': OR, 'OR_lo': lo, 'OR_hi': hi, 'p_or': p_or, 'auc': A})

with open(OUT + 'phase1_results.json', 'w') as f:
    json.dump(results, f, ensure_ascii=False, indent=1)

# ---------- print summary ----------
print('\n=== Table 1 (连续变量, mean±SD, P 为 Mann-Whitney) ===')
for r in results['table1']:
    if 'ev_med' in r:
        print(f"{r['var']:24s} 事件{n_ev if False else r['n_ev']}: {r['ev']:>16s} | 对照{r['n_ne']}: {r['ne']:>16s} | P_t={r['p_t']:.4f} P_MW={r['p_mw']:.4f}")
    else:
        print(f"{r['var']:24s} 事件: {r['ev']:>10s} | 对照: {r['ne']:>10s} | P={r['p_cat']:.4f} ({r['method']})")
print('\n=== 剂量学单因素 ===')
for r in results['dose_univar']:
    print(f"{r['var']:16s} {r['ev']:>14s} vs {r['ne']:>14s} | P_MW={r['p_mw']:.4f} | OR={r['OR']:.4f} ({r['OR_lo']:.4f}-{r['OR_hi']:.4f}) P={r['p_or']:.4f} | AUC={r['auc']:.3f} ({r['auc_lo']:.3f}-{r['auc_hi']:.3f})")
print('\n=== 临床/基线单因素 ===')
for r in results['clin_univar']:
    print(f"{r['var']:12s} n={r['n']} OR={r['OR']:.4f} ({r['OR_lo']:.4f}-{r['OR_hi']:.4f}) P={r['p_or']:.4f} AUC={r['auc']:.3f}")
