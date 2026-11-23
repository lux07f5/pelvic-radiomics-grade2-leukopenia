# -*- coding: utf-8 -*-
"""Phase 6: pelvis-volume-first (main line) analysis.
- pelvis volume alone: OR (per cm3 and per 100 cm3), AUC, optimism, OOF calibration
- pelvis volume + Rad-score: incremental value, DeLong, IDI, continuous NRI
- pelvis volume cut-off (Youden), tertile event rates, operating characteristics
- DCA for pelvis / Rad-score / combined
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, csv, math, random, sys
from lib.stats_lib import *

OUT = config.OUT
p2 = json.load(open(OUT + 'phase2_results.json'))
D = json.load(open(OUT + 'radiomics_matrix.json'))
rows = list(csv.DictReader(open(OUT + 'merged_dataset.csv', encoding='utf-8-sig')))
y = p2['y']; n = len(y); names = D['names']
rowmap = {r['patient_id']: r for r in rows}
prev = sum(y) / n

pv = [float(rowmap[nm]['pelvis_volume']) for nm in names]
_mv = next(j for j in range(len(D['feat_names']))
           if D['feat_class'][j] == 'original:shape' and D['feat_names'][j] == 'MeshVolume')
mv = [D['X'][i][_mv] for i in range(n)]                       # original:shape MeshVolume
rad_lin = [p2['rad_intercept'] + sum(rf['coef'] * D['X'][i][rf['matrix_row']] for rf in p2['rad_features']) for i in range(n)]

def probs(lin): return [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in lin]

# ---------------- DeLong ----------------
def ty_sum_le(x, neg):
    # number of negatives < x plus 0.5 * equal
    c = 0.0
    for t in neg:
        if t < x: c += 1.0
        elif t == x: c += 0.5
    return c

def pos_sum_le(x, pos):
    c = 0.0
    for t in pos:
        if t < x: c += 1.0
        elif t == x: c += 0.5
    return c

def delong_paired(sa, sb, yv):
    """Two-sided DeLong test for paired ROC curves."""
    pa = [s for s, t in zip(sa, yv) if t == 1]; na = [s for s, t in zip(sa, yv) if t == 0]
    pb = [s for s, t in zip(sb, yv) if t == 1]; nb = [s for s, t in zip(sb, yv) if t == 0]
    m, k = len(pa), len(na)
    def place(pos, neg):
        V = [ty_sum_le(p, neg) / len(neg) for p in pos]
        H = [pos_sum_le(q, pos) / len(pos) for q in neg]
        return V, H
    V1, H1 = place(pa, na)
    V2, H2 = place(pb, nb)
    A1 = sum(V1) / m; A2 = sum(V2) / m
    var1 = sum((v - A1) ** 2 for v in V1) / (m * (m - 1)) + sum((h - A1) ** 2 for h in H1) / (k * (k - 1))
    var2 = (sum((v - A2) ** 2 for v in V2) / (m * (m - 1)) + sum((h - A2) ** 2 for h in H2) / (k * (k - 1)))
    cov = (sum((V1[i] - A1) * (V2[i] - A2) for i in range(m)) / (m * (m - 1))
           + sum((H1[i] - A1) * (H2[i] - A2) for i in range(k)) / (k * (k - 1)))
    var_d = var1 + var2 - 2 * cov
    if var_d <= 0: return A1, A2, float('nan'), None
    z = (A1 - A2) / math.sqrt(var_d)
    p = 2.0 * (1.0 - norm_cdf(abs(z)))
    return A1, A2, p, (A1 - A2, math.sqrt(var_d))

# ---------------- models ----------------
def fit(X):
    b, se, cv = logistic_fit(X, y)
    pr = logit_predict(b, X)
    return b, se, pr

def optimism(X, B=1000, seed=7):
    rnd = random.Random(seed); opts = []
    for _ in range(B):
        idx = [rnd.randrange(n) for _ in range(n)]
        yy = [y[i] for i in idx]
        if sum(yy) in (0, n): continue
        Xb = [X[i] for i in idx]
        bb, _, cvv = logistic_fit(Xb, yy, max_iter=60)
        if not cvv: continue
        opts.append(auc(logit_predict(bb, Xb), yy) - auc(logit_predict(bb, X), y))
    return mean(opts)

def oof_lin(X, nfolds=10, seed=42):
    rnd = random.Random(seed)
    pos = [i for i in range(n) if y[i] == 1]; neg = [i for i in range(n) if y[i] == 0]
    rnd.shuffle(pos); rnd.shuffle(neg)
    folds = [[] for _ in range(nfolds)]
    for k, i in enumerate(pos): folds[k % nfolds].append(i)
    for k, i in enumerate(neg): folds[k % nfolds].append(i)
    lin = [None] * n
    for f in range(nfolds):
        test = folds[f]
        Xtr = [X[i] for i in range(n) if i not in test]; ytr = [y[i] for i in range(n) if i not in test]
        b, _, _ = logistic_fit(Xtr, ytr)
        for i in test:
            lin[i] = b[0] + sum(b[j + 1] * X[i][j] for j in range(len(X[0])))
    return lin

def calibration(lin, yv, g=5):
    pr = [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in lin]
    b, se, cv = logistic_fit([[l] for l in lin], yv); slope = b[1]
    a = 0.0
    for _ in range(200):
        mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, a + l)))) for l in lin]
        gr = sum(yv[i] - mu[i] for i in range(len(yv)))
        he = sum(max(m * (1 - m), 1e-9) for m in mu)
        da = gr / he; a += da
        if abs(da) < 1e-9: break
    order = sorted(range(len(pr)), key=lambda i: pr[i])
    groups = [order[int(i * len(pr) / g):int((i + 1) * len(pr) / g)] for i in range(g)]
    HL = 0.0; bins = []
    for grp in groups:
        O = sum(yv[i] for i in grp); E = sum(pr[i] for i in grp); ng = len(grp)
        if E > 0.5 and ng - E > 0.5:
            HL += (O - E) ** 2 / (E * (1 - E / ng))
        bins.append({'pred': E / ng, 'obs': O / ng, 'n': ng})
    return {'slope': slope, 'intercept': a, 'HL': HL, 'HL_p': chi2_sf(HL, g - 2), 'bins': bins}

X_pv = [[v] for v in pv]
X_rad = [[l] for l in rad_lin]
X_both = [[pv[i], rad_lin[i]] for i in range(n)]

res = {}
b_pv, se_pv, pr_pv = fit(X_pv)
b_rad, se_rad, pr_rad = fit(X_rad)
b_both, se_both, pr_both = fit(X_both)

for tag, b, se, pr, X in (('pelvis', b_pv, se_pv, pr_pv, X_pv),
                          ('rad', b_rad, se_rad, pr_rad, X_rad),
                          ('both', b_both, se_both, pr_both, X_both)):
    A = auc(pr, y); lo, hi = bootstrap_auc_ci(pr, y, B=2000)
    opt = optimism(X)
    cal = calibration(oof_lin(X), y)
    res[tag] = {'beta': b, 'se': se, 'auc': A, 'lo': lo, 'hi': hi,
                'auc_adj': A - opt, 'optimism': opt, 'cal': cal,
                'roc': roc_curve(pr, y), 'probs': pr}
    print(f"[{tag}] AUC={A:.4f} ({lo:.3f}-{hi:.3f}) adj={A-opt:.4f} | OOF slope={cal['slope']:.3f} int={cal['intercept']:+.3f} HL P={cal['HL_p']:.3f}")
    for j, nm in enumerate(['(Intercept)'] + [f'x{j}' for j in range(len(X[0]))]):
        print(f"    {nm}: beta={b[j]:+.6f} se={se[j]:.6f} OR={math.exp(b[j]):.6f} P={2*(1-norm_cdf(abs(b[j]/se[j]))):.4g}")

# per-100 cm3 OR
print('\n骨盆体积 per 100 cm3: OR=%.3f (%.3f-%.3f) P=%.4g' % (
    math.exp(b_pv[1] * 100), math.exp((b_pv[1] - 1.96 * se_pv[1]) * 100), math.exp((b_pv[1] + 1.96 * se_pv[1]) * 100),
    2 * (1 - norm_cdf(abs(b_pv[1] / se_pv[1])))))

# ---------------- DeLong ----------------
dl = {}
for a, b, la, lb in ((pr_both, pr_pv, 'both', 'pelvis'), (pr_both, pr_rad, 'both', 'rad'), (pr_rad, pr_pv, 'rad', 'pelvis')):
    A1, A2, p, dd = delong_paired(a, b, y)
    dl[f'{la}_vs_{lb}'] = {'A1': A1, 'A2': A2, 'diff': dd[0], 'se': dd[1], 'p': p}
    print(f'DeLong {la} vs {lb}: dAUC={dd[0]:+.4f} (SE {dd[1]:.4f}) P={p:.3f}')

# ---------------- IDI / continuous NRI ----------------
def idi_nri(p_old, p_new, yv):
    e = [i for i in range(len(yv)) if yv[i] == 1]; ne = [i for i in range(len(yv)) if yv[i] == 0]
    IS_old = mean([p_old[i] for i in e]) - mean([p_old[i] for i in ne])
    IS_new = mean([p_new[i] for i in e]) - mean([p_new[i] for i in ne])
    up_e = mean([1.0 if p_new[i] > p_old[i] else 0.0 for i in e])
    dn_e = mean([1.0 if p_new[i] < p_old[i] else 0.0 for i in e])
    up_n = mean([1.0 if p_new[i] > p_old[i] else 0.0 for i in ne])
    dn_n = mean([1.0 if p_new[i] < p_old[i] else 0.0 for i in ne])
    return IS_new - IS_old, (up_e - dn_e) + (dn_n - up_n)

idi, nri = idi_nri(pr_pv, pr_both, y)
print(f'\nIDI={idi:+.4f}  continuous NRI={nri:+.4f}')
rnd = random.Random(11); bs_i, bs_n = [], []
for _ in range(1000):
    idx = [rnd.randrange(n) for _ in range(n)]
    yy = [y[i] for i in idx]
    if sum(yy) in (0, n): continue
    po = [pr_pv[i] for i in idx]; pn = [pr_both[i] for i in idx]
    a1, a2 = idi_nri(po, pn, yy)
    bs_i.append(a1); bs_n.append(a2)
bs_i.sort(); bs_n.sort()
ci_i = (bs_i[int(0.025 * len(bs_i))], bs_i[int(0.975 * len(bs_i))])
ci_n = (bs_n[int(0.025 * len(bs_n))], bs_n[int(0.975 * len(bs_n))])
p_i = 2 * min(sum(1 for v in bs_i if v <= 0) / len(bs_i), sum(1 for v in bs_i if v >= 0) / len(bs_i))
p_n = 2 * min(sum(1 for v in bs_n if v <= 0) / len(bs_n), sum(1 for v in bs_n if v >= 0) / len(bs_n))
print(f'IDI 95%CI {ci_i[0]:+.4f}~{ci_i[1]:+.4f} P={p_i:.3f}; NRI 95%CI {ci_n[0]:+.4f}~{ci_n[1]:+.4f} P={p_n:.3f}')

# ---------------- pelvis volume operating characteristics ----------------
pv_e = [pv[i] for i in range(n) if y[i] == 1]
pv_n = [pv[i] for i in range(n) if y[i] == 0]
print('\n骨盆体积 事件组 %.1f±%.1f 中位 %.1f (IQR %.1f-%.1f)' % (mean(pv_e), sd(pv_e), median(pv_e), iqr(pv_e)[0], iqr(pv_e)[1]))
print('骨盆体积 非事件 %.1f±%.1f 中位 %.1f (IQR %.1f-%.1f)' % (mean(pv_n), sd(pv_n), median(pv_n), iqr(pv_n)[0], iqr(pv_n)[1]))

best = None
for cut in sorted(set(pv)):
    tp = sum(1 for i in range(n) if pv[i] <= cut and y[i] == 1)
    fp = sum(1 for i in range(n) if pv[i] <= cut and y[i] == 0)
    fn = len(pv_e) - tp; tn = len(pv_n) - fp
    sens = tp / (tp + fn); spec = tn / (tn + fp)
    j = sens + spec - 1
    if best is None or j > best['youden']:
        best = {'cut': cut, 'sens': sens, 'spec': spec, 'youden': j, 'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
                'ppv': tp / max(tp + fp, 1), 'npv': tn / max(tn + fn, 1), 'acc': (tp + tn) / n}
print('Youden 最佳截点 %.1f cm3: 敏感度 %.3f 特异度 %.3f PPV %.3f NPV %.3f 准确率 %.3f' % (
    best['cut'], best['sens'], best['spec'], best['ppv'], best['npv'], best['acc']))

qs = sorted(pv)
t1 = qs[int(n / 3)]; t2 = qs[int(2 * n / 3)]
def tert(v):
    return 0 if v <= t1 else (1 if v <= t2 else 2)
ter = [[0, 0] for _ in range(3)]
for i in range(n):
    k = tert(pv[i]); ter[k][0] += 1; ter[k][1] += y[i]
print('三分位: T1<=%.1f: %d/%d (%.1f%%), T2 %.1f-%.1f: %d/%d (%.1f%%), T3>%.1f: %d/%d (%.1f%%)' % (
    t1, ter[0][1], ter[0][0], 100 * ter[0][1] / ter[0][0], t1, t2, ter[1][1], ter[1][0], 100 * ter[1][1] / ter[1][0],
    t2, ter[2][1], ter[2][0], 100 * ter[2][1] / ter[2][0]))

# Cochran-Armitage trend test across tertiles
N = sum(t[0] for t in ter); R = sum(t[1] for t in ter); pbar = R / N
scores_t = [0, 1, 2]
num_t = sum(scores_t[k] * (ter[k][1] - ter[k][0] * pbar) for k in range(3))
var_t = pbar * (1 - pbar) * (sum(ter[k][0] * scores_t[k] ** 2 for k in range(3))
                             - (sum(ter[k][0] * scores_t[k] for k in range(3)) ** 2) / N)
z_t = num_t / math.sqrt(var_t)
p_trend = 2 * (1 - norm_cdf(abs(z_t)))
print('Cochran-Armitage 趋势 z=%.3f P=%.4g' % (z_t, p_trend))

# correlation with MeshVolume
r = 0.0
mp, mm = mean(pv), mean(mv)
num = sum((pv[i] - mp) * (mv[i] - mm) for i in range(n))
den = math.sqrt(sum((v - mp) ** 2 for v in pv) * sum((v - mm) ** 2 for v in mv))
r = num / den
print('骨盆体积 vs MeshVolume r=%.4f' % r)

# ---------------- DCA ----------------
ths = [t / 100 for t in range(1, 61)]
dca = {'thresholds': ths, 'curves': {'全部治疗': [prev - (1 - prev) * t / (1 - t) for t in ths], '不治疗': [0.0] * len(ths)}}
for lab, pr in (('骨盆体积', pr_pv), ('Rad-score', pr_rad), ('联合', pr_both)):
    nb = []
    for t in ths:
        tp = sum(1 for i in range(n) if pr[i] >= t and y[i] == 1)
        fp = sum(1 for i in range(n) if pr[i] >= t and y[i] == 0)
        nb.append(tp / n - fp / n * t / (1 - t))
    dca['curves'][lab] = nb
for t in (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40):
    k = ths.index(t)
    print('DCA t=%.2f 骨盆 %+.3f Rad %+.3f 联合 %+.3f 全治 %+.3f' % (
        t, dca['curves']['骨盆体积'][k], dca['curves']['Rad-score'][k], dca['curves']['联合'][k], dca['curves']['全部治疗'][k]))

rad_e = [rad_lin[i] for i in range(n) if y[i] == 1]
rad_n = [rad_lin[i] for i in range(n) if y[i] == 0]
print('\nRad-score(线性预测子) 事件 %.3f±%.3f 非事件 %.3f±%.3f MW P=%.4g' % (
    mean(rad_e), sd(rad_e), mean(rad_n), sd(rad_n), mann_whitney(rad_e, rad_n)[1]))

json.dump({'models': {k: {kk: vv for kk, vv in v.items() if kk != 'probs'} for k, v in res.items()},
           'per100': {'OR': math.exp(b_pv[1] * 100), 'lo': math.exp((b_pv[1] - 1.96 * se_pv[1]) * 100),
                      'hi': math.exp((b_pv[1] + 1.96 * se_pv[1]) * 100),
                      'p': 2 * (1 - norm_cdf(abs(b_pv[1] / se_pv[1])))},
           'delong': dl, 'idi': {'value': idi, 'ci': ci_i, 'p': p_i},
           'nri': {'value': nri, 'ci': ci_n, 'p': p_n},
           'cutoff': best, 'tertiles': {'t1': t1, 't2': t2, 'counts': ter, 'z': z_t, 'p_trend': p_trend},
           'corr_meshvolume': r, 'dca': dca,
           'desc': {'pv_event': [mean(pv_e), sd(pv_e), median(pv_e), iqr(pv_e)],
                    'pv_nonevent': [mean(pv_n), sd(pv_n), median(pv_n), iqr(pv_n)]}},
          open(OUT + 'phase6_mainline.json', 'w'), ensure_ascii=False, indent=1)
print('\n保存 phase6_mainline.json')
