# -*- coding: utf-8 -*-
"""Phase 4: calibration (intercept/slope, HL test, curve), DCA, Firth refit."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, csv, math, sys
from lib.stats_lib import *

OUT = config.OUT
p2 = json.load(open(OUT + 'phase2_results.json'))
D = json.load(open(OUT + 'radiomics_matrix.json'))
rows = list(csv.DictReader(open(OUT + 'merged_dataset.csv', encoding='utf-8-sig')))
y = p2['y']; n = len(y); names = D['names']
rowmap = {r['patient_id']: r for r in rows}
prev = sum(y) / n

# rad linear predictor (raw scale, from saved coefficients)
rad_lin = [p2['rad_intercept'] + sum(rf['coef'] * D['X'][i][rf['matrix_row']] for rf in p2['rad_features']) for i in range(n)]
v5 = [float(rowmap[nm]['V5']) for nm in names]
pv = [float(rowmap[nm]['pelvis_volume']) for nm in names]
v5z = [(v - mean(v5)) / sd(v5) for v in v5]

# combined model (standard logistic, for probs)
betaC, seC, _ = logistic_fit([[rad_lin[i], v5z[i]] for i in range(n)], y)
comb_lin = [betaC[0] + betaC[1] * rad_lin[i] + betaC[2] * v5z[i] for i in range(n)]

# pelvis volume model
betaB, seB, _ = logistic_fit([[v] for v in pv], y)
pv_lin = [betaB[0] + betaB[1] * v for v in pv]

def probs(lin): return [1.0 / (1.0 + math.exp(-max(-30, min(30, e)))) for e in lin]

models = {
    '骨盆体积': pv_lin,
    'Rad-score': rad_lin,
    '联合模型': comb_lin,
}

# ---------- calibration ----------
def calibration(lin, yv, g=5):
    pr = probs(lin)
    # calibration slope: y ~ LP
    b, se, cv = logistic_fit([[l] for l in lin], yv)
    slope = b[1]
    # calibration intercept: y ~ 1 with offset LP (single-param Newton)
    a = 0.0
    for _ in range(100):
        mu = [1.0 / (1.0 + math.exp(-max(-30, min(30, a + l)))) for l in lin]
        gr = sum(yv[i] - mu[i] for i in range(n))
        he = sum(max(m * (1 - m), 1e-9) for m in mu)
        da = gr / he
        a += da
        if abs(da) < 1e-8: break
    # HL test with g groups
    order = sorted(range(len(pr)), key=lambda i: pr[i])
    groups = [order[int(i * len(pr) / g):int((i + 1) * len(pr) / g)] for i in range(g)]
    HL = 0.0
    bins = []
    for gi, grp in enumerate(groups):
        O = sum(yv[i] for i in grp)
        E = sum(pr[i] for i in grp)
        ng = len(grp)
        if E > 0.5 and ng - E > 0.5:
            HL += (O - E) ** 2 / (E * (1 - E / ng))
        bins.append({'lo': pr[grp[0]], 'hi': pr[grp[-1]], 'pred': E / ng, 'obs': O / ng, 'n': ng})
    p_hl = chi2_sf(HL, g - 2)
    return {'slope': slope, 'intercept': a, 'HL': HL, 'HL_p': p_hl, 'groups': g, 'bins': bins}

def oof_lin(X, yv, nfolds=10, seed=42):
    """10-fold out-of-fold linear predictors for the FINAL logistic step (feature set fixed)."""
    rnd = random.Random(seed)
    pos = [i for i in range(n) if yv[i] == 1]
    neg = [i for i in range(n) if yv[i] == 0]
    rnd.shuffle(pos); rnd.shuffle(neg)
    folds = [[] for _ in range(nfolds)]
    for k, i in enumerate(pos): folds[k % nfolds].append(i)
    for k, i in enumerate(neg): folds[k % nfolds].append(i)
    lin = [None] * n
    for f in range(nfolds):
        test = folds[f]
        Xtr = [X[i] for i in range(n) if i not in test]
        ytr = [yv[i] for i in range(n) if i not in test]
        b, _, _ = logistic_fit(Xtr, ytr)
        for i in test:
            lin[i] = b[0] + sum(b[j + 1] * X[i][j] for j in range(len(X[0])))
    return lin

# out-of-fold calibration (honest, non-degenerate)
oof_inputs = {'骨盆体积': [[pv[i]] for i in range(n)],
              'Rad-score': [[rad_lin[i]] for i in range(n)],
              '联合模型': [[rad_lin[i], v5z[i]] for i in range(n)]}
calib_oof = {}
for lab, X_ in oof_inputs.items():
    lin_oof = oof_lin(X_, y)
    calib_oof[lab] = calibration(lin_oof, y)
    print(f"OOF 校准 [{lab}]: slope={calib_oof[lab]['slope']:.3f}, intercept={calib_oof[lab]['intercept']:.3f}, HL={calib_oof[lab]['HL']:.3f}, P={calib_oof[lab]['HL_p']:.3f}")

calib = {}
for lab, lin in models.items():
    calib[lab] = calibration(lin, y)
    print(f"校准 [{lab}]: slope={calib[lab]['slope']:.3f}, intercept={calib[lab]['intercept']:.3f}, HL(4df)={calib[lab]['HL']:.3f}, P={calib[lab]['HL_p']:.3f}")

# ---------- DCA ----------
ths = [t / 100 for t in range(1, 61)]
dca = {'thresholds': ths, 'curves': {}}
dca['curves']['全部治疗'] = [prev - (1 - prev) * t / (1 - t) for t in ths]
dca['curves']['不治疗'] = [0.0] * len(ths)
for lab, lin in models.items():
    pr = probs(lin)
    nb = []
    for t in ths:
        tp = sum(1 for i in range(n) if pr[i] >= t and y[i] == 1)
        fp = sum(1 for i in range(n) if pr[i] >= t and y[i] == 0)
        nb.append(tp / n - fp / n * t / (1 - t))
    dca['curves'][lab] = nb
print('DCA 完成')

# ---------- Firth refit ----------
fb, fse, fcv = firth_logistic_fit([[rad_lin[i], v5z[i]] for i in range(n)], y)
firth_combined = {'beta_rad': fb[1], 'se_rad': fse[1], 'OR_rad': math.exp(fb[1]),
                  'beta_v5': fb[2], 'se_v5': fse[2], 'OR_v5': math.exp(fb[2]),
                  'intercept': fb[0], 'converged': fcv}
print(f"Firth 联合模型: rad_lin OR={firth_combined['OR_rad']:.3f} (beta={fb[1]:.4f}±{fse[1]:.4f}), V5z OR={firth_combined['OR_v5']:.3f} (beta={fb[2]:.4f}±{fse[2]:.4f}), conv={fcv}")
fbB, fseB, _ = firth_logistic_fit([[v] for v in pv], y)
print(f"Firth 骨盆体积: OR/cm³={math.exp(fbB[1]):.5f} (beta={fbB[1]:.6f}±{fseB[1]:.6f})")

json.dump({'calibration': calib, 'calibration_oof': calib_oof, 'dca': dca,
           'firth_combined': firth_combined,
           'firth_pelvis': {'beta': fbB[1], 'se': fseB[1], 'OR': math.exp(fbB[1])}},
          open(OUT + 'phase4_results.json', 'w'), ensure_ascii=False, indent=1)
print('保存 phase4_results.json')
