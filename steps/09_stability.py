# -*- coding: utf-8 -*-
"""Repeated 10-fold CV LASSO stability: 100 repeats, selection frequency per feature."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, sys, time, os
from lib.stats_lib import *

OUT = config.OUT
p2 = json.load(open(OUT + 'phase2_results.json'))
D = json.load(open(OUT + 'radiomics_matrix.json'))
y = p2['y']; n = len(y)

# rebuild the same top-30 filtered set as phase2
sel_feats = p2['after_dedup'][:30]
Xs = [[D['X'][i][f['matrix_row']] for f in sel_feats] for i in range(n)]
mu = [mean([r[j] for r in Xs]) for j in range(len(sel_feats))]
sdd = [sd([r[j] for r in Xs]) or 1.0 for j in range(len(sel_feats))]
Xs_std = [[(r[j] - mu[j]) / sdd[j] for j in range(len(sel_feats))] for r in Xs]

R = int(os.environ.get('RAD_STABILITY_REPEATS', '100'))
cnt_min = [0] * len(sel_feats)
cnt_1se = [0] * len(sel_feats)
PROG = OUT + 'stability_progress.jsonl'
# resume: load completed repeats
done = 0
if os.path.exists(PROG):
    for line in open(PROG):
        line = line.strip()
        if not line: continue
        rec = json.loads(line)
        for j in rec['min']: cnt_min[j] += 1
        for j in rec['1se']: cnt_1se[j] += 1
        done += 1
    print(f'续跑：已完成 {done}/{R} 次', flush=True)
t0 = time.time()
for r in range(done, R):
    s_min, s_1se = lasso_select(Xs_std, y, nfolds=10, nlam=50, seed=1000 + r)
    for j in s_min: cnt_min[j] += 1
    for j in s_1se: cnt_1se[j] += 1
    with open(PROG, 'a') as fh:
        fh.write(json.dumps({'min': sorted(s_min), '1se': sorted(s_1se)}) + '\n')
    if (r + 1) % 5 == 0:
        el = time.time() - t0
        print(f'[{r+1}/{R}] 本轮用时 {el:.0f}s, 预计剩余 {el/(r+1-done)*(R-r-1):.0f}s', flush=True)

res = []
for j, f in enumerate(sel_feats):
    res.append({'class': f['class'], 'name': f['name'], 'matrix_row': f['matrix_row'],
                'freq_min': cnt_min[j] / R, 'freq_1se': cnt_1se[j] / R})
res.sort(key=lambda d: -d['freq_min'])
print('\n=== 入选频率 (λ.min, 100 次重复 CV) ===')
for d in res:
    if d['freq_min'] > 0:
        print(f"  {d['freq_min']*100:5.0f}%  (1se {d['freq_1se']*100:3.0f}%)  {d['class']}_{d['name']}")
json.dump({'repeats': R, 'features': res}, open(OUT + 'stability_results.json', 'w'), ensure_ascii=False, indent=1)
print('保存 stability_results.json')
