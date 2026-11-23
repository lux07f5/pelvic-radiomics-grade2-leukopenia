# -*- coding: utf-8 -*-
"""Parallel worker for repeated-CV LASSO stability. Usage: worker.py START END OUTFILE"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, sys, time, os
from lib.stats_lib import *

OUT = config.OUT
r0, r1, outf = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]

p2 = json.load(open(OUT + 'phase2_results.json'))
D = json.load(open(OUT + 'radiomics_matrix.json'))
y = p2['y']; n = len(y)
sel_feats = p2['after_dedup'][:30]
Xs = [[D['X'][i][f['matrix_row']] for f in sel_feats] for i in range(n)]
mu = [mean([r[j] for r in Xs]) for j in range(len(sel_feats))]
sdd = [sd([r[j] for r in Xs]) or 1.0 for j in range(len(sel_feats))]
Xs_std = [[(r[j] - mu[j]) / sdd[j] for j in range(len(sel_feats))] for r in Xs]

done = set()
if os.path.exists(outf):
    for line in open(outf):
        line = line.strip()
        if line: done.add(json.loads(line)['r'])
    print(f'[{outf}] 续跑，已完成 {len(done)} 次', flush=True)

t0 = time.time()
for r in range(r0, r1):
    if r in done: continue
    s_min, s_1se = lasso_select(Xs_std, y, nfolds=10, nlam=50, seed=1000 + r)
    with open(outf, 'a') as fh:
        fh.write(json.dumps({'r': r, 'min': sorted(s_min), '1se': sorted(s_1se)}) + '\n')
print(f'[{outf}] 完成 {r1 - r0} 次, 用时 {time.time() - t0:.0f}s', flush=True)
