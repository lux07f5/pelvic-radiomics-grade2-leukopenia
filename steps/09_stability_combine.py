# -*- coding: utf-8 -*-
"""Combine worker outputs + legacy progress into stability_results.json."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import json, sys, glob

OUT = config.OUT
p2 = json.load(open(OUT + 'phase2_results.json'))
sel_feats = p2['after_dedup'][:30]
P = len(sel_feats)
R = int(os.environ.get('RAD_STABILITY_REPEATS', '100'))

cnt_min = [0] * P
cnt_1se = [0] * P
seen = set()
worker_files = sorted(glob.glob(OUT + 'stability_w*.jsonl'))
# when the run was parallelised, the worker files are the whole truth; the
# single-process progress file would otherwise be counted twice
files = worker_files if worker_files else [OUT + 'stability_progress.jsonl']
nrec = 0
for f in files:
    try:
        for line in open(f):
            line = line.strip()
            if not line: continue
            rec = json.loads(line)
            r = rec.get('r', -1)
            if r >= 0:
                if r in seen: continue
                seen.add(r)
            for j in rec['min']: cnt_min[j] += 1
            for j in rec['1se']: cnt_1se[j] += 1
            nrec += 1
    except FileNotFoundError:
        pass

print(f'共读取 {nrec} 次重复（唯一 seed 数 {len(seen) if seen else "legacy"}）')
res = []
for j, f in enumerate(sel_feats):
    res.append({'class': f['class'], 'name': f['name'], 'matrix_row': f['matrix_row'],
                'freq_min': cnt_min[j] / R, 'freq_1se': cnt_1se[j] / R})
res.sort(key=lambda d: -d['freq_min'])
print('\n=== 入选频率 (λ.min, %d 次重复 10 折 CV) ===' % R)
for d in res:
    if d['freq_min'] > 0:
        print(f"  {d['freq_min']*100:5.0f}%  (1se {d['freq_1se']*100:3.0f}%)  {d['class']}_{d['name']} [row {d['matrix_row']+2}]")
json.dump({'repeats': R, 'completed': nrec, 'features': res},
          open(OUT + 'stability_results.json', 'w'), ensure_ascii=False, indent=1)
print('\n保存 stability_results.json')
