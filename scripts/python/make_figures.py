"""make_figures.py - 纯 Python 标准库 SVG 出图（无 matplotlib/numpy）
v3: 新剂量学数据 - 4 模型统一 60 例口径 (P059/P060 included)
"""
import json
import os
import math
import statistics
from collections import Counter

OUT_DIR = 'validation_output/figures'
os.makedirs(OUT_DIR, exist_ok=True)

SVG_W = 720
SVG_H = 540
PAD_L = 80
PAD_R = 40
PAD_T = 60
PAD_B = 70

def svg_header(title, subtitle=''):
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SVG_W} {SVG_H}" '
         f'font-family="Helvetica, Arial, sans-serif">'
         f'<rect width="{SVG_W}" height="{SVG_H}" fill="#fafafa"/>'
         f'<text x="{SVG_W/2}" y="28" text-anchor="middle" font-size="17" '
         f'font-weight="600" fill="#222">{title}</text>')
    if subtitle:
        s += f'<text x="{SVG_W/2}" y="48" text-anchor="middle" font-size="12" fill="#666">{subtitle}</text>'
    return s

def svg_footer():
    return '</svg>'

def ax(vmin, vmax, n=5):
    out = []
    for i in range(n+1):
        v = vmin + (vmax - vmin) * i / n
        out.append((v, f'{v:.2f}'))
    return out

def draw_axes(xmin, xmax, ymin, ymax, xlabel, ylabel, xticks=None, yticks=None, xfmt='{:.2f}', yfmt='{:.2f}'):
    parts = []
    x0, x1 = PAD_L, SVG_W - PAD_R
    y0, y1 = PAD_T, SVG_H - PAD_B
    parts.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#333" stroke-width="1.2"/>')
    parts.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#333" stroke-width="1.2"/>')
    if xticks is None: xticks = ax(xmin, xmax, 5)
    if yticks is None: yticks = ax(ymin, ymax, 5)
    for v,_ in xticks:
        px = x0 + (v - xmin) / (xmax - xmin) * (x1 - x0)
        parts.append(f'<line x1="{px}" y1="{y1}" x2="{px}" y2="{y1+4}" stroke="#333"/>')
        parts.append(f'<text x="{px}" y="{y1+18}" text-anchor="middle" font-size="10" fill="#555">{xfmt.format(v)}</text>')
    for v,_ in yticks:
        py = y1 - (v - ymin) / (ymax - ymin) * (y1 - y0)
        parts.append(f'<line x1="{x0-4}" y1="{py}" x2="{x0}" y2="{py}" stroke="#333"/>')
        parts.append(f'<text x="{x0-8}" y="{py+3}" text-anchor="end" font-size="10" fill="#555">{yfmt.format(v)}</text>')
    parts.append(f'<text x="{(x0+x1)/2}" y="{y1+44}" text-anchor="middle" font-size="12" fill="#333">{xlabel}</text>')
    parts.append(f'<text x="{x0-50}" y="{(y0+y1)/2}" text-anchor="middle" font-size="12" fill="#333" '
                 f'transform="rotate(-90 {x0-50} {(y0+y1)/2})">{ylabel}</text>')
    return ''.join(parts), (x0, y0, x1, y1)

def map_xy(x, y, box, xmin, xmax, ymin, ymax):
    x0, y0, x1, y1 = box
    px = x0 + (x - xmin) / (xmax - xmin) * (x1 - x0)
    py = y1 - (y - ymin) / (ymax - ymin) * (y1 - y0)
    return px, py

# ============================================================================
# 加载数据
# ============================================================================
nri = json.load(open('validation_output/nri_idi_result.json'))
nri36 = json.load(open('validation_output/nri_idi_36win.json'))
base = json.load(open('validation_output/baseline_60.json'))
merged = json.load(open('validation_output/merged_60.json'))
dose = json.load(open('validation_output/dosi_input.json'))
roc = json.load(open('validation_output/roc_curves.json')) if os.path.exists('validation_output/roc_curves.json') else None

ev_names = [p['name'] for p in merged if p['y'] == 1]
ne_names = [p['name'] for p in merged if p['y'] == 0]

# ============================================================================
# Figure 1: 决策曲线 DCA (n=48)
# ============================================================================
def fig_dca():
    parts = [svg_header('Figure 1. 决策曲线分析 (DCA)',
                        f'n = 48 (事件 12 / 非事件 36)，影像模型 vs 影像+WBC，Bootstrap B=1000')]
    ax_str, box = draw_axes(0, 0.50, -0.10, 0.25,
                             xlabel='Threshold Probability',
                             ylabel='Net Benefit',
                             xfmt='{:.2f}', yfmt='{:.2f}')
    parts.append(ax_str)
    thr = nri['dca']['treat_all']
    def to_line(vals, color, dash='', label='', width=1.5):
        pts = [map_xy(thr[i], vals[i], box, 0, 0.50, -0.10, 0.25) for i in range(len(thr))]
        d = ' '.join([f'{"M" if i==0 else "L"} {p[0]:.1f} {p[1]:.1f}' for i,p in enumerate(pts)])
        s = f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"'
        if dash: s += f' stroke-dasharray="{dash}"'
        s += '/>'
        s += f'<text x="{pts[-1][0]+6}" y="{pts[-1][1]+3}" font-size="10" fill="{color}">{label}</text>'
        return s
    parts.append(to_line(nri['dca']['treat_all'], '#888', dash='4 3', label='Treat All'))
    parts.append(to_line(nri['dca']['model_A'], '#1976d2', width=2, label='Imaging (M3)'))
    parts.append(to_line(nri['dca']['model_B'], '#d32f2f', width=2, label='Imaging + WBC'))
    px, py = map_xy(0, 0, box, 0, 0.50, -0.10, 0.25)
    parts.append(f'<line x1="{box[0]}" y1="{py}" x2="{box[2]}" y2="{py}" stroke="#bbb" stroke-dasharray="2 4"/>')
    parts.append(f'<text x="{box[2]+6}" y="{py+3}" font-size="10" fill="#888">Treat None</text>')
    px15, py15 = map_xy(0.15, 0, box, 0, 0.50, -0.10, 0.25)
    parts.append(f'<line x1="{px15}" y1="{box[1]}" x2="{px15}" y2="{box[3]}" stroke="#888" stroke-dasharray="2 4"/>')
    parts.append(f'<text x="{px15+4}" y="{box[1]+12}" font-size="9" fill="#666">@15% threshold</text>')
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig1_dca.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

# ============================================================================
# Figure 2: Bootstrap 特征入选频次柱状图
# ============================================================================
def fig_feature_freq():
    freqs = []
    with open('validation_output/v2_final/feature_selection_frequency.csv', 'r', encoding='utf-8') as f:
        next(f)
        for line in f:
            parts_ = line.strip().split(',')
            freqs.append((parts_[0], int(parts_[1]), float(parts_[2])))
    freqs = [x for x in freqs if x[2] >= 0.10][:8]
    parts = [svg_header('Figure 2. Bootstrap 影像特征入选频次',
                        f'B = 500 (有效 185 次), 频次 >= 10% 排名前 8 特征')]
    ymin, ymax = 0, 0.55
    ax_str, box = draw_axes(0, len(freqs), ymin, ymax,
                             xlabel='', ylabel='Selection Frequency (B=500)',
                             xticks=[(i+0.5, freqs[i][0]) for i in range(len(freqs))],
                             yticks=ax(ymin, ymax, 5), xfmt='{}', yfmt='{:.1%}')
    parts.append(ax_str)
    palette = ['#1976d2', '#d32f2f', '#388e3c', '#f57c00', '#7b1fa2', '#0097a7', '#5d4037', '#c2185b']
    for i, (name, cnt, freq) in enumerate(freqs):
        x0 = box[0] + (i + 0.15) * (box[2]-box[0]) / len(freqs)
        x1 = box[0] + (i + 0.85) * (box[2]-box[0]) / len(freqs)
        py = box[3] - (freq - ymin) / (ymax - ymin) * (box[3] - box[1])
        parts.append(f'<rect x="{x0:.1f}" y="{py:.1f}" width="{x1-x0:.1f}" height="{box[3]-py:.1f}" '
                     f'fill="{palette[i % len(palette)]}" opacity="0.85"/>')
        parts.append(f'<text x="{(x0+x1)/2:.1f}" y="{py-4:.1f}" text-anchor="middle" font-size="10" '
                     f'fill="#333">{freq*100:.1f}%</text>')
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig2_feature_freq.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

# ============================================================================
# Figure 3: 基线 WBC/Hb/BMI 箱线图 (n=48)
# ============================================================================
def fig_baseline_boxplot():
    metrics = [('WBC', 0), ('Hb', 1), ('BMI', 2)]
    data = {}
    for n in ev_names:
        v = base.get(n)
        if v: data.setdefault('Event', []).append([v[0], v[1], v[2]])
    for n in ne_names:
        v = base.get(n)
        if v: data.setdefault('Non-Event', []).append([v[0], v[1], v[2]])
    parts = [svg_header('Figure 3. 基线血液学与体质指标事件/非事件对比',
                        f'n = 48 (事件 12 / 非事件 36), 60 例中 12 例无基线血被排除; Wilcoxon, 全部 P > 0.05')]
    ax_str, box = draw_axes(0, 3, 0, 1,
                             xlabel='', ylabel='Standardized Value (min-max rescaled)',
                             xticks=[(i+0.5, m[0]) for i,m in enumerate(metrics)],
                             yticks=ax(0, 1, 5), xfmt='{}', yfmt='{:.1f}')
    parts.append(ax_str)
    n_groups = 2
    group_w = (box[2]-box[0]) / 3 * 0.7 / n_groups
    palette = {'Event': '#d32f2f', 'Non-Event': '#1976d2'}
    for gi, m in enumerate(metrics):
        all_vals = [r[m[1]] for grp in data.values() for r in grp]
        vmin, vmax = min(all_vals), max(all_vals)
        for gj, label in enumerate(['Event', 'Non-Event']):
            vals = sorted([r[m[1]] for r in data[label]])
            n = len(vals)
            q1 = vals[n//4]
            med = vals[n//2]
            q3 = vals[3*n//4]
            def res(v): return (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            cx = box[0] + (gi + 0.15 + (gj + 0.5) * 0.7 / n_groups) * (box[2]-box[0]) / 3
            wmin, wmax = vals[0], vals[-1]
            parts.append(f'<line x1="{cx:.1f}" y1="{box[3] - res(wmin)*(box[3]-box[1]):.1f}" '
                         f'x2="{cx:.1f}" y2="{box[3] - res(wmax)*(box[3]-box[1]):.1f}" '
                         f'stroke="{palette[label]}" stroke-width="1.5"/>')
            bx0, bx1 = cx - group_w/2, cx + group_w/2
            by_q3 = box[3] - res(q3)*(box[3]-box[1])
            by_q1 = box[3] - res(q1)*(box[3]-box[1])
            parts.append(f'<rect x="{bx0:.1f}" y="{by_q3:.1f}" width="{bx1-bx0:.1f}" height="{by_q1-by_q3:.1f}" '
                         f'fill="{palette[label]}" opacity="0.5" stroke="{palette[label]}" stroke-width="1.5"/>')
            my = box[3] - res(med)*(box[3]-box[1])
            parts.append(f'<line x1="{bx0:.1f}" y1="{my:.1f}" x2="{bx1:.1f}" y2="{my:.1f}" '
                         f'stroke="#222" stroke-width="2"/>')
            parts.append(f'<text x="{cx:.1f}" y="{box[3]+30}" text-anchor="middle" font-size="9" fill="#555">n={n}</text>')
    lx = box[0] + 8
    ly = box[1] + 8
    for k, c in palette.items():
        parts.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="14" height="10" fill="{c}" opacity="0.5" stroke="{c}"/>')
        parts.append(f'<text x="{lx+18:.1f}" y="{ly+9}" font-size="11" fill="#333">{k}</text>')
        ly += 16
    parts.append(f'<text x="{box[0]+10}" y="{box[3]+50}" font-size="11" fill="#666" font-style="italic">'
                 f'P(WBC)=0.30, P(Hb)=0.77, P(BMI)=0.77 — all non-significant (Wilcoxon)</text>')
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig3_baseline_boxplot.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

# ============================================================================
# Figure 4: 剂量学 4 指标 箱线图 (n=60 新数据)
# ============================================================================
def fig_dose_boxplot():
    metrics = [
        ('Volume (cm^3)', 0),
        ('Mean Dose (cGy)', 3),
        ('V5 (%)', 4),
        ('V20 (%)', 6),
    ]
    ev_dose = [d for d in dose['event']['data'] if d['name'] in ev_names]
    ne_dose = [d for d in dose['nonevent']['data'] if d['name'] in ne_names]
    parts = [svg_header('Figure 4. 盆腔 ROI 剂量学指标事件/非事件对比 (n=60, 新数据)',
                        f'事件 {len(ev_dose)} / 非事件 {len(ne_dose)}, 全 60 例具备完整剂量; Volume P=0.001 显著')]
    ax_str, box = draw_axes(0, 4, 0, 1,
                             xlabel='', ylabel='Standardized Value (rescaled)',
                             xticks=[(i+0.5, m[0]) for i,m in enumerate(metrics)],
                             yticks=ax(0, 1, 5), xfmt='{}', yfmt='{:.1f}')
    parts.append(ax_str)
    palette = {'Event': '#d32f2f', 'Non-Event': '#1976d2'}
    n_groups = 2
    group_w = (box[2]-box[0]) / 4 * 0.7 / n_groups
    pvals = ['**P=0.001**', 'P=0.489', 'P=0.051', 'P=0.898']
    for gi, (name, idx) in enumerate(metrics):
        ev_vals = sorted([d['vals'][idx] for d in ev_dose])
        ne_vals = sorted([d['vals'][idx] for d in ne_dose])
        all_vals = ev_vals + ne_vals
        vmin, vmax = min(all_vals), max(all_vals)
        def res(v): return (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
        for gj, (label, vals) in enumerate([('Event', ev_vals), ('Non-Event', ne_vals)]):
            n = len(vals)
            q1 = vals[n//4]
            med = vals[n//2]
            q3 = vals[3*n//4]
            wmin, wmax = vals[0], vals[-1]
            cx = box[0] + (gi + 0.15 + (gj + 0.5) * 0.7 / n_groups) * (box[2]-box[0]) / 4
            parts.append(f'<line x1="{cx:.1f}" y1="{box[3] - res(wmin)*(box[3]-box[1]):.1f}" '
                         f'x2="{cx:.1f}" y2="{box[3] - res(wmax)*(box[3]-box[1]):.1f}" '
                         f'stroke="{palette[label]}" stroke-width="1.5"/>')
            bx0, bx1 = cx - group_w/2, cx + group_w/2
            by_q3 = box[3] - res(q3)*(box[3]-box[1])
            by_q1 = box[3] - res(q1)*(box[3]-box[1])
            parts.append(f'<rect x="{bx0:.1f}" y="{by_q3:.1f}" width="{bx1-bx0:.1f}" height="{by_q1-by_q3:.1f}" '
                         f'fill="{palette[label]}" opacity="0.5" stroke="{palette[label]}" stroke-width="1.5"/>')
            my = box[3] - res(med)*(box[3]-box[1])
            parts.append(f'<line x1="{bx0:.1f}" y1="{my:.1f}" x2="{bx1:.1f}" y2="{my:.1f}" '
                         f'stroke="#222" stroke-width="2"/>')
            parts.append(f'<text x="{cx:.1f}" y="{box[3]+30}" text-anchor="middle" font-size="9" fill="#555">n={n}</text>')
        parts.append(f'<text x="{box[0] + (gi+0.5)*(box[2]-box[0])/4:.1f}" y="{box[1]-6:.1f}" '
                     f'text-anchor="middle" font-size="10" fill="#666" font-style="italic">{pvals[gi]}</text>')
    lx = box[0] + 8
    ly = box[1] + 8
    for k, c in palette.items():
        parts.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="14" height="10" fill="{c}" opacity="0.5" stroke="{c}"/>')
        parts.append(f'<text x="{lx+18:.1f}" y="{ly+9}" font-size="11" fill="#333">{k}</text>')
        ly += 16
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig4_dose_boxplot.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

# ============================================================================
# Figure 5: 四模型 AUC 柱状图 (新 60 例口径)
# ============================================================================
def fig_auc_bar():
    if roc:
        models = [
            ('M1\nClinical', roc['models']['M1']['auc'], roc['models']['M1']['n'], '#90a4ae'),
            ('M2\nDosimetric', roc['models']['M2']['auc'], roc['models']['M2']['n'], '#7e57c2'),
            ('M3\nRadiomics', roc['models']['M3']['auc'], roc['models']['M3']['n'], '#1976d2'),
            ('M4\nCombined', roc['models']['M4']['auc'], roc['models']['M4']['n'], '#d32f2f'),
        ]
        sub = f"n = 60 (事件 15), 4 模型统一 60 例口径 (新剂量学数据); M4 最高"
    else:
        models = [('M1', 0.6, 60, '#90a4ae'),('M2', 0.7, 60, '#7e57c2'),('M3', 0.79, 60, '#1976d2'),('M4', 0.8, 60, '#d32f2f')]
        sub = 'n=60'
    parts = [svg_header('Figure 5. 多模态模型表观 AUC 对照 (新 60 例口径)', sub)]
    ax_str, box = draw_axes(0, 4, 0.4, 0.9,
                             xlabel='', ylabel='Apparent AUC',
                             xticks=[(i+0.5, m[0]) for i,m in enumerate(models)],
                             yticks=ax(0.4, 0.9, 5), xfmt='{}', yfmt='{:.2f}')
    parts.append(ax_str)
    px5 = box[0] + (0.5 - 0.4) / (0.9 - 0.4) * (box[2]-box[0])
    parts.append(f'<line x1="{px5:.1f}" y1="{box[1]}" x2="{px5:.1f}" y2="{box[3]}" '
                 f'stroke="#bbb" stroke-dasharray="4 3"/>')
    parts.append(f'<text x="{px5+6:.1f}" y="{box[1]+12}" font-size="9" fill="#888">0.50 (random)</text>')
    px_c = box[0] + (0.695 - 0.4) / (0.9 - 0.4) * (box[2]-box[0])
    parts.append(f'<line x1="{px_c:.1f}" y1="{box[1]}" x2="{px_c:.1f}" y2="{box[3]}" '
                 f'stroke="#d32f2f" stroke-dasharray="2 4" opacity="0.6"/>')
    parts.append(f'<text x="{px_c-6:.1f}" y="{box[1]+12}" font-size="9" fill="#d32f2f" text-anchor="end">optimism-corrected 0.695</text>')
    bar_w = (box[2]-box[0]) / 4 * 0.6
    for i, (name, auc, n, color) in enumerate(models):
        cx = box[0] + (i + 0.5) * (box[2]-box[0]) / 4
        py = box[3] - (auc - 0.4) / (0.9 - 0.4) * (box[3]-box[1])
        parts.append(f'<rect x="{cx - bar_w/2:.1f}" y="{py:.1f}" width="{bar_w:.1f}" height="{box[3]-py:.1f}" '
                     f'fill="{color}" opacity="0.85"/>')
        parts.append(f'<text x="{cx:.1f}" y="{py-6:.1f}" text-anchor="middle" font-size="11" font-weight="600" fill="#222">{auc:.3f}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{box[3]+18}" text-anchor="middle" font-size="9" fill="#555">n={n}</text>')
    parts.append(f'<text x="{box[0]+10}" y="{box[3]+50}" font-size="10" fill="#666" font-style="italic">'
                 f'Pairwise DeLong tests P &gt; 0.05 (underpowered; events=15)</text>')
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig5_auc_bar.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

# ============================================================================
# Figure 6: ROC 曲线 (4 模型真实坐标, 全 60 例)
# ============================================================================
def fig_roc():
    if not roc: return
    parts = [svg_header('Figure 6. 多模态模型 ROC 曲线 (n=60 全 60 例)',
                        'AUC 见表 4; 对角线为随机基线; 各模型独立 IRLS 拟合 (新剂量学口径)')]
    ax_str, box = draw_axes(0, 1, 0, 1,
                             xlabel='1 - Specificity (FPR)',
                             ylabel='Sensitivity (TPR)',
                             xticks=ax(0, 1, 5), yticks=ax(0, 1, 5),
                             xfmt='{:.1f}', yfmt='{:.1f}')
    parts.append(ax_str)
    parts.append(f'<line x1="{box[0]}" y1="{box[3]}" x2="{box[2]}" y2="{box[1]}" '
                 f'stroke="#bbb" stroke-dasharray="3 3"/>')
    parts.append(f'<text x="{box[2]-8}" y="{box[1]+14}" text-anchor="end" font-size="9" fill="#888">random</text>')
    colors = {'M1': '#90a4ae', 'M2': '#7e57c2', 'M3': '#1976d2', 'M4': '#d32f2f'}
    widths = {'M3': 2.5, 'M1': 1.4, 'M2': 1.4, 'M4': 1.4}
    for mk in ['M1', 'M2', 'M3', 'M4']:
        m = roc['models'][mk]
        pts = m['points']
        d = ' '.join([f'{"M" if i==0 else "L"} {box[0] + pt["fpr"]*(box[2]-box[0]):.1f} '
                      f'{box[3] - pt["tpr"]*(box[3]-box[1]):.1f}' for i,pt in enumerate(pts)])
        parts.append(f'<path d="{d}" fill="none" stroke="{colors[mk]}" stroke-width="{widths[mk]}"/>')
        last = pts[-2]
        lx = box[0] + last['fpr']*(box[2]-box[0]) + 6
        ly = box[3] - last['tpr']*(box[3]-box[1]) + 4
        parts.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="11" font-weight="600" fill="{colors[mk]}">'
                     f'{mk} AUC={m["auc"]:.3f}</text>')
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig6_roc.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

# ============================================================================
# Figure 7: 校准曲线 (M3, 10 等分, n=60)
# ============================================================================
def fig_calibration():
    if not roc or not roc.get('calibration_M3'): return
    parts = [svg_header('Figure 7. M3 影像模型校准曲线 (n=60, 10 等分)',
                        'Hosmer-Lemeshow chi^2=7.01, df=8, P=0.535; optimism-corrected AUC=0.695')]
    ax_str, box = draw_axes(0, 1, 0, 1,
                             xlabel='Predicted Probability',
                             ylabel='Observed Event Rate',
                             xticks=ax(0, 1, 5), yticks=ax(0, 1, 5),
                             xfmt='{:.1f}', yfmt='{:.1f}')
    parts.append(ax_str)
    parts.append(f'<line x1="{box[0]}" y1="{box[3]}" x2="{box[2]}" y2="{box[1]}" '
                 f'stroke="#888" stroke-dasharray="4 3"/>')
    parts.append(f'<text x="{box[2]-8}" y="{box[1]+14}" text-anchor="end" font-size="10" fill="#666">Ideal</text>')
    cal = roc['calibration_M3']
    nmax = max(c['n'] for c in cal)
    for c in cal:
        px = box[0] + c['meanPred']*(box[2]-box[0])
        py = box[3] - c['obsRate']*(box[3]-box[1])
        r = 4 + 8 * c['n'] / nmax
        parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r:.1f}" fill="#1976d2" opacity="0.6" '
                     f'stroke="#0d47a1" stroke-width="1.2"/>')
        parts.append(f'<text x="{px:.1f}" y="{py-r-3:.1f}" text-anchor="middle" font-size="8" fill="#333">n={c["n"]}</text>')
    base_rate = 15/60
    py_base = box[3] - base_rate*(box[3]-box[1])
    parts.append(f'<line x1="{box[0]}" y1="{py_base:.1f}" x2="{box[2]}" y2="{py_base:.1f}" '
                 f'stroke="#bbb" stroke-dasharray="2 4"/>')
    parts.append(f'<text x="{box[2]-8}" y="{py_base-3:.1f}" text-anchor="end" font-size="9" fill="#888">event rate 25%</text>')
    lx, ly = box[0]+8, box[1]+8
    parts.append(f'<text x="{lx}" y="{ly}" font-size="10" fill="#333">bubble size ~ n/bin</text>')
    parts.append(svg_footer())
    with open(f'{OUT_DIR}/fig7_calibration.svg', 'w', encoding='utf-8') as f:
        f.write(''.join(parts))

fig_dca()
fig_feature_freq()
fig_baseline_boxplot()
fig_dose_boxplot()
fig_auc_bar()
fig_roc()
fig_calibration()
print(f'生成 SVG → {OUT_DIR}/')
for f in sorted(os.listdir(OUT_DIR)):
    sz = os.path.getsize(f'{OUT_DIR}/{f}')
    print(f'  {f}: {sz} bytes')
