"""
nomogram.py - 纯标准库 SVG nomogram
M3 影像组学模型 (DepVar + MeshVolume, β0=-1.555, βDV=-0.468, βMV=-1.167)
"""
import json, statistics, math

# === 真实数据 (60 例完整队列) ===
v2 = json.load(open('validation_output/input_data_v2.json'))
fn = v2['feature_names']
iDV = fn.index('DependenceVariance')
iMV = fn.index('MeshVolume')
X = v2['X']
y = v2['y']
DVs = [X[i][iDV] for i in range(60)]
MVs = [X[i][iMV] for i in range(60)]

# z-score 归一化
DV_mu, DV_sd = statistics.mean(DVs), statistics.stdev(DVs)
MV_mu, MV_sd = statistics.mean(MVs), statistics.stdev(MVs)

# v2_final 标准 β 系数
b0, bDV, bMV = -1.555, -0.468, -1.167

# LP & Prob (60 例)
LPs = [b0 + bDV*(DVs[i]-DV_mu)/DV_sd + bMV*(MVs[i]-MV_mu)/MV_sd for i in range(60)]
Probs = [1/(1+math.exp(-lp)) for lp in LPs]

# Nomogram 范围
DV_min, DV_max = min(DVs), max(DVs)
MV_min, MV_max = min(MVs), max(MVs)
LP_min_data, LP_max_data = min(LPs), max(LPs)

# === SVG 绘制 ===
W, H = 720, 700
parts = []
parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#333">')
parts.append(f'<rect width="{W}" height="{H}" fill="#fafafa"/>')
parts.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-size="18" font-weight="700" fill="#1a237e">盆腔放疗 Grade≥2 白细胞减少影像组学诺莫图 (Nomogram)</text>')
parts.append(f'<text x="{W/2}" y="50" text-anchor="middle" font-size="11" fill="#555">M3 模型: DependenceVariance (GLCM) + MeshVolume (Shape) | n=60, 事件=15 | AUC=0.793, 校正 AUC=0.695</text>')

# 标尺绘制函数
def axis(x1, y1, x2, y2, label, ticks, label_pos='top', tick_fmt=None, color='#444'):
    s = []
    s.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="2"/>')
    for tv, tl in ticks:
        # tv 是归一化位置 (0-1), tl 是标签
        tx = x1 + tv * (x2 - x1)
        s.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" stroke-width="1.5"/>'.format(tx, y2-5, tx, y2+5, color))
        txt = tick_fmt(tl) if tick_fmt else str(tl)
        anchor = 'middle'
        label_y = y2 + 18 if label_pos == 'top' else y1 - 8
        s.append('<text x="{:.1f}" y="{:.1f}" text-anchor="{}" font-size="10" fill="#333">{}</text>'.format(tx, label_y, anchor, txt))
    # 轴标签
    lbl_y = y1 - 12 if label_pos == 'top' else y2 + 30
    s.append(f'<text x="{(x1+x2)/2}" y="{lbl_y}" text-anchor="middle" font-size="12" font-weight="600" fill="#222">{label}</text>')
    return s

# Nomogram 设计：5 行从上到下
# Row 1 (y=100): Points (0-100) - 顶轴
# Row 2 (y=180): DepVar (GLCM)
# Row 3 (y=260): MeshVolume (shape)
# Row 4 (y=340): Total Points (0-200) - 中轴
# Row 5 (y=440): Predicted Probability (0-1)

# x 范围
ax_x1, ax_x2 = 180, 660
n_ticks = 11

# 1. Points 顶轴
parts += axis(ax_x1, 90, ax_x2, 90, 'Points', [(i/10, i*10) for i in range(11)], tick_fmt=lambda v: f'{int(v)}')
parts.append(f'<text x="{ax_x1-10}" y="94" text-anchor="end" font-size="10" fill="#888">0</text>')
parts.append(f'<text x="{ax_x2+10}" y="94" text-anchor="start" font-size="10" fill="#888">100</text>')

# 2. DepVar 变量轴
# βDV < 0, 高值=低分
# x=16.29 → 0 pts (右边), x=6.10 → 100 pts (左边)
# 在 axis 上, x1=左 0, x2=右 100
# 但因为 β<0, 我们把 x_max 放在 0 分端 (左侧), x_min 放在 100 分端 (右侧)
# 实际数据: x_max=16.29 (高分=低风险) → 0 分, x_min=6.10 (低分=高风险) → 100 分
# axis 上左=0, 右=100, 所以 x_max(16.29) → 左 0 pts, x_min(6.10) → 右 100 pts
# ticks: 6.10, 7.0, 8.0, ..., 16.29 → 从右到左
def fmt_dv(v):
    if v == int(v): return f'{int(v)}'
    return f'{v:.1f}'

# 反向映射: x_value → points (0-100)
# x=16.29 → 0, x=6.10 → 100
# points(x) = 100 * (16.29 - x) / (16.29 - 6.10)
# x_value at pts t: 16.29 - t/100 * 10.19
DV_ticks_vals = [6.10, 7.5, 9.0, 10.5, 12.0, 13.5, 15.0, 16.29]
DV_ticks_pts = [100 * (16.29 - v) / 10.19 for v in DV_ticks_vals]
DV_ticks = [(p/100, fmt_dv(v)) for v, p in zip(DV_ticks_vals, DV_ticks_pts)]
parts += axis(ax_x1, 170, ax_x2, 170, 'DependenceVariance (GLCM)', DV_ticks, tick_fmt=lambda v: v, color='#1976d2')
# 标注方向
parts.append(f'<text x="{ax_x1-10}" y="174" text-anchor="end" font-size="9" fill="#1976d2">↓高风险</text>')
parts.append(f'<text x="{ax_x2+10}" y="174" text-anchor="start" font-size="9" fill="#1976d2">低风险↑</text>')

# 3. MeshVolume 变量轴
# βMV < 0, x_max=1.838M → 0 pts, x_min=1.009M → 100 pts
MV_ticks_vals = [1.01, 1.15, 1.30, 1.45, 1.60, 1.75, 1.84]  # 单位: 百万
MV_ticks_vals_full = [v*1e6 for v in MV_ticks_vals]
MV_ticks_pts = [100 * (1.838e6 - v) / (1.838e6 - 1.009e6) for v in MV_ticks_vals_full]
MV_ticks = [(p/100, f'{v:.2f}') for v, p in zip(MV_ticks_vals, MV_ticks_pts)]
parts += axis(ax_x1, 250, ax_x2, 250, 'MeshVolume (×10⁶)', MV_ticks, tick_fmt=lambda v: v, color='#d32f2f')
parts.append(f'<text x="{ax_x1-10}" y="254" text-anchor="end" font-size="9" fill="#d32f2f">↓高风险</text>')
parts.append(f'<text x="{ax_x2+10}" y="254" text-anchor="start" font-size="9" fill="#d32f2f">低风险↑</text>')

# 4. Total Points 中轴
parts += axis(ax_x1, 330, ax_x2, 330, 'Total Points', [(i/10, i*20) for i in range(11)], tick_fmt=lambda v: f'{int(v)}')
parts.append(f'<text x="{ax_x1-10}" y="334" text-anchor="end" font-size="10" fill="#888">0</text>')
parts.append(f'<text x="{ax_x2+10}" y="334" text-anchor="start" font-size="10" fill="#888">200</text>')

# 5. Predicted Probability 底轴
# TP → LP → Prob: LP = LP_min + TP/200 * (LP_max - LP_min)
# 用关键概率点: 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7
def inv_sigmoid(p):
    return math.log(p/(1-p))
prob_ticks = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
prob_lps = [inv_sigmoid(p) for p in prob_ticks]
# LP → TP
prob_tps = [(lp - LP_min_data) / (LP_max_data - LP_min_data) * 200 for lp in prob_lps]
prob_ticks_pos = [(max(0, min(1, tp/200)), f'{p*100:.0f}%' if p<0.1 else f'{p*100:.0f}') for tp, p in zip(prob_tps, prob_ticks)]
# filter 0-1
prob_ticks_pos = [(t, l) for t, l in prob_ticks_pos if 0 <= t <= 1]
parts += axis(ax_x1, 410, ax_x2, 410, 'Predicted Probability', prob_ticks_pos, tick_fmt=lambda v: v, color='#388e3c')

# 标尺间连接线
parts.append(f'<line x1="{ax_x1-30}" y1="90" x2="{ax_x1-30}" y2="410" stroke="#bbb" stroke-dasharray="2 3"/>')
parts.append(f'<line x1="{ax_x2+30}" y1="90" x2="{ax_x2+30}" y2="410" stroke="#bbb" stroke-dasharray="2 3"/>')

# 临床示例 (一个事件组典型病人, 已去标识化)
# 找一个事件组典型: LP 接近 0
# 找 60 例里 LP 接近 0 的事件
typical_idx = sorted(range(60), key=lambda i: abs(LPs[i]))[5]  # 接近中位的
typical_name = v2['patients'][typical_idx]
typical_LP = LPs[typical_idx]
typical_P = Probs[typical_idx]
typical_y = y[typical_idx]
y_label = '事件' if typical_y == 1 else '非事件'
# 算 points
tp = (typical_LP - LP_min_data) / (LP_max_data - LP_min_data) * 200
dv_pts = 100 * (16.29 - DVs[typical_idx]) / 10.19
mv_pts = 100 * (1.838e6 - MVs[typical_idx]) / (1.838e6 - 1.009e6)
parts.append(f'<text x="{W/2}" y="455" text-anchor="middle" font-size="12" font-weight="600" fill="#1a237e">使用示例 (Typical case)</text>')
parts.append(f'<text x="{W/2}" y="478" text-anchor="middle" font-size="11" fill="#444">'
             f'假设患者: {typical_name} ({y_label}) | DV={DVs[typical_idx]:.2f} → {dv_pts:.0f} pts | '
             f'MV={MVs[typical_idx]/1e6:.3f}×10⁶ → {mv_pts:.0f} pts | '
             f'Total={dv_pts+mv_pts:.0f} pts → P(Grade≥2 WBC↓) ≈ {typical_P*100:.1f}%</text>')

# β 系数标注
parts.append(f'<text x="{W/2}" y="510" text-anchor="middle" font-size="11" fill="#666">'
             f'β 系数 (z-score 标准化后): β0 = {b0:.3f}, βDV = {bDV:.3f}, βMV = {bMV:.3f}</text>')
parts.append(f'<text x="{W/2}" y="528" text-anchor="middle" font-size="10" fill="#888">'
             f'标准化参数: DV μ={DV_mu:.3f}, σ={DV_sd:.3f} | MV μ={MV_mu:.0f}, σ={MV_sd:.0f}</text>')

# 局限性
parts.append(f'<rect x="40" y="550" width="{W-80}" height="120" fill="#fff3e0" stroke="#ffb74d" stroke-width="1" rx="4"/>')
parts.append(f'<text x="55" y="572" font-size="12" font-weight="700" fill="#e65100">⚠ 局限性 (Limitations)</text>')
limits = [
    '• 样本量小 (n=60, 事件 15), LASSO 约束使模型骨架仅 2 特征, 但 EPV=7.5 仍偏紧',
    '• 内部 Bootstrap 乐观校正 (B=500): 表观 AUC 0.793 → 校正 AUC 0.695 (optimism=0.098)',
    '• 缺乏外部验证: UCLA 同题研究 (152 例, AUC 0.71-0.77) 可作锚点但人群/方案不同',
    '• 影像组学特征依赖 PyRadiomics 提取参数 (bin width/归一化), 跨中心可复现性待证',
    '• 仅做"是否发生"二分类, 未区分严重程度 (Grade 2/3) 与时间-事件终点 (kinetics)',
    '• Nomogram 的概率刻度基于 60 例样本 LP 实际范围, 外推时需谨慎'
]
for i, l in enumerate(limits):
    parts.append(f'<text x="55" y="{595 + i*14}" font-size="10" fill="#333">{l}</text>')

parts.append('</svg>')

with open('validation_output/figures/fig8_nomogram.svg', 'w', encoding='utf-8') as f:
    f.write('\n'.join(parts))
print('Nomogram SVG: validation_output/figures/fig8_nomogram.svg')

# === 同时输出 nomogram 参数 JSON, 给 HTML 交互式用 ===
nomo_params = {
    'b0': b0, 'bDV': bDV, 'bMV': bMV,
    'DV_mu': DV_mu, 'DV_sd': DV_sd, 'MV_mu': MV_mu, 'MV_sd': MV_sd,
    'DV_min': DV_min, 'DV_max': DV_max,
    'MV_min': MV_min, 'MV_max': MV_max,
    'LP_min': LP_min_data, 'LP_max': LP_max_data,
    'apparent_auc': 0.793, 'corrected_auc': 0.695,
    'N': 60, 'events': 15
}
with open('validation_output/nomogram_params.json', 'w', encoding='utf-8') as f:
    json.dump(nomo_params, f, ensure_ascii=False, indent=2)
print('参数 JSON: validation_output/nomogram_params.json')
