"""Supplementary Figures S1 and S2, built at final physical size.

Same point-based canvas as the main figures: the SVG declares a millimetre
width and the viewBox is in points, so every font-size is the printed size.
S1 is a three-panel page and takes the double column (190 mm); S2 is a single
panel with a side legend and takes the 1.5 column width (140 mm).

Typography and stroke weights follow the Elsevier floor: 7 pt normal text,
nothing below 6 pt, lines between 0.25 and 1.5 pt.

Inputs (read-only): outputs/phase6_mainline.json
Outputs: Supplementary_Figure_S1_calibration.svg
         Supplementary_Figure_S2_dca.svg
         Supplementary_Figure_preview.html
"""

import json
import math
import os

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
BASE = config.OUT_DIR
SRC = os.path.join(BASE, 'phase6_mainline.json')
OUT1 = config.fig('Supplementary_Figure_S1_calibration.svg')
OUT2 = config.fig('Supplementary_Figure_S2_dca.svg')
OUTP = config.fig('Supplementary_Figure_preview.html')

PT_MM = 25.4 / 72.0
W = 190.0 / PT_MM           # 538.58 pt
W2 = 140.0 / PT_MM          # 396.85 pt, the 1.5 column width

INK = '#1a1a1a'
GREY = '#5f5f5f'
MID = '#9a9a9a'
RULE = '#c8c8c8'
COL = {'pelvis': '#1f4e79', 'rad': '#a8322d', 'both': '#2f6b3a'}
LABEL = {
    'pelvis': 'Pelvic bone volume',
    'rad': 'Rad-score (4 features)',
    'both': 'Bone volume + Rad-score',
}

FS_AXIS = 8.0
FS_TICK = 7.0
FS_NOTE = 7.0
FS_TITLE = 7.5


def wilson(k, n, z=1.96):
    """95% Wilson score interval for a proportion."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def head(w_pt, h_pt):
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%.2fmm" height="%.2fmm" '
           'viewBox="0 0 %.2f %.2f" font-family="Arial, Helvetica, sans-serif">'
           % (w_pt * PT_MM, h_pt * PT_MM, w_pt, h_pt)]
    out.append('<rect width="%.2f" height="%.2f" fill="#ffffff"/>' % (w_pt, h_pt))
    return out


# ----------------------------------------------------------------------------
# Figure S1: calibration, three panels
# ----------------------------------------------------------------------------

def figure_s1(path):
    data = json.load(open(SRC, encoding='utf-8'))
    models = data['models']
    order = ['pelvis', 'rad', 'both']

    pw, ph = 138.0, 148.0
    gap = 34.0
    left, top = 40.0, 30.0
    H = top + ph + 62.0

    out = head(W, H)

    for i, key in enumerate(order):
        x0 = left + i * (pw + gap)
        y0 = top
        out.append('<g>')
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
                   'stroke="%s" stroke-width="0.5"/>' % (x0, y0, pw, ph, MID))

        def X(v, x0=x0):
            return x0 + v * pw

        def Y(v, y0=y0):
            return y0 + ph - v * ph

        # identity line
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.5" stroke-dasharray="3 2.5"/>'
                   % (X(0), Y(0), X(1), Y(1), MID))

        for t in (0.0, 0.25, 0.5, 0.75, 1.0):
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="0.5"/>' % (X(t), Y(0), X(t), Y(0) + 3, GREY))
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="0.5"/>' % (X(0), Y(t), X(0) - 3, Y(t), GREY))
            out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                       'text-anchor="middle">%.2g</text>'
                       % (X(t), Y(0) + 11, FS_TICK, GREY, t))
            out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                       'text-anchor="end">%.2g</text>'
                       % (X(0) - 4, Y(t) + 2.4, FS_TICK, GREY, t))

        bins = models[key]['cal']['bins']
        cal = models[key]['cal']
        # error bars first so the markers sit on top
        for b in bins:
            k = round(b['obs'] * b['n'])
            lo, hi = wilson(k, b['n'])
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="0.8"/>'
                       % (X(b['pred']), Y(lo), X(b['pred']), Y(hi), COL[key]))
        for b in bins:
            out.append('<circle cx="%.1f" cy="%.1f" r="2.6" fill="%s"/>'
                       % (X(b['pred']), Y(b['obs']), COL[key]))

        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle" font-weight="bold">%s</text>'
                   % (x0 + pw / 2, y0 - 8, FS_TITLE, INK, esc(LABEL[key])))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">slope %.2f  HL P = %.3f</text>'
                   % (x0 + pw / 2, y0 + ph + 14, FS_NOTE, GREY,
                      cal['slope'], cal['HL_p']))
        out.append('</g>')

    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">Predicted probability</text>'
               % (left + (3 * pw + 2 * gap) / 2, H - 10, FS_AXIS, INK))
    out.append('<text x="11" y="%.1f" font-size="%.1f" fill="%s" text-anchor="middle" '
               'transform="rotate(-90 11 %.1f)">Observed proportion</text>'
               % (top + ph / 2, FS_AXIS, INK, top + ph / 2))
    out.append('</svg>')
    open(path, 'w', encoding='utf-8').write('\n'.join(out))


# ----------------------------------------------------------------------------
# Figure S2: decision curve analysis
# ----------------------------------------------------------------------------

def figure_s2(path):
    data = json.load(open(SRC, encoding='utf-8'))
    th = data['dca']['thresholds']
    cv = data['dca']['curves']

    H = 190.0
    left, top = 42.0, 26.0
    pw = 230.0
    ph = H - top - 44.0

    ymin, ymax = -0.10, 0.25

    def X(t):
        return left + (t - th[0]) / (th[-1] - th[0]) * pw

    def Y(v):
        v = max(ymin, min(ymax, v))
        return top + ph - (v - ymin) / (ymax - ymin) * ph

    series = [
        ('全部治疗', 'Treat all', '#b0b0b0', '4 3', 0.8),
        ('骨盆体积', 'Pelvic bone volume', COL['pelvis'], None, 1.2),
        ('Rad-score', 'Rad-score', COL['rad'], '5 2.6', 1.2),
        ('联合', 'Bone volume + Rad-score', COL['both'], '7 2 1.6 2', 1.2),
    ]

    out = head(W2, H)
    out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
               'stroke="%s" stroke-width="0.5"/>' % (left, top, pw, ph, MID))

    for v in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25):
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="%.2f"/>'
                   % (left, Y(v), left + pw, Y(v), MID if v == 0 else RULE,
                      0.5 if v == 0 else 0.25))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end">%.2f</text>'
                   % (left - 5, Y(v) + 2.4, FS_TICK, GREY, v))
    if ymin < 0:
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.5"/>' % (left, Y(0), left + pw, Y(0), INK))

    for t in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.5"/>' % (X(t), top + ph, X(t), top + ph + 3, GREY))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%.1f</text>'
                   % (X(t), top + ph + 12, FS_TICK, GREY, t))

    for zh, en, col, dash, width in series:
        pts = ' '.join('%.1f,%.1f' % (X(t), Y(v)) for t, v in zip(th, cv[zh]))
        d = ' stroke-dasharray="%s"' % dash if dash else ''
        out.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="%.1f"%s '
                   'stroke-linejoin="round"/>' % (pts, col, width, d))

    # net benefit at 0.20, the value quoted in the text
    i20 = min(range(len(th)), key=lambda k: abs(th[k] - 0.20))
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
               'stroke-width="0.5" stroke-dasharray="2.5 2"/>'
               % (X(0.20), Y(0), X(0.20), Y(0.20), MID))
    out.append('<circle cx="%.1f" cy="%.1f" r="2.2" fill="%s"/>'
               % (X(0.20), Y(cv['骨盆体积'][i20]), COL['pelvis']))

    lx, ly = left + pw + 10, top + 14
    for i, (zh, en, col, dash, width) in enumerate(series):
        yy = ly + i * 14
        d = ' stroke-dasharray="%s"' % dash if dash else ''
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="%.1f"%s/>' % (lx, yy, lx + 16, yy, col, width, d))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">%s</text>'
                   % (lx + 20, yy + 2.4, FS_NOTE, INK, esc(en)))
    yy = ly + len(series) * 14
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
               'stroke-width="0.8"/>' % (lx, yy, lx + 16, yy, GREY))
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">Treat none</text>'
               % (lx + 20, yy + 2.4, FS_NOTE, GREY))

    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">Threshold probability</text>'
               % (left + pw / 2, H - 10, FS_AXIS, INK))
    out.append('<text x="11" y="%.1f" font-size="%.1f" fill="%s" text-anchor="middle" '
               'transform="rotate(-90 11 %.1f)">Net benefit</text>'
               % (top + ph / 2, FS_AXIS, INK, top + ph / 2))
    out.append('</svg>')
    open(path, 'w', encoding='utf-8').write('\n'.join(out))


# ----------------------------------------------------------------------------
# Preview page
# ----------------------------------------------------------------------------

def preview(path):
    def embed(p):
        s = open(p, encoding='utf-8').read()
        return s.split('?>')[-1] if s.startswith('<?xml') else s

    f1, f2 = embed(OUT1), embed(OUT2)
    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Supplementary figures S1-S2</title>
<style>
body{margin:0;padding:26px 30px 46px;background:#f4f5f7;
     font-family:Arial,Helvetica,"PingFang SC",sans-serif;color:#1a1a1a}
h1{font-size:15px;font-weight:700;margin:0 0 3px}
p.sub{font-size:12px;color:#5f5f5f;margin:0 0 20px;line-height:1.6}
.card{background:#fff;border:1px solid #dcdee2;border-radius:3px;
      padding:16px 18px 10px;margin-bottom:18px;max-width:760px}
.card svg{width:100%;height:auto;display:block}
.cap{font-size:11.5px;color:#3d3d3d;line-height:1.6;margin-top:10px;
     border-top:1px solid #eceef1;padding-top:8px}
.tag{display:inline-block;font-size:10px;font-weight:700;letter-spacing:.05em;
     color:#1f4e79;background:#e9eff6;border-radius:2px;padding:2px 6px;margin-bottom:12px}
</style></head><body>
<h1>Supplementary figures, rebuilt at final size</h1>
<p class="sub">S1 at the 190 mm double column, S2 at the 140 mm width. Same
point-based typography as the main figures.</p>
<div class="card"><span class="tag">FIGURE S1</span>__F1__
<div class="cap"><b>Figure S1.</b> Calibration of the three models on out-of-fold
predictions. Points are observed proportions within bins of predicted probability,
bars the 95% Wilson interval; the dashed line is the identity. Slope and the
Hosmer-Lemeshow P value are given beneath each panel.</div></div>
<div class="card"><span class="tag">FIGURE S2</span>__F2__
<div class="cap"><b>Figure S2.</b> Decision curve analysis for the three models,
against treating all and treating none. The marker on the bone volume curve is the
net benefit at a threshold of 0.20, the value quoted in the text.</div></div>
</body></html>"""
    html = html.replace('__F1__', f1).replace('__F2__', f2)
    open(path, 'w', encoding='utf-8').write(html)


if __name__ == '__main__':
    figure_s1(OUT1)
    figure_s2(OUT2)
    preview(OUTP)
    for p in (OUT1, OUT2, OUTP):
        print('written', os.path.basename(p), os.path.getsize(p), 'bytes')
