"""Main-text Figures 1 to 3, built at final physical size for Radiotherapy and Oncology.

Canvas is in points with an explicit millimetre width, so that a font-size of 7
is 7 pt after printing and not merely 7 units on some arbitrary canvas. Target
width is the Elsevier double column, 190 mm. Everything here is pure line art,
which the publisher asks at 1000 dpi: render at 190 mm / 1000 dpi (7480 px) and
no resampling is needed.

Typography follows the Elsevier floor: 7 pt for normal text, nothing below 6 pt,
axis and data lines between 0.5 and 1.5 pt. Series are separated by line style
and fill depth as well as by hue, so the figure survives a grayscale proof.

Inputs (read-only): outputs/merged_dataset.csv, outputs/phase6_mainline.json,
                    outputs/stability_results.json
Outputs: Figure_1_volume.svg, Figure_2_roc.svg, Figure_3_stability.svg,
         Figure_preview.html
"""

import csv
import json
import math
import os

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
BASE = config.OUT_DIR
CSV = os.path.join(BASE, 'merged_dataset.csv')
MAIN = os.path.join(BASE, 'phase6_mainline.json')
STAB = os.path.join(BASE, 'stability_results.json')

OUT1 = config.fig('Figure_1_volume.svg')
OUT2 = config.fig('Figure_2_roc.svg')
OUT3 = config.fig('Figure_3_stability.svg')
OUTP = config.fig('Figure_preview.html')

PT_MM = 25.4 / 72.0          # one point in millimetres
W = 190.0 / PT_MM            # 538.58 pt, the double column width

INK = '#1a1a1a'
GREY = '#5f5f5f'
MID = '#9a9a9a'
RULE = '#c8c8c8'
COL = {'pelvis': '#1f4e79', 'rad': '#a8322d', 'both': '#2f6b3a'}
DASH = {'pelvis': None, 'rad': '5 2.6', 'both': '7 2 1.6 2'}

# font sizes, in points on the printed page
FS_PANEL = 9.0
FS_AXIS = 8.0
FS_TICK = 7.0
FS_NOTE = 7.0
FS_SMALL = 6.5
FS_GROUP = 7.5

LABEL = {
    'pelvis': 'Pelvic bone volume',
    'rad': 'Rad-score',
    'both': 'Bone volume + Rad-score',
}
AUC = {'pelvis': '0.784', 'rad': '0.810', 'both': '0.821'}


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def quantile(vals, q):
    """Linear interpolation, matching numpy's default (type 7)."""
    v = sorted(vals)
    if len(v) == 1:
        return v[0]
    pos = (len(v) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return v[lo]
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


def seedf(s):
    """Deterministic pseudo-random in [0,1) from a patient name."""
    h = 0
    for ch in s:
        h = (h * 131 + ord(ch)) % 1000003
    return (h % 997) / 997.0


def head(h_pt):
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%.2fmm" height="%.2fmm" '
           'viewBox="0 0 %.2f %.2f" font-family="Arial, Helvetica, sans-serif">'
           % (W * PT_MM, h_pt * PT_MM, W, h_pt)]
    out.append('<rect width="%.2f" height="%.2f" fill="#ffffff"/>' % (W, h_pt))
    return out


def panel_letter(out, x, y, letter):
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'font-weight="bold">%s</text>' % (x, y, FS_PANEL, INK, letter))


def ylab(out, y, text):
    out.append('<text x="11" y="%.1f" font-size="%.1f" fill="%s" text-anchor="middle" '
               'transform="rotate(-90 11 %.1f)">%s</text>'
               % (y, FS_AXIS, INK, y, esc(text)))


# ----------------------------------------------------------------------------
# Figure 1: bone volume distribution (A) and tertile event rate (B)
# ----------------------------------------------------------------------------

def figure1(path):
    rows = list(csv.DictReader(open(CSV, encoding='utf-8-sig')))
    ev = [(r['patient_id'], float(r['pelvis_volume'])) for r in rows if r['group'] == '1']
    no = [(r['patient_id'], float(r['pelvis_volume'])) for r in rows if r['group'] == '0']
    main = json.load(open(MAIN, encoding='utf-8'))

    H = 250.0
    left, top, ph = 44.0, 26.0, 168.0
    pwA, pwB = 208.0, 168.0
    xB = left + pwA + 74.0
    ylo, yhi = 820.0, 1820.0

    out = head(H)

    def YA(v):
        return top + ph - (v - ylo) / (yhi - ylo) * ph

    # ---------------- panel A ----------------
    out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
               'stroke="%s" stroke-width="0.5"/>' % (left, top, pwA, ph, MID))
    for v in range(900, 1801, 200):
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.25" stroke-dasharray="2 3"/>'
                   % (left, YA(v), left + pwA, YA(v), RULE))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end">%d</text>'
                   % (left - 5, YA(v) + 2.4, FS_TICK, GREY, v))

    groups = [
        (no, left + pwA * 0.30, '#cfdae6', COL['pelvis'], COL['pelvis'],
         'No grade \u22652', 'n = 45'),
        (ev, left + pwA * 0.72, COL['rad'], '#7d221e', '#4a1512',
         'Grade \u22652', 'n = 15'),
    ]
    bw = 44.0

    for vals, cx, fill, stroke, dotc, lab, nlab in groups:
        ys = [v for _, v in vals]
        q1, q2, q3 = quantile(ys, 0.25), quantile(ys, 0.5), quantile(ys, 0.75)
        iqr = q3 - q1
        lo = min(v for v in ys if v >= q1 - 1.5 * iqr)
        hi = max(v for v in ys if v <= q3 + 1.5 * iqr)
        # whiskers
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.8"/>' % (cx, YA(hi), cx, YA(q3), stroke))
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.8"/>' % (cx, YA(q1), cx, YA(lo), stroke))
        for v in (hi, lo):
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="0.8"/>' % (cx - 9, YA(v), cx + 9, YA(v), stroke))
        # box
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" '
                   'stroke="%s" stroke-width="0.8"/>'
                   % (cx - bw / 2, YA(q3), bw, YA(q1) - YA(q3), fill, stroke))
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="1.2"/>'
                   % (cx - bw / 2, YA(q2), cx + bw / 2, YA(q2), stroke))
        # individual patients
        for name, v in vals:
            jx = cx + (seedf(name) - 0.5) * 1.5 * bw
            out.append('<circle cx="%.1f" cy="%.1f" r="1.9" fill="%s" '
                       'fill-opacity="0.6"/>' % (jx, YA(v), dotc))
        # mean and one standard deviation
        m = sum(ys) / len(ys)
        s = math.sqrt(sum((v - m) ** 2 for v in ys) / (len(ys) - 1))
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="1.4" stroke-linecap="round"/>'
                   % (cx - 17, YA(m), cx + 17, YA(m), INK))
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.5" stroke-dasharray="2.5 2"/>'
                   % (cx, YA(m - s), cx, YA(m + s), INK))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle" font-weight="bold">%.1f</text>'
                   % (cx, YA(m) - 5, FS_NOTE, INK, m))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%s</text>'
                   % (cx, top + ph + 12, FS_GROUP, INK, esc(lab)))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%s</text>'
                   % (cx, top + ph + 23, FS_SMALL, GREY, nlab))

    ylab(out, top + ph / 2, 'Pelvic bone volume (cm\u00b3)')
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">P = 0.001</text>'
               % (left + pwA / 2, top - 8, FS_NOTE, INK))
    panel_letter(out, left - 33, top - 7, 'A')

    # ---------------- panel B ----------------
    t = main['tertiles']['counts']
    ymax = 55.0

    def YB(v):
        return top + ph - v / ymax * ph

    out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
               'stroke="%s" stroke-width="0.5"/>' % (xB, top, pwB, ph, MID))
    for v in (0, 20, 40):
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="%.2f"/>'
                   % (xB, YB(v), xB + pwB, YB(v),
                      MID if v == 0 else RULE, 0.5 if v == 0 else 0.25))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end">%d</text>'
                   % (xB - 5, YB(v) + 2.4, FS_TICK, GREY, v))

    barw = 28.0
    shades = ['#1f4e79', '#5b87b0', '#b7cadd']
    names = ['T1', 'T2', 'T3']
    subs = ['\u22641158.7', '1158.7-1356.1', '>1356.1']
    for i, (nn, k) in enumerate(zip(names, t)):
        cx = xB + pwB * (0.19 + 0.31 * i)
        rate = 100.0 * k[1] / k[0]
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>'
                   % (cx - barw / 2, YB(rate), barw, YB(0) - YB(rate), shades[i]))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle" font-weight="bold">%.1f</text>'
                   % (cx, YB(rate) - 5, FS_NOTE, INK, rate))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%d/%d</text>'
                   % (cx, YB(rate) - 14, FS_SMALL, GREY, k[1], k[0]))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%s</text>'
                   % (cx, top + ph + 12, FS_GROUP, INK, nn))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%s</text>'
                   % (cx, top + ph + 23, FS_SMALL, GREY, subs[i]))

    ylab(out, top + ph / 2, 'Grade \u22652 leukopenia (%)')
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">P for trend = 0.018</text>'
               % (xB + pwB / 2, top - 8, FS_NOTE, INK))
    panel_letter(out, xB - 20, top - 7, 'B')

    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">Pelvic bone volume tertile (cm\u00b3)</text>'
               % (xB + pwB / 2, H - 8, FS_AXIS, INK))
    out.append('</svg>')
    open(path, 'w', encoding='utf-8').write('\n'.join(out))


# ----------------------------------------------------------------------------
# Figure 2: ROC for the three models (A), operating point of bone volume (B)
# ----------------------------------------------------------------------------

def roc_panel(out, x0, y0, pw, ph, keys, data, marker_cut=False):
    def X(v):
        return x0 + v * pw

    def Y(v):
        return y0 + ph - v * ph

    out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
               'stroke="%s" stroke-width="0.5"/>' % (x0, y0, pw, ph, MID))
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
                   % (X(t), Y(0) + 12, FS_TICK, GREY, t))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end">%.2g</text>'
                   % (X(0) - 5, Y(t) + 2.4, FS_TICK, GREY, t))

    for k in keys:
        pts = ' '.join('%.1f,%.1f' % (X(f), Y(t)) for f, t in data['models'][k]['roc'])
        d = ' stroke-dasharray="%s"' % DASH[k] if DASH[k] else ''
        out.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.2"%s '
                   'stroke-linejoin="round"/>' % (pts, COL[k], d))

    if marker_cut:
        cut = data['cutoff']
        fpr, tpr = 1 - cut['spec'], cut['sens']
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.5" stroke-dasharray="2.5 2"/>'
                   % (X(0), Y(tpr), X(fpr), Y(tpr), INK))
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.5" stroke-dasharray="2.5 2"/>'
                   % (X(fpr), Y(0), X(fpr), Y(tpr), INK))
        out.append('<circle cx="%.1f" cy="%.1f" r="3" fill="#ffffff" stroke="%s" '
                   'stroke-width="1.1"/>' % (X(fpr), Y(tpr), COL['pelvis']))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end">1170.4 cm\u00b3</text>'
                   % (X(fpr) - 5, Y(tpr) - 6, FS_NOTE, INK))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">'
                   'sens 73.3%%</text>' % (X(fpr) + 6, Y(tpr) + 3, FS_SMALL, GREY))


def figure2(path):
    data = json.load(open(MAIN, encoding='utf-8'))
    H = 250.0
    left, top, ph = 44.0, 26.0, 160.0
    pwA = 196.0
    xB = left + pwA + 58.0
    pwB = 196.0

    out = head(H)
    roc_panel(out, left, top, pwA, ph, ['pelvis', 'rad', 'both'], data)
    roc_panel(out, xB, top, pwB, ph, ['pelvis'], data, marker_cut=True)

    # legend under panel A; line style carries the same distinction as colour
    ly = top + ph + 20
    for i, k in enumerate(['pelvis', 'rad', 'both']):
        yy = ly + i * 13
        d = ' stroke-dasharray="%s"' % DASH[k] if DASH[k] else ''
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="1.2"%s/>' % (left + 2, yy, left + 22, yy, COL[k], d))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">%s</text>'
                   % (left + 26, yy + 2.4, FS_NOTE, INK, esc(LABEL[k])))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end">AUC %s</text>'
                   % (left + pwA - 1, yy + 2.4, FS_NOTE, GREY, AUC[k]))

    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">1 - specificity</text>'
               % (left + pwA / 2, H - 8, FS_AXIS, INK))
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">1 - specificity</text>'
               % (xB + pwB / 2, H - 8, FS_AXIS, INK))
    ylab(out, top + ph / 2, 'Sensitivity')
    panel_letter(out, left - 33, top - 7, 'A')
    panel_letter(out, xB - 20, top - 7, 'B')
    out.append('</svg>')
    open(path, 'w', encoding='utf-8').write('\n'.join(out))


# ----------------------------------------------------------------------------
# Figure 3: selection frequency across 100 repeated cross-validations
# ----------------------------------------------------------------------------

def figure3(path):
    feats = json.load(open(STAB, encoding='utf-8'))['features']
    feats = sorted(feats, key=lambda f: -f['freq_min'])[:24]

    rowh = 10.6
    left, top = 182.0, 30.0
    pw = W - left - 42.0
    ph = rowh * len(feats)
    H = top + ph + 44.0

    out = head(H)

    def X(v):
        return left + v * pw

    out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
               'stroke="%s" stroke-width="0.5"/>' % (left, top, pw, ph, MID))
    for v in (0, 0.2, 0.4, 0.6, 0.8):
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="0.25" stroke-dasharray="2 3"/>'
                   % (X(v), top, X(v), top + ph, RULE))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="middle">%d</text>'
                   % (X(v), top + ph + 12, FS_TICK, GREY, int(v * 100)))

    primary = {
        'wavelet-HHL:ngtdm Busyness',
        'log-sigma-1-mm-3D:ngtdm Busyness',
        'original:shape MeshVolume',
        'wavelet-HHH:glrlm RunVariance',
    }

    for i, f in enumerate(feats):
        y = top + i * rowh
        key = '%s %s' % (f['class'], f['name'])
        is_p = key in primary
        col = COL['pelvis'] if is_p else '#b7cadd'
        if f['freq_min'] > 0:
            out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>'
                       % (left, y + 1.8, f['freq_min'] * pw, rowh - 3.4, col))
        else:
            # a candidate that was never selected still gets a mark at zero,
            # so that an empty row reads as a result and not as a rendering fault
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="0.8"/>'
                       % (left, y + rowh / 2, left + 2.5, y + rowh / 2, '#8fa8c0'))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
                   'text-anchor="end"%s>%s</text>'
                   % (left - 5, y + rowh - 3.1, FS_SMALL, INK if is_p else GREY,
                      ' font-weight="bold"' if is_p else '', esc(key)))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">%d</text>'
                   % (f['freq_min'] * pw + left + 4, y + rowh - 3.1,
                      FS_SMALL, INK if is_p else GREY, round(f['freq_min'] * 100)))

    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" '
               'text-anchor="middle">Selected in (%%) of 100 repeats</text>'
               % (left + pw / 2, H - 10, FS_AXIS, INK))

    lx, ly = left + 2, top - 11
    out.append('<rect x="%.1f" y="%.1f" width="7" height="7" fill="%s"/>'
               % (lx, ly - 6, COL['pelvis']))
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">'
               'Four primary features</text>' % (lx + 11, ly, FS_SMALL, INK))
    out.append('<rect x="%.1f" y="%.1f" width="7" height="7" fill="#b7cadd"/>'
               % (lx + 118, ly - 6))
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">'
               'Other candidates</text>' % (lx + 129, ly, FS_SMALL, GREY))
    out.append('</svg>')
    open(path, 'w', encoding='utf-8').write('\n'.join(out))


# ----------------------------------------------------------------------------
# Preview page
# ----------------------------------------------------------------------------

def preview(path):
    def embed(p):
        s = open(p, encoding='utf-8').read()
        return s.split('?>')[-1] if s.startswith('<?xml') else s

    f1, f2, f3 = embed(OUT1), embed(OUT2), embed(OUT3)
    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Main figures 1-3</title>
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
<h1>Main-text display items, built at 190 mm</h1>
<p class="sub">Radiotherapy and Oncology. The canvas is 190 mm wide with
point-based typography, so 7 pt here is 7 pt in print. Rendered at browser
width; the printed figure is about twice this size.</p>
<div class="card"><span class="tag">FIGURE 1</span>__F1__
<div class="cap"><b>Figure 1.</b> Pelvic bone volume and grade 2 or worse leukopenia.
(A) Distribution in the two groups; the box spans the interquartile range with the
median as a heavy line, whiskers extend to 1.5 times the interquartile range, the
horizontal bar marks the mean and the dashed line one standard deviation.
(B) Event rate by tertile of pelvic bone volume, with the number of events per tertile
above each bar.</div></div>
<div class="card"><span class="tag">FIGURE 2</span>__F2__
<div class="cap"><b>Figure 2.</b> Discrimination and operating point.
(A) Receiver operating characteristic curves for the three models.
(B) Pelvic bone volume alone, with the Youden-optimal cut-off of 1170.4 cm&#179;
marked. Point estimates of the area under the curve are in the legend;
confidence intervals are in Table 3.</div></div>
<div class="card"><span class="tag">FIGURE 3</span>__F3__
<div class="cap"><b>Figure 3.</b> Selection stability across 100 repeated 10-fold
cross-validations. Bars show the proportion of repeats in which each candidate was
retained at the minimum cross-validated lambda. The 24 most frequently selected
candidates are shown, in descending order.</div></div>
</body></html>"""
    html = html.replace('__F1__', f1).replace('__F2__', f2).replace('__F3__', f3)
    open(path, 'w', encoding='utf-8').write(html)


if __name__ == '__main__':
    figure1(OUT1)
    figure2(OUT2)
    figure3(OUT3)
    preview(OUTP)
    for p in (OUT1, OUT2, OUT3, OUTP):
        print('written', os.path.basename(p), os.path.getsize(p), 'bytes')
