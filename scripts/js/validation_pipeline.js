// =====================================================================
// 全流程内部验证计算内核 v2（Node.js 纯标准库）
// ---------------------------------------------------------------------
// 管线（与 R 版 radiomics_bm_pipeline_validation.R 对应）：
//   单因素 Wilcoxon p<0.05 筛选 -> LASSO(logistic, 坐标下降, 10折CV,
//   deviance, lambda.1se) -> 无惩罚 logistic(IRLS) 重拟合 -> AUC
//   Bootstrap 重抽样内完整重跑上述管线，乐观度=AUC(样本内)-AUC(原样本)
// v2 提速：glmnet 同款外层二次逼近(IRLS式)+内层加权最小二乘坐标下降
//          （内层无 exp），活动集循环，收敛容差自适应
// 用法: node validation_pipeline.js <input.json> <outdir> <B> [seed]
// =====================================================================
"use strict";
const fs = require("fs");
const path = require("path");

// ---------------- 数学工具 ----------------
function sigmoid(x) { x = Math.max(-30, Math.min(30, x)); return 1 / (1 + Math.exp(-x)); }
function soft(x, t) { return x > 0 ? Math.max(x - t, 0) : Math.min(x + t, 0); }

function erf(x) { // Abramowitz-Stegun 7.1.26
  const s = x < 0 ? -1 : 1; x = Math.abs(x);
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741,
        a4 = -1.453152027, a5 = 1.061405429, pp = 0.3275911;
  const t = 1 / (1 + pp * x);
  const y = 1 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return s * y;
}
const phi = z => 0.5 * (1 + erf(z / Math.SQRT2));

function gammaln(x) {
  const g = [76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let xx = x, tmp = xx + 5.5, ser = 1.000000000190015;
  tmp -= (xx + 0.5) * Math.log(tmp);
  for (let j = 0; j < 6; j++) { xx += 1; ser += g[j] / xx; }
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}
function gser(a, x) {
  const IT = 600, EPS = 3e-14;
  let ap = a, sum = 1 / a, del = sum;
  for (let n = 1; n <= IT; n++) { ap++; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * EPS) break; }
  return sum * Math.exp(-x + a * Math.log(x) - gammaln(a));
}
function gcf(a, x) {
  const IT = 600, EPS = 3e-14, FPMIN = 1e-300;
  let b = x + 1 - a, c = 1 / FPMIN, d = 1 / b, h = d;
  for (let i = 1; i <= IT; i++) {
    const an = -i * (i - a);
    b += 2; d = an * d + b; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; const del = d * c; h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return Math.exp(-x + a * Math.log(x) - gammaln(a)) * h;
}
function chi2Sf(df, x) { if (!(x > 0)) return 1; const a = df / 2, z = x / 2; return (z < a + 1) ? (1 - gser(a, z)) : gcf(a, z); }

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ---------------- 统计检验 ----------------
function wilcoxP(x0, x1) {
  const n0 = x0.length, n1 = x1.length, N = n0 + n1;
  if (!n0 || !n1 || N < 4) return null;
  const all = [];
  for (const v of x0) all.push({ v, g: 0 });
  for (const v of x1) all.push({ v, g: 1 });
  all.sort((a, b) => a.v - b.v);
  const ranks = new Array(N);
  let i = 0, tieSum = 0;
  while (i < N) {
    let j = i;
    while (j + 1 < N && all[j + 1].v === all[i].v) j++;
    const r = (i + j) / 2 + 1;
    const t = j - i + 1;
    if (t > 1) tieSum += t * t * t - t;
    for (let k = i; k <= j; k++) ranks[k] = r;
    i = j + 1;
  }
  let R0 = 0;
  for (let k = 0; k < N; k++) if (all[k].g === 0) R0 += ranks[k];
  const U1 = R0 - n0 * (n0 + 1) / 2, mu = n0 * n1 / 2;
  const varu = n0 * n1 / 12 * ((N + 1) - tieSum / (N * (N - 1)));
  if (!(varu > 0)) return 1.0;
  const z = (Math.abs(U1 - mu) - 0.5) / Math.sqrt(varu);
  return 2 * (1 - phi(z));
}

function rankdataAvg(arr) {
  const idx = arr.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
  const ranks = new Array(arr.length);
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const r = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++) ranks[idx[k][1]] = r;
    i = j + 1;
  }
  return ranks;
}

function aucOf(pred, y) {
  if (new Set(pred).size < 2) return NaN;
  const r = rankdataAvg(pred);
  let n1 = 0, R1 = 0;
  for (let i = 0; i < y.length; i++) if (y[i] === 1) { n1++; R1 += r[i]; }
  const n0 = y.length - n1;
  if (!n0 || !n1) return NaN;
  return (R1 - n1 * (n1 + 1) / 2) / (n0 * n1);
}

// ---------------- LASSO logistic 路径（外层二次逼近+内层WLS坐标下降） ----------------
function lassoPath(X, y, lamsFixed, nlam, ratio) {
  const n = X.length, p = X[0].length;
  const cols = [], mean = [], sd = [], keep = [];
  for (let j = 0; j < p; j++) {
    let m = 0; for (let i = 0; i < n; i++) m += X[i][j]; m /= n;
    let s = 0; for (let i = 0; i < n; i++) { const d = X[i][j] - m; s += d * d; }
    s = Math.sqrt(s / n);
    if (s > 1e-12) {
      keep.push(j); mean.push(m); sd.push(s);
      const c = new Float64Array(n);
      for (let i = 0; i < n; i++) c[i] = (X[i][j] - m) / s;
      cols.push(c);
    }
  }
  const pk = keep.length;
  if (pk === 0) return null;
  let lams;
  if (lamsFixed) { lams = lamsFixed; }
  else {
    const ybar = y.reduce((a, b) => a + b, 0) / n;
    let gm = 0;
    for (let k = 0; k < pk; k++) {
      const c = cols[k]; let g = 0;
      for (let i = 0; i < n; i++) g += c[i] * (y[i] - ybar);
      gm = Math.max(gm, Math.abs(g));
    }
    const lmax = gm / n;
    if (!(lmax > 0)) return null;
    lams = [];
    for (let t = 0; t < nlam; t++) lams.push(lmax * Math.pow(ratio, t / Math.max(nlam - 1, 1)));
  }
  const nl = lams.length;
  const ybar0 = y.reduce((a, b) => a + b, 0) / n;
  let b0 = Math.log(Math.max(ybar0, 1e-6) / Math.max(1 - ybar0, 1e-6));
  const b = new Float64Array(pk);
  const eta = new Float64Array(n).fill(b0);
  const pr = new Float64Array(n);
  for (let i = 0; i < n; i++) pr[i] = sigmoid(b0);
  const B0raw = new Array(nl), Braw = new Array(nl);

  for (let li = 0; li < nl; li++) {
    const lamEff = n * lams[li];
    for (let cyc = 0; cyc < 200; cyc++) {
      // 截距牛顿步
      let sw = 0, sr = 0;
      for (let i = 0; i < n; i++) { const w = pr[i] * (1 - pr[i]); sw += w; sr += y[i] - pr[i]; }
      if (sw > 1e-12) {
        const d0 = sr / sw; b0 += d0;
        for (let i = 0; i < n; i++) eta[i] += d0;
        for (let i = 0; i < n; i++) pr[i] = sigmoid(eta[i]);
      }
      let maxd = 0;
      const fullSweep = (cyc % 8 === 0);
      for (let k = 0; k < pk; k++) {
        if (!fullSweep && b[k] === 0) continue;          // 活动集：仅迭代非零坐标
        const c = cols[k];
        let g = 0, H = 0;
        for (let i = 0; i < n; i++) {
          const w = pr[i] * (1 - pr[i]);
          const ci = c[i];
          g += ci * (y[i] - pr[i]);
          H += ci * ci * w;
        }
        if (H < 1e-12) continue;
        const old = b[k];
        const nb = soft(g + H * old, lamEff) / H;
        const d = nb - old;
        if (d !== 0) {
          b[k] = nb;
          for (let i = 0; i < n; i++) eta[i] += d * c[i];
          for (let i = 0; i < n; i++) pr[i] = sigmoid(eta[i]);
        }
        const ad = Math.abs(d);
        if (ad > maxd) maxd = ad;
      }
      let maxb = 0;
      for (let k = 0; k < pk; k++) { const ab = Math.abs(b[k]); if (ab > maxb) maxb = ab; }
      if (maxd <= 1e-6 * (1 + maxb)) break;
    }
    const cr = new Array(pk);
    let b0r = b0;
    for (let k = 0; k < pk; k++) { cr[k] = b[k] / sd[k]; b0r -= b[k] * mean[k] / sd[k]; }
    B0raw[li] = b0r; Braw[li] = cr;
  }
  return { lams, B0raw, Braw, keep };
}

function predAt(model, m, X) {
  const b0 = model.B0raw[m], cr = model.Braw[m], keep = model.keep;
  const out = new Array(X.length);
  for (let i = 0; i < X.length; i++) {
    let e = b0; const xi = X[i];
    for (let k = 0; k < keep.length; k++) e += cr[k] * xi[keep[k]];
    out[i] = sigmoid(e);
  }
  return out;
}

// ---------------- CV 选 lambda.1se ----------------
function stratFolds(y, nfold, rng) {
  const fold = new Array(y.length).fill(0);
  for (const cls of [0, 1]) {
    const idx = [];
    for (let i = 0; i < y.length; i++) if (y[i] === cls) idx.push(i);
    for (let i = idx.length - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      const t = idx[i]; idx[i] = idx[j]; idx[j] = t;
    }
    idx.forEach((pi, k) => { fold[pi] = k % nfold; });
  }
  return fold;
}

function cvFit(X, y, rng, nfold, nlam, ratio) {
  const full = lassoPath(X, y, null, nlam, ratio);
  if (!full) return null;
  const lams = full.lams, nl = lams.length, n = X.length;
  const fold = stratFolds(y, nfold, rng);
  const foldDev = Array.from({ length: nfold }, () => new Float64Array(nl));
  const foldN = new Array(nfold).fill(0);
  for (let k = 0; k < nfold; k++) {
    const Xtr = [], ytr = [], Xte = [], yte = [];
    for (let i = 0; i < n; i++) {
      if (fold[i] === k) { Xte.push(X[i]); yte.push(y[i]); foldN[k]++; }
      else { Xtr.push(X[i]); ytr.push(y[i]); }
    }
    if (!yte.length || !ytr.length) continue;
    const f = lassoPath(Xtr, ytr, lams, nlam, ratio);
    if (!f) return null;
    for (let m = 0; m < nl; m++) {
      const pp = predAt(f, m, Xte);
      let d = 0;
      for (let i = 0; i < yte.length; i++) {
        const q = Math.min(Math.max(pp[i], 1e-7), 1 - 1e-7);
        d += -2 * (yte[i] * Math.log(q) + (1 - yte[i]) * Math.log(1 - q));
      }
      foldDev[k][m] = d;
    }
  }
  const cvm = new Array(nl), se = new Array(nl);
  for (let m = 0; m < nl; m++) {
    const fm = []; let tot = 0;
    for (let k = 0; k < nfold; k++) {
      if (!foldN[k]) continue;
      const v = foldDev[k][m] / foldN[k];
      fm.push(v); tot += foldDev[k][m];
    }
    cvm[m] = tot / n;
    const mu = fm.reduce((a, b) => a + b, 0) / fm.length;
    const s2 = fm.reduce((a, b) => a + (b - mu) * (b - mu), 0) / (fm.length - 1);
    se[m] = Math.sqrt(s2 / fm.length);
  }
  let iMin = 0;
  for (let m = 1; m < nl; m++) if (cvm[m] < cvm[iMin]) iMin = m;
  const thr = cvm[iMin] + se[iMin];
  let i1se = iMin;
  for (let m = 0; m < nl; m++) { if (cvm[m] <= thr) { i1se = m; break; } }
  return { i1se, full, cvm, se, lams };
}

// ---------------- 无惩罚 logistic (IRLS) ----------------
function solveLin(A, rhs) {
  const n = rhs.length;
  const M = A.map((row, i) => [...row, rhs[i]]);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    if (Math.abs(M[piv][c]) < 1e-12) return null;
    const t = M[c]; M[c] = M[piv]; M[piv] = t;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = M[r][c] / M[c][c];
      for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
    }
  }
  const x = new Array(n);
  for (let i = 0; i < n; i++) x[i] = M[i][n] / M[i][i];
  return x;
}

function irlsLogistic(Xsel, y) {
  const n = Xsel.length, k = Xsel[0].length, pp = k + 1;
  const Xd = Xsel.map(r => [1, ...r]);
  let beta = new Array(pp).fill(0);
  for (let it = 0; it < 60; it++) {
    const eta = new Array(n), pr = new Array(n), w = new Array(n), z = new Array(n);
    for (let i = 0; i < n; i++) {
      let e = 0;
      for (let a = 0; a < pp; a++) e += Xd[i][a] * beta[a];
      eta[i] = e; pr[i] = sigmoid(e);
      w[i] = Math.max(pr[i] * (1 - pr[i]), 1e-8);
      z[i] = e + (y[i] - pr[i]) / w[i];
    }
    const A = Array.from({ length: pp }, () => new Array(pp).fill(0));
    const rhs = new Array(pp).fill(0);
    for (let i = 0; i < n; i++) {
      const xi = Xd[i], wi = w[i], zi = z[i];
      for (let a = 0; a < pp; a++) {
        rhs[a] += wi * xi[a] * zi;
        for (let b = 0; b < pp; b++) A[a][b] += wi * xi[a] * xi[b];
      }
    }
    for (let a = 0; a < pp; a++) A[a][a] += 1e-9;
    const bn = solveLin(A, rhs);
    if (!bn) return null;
    let dm = 0;
    for (let a = 0; a < pp; a++) dm = Math.max(dm, Math.abs(bn[a] - beta[a]));
    beta = bn;
    if (dm < 1e-9) break;
  }
  for (let a = 1; a < pp; a++) if (Math.abs(beta[a]) > 30) return null;
  return beta;
}

function predictGlm(beta, Xsel) {
  return Xsel.map(r => {
    let e = beta[0];
    for (let j = 0; j < r.length; j++) e += beta[j + 1] * r[j];
    return sigmoid(e);
  });
}

// ---------------- 完整管线（带失败原因） ----------------
function fitPipeline(X, y, names, rng, nlam, ratio, pScreen) {
  const idx0 = [], idx1 = [];
  for (let i = 0; i < y.length; i++) (y[i] === 1 ? idx1 : idx0).push(i);
  if (idx0.length < 3 || idx1.length < 3) return { ok: false, reason: "class_too_small" };
  const cand = [];
  for (let j = 0; j < names.length; j++) {
    const col = X.map(r => r[j]);
    if (new Set(col).size < 3) continue;
    const x0 = idx0.map(i => col[i]), x1 = idx1.map(i => col[i]);
    const pv = wilcoxP(x0, x1);
    if (pv !== null && pv < pScreen) cand.push(j);
  }
  if (cand.length < 2) return { ok: false, reason: "screen_lt2" };
  const Xc = X.map(r => cand.map(j => r[j]));
  const cv = cvFit(Xc, y, rng, 10, nlam, ratio);
  if (!cv) return { ok: false, reason: "cv_fail" };
  const m = cv.i1se, cr = cv.full.Braw[m], keep = cv.full.keep;
  const selLocal = [];
  for (let k = 0; k < keep.length; k++) if (Math.abs(cr[k]) > 1e-12) selLocal.push(k);
  if (!selLocal.length) return { ok: false, reason: "lasso_zero" };
  const selOrig = selLocal.map(k => cand[keep[k]]);
  const Xsel = X.map(r => selOrig.map(j => r[j]));
  const beta = irlsLogistic(Xsel, y);
  if (!beta) return { ok: false, reason: "irls_fail" };
  return { ok: true, selOrig, beta, selNames: selOrig.map(j => names[j]), lambda: cv.lams[m], nCand: cand.length };
}

// ---------------- HL 检验 ----------------
function hlTest(y, pp) {
  const n = y.length, g = Math.min(10, Math.max(2, Math.floor(n / 6)));
  const idx = [...pp.keys()].sort((a, b) => pp[a] - pp[b]);
  const groups = Array.from({ length: g }, () => ({ o: 0, e: 0, m: 0 }));
  idx.forEach((pi, k) => {
    const gi = Math.min(g - 1, Math.floor(k * g / n));
    groups[gi].o += y[pi]; groups[gi].e += pp[pi]; groups[gi].m++;
  });
  let X2 = 0, used = 0;
  for (const gr of groups) {
    if (!gr.m) continue;
    const denom = gr.e * (1 - gr.e / gr.m);
    if (denom <= 1e-12) continue;
    X2 += (gr.o - gr.e) * (gr.o - gr.e) / denom;
    used++;
  }
  const df = Math.max(1, used - 2);
  return { X2, df, p: chi2Sf(df, X2) };
}

function loglik(y, pp) {
  let s = 0;
  for (let i = 0; i < y.length; i++) {
    const q = Math.min(Math.max(pp[i], 1e-12), 1 - 1e-12);
    s += y[i] * Math.log(q) + (1 - y[i]) * Math.log(1 - q);
  }
  return s;
}

// ---------------- 主流程 ----------------
async function main() {
  const [, , inPath, outDir, Bstr, seedStr] = process.argv;
  const B = parseInt(Bstr || "500", 10);
  const SEED = parseInt(seedStr || "20251207", 10);
  const NLAM = 100, RATIO = 1e-4, P_SCREEN = 0.05;
  const rng = mulberry32(SEED);
  fs.mkdirSync(outDir, { recursive: true });

  const inp = JSON.parse(fs.readFileSync(inPath, "utf8"));
  const names = inp.feature_names, y = inp.y, X = inp.X;
  const n = y.length;
  console.log(`样本 n=${n} | 事件 ${y.reduce((a, b) => a + b, 0)} (${(100 * y.reduce((a, b) => a + b, 0) / n).toFixed(1)}%) | 特征数 ${names.length}`);

  const t0 = Date.now();
  const full = fitPipeline(X, y, names, rng, NLAM, RATIO, P_SCREEN);
  if (!full.ok) { console.error("全样本管线失败: " + full.reason); process.exit(1); }
  const XselFull = X.map(r => full.selOrig.map(j => r[j]));
  const predFull = predictGlm(full.beta, XselFull);
  const aucApp = aucOf(predFull, y);
  console.log(`全样本候选特征数(Wilcoxon p<0.05): ${full.nCand}`);
  console.log(`全样本入选特征(${full.selNames.length}): ${full.selNames.join(", ")}`);
  console.log(`lambda.1se = ${full.lambda.toPrecision(4)}`);
  console.log(`表观 AUC = ${aucApp.toFixed(4)}`);

  const ybar = y.reduce((a, b) => a + b, 0) / n;
  const devNull = -2 * loglik(y, new Array(n).fill(ybar));
  const devFull = -2 * loglik(y, predFull);
  const lrtX2 = devNull - devFull, lrtDf = full.selOrig.length;
  const lrtP = chi2Sf(lrtDf, lrtX2);
  console.log(`整体模型 LRT: X2=${lrtX2.toFixed(2)}, df=${lrtDf}, p=${lrtP.toPrecision(4)}`);
  const hl = hlTest(y, predFull);
  console.log(`Hosmer-Lemeshow: X2=${hl.X2.toFixed(2)}, df=${hl.df}, p=${hl.p.toFixed(4)}`);

  console.log(`\n开始 Bootstrap 全流程验证 (B=${B}, 每次重抽样内完整重跑 筛选+LASSO+logistic)...`);
  const optimism = [], freq = {}, selSize = [], candN = [];
  const reasons = {};
  let skipped = 0;
  const t1 = Date.now();
  for (let b = 1; b <= B; b++) {
    const idx = [];
    for (let i = 0; i < n; i++) idx.push(Math.floor(rng() * n));
    const Xb = idx.map(i => X[i]), yb = idx.map(i => y[i]);
    const res = fitPipeline(Xb, yb, names, rng, NLAM, RATIO, P_SCREEN);
    if (!res.ok) { skipped++; reasons[res.reason] = (reasons[res.reason] || 0) + 1; continue; }
    const predB = predictGlm(res.beta, Xb.map(r => res.selOrig.map(j => r[j])));
    const predO = predictGlm(res.beta, X.map(r => res.selOrig.map(j => r[j])));
    const aB = aucOf(predB, yb), aO = aucOf(predO, y);
    if (!isFinite(aB) || !isFinite(aO)) { skipped++; reasons["auc_nan"] = (reasons["auc_nan"] || 0) + 1; continue; }
    optimism.push(aB - aO);
    selSize.push(res.selNames.length);
    candN.push(res.nCand);
    for (const nm of res.selNames) freq[nm] = (freq[nm] || 0) + 1;
    if (b % 50 === 0) {
      const el = (Date.now() - t1) / 1000;
      console.log(`  bootstrap ${b}/${B} (有效 ${optimism.length}, 跳过 ${skipped}, 用时 ${el.toFixed(0)}s)`);
    }
  }
  const valid = optimism.length;
  const meanO = optimism.reduce((a, b) => a + b, 0) / valid;
  const sortedO = [...optimism].sort((a, b) => a - b);
  const q = p => sortedO[Math.min(valid - 1, Math.max(0, Math.floor(p * valid)))];
  const aucCorr = aucApp - meanO;
  console.log(`\n================ 验证结果汇总 ================`);
  console.log(`表观 AUC          : ${aucApp.toFixed(4)}`);
  console.log(`乐观度均值        : ${meanO.toFixed(4)} (95%CI: ${q(0.025).toFixed(4)} ~ ${q(0.975).toFixed(4)})`);
  console.log(`乐观校正后 AUC    : ${aucCorr.toFixed(4)}`);
  console.log(`有效 bootstrap    : ${valid}/${B} (跳过 ${skipped}: ${JSON.stringify(reasons)})`);
  const meanSel = selSize.reduce((a, b) => a + b, 0) / (selSize.length || 1);
  const meanCand = candN.reduce((a, b) => a + b, 0) / (candN.length || 1);
  console.log(`重抽样内平均候选特征数: ${meanCand.toFixed(1)} | 平均入选特征数: ${meanSel.toFixed(2)}`);

  const freqArr = Object.entries(freq).sort((a, b) => b[1] - a[1]);
  console.log(`\n特征入选频次(前15):`);
  for (const [nm, c] of freqArr.slice(0, 15)) console.log(`  ${nm}: ${c} (${(100 * c / valid).toFixed(1)}%)`);

  const summary = [
    ["metric", "value"],
    ["apparent_AUC", aucApp.toFixed(4)],
    ["optimism_mean", meanO.toFixed(4)],
    ["optimism_2.5pct", q(0.025).toFixed(4)],
    ["optimism_97.5pct", q(0.975).toFixed(4)],
    ["optimism_corrected_AUC", aucCorr.toFixed(4)],
    ["valid_bootstrap_n", String(valid)],
    ["skipped_bootstrap_n", String(skipped)],
    ["mean_selected_features", meanSel.toFixed(2)],
    ["LRT_X2", lrtX2.toFixed(4)],
    ["LRT_df", String(lrtDf)],
    ["LRT_p", lrtP.toPrecision(6)],
    ["HL_X2", hl.X2.toFixed(4)],
    ["HL_df", String(hl.df)],
    ["HL_p", hl.p.toFixed(6)],
    ["apparent_features", full.selNames.join(";")],
    ["lambda_1se", full.lambda.toPrecision(6)],
  ];
  fs.writeFileSync(path.join(outDir, "full_pipeline_validation_summary.csv"),
    summary.map(r => r.map(v => /[",]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : v).join(",")).join("\n"));
  fs.writeFileSync(path.join(outDir, "feature_selection_frequency.csv"),
    ["feature,count,frequency"].concat(freqArr.map(([nm, c]) => `${nm},${c},${(c / valid).toFixed(3)}`)).join("\n"));
  fs.writeFileSync(path.join(outDir, "results.json"), JSON.stringify({
    apparent_auc: aucApp, optimism_mean: meanO, optimism_ci: [q(0.025), q(0.975)],
    corrected_auc: aucCorr, valid: valid, skipped: skipped, skip_reasons: reasons,
    apparent_features: full.selNames, lambda_1se: full.lambda,
    lrt: { X2: lrtX2, df: lrtDf, p: lrtP }, hl,
    freq: freqArr, optimism_values: optimism, runtime_sec: (Date.now() - t0) / 1000,
  }, null, 2));
  console.log(`\n结果已保存至 ${outDir}`);
}
main().catch(e => { console.error(e); process.exit(1); });
