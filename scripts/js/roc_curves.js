// roc_curves.js v3 - 4-model ROC curves, unified 60-case cohort (P059/P060 included)
const fs = require('fs');

function sigmoid(z){ return 1/(1+Math.exp(-z)); }
Math.erfc = function(x){
  const t = 1/(1+0.5*Math.abs(x));
  const tau = t*Math.exp(-x*x - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
    t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*0.17087277)))))))));
  return x >= 0 ? tau : 2 - tau;
};
function solveLin(A, b){
  const n = A.length;
  for (let i=0;i<n;i++) A[i][i] += 1e-9;
  const M = A.map((r,i)=>[...r, b[i]]);
  for (let c=0;c<n;c++){
    let p = c;
    for (let r=c+1;r<n;r++) if (Math.abs(M[r][c])>Math.abs(M[p][c])) p = r;
    [M[c],M[p]] = [M[p],M[c]];
    const d = M[c][c] || 1e-12;
    for (let j=c;j<=n;j++) M[c][j] /= d;
    for (let r=0;r<n;r++) if (r!==c && M[r][c]!==0){
      const f = M[r][c];
      for (let j=c;j<=n;j++) M[r][j] -= f*M[c][j];
    }
  }
  return M.map(r=>r[n]);
}
function irlsLogistic(X, y){
  const n = X.length, p = X[0].length;
  let beta = new Array(p).fill(0);
  for (let it=0; it<200; it++){
    const W = X.map(r=>{ const z=r.reduce((s,v,j)=>s+v*beta[j],0); return {w: Math.max(sigmoid(z)*(1-sigmoid(z)),1e-8), pr: sigmoid(z)}; });
    const A = Array.from({length:p},()=>new Array(p).fill(0));
    const g = new Array(p).fill(0);
    for (let i=0;i<n;i++){
      const gr = y[i]-W[i].pr;
      for (let a=0;a<p;a++){ g[a] += X[i][a]*gr; for (let b=0;b<p;b++) A[a][b] += X[i][a]*X[i][b]*W[i].w; }
    }
    const step = solveLin(A, g);
    beta = beta.map((v,j)=>v+step[j]);
    if (Math.max(...step.map(Math.abs)) < 1e-9) break;
    if (Math.max(...beta.map(Math.abs)) > 30) return null;
  }
  return beta;
}
function aucOf(score, y){
  const c1=[], c0=[];
  score.forEach((s,i)=> (y[i]?c1:c0).push(s));
  let U=0, T=0;
  for (const a of c1) for (const b of c0){ if (a>b) U++; else if (a===b) T++; }
  return (U+0.5*T)/(c1.length*c0.length);
}

const data = JSON.parse(fs.readFileSync('validation_output/input_data_v2.json', 'utf-8'));
const merged = JSON.parse(fs.readFileSync('validation_output/merged_60.json', 'utf-8'));
const clMap = new Map(merged.map(p=>[p.name, p]));

const f = name => data.feature_names.indexOf(name);
const fId = {DV: f('DependenceVariance'), MV: f('MeshVolume')};
console.log('Feature indices:', fId);

const order = data.patients;
const idxAll = order.map((_,i)=>i);
const y = data.y;
const X = data.X;
const N_all = y.length, E_all = y.filter(v=>v===1).length;
console.log(`全样本: N=${N_all}, E=${E_all}`);

// 4 模型统一 60 例
const dose = JSON.parse(fs.readFileSync('validation_output/dosi_input.json', 'utf-8'));
const doseMap = new Map();
for (const grp of ['event', 'nonevent']){
  for (const d of dose[grp].data) doseMap.set(d.name, d);
}

const models = {};

// M1: 临床 (60例) - 年龄+宫颈癌
const Xm1 = idxAll.map(i => {
  const cl = clMap.get(order[i]);
  const ages = idxAll.map(j => clMap.get(order[j]).age);
  const mu = ages.reduce((s,v)=>s+v,0)/ages.length;
  const sd = Math.sqrt(ages.reduce((s,v)=>s+(v-mu)*(v-mu),0)/(ages.length-1)) || 1;
  return [1, (cl.age-mu)/sd, cl.cervical ? 1 : 0];
});
let b1 = irlsLogistic(Xm1, y);
if (!b1) b1 = [0,0,0];
const score1 = Xm1.map(r=>sigmoid(r.reduce((s,v,j)=>s+v*b1[j],0)));
models.M1 = {beta: b1, score: score1, y: y, n: N_all, e: E_all, auc: aucOf(score1, y), name: 'M1 Clinical (Age + Cervical)'};
console.log(`M1 n=${N_all} β=${b1.map(v=>v.toFixed(3))} AUC=${models.M1.auc.toFixed(4)}`);

// M2: 剂量 (60例) - Volume + V20
const Xm2 = idxAll.map(i => {
  const d = doseMap.get(order[i]);
  return [1, d ? d['vals'][0] : 0, d ? d['vals'][6] : 0];
});
let b2 = irlsLogistic(Xm2, y);
if (!b2) b2 = [0,0,0];
const score2 = Xm2.map(r=>sigmoid(r.reduce((s,v,j)=>s+v*b2[j],0)));
models.M2 = {beta: b2, score: score2, y: y, n: N_all, e: E_all, auc: aucOf(score2, y), name: 'M2 Dosimetric (Volume + V20)'};
console.log(`M2 n=${N_all} β=${b2.map(v=>v.toFixed(3))} AUC=${models.M2.auc.toFixed(4)}`);

// M3: 影像 (60例) - DV + MV
const Xm3 = X.map(r => [1, r[fId.DV], r[fId.MV]]);
let b3 = irlsLogistic(Xm3, y);
if (!b3) b3 = [0,0,0];
const score3 = Xm3.map(r=>sigmoid(r.reduce((s,v,j)=>s+v*b3[j],0)));
models.M3 = {beta: b3, score: score3, y: y, n: N_all, e: E_all, auc: aucOf(score3, y), name: 'M3 Radiomics (DV + MV)'};
console.log(`M3 n=${N_all} β=${b3.map(v=>v.toFixed(3))} AUC=${models.M3.auc.toFixed(4)}`);

// M4: 联合 (60例) - 临床+剂量+影像
const Xm4 = idxAll.map(i => {
  const cl = clMap.get(order[i]);
  const r = X[i];
  const d = doseMap.get(order[i]);
  const ages = idxAll.map(j => clMap.get(order[j]).age);
  const ageMu = ages.reduce((s,v)=>s+v,0)/ages.length;
  const ageSd = Math.sqrt(ages.reduce((s,v)=>s+(v-ageMu)*(v-ageMu),0)/(ages.length-1)) || 1;
  return [1, (cl.age-ageMu)/ageSd, cl.cervical?1:0, d?d['vals'][0]:0, d?d['vals'][6]:0, r[fId.DV], r[fId.MV]];
});
let b4 = irlsLogistic(Xm4, y);
if (!b4) b4 = new Array(Xm4[0].length).fill(0);
const score4 = Xm4.map(r=>sigmoid(r.reduce((s,v,j)=>s+v*b4[j],0)));
models.M4 = {beta: b4, score: score4, y: y, n: N_all, e: E_all, auc: aucOf(score4, y), name: 'M4 Combined (M1+M2+M3)'};
console.log(`M4 n=${N_all} β=${b4.map(v=>v.toFixed(3))} AUC=${models.M4.auc.toFixed(4)}`);

// === ROC 坐标 ===
function rocCurve(score, y){
  const sorted = score.map((s,i)=>({s, y: y[i]})).sort((a,b)=>b.s - a.s);
  let tp = 0, fp = 0;
  const totalPos = y.filter(v=>v===1).length;
  const totalNeg = y.length - totalPos;
  const points = [{fpr:0, tpr:0}];
  let prevScore = Infinity;
  for (const p of sorted){
    if (p.s !== prevScore){
      points.push({fpr: fp/totalNeg, tpr: tp/totalPos});
      prevScore = p.s;
    }
    if (p.y === 1) tp++; else fp++;
  }
  points.push({fpr:1, tpr:1});
  return points;
}

// === 校准 (M3 60 例) ===
function calibration(score, y, nBins=10){
  const sorted = score.map((s,i)=>({s, y: y[i]})).sort((a,b)=>a.s - b.s);
  const bins = [];
  for (let i=0; i<nBins; i++){
    const lo = Math.floor(i * score.length / nBins);
    const hi = Math.floor((i+1) * score.length / nBins);
    const sub = sorted.slice(lo, hi);
    if (sub.length === 0) continue;
    const meanPred = sub.reduce((s,p)=>s+p.s, 0) / sub.length;
    const obsRate = sub.reduce((s,p)=>s+p.y, 0) / sub.length;
    bins.push({meanPred, obsRate, n: sub.length});
  }
  return bins;
}

const rocOut = {};
for (const [k, m] of Object.entries(models)){
  rocOut[k] = {name: m.name, auc: m.auc, n: m.n, e: m.e, beta: m.beta, points: rocCurve(m.score, m.y)};
}

const calM3 = calibration(score3, y, 10);

const out = {
  n_all: N_all, e_all: E_all,
  models: rocOut,
  calibration_M3: calM3
};
fs.writeFileSync('validation_output/roc_curves.json', JSON.stringify(out));
console.log('Wrote validation_output/roc_curves.json');
