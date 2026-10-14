// ============================================================
// ⚠️  DE-IDENTIFIED VERSION (PHI removed for public release)
// Original patient names replaced with P059/P060 codes.
// De-identification map is NOT included in this repository.
// ============================================================

// 新队列(汇总表口径, 60例/15事件)四模型对照 v2
const fs = require('fs');
function sigmoid(z){ return 1/(1+Math.exp(-z)); }
Math.erfc = function(x){
  const t = 1/(1+0.5*Math.abs(x));
  const tau = t*Math.exp(-x*x - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
    t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*0.17087277)))))))));
  return x >= 0 ? tau : 2 - tau;
};
function normSf(x){ return Math.erfc(Math.abs(x)/Math.SQRT2); }
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
  for (let it=0; it<100; it++){
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
function delongParts(score, y){
  const idx1=[], idx0=[];
  y.forEach((v,i)=> (v?idx1:idx0).push(i));
  const n1=idx1.length, n0=idx0.length;
  const auc = aucOf(score,y);
  const V10 = idx1.map(i=>{ let s=0; for (const j of idx0){ if (score[i]>score[j]) s++; else if (score[i]===score[j]) s+=0.5; } return s/n0-auc; });
  const V01 = idx0.map(j=>{ let s=0; for (const i of idx1){ if (score[i]>score[j]) s++; else if (score[i]===score[j]) s+=0.5; } return s/n1-auc; });
  return {auc, V10, V01, n1, n0};
}
function delongPaired(m1, m2){
  const n1=m1.n1, n0=m1.n0;
  const var1 = m1.V10.reduce((s,v)=>s+v*v,0)/(n1*(n1-1)) + m1.V01.reduce((s,v)=>s+v*v,0)/(n0*(n0-1));
  const var2 = m2.V10.reduce((s,v)=>s+v*v,0)/(n1*(n1-1)) + m2.V01.reduce((s,v)=>s+v*v,0)/(n0*(n0-1));
  let cov=0; for (let i=0;i<n1;i++) cov += m1.V10[i]*m2.V10[i]; cov /= n1*(n1-1);
  let cov0=0; for (let j=0;j<n0;j++) cov0 += m1.V01[j]*m2.V01[j]; cov0 /= n0*(n0-1);
  const z = (m1.auc-m2.auc)/Math.sqrt(Math.max(var1+var2-2*(cov+cov0),1e-15));
  return {delta: m1.auc-m2.auc, z, p: normSf(z)};
}
function zscore(col){
  const mu = col.reduce((a,b)=>a+b,0)/col.length;
  const sd = Math.sqrt(col.reduce((s,v)=>s+(v-mu)**2,0)/col.length);
  return col.map(v=>(v-mu)/sd);
}

// ---------- 数据 ----------
const v2 = JSON.parse(fs.readFileSync('validation_output/input_data_v2.json','utf8'));
const dosiOld = JSON.parse(fs.readFileSync('validation_output/dosi/dosi_merged_v2_60.json','utf8'));
const doseBy = {}; dosiOld.data.forEach(d=>doseBy[d.name]=d.dosi);
const clinAdd = {"P059":{age:45,sex_F:1,cervical:1}, "P060":{age:72,sex_F:1,cervical:1}};
const oldClin = JSON.parse(fs.readFileSync('validation_output/merged_60.json','utf8'));
const clinBy = {}; oldClin.forEach(m=>clinBy[m.name]=m);

const pts = v2.patients, y = v2.y, X = v2.X, fn = v2.feature_names;
const iDV = fn.indexOf('DependenceVariance'), iMV = fn.indexOf('MeshVolume');
const age = pts.map(p=>(clinBy[p]||clinAdd[p]).age);
const cer = pts.map(p=>(clinBy[p]||clinAdd[p]).cervical);
const d58 = pts.map((_,i)=>i).filter(i=>pts[i] in doseBy);

// M3 影像lp: DependenceVariance + MeshVolume (全60例, 标准化 IRLS)
const zDV = zscore(X.map(r=>r[iDV]));
const zMV = zscore(X.map(r=>r[iMV]));
const Xr = pts.map((_,i)=>[1, zDV[i], zMV[i]]);
const betaR = irlsLogistic(Xr, y);
const radlpAll = Xr.map(r=>r.reduce((s,v,j)=>s+v*betaR[j],0));

// 剂量变量(58例内标准化)
const dv58 = d58.map(i=>doseBy[pts[i]][0]);
const v558 = d58.map(i=>doseBy[pts[i]][4]);
const zv = zscore(dv58), z5 = zscore(v558);
const lp58 = d58.map(i=>radlpAll[i]);
const zlp58 = zscore(lp58);
const y58 = d58.map(i=>y[i]);
const zAge58 = zscore(d58.map(i=>age[i]));
const cer58 = d58.map(i=>cer[i]);

function fitModel(name, cols, yy){
  const Xf = yy.map((_,i)=>[1, ...cols.map(c=>c[i])]);
  const beta = irlsLogistic(Xf, yy);
  if (!beta){ console.log(name+': IRLS失败'); return null; }
  const score = Xf.map(r=>r.reduce((s,v,j)=>s+v*beta[j],0));
  console.log(`${name}: AUC = ${aucOf(score, yy).toFixed(4)}`);
  return delongParts(score, yy);
}
console.log('===== 新队列四模型(表观AUC, 60例/15事件) =====');
const M1 = fitModel('M1 临床(年龄,宫颈癌) n=60', [zscore(age), cer], y);
const M3 = delongParts(radlpAll, y);
console.log(`M3 影像(DepVar+MeshVol) n=60: AUC = ${M3.auc.toFixed(4)}`);
const M2 = fitModel('M2 剂量(Volume,V5) n=58 ', [zv, z5], y58);
const M4 = fitModel('M4 联合(Volume+影像lp) n=58', [zv, zlp58], y58);
// M3 在58例子集上的版本(供与M2/M4配对)
const M3_58 = delongParts(lp58, y58);
const M1_58 = (()=>{ const Xf = y58.map((_,i)=>[1, zAge58[i], cer58[i]]); const b = irlsLogistic(Xf, y58);
  const sc = Xf.map(r=>r.reduce((s,v,j)=>s+v*b[j],0)); return delongParts(sc, y58); })();

console.log('\n===== DeLong 配对检验 =====');
const pairs = [
  ['M3影像 vs M1临床 (n=60)', M3, M1],
  ['M3影像 vs M2剂量 (n=58)', M3_58, M2],
  ['M4联合 vs M1临床 (n=58)', M4, M1_58],
  ['M4联合 vs M2剂量 (n=58)', M4, M2],
  ['M4联合 vs M3影像 (n=58)', M4, M3_58],
];
const out = [];
for (const [lab,a,b] of pairs){
  if (!a||!b) continue;
  const r = delongPaired(a, b);
  console.log(`${lab}: ΔAUC = ${r.delta>=0?'+':''}${r.delta.toFixed(4)}, z = ${r.z.toFixed(3)}, p = ${r.p<0.001?'<0.001':r.p.toFixed(3)}`);
  out.push({lab, ...r});
}
console.log('\n===== 单变量AUC =====');
for (const [lab, arr, yy] of [['MeshVolume(反向,60)', X.map(r=>-r[iMV]), y], ['DependenceVariance(60)', X.map(r=>r[iDV]), y],
  ['剂量Volume(反向,58)', dv58.map(v=>-v), y58], ['V5(58)', v558, y58], ['影像lp(60)', radlpAll, y]]){
  console.log(`  ${lab}: AUC = ${aucOf(arr, yy).toFixed(4)}`);
}
fs.writeFileSync('../data/model_comparison_v2_results.json'  // DEID, JSON.stringify({
  M1: M1.auc, M2: M2.auc, M3: M3.auc, M4: M4.auc, delong: out
}, null, 2));
