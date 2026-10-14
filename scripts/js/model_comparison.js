// 四模型对照: 临床 vs 剂量 vs 影像组学 vs 联合 (60例, 15事件)
// 用法: node model_comparison.js
const fs = require('fs');

// ---------- 数学工具 ----------
function sigmoid(z){ return 1/(1+Math.exp(-z)); }
function normSf(x){ // 双侧p
  const z = Math.abs(x)/Math.SQRT2;
  return Math.erfc(z);
}
Math.erfc = function(x){
  // Lentz连分式或级数; 用互补误差函数数值实现
  const t = 1/(1+0.5*Math.abs(x));
  const tau = t*Math.exp(-x*x - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
    t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*0.17087277)))))))));
  return x >= 0 ? tau : 2 - tau;
};

function solveLin(A, b){ // 高斯消元 + 微小岭
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

function irlsLogistic(X, y){ // X: n×p (已含截距列), 返回 beta 或 null
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
  const n1 = y.reduce((s,v)=>s+v,0), n0 = y.length-n1;
  let c1=[],c0=[];
  score.forEach((s,i)=> (y[i]?c1:c0).push(s));
  let U=0, T=0;
  for (const a of c1) for (const b of c0){ if (a>b) U++; else if (a===b) T++; }
  return (U+0.5*T)/(n1*n0);
}

// DeLong: 返回 {auc, V10[], V01[]} 供配对检验
function delongParts(score, y){
  const idx1=[], idx0=[];
  y.forEach((v,i)=> (v?idx1:idx0).push(i));
  const n1=idx1.length, n0=idx0.length;
  const auc = aucOf(score,y);
  const V10 = idx1.map(i=>{
    let s=0; for (const j of idx0){ if (score[i]>score[j]) s++; else if (score[i]===score[j]) s+=0.5; }
    return s/n0 - auc;
  });
  const V01 = idx0.map(j=>{
    let s=0; for (const i of idx1){ if (score[i]>score[j]) s++; else if (score[i]===score[j]) s+=0.5; }
    return s/n1 - auc;
  });
  return {auc, V10, V01, n1, n0};
}
function delongPaired(m1, m2){ // 同一批患者的两个模型
  const n1=m1.n1, n0=m1.n0;
  const var1 = m1.V10.reduce((s,v)=>s+v*v,0)/(n1*(n1-1)) + m1.V01.reduce((s,v)=>s+v*v,0)/(n0*(n0-1));
  const var2 = m2.V10.reduce((s,v)=>s+v*v,0)/(n1*(n1-1)) + m2.V01.reduce((s,v)=>s+v*v,0)/(n0*(n0-1));
  let cov = 0; for (let i=0;i<n1;i++) cov += m1.V10[i]*m2.V10[i]; cov /= n1*(n1-1);
  let cov0 = 0; for (let j=0;j<n0;j++) cov0 += m1.V01[j]*m2.V01[j]; cov0 /= n0*(n0-1);
  const v = var1 + var2 - 2*(cov+cov0);
  const z = (m1.auc-m2.auc)/Math.sqrt(Math.max(v,1e-15));
  return {delta: m1.auc-m2.auc, z, p: normSf(z)};
}

// ---------- 数据 ----------
const merged = JSON.parse(fs.readFileSync('validation_output/merged_60.json','utf8'));
const inp = JSON.parse(fs.readFileSync('validation_output/input_data.json','utf8'));
const fn = inp.feature_names, Xrad = inp.X;
const iB = fn.indexOf('Busyness'), iC = fn.indexOf('Contrast'), iS = fn.indexOf('SurfaceVolumeRatio');
const byName = {};
inp.patients.forEach((n,i)=> byName[n]=i);

const n = merged.length;
const y = merged.map(m=>m.y);

function zscore(col){
  const mu = col.reduce((a,b)=>a+b,0)/col.length;
  const sd = Math.sqrt(col.reduce((s,v)=>s+(v-mu)**2,0)/col.length);
  return col.map(v=>(v-mu)/sd);
}
const z = {
  age: zscore(merged.map(m=>m.age)),
  sex: merged.map(m=>m.sex_F),
  cer: merged.map(m=>m.cervical),
  vol: zscore(merged.map(m=>m.vol)),
  V5:  zscore(merged.map(m=>m.V5)),
};
// 影像lp: 3特征 IRLS (标准化)
const zb = zscore(merged.map(m=>Xrad[byName[m.name]][iB]));
const zc = zscore(merged.map(m=>Xrad[byName[m.name]][iC]));
const zs = zscore(merged.map(m=>Xrad[byName[m.name]][iS]));
const X3 = merged.map((_,i)=>[1, zb[i], zc[i], zs[i]]);
const beta3 = irlsLogistic(X3, y);
const radlp = merged.map((_,i)=> X3[i].reduce((s,v,j)=>s+v*beta3[j],0));
const zrad = zscore(radlp);

// ---------- 模型 ----------
function fitModel(name, cols){
  const X = merged.map((_,i)=>[1, ...cols.map(c=>c[i])]);
  const beta = irlsLogistic(X, y);
  if (!beta){ console.log(`${name}: IRLS失败`); return null; }
  const score = X.map(r=>r.reduce((s,v,j)=>s+v*beta[j],0));
  const dp = delongParts(score, y);
  console.log(`${name}: AUC = ${dp.auc.toFixed(4)}  β = [${beta.map(b=>b.toFixed(3)).join(', ')}]`);
  return {name, dp, beta};
}
console.log('===== 四模型(表观AUC, 60例/15事件) =====');
const M_clin = fitModel('M1 临床(年龄,性别,宫颈)', [z.age, z.sex, z.cer]);
const M_dose = fitModel('M2 剂量(Volume,V5)    ', [z.vol, z.V5]);
const M_rad  = {name:'M3 影像(3特征glm)      ', dp: delongParts(radlp, y)};
console.log(`${M_rad.name}: AUC = ${M_rad.dp.auc.toFixed(4)}`);
const M_comb = fitModel('M4 联合(Volume+影像lp)', [z.vol, zrad]);

console.log('\n===== DeLong 配对检验 =====');
const pairs = [
  ['M3影像 vs M1临床', M_rad, M_clin],
  ['M3影像 vs M2剂量', M_rad, M_dose],
  ['M4联合 vs M1临床', M_comb, M_clin],
  ['M4联合 vs M2剂量', M_comb, M_dose],
  ['M4联合 vs M3影像', M_comb, M_rad],
];
for (const [lab,a,b] of pairs){
  if (!a||!b) continue;
  const r = delongPaired(a.dp, b.dp);
  console.log(`${lab}: ΔAUC = ${r.delta>=0?'+':''}${r.delta.toFixed(4)}, z = ${r.z.toFixed(3)}, p = ${r.p<0.001?'<0.001':r.p.toFixed(3)}`);
}

// ---------- 单变量 AUC ----------
console.log('\n===== 单变量AUC =====');
for (const [lab, arr] of [['Volume', merged.map(m=>-m.vol)], ['V5', merged.map(m=>m.V5)], ['影像lp', radlp], ['年龄', merged.map(m=>-m.age)], ['性别(女)', merged.map(m=>m.sex_F)], ['宫颈癌', merged.map(m=>m.cervical)]]){
  console.log(`  ${lab}: AUC = ${aucOf(arr, y).toFixed(4)}`);
}

// 保存结果
fs.writeFileSync('validation_output/model_comparison_results.json', JSON.stringify({
  M_clin: {auc:M_clin.auc||M_clin.dp.auc}, n:60, events:15,
  models: {clin: M_clin&&M_clin.dp.auc, dose: M_dose&&M_dose.dp.auc, rad: M_rad.dp.auc, comb: M_comb&&M_comb.dp.auc},
  delong: pairs.filter(([,a,b])=>a&&b).map(([lab,a,b])=>({lab, ...delongPaired(a.dp,b.dp)}))
}, null, 2));
