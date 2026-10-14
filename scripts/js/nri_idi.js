// 影像模型 ± 基线WBC 的 NRI/IDI bootstrap 敏感性分析
// 输入: baseline_60.json (基线) + input_data_v2.json (影像特征) + 60例y
// 模型A: 仅影像(DepVar+MeshVolume), IRLS+高斯消元
// 模型B: 影像 + WBC (按Firth式小样本惩罚, 用L2=0.1近似)
// 输出: NRI_continuous, IDI, ΔAUC, bootstrap 95%CI

const fs = require('fs');

const v2 = JSON.parse(fs.readFileSync('validation_output/input_data_v2.json', 'utf8'));
const bl = JSON.parse(fs.readFileSync('validation_output/baseline_60.json', 'utf8'));
const pts = v2.patients, X = v2.X, y = v2.y;
const fn = v2.feature_names;
const iDV = fn.indexOf('DependenceVariance');
const iMV = fn.indexOf('MeshVolume');

function sigmoid(z){ return 1/(1+Math.exp(-z)); }

// 标准化
function stdCols(arr){
  const n = arr.length;
  const mu = arr.reduce((a,b)=>a+b,0)/n;
  const sd = Math.sqrt(arr.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n);
  return { mu, sd, z: arr.map(v=>(v-mu)/(sd||1)) };
}

const zDV = stdCols(pts.map((_,i)=>X[i][iDV]));
const zMV = stdCols(pts.map((_,i)=>X[i][iMV]));

// 筛选有基线WBC且 gap ≤ 14d (窗内) 的样本
// 字段: [WBC, Hb, BMI, gap_days, RT, baseline_date, surgery_date, note]
const validIdx14 = [];
for (let i=0;i<60;i++) {
  const r = bl[pts[i]];
  if (r && r[0] != null && r[3] != null && r[3] <= 14) validIdx14.push(i);
}
const validIdxAny = [];
for (let i=0;i<60;i++) {
  const r = bl[pts[i]];
  if (r && r[0] != null) validIdxAny.push(i);
}
console.log(`基线WBC窗内(gap≤14): ${validIdx14.length} 例, 事件 ${validIdx14.filter(i=>y[i]==1).length}`);
console.log(`基线WBC任何有效: ${validIdxAny.length} 例, 事件 ${validIdxAny.filter(i=>y[i]==1).length}`);

// 改用窗内 36 例作为主要分析
const validIdx = validIdx14;
const sub = validIdx.map(i=>({
  y: y[i], zDV: zDV.z[i], zMV: zMV.z[i],
  wbc: bl[pts[i]][0]
}));
// WBC标准化
const wbcArr = sub.map(s=>s.wbc);
const wbcStats = stdCols(wbcArr);
sub.forEach(s=>s.zWBC = wbcStats.z[sub.indexOf(s)]);

function fitLogit(design, yy, lambda=0){
  // design: n×p 矩阵(含截距列), yy: 0/1
  const n = yy.length;
  const p = design[0].length;
  let beta = new Array(p).fill(0);
  for (let it=0; it<120; it++) {
    const W = [], pr = [];
    for (let i=0; i<n; i++) {
      let z = 0; for (let j=0; j<p; j++) z += beta[j]*design[i][j];
      const s = sigmoid(z); pr.push(s);
      W.push(Math.max(s*(1-s), 1e-8));
    }
    const A = Array.from({length:p},()=>new Array(p).fill(0));
    const g = new Array(p).fill(0);
    for (let i=0; i<n; i++) {
      const r = yy[i]-pr[i];
      for (let a=0; a<p; a++) {
        g[a] += design[i][a]*r;
        for (let b=0; b<p; b++) A[a][b] += design[i][a]*design[i][b]*W[i];
      }
    }
    // L2惩罚(对非截距)
    for (let j=1; j<p; j++) A[j][j] += lambda;
    // 解 A*step = g
    const M = A.map((row,i)=>row.concat([g[i]]));
    for (let c=0;c<p;c++){
      let piv=c;
      for (let r=c+1;r<p;r++) if (Math.abs(M[r][c])>Math.abs(M[piv][c])) piv=r;
      [M[c],M[piv]]=[M[piv],M[c]];
      const d=M[c][c]||1e-12;
      for (let j=c;j<=p;j++) M[c][j]/=d;
      for (let r=0;r<p;r++){ if(r!==c&&M[r][c]!==0){const f=M[r][c]; for(let j=c;j<=p;j++) M[r][j]-=f*M[c][j];} }
    }
    const step = M.map(row=>row[p]);
    if (Math.max(...step.map(Math.abs)) < 1e-8) break;
    for (let j=0;j<p;j++) beta[j] += step[j];
  }
  return beta;
}

function predict(design, beta){
  return design.map(row => {
    let z=0; for (let j=0;j<beta.length;j++) z += beta[j]*row[j];
    return sigmoid(z);
  });
}

function aucOf(sc, yy){
  const c1 = sc.map((s,i)=>({s,y:yy[i]})).filter(o=>o.y===1).map(o=>o.s);
  const c0 = sc.map((s,i)=>({s,y:yy[i]})).filter(o=>o.y===0).map(o=>o.s);
  let U=0; for (const a of c1) for (const b of c0) if (a>b) U+=1; else if (a===b) U+=0.5;
  return U/(c1.length*c0.length);
}

// === 拟合基线影像模型 M_imaging ===
const dA = sub.map(s=>[1, s.zDV, s.zMV]);
const betaA = fitLogit(dA, sub.map(s=>s.y), 0.001);  // 微弱L2稳定
const pA = predict(dA, betaA);
const aucA = aucOf(pA, sub.map(s=>s.y));
console.log(`M_imaging (影像 only):  β=[${betaA.map(b=>b.toFixed(3)).join(', ')}],  AUC=${aucA.toFixed(4)}`);

// === 拟合影像+WBC 模型 M_imaging_WBC (含WBC主效应, 不含交互以保EPV) ===
const dB = sub.map(s=>[1, s.zDV, s.zMV, s.zWBC]);
const betaB = fitLogit(dB, sub.map(s=>s.y), 0.001);
const pB = predict(dB, betaB);
const aucB = aucOf(pB, sub.map(s=>s.y));
console.log(`M_imaging_WBC:           β=[${betaB.map(b=>b.toFixed(3)).join(', ')}],  AUC=${aucB.toFixed(4)}`);

// === NRI_continuous 和 IDI (apparent) ===
function nriIdi(pOld, pNew, yy){
  const idx = [...Array(yy.length).keys()];
  const ev = idx.filter(i=>yy[i]===1);
  const ne = idx.filter(i=>yy[i]===0);
  // NRI_cont: 对事件, p_new-p_old>0 为"上移"=+1; 对非事件, <0 为"下移"=+1
  let up_ev=0, dn_ne=0;
  for (const i of ev) if (pNew[i] > pOld[i]) up_ev++;
  for (const i of ne) if (pNew[i] < pOld[i]) dn_ne++;
  const nri = up_ev/ev.length - (ne.length - dn_ne)/ne.length;
  // IDI = mean(pNew - pOld)
  const idi = (pNew.reduce((a,b)=>a+b,0) - pOld.reduce((a,b)=>a+b,0)) / yy.length;
  return { nri, idi, up_ev, dn_ne, n_ev:ev.length, n_ne:ne.length };
}
const m1 = nriIdi(pA, pB, sub.map(s=>s.y));
console.log(`NRI_continuous = ${m1.nri.toFixed(3)} (事件上移${m1.up_ev}/${m1.n_ev} - 非事件下移${m1.dn_ne}/${m1.n_ne})`);
console.log(`IDI = ${m1.idi.toFixed(4)} (${(m1.idi*100).toFixed(2)}%净概率提升)`);

// === Bootstrap 95%CI (B=1000) ===
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
const rng = mulberry32(20261009);
const B = 1000;
const nriBoot = [], idiBoot = [], dAucBoot = [];
for (let b=0; b<B; b++) {
  // 分层bootstrap with-replacement
  const ev = sub.map((s,i)=>({s,i})).filter(o=>o.s.y===1);
  const ne = sub.map((s,i)=>({s,i})).filter(o=>o.s.y===0);
  function sampleWR(arr) {
    return Array.from({length: arr.length}, () => arr[Math.floor(rng()*arr.length)].i);
  }
  const evB = sampleWR(ev);
  const neB = sampleWR(ne);
  const boot = [...evB, ...neB];
  const ss = boot.map(i => sub[i]);
  const ddA = ss.map(s=>[1, s.zDV, s.zMV]);
  const ddB = ss.map(s=>[1, s.zDV, s.zMV, s.zWBC]);
  const yyB = ss.map(s=>s.y);
  try {
    const bA = fitLogit(ddA, yyB, 0.001);
    const bB = fitLogit(ddB, yyB, 0.001);
    const pAb = predict(ddA, bA);
    const pBb = predict(ddB, bB);
    const mi = nriIdi(pAb, pBb, yyB);
    nriBoot.push(mi.nri);
    idiBoot.push(mi.idi);
    dAucBoot.push(aucOf(pBb, yyB) - aucOf(pAb, yyB));
  } catch(e) { /* 跳过拟合失败 */ }
}
function ci(arr){
  const s = arr.slice().sort((a,b)=>a-b);
  return [s[Math.floor(s.length*0.025)], s[Math.floor(s.length*0.975)]];
}
const nriCI = ci(nriBoot);
const idiCI = ci(idiBoot);
const dAucCI = ci(dAucBoot);
const pct = (x) => `${(x*100).toFixed(1)}%`;
console.log(`\n--- Bootstrap (B=${nriBoot.length}, 95%CI) ---`);
console.log(`NRI_continuous: ${m1.nri.toFixed(3)} [${nriCI[0].toFixed(3)}, ${nriCI[1].toFixed(3)}]`);
console.log(`IDI:            ${m1.idi.toFixed(4)} [${idiCI[0].toFixed(4)}, ${idiCI[1].toFixed(4)}]`);
console.log(`ΔAUC:           ${(aucB-aucA).toFixed(4)} [${dAucCI[0].toFixed(4)}, ${dAucCI[1].toFixed(4)}]`);

// === 决策曲线分析(三条): treat-all, treat-none, 影像模型, 影像+WBC ===
// net benefit = TP/N - FP/N * (pt/(1-pt))
function dca(probs, yy, thresholds){
  const N = yy.length;
  return thresholds.map(pt => {
    let tp=0, fp=0;
    for (let i=0; i<N; i++) {
      const pred = probs[i] >= pt ? 1 : 0;
      if (pred===1 && yy[i]===1) tp++;
      if (pred===1 && yy[i]===0) fp++;
    }
    return tp/N - fp/N * (pt/(1-pt));
  });
}
const thr = Array.from({length:41}, (_,k)=>0.05 + k*0.01);
const nbAll = thr.map(pt => { const p=yy_avg(sub.map(s=>s.y)); return p - (1-p)*(pt/(1-pt)); });
function yy_avg(arr){ return arr.reduce((a,b)=>a+b,0)/arr.length; }
const pBase = yy_avg(sub.map(s=>s.y));
const nbAll2 = thr.map(pt => pBase - (1-pBase)*(pt/(1-pt)));
const nbA = dca(pA, sub.map(s=>s.y), thr);
const nbB = dca(pB, sub.map(s=>s.y), thr);

// 输出关键数据
const out = {
  n_valid: sub.length, n_event: sub.filter(s=>s.y===1).length, n_noevent: sub.filter(s=>s.y===0).length,
  model_A_imaging: { beta: betaA, auc: aucA, n: sub.length },
  model_B_imaging_WBC: { beta: betaB, auc: aucB, n: sub.length },
  apparent: { nri: m1.nri, idi: m1.idi, delta_auc: aucB-aucA,
              nri_events_up: m1.up_ev, nri_nonevents_down: m1.dn_ne },
  bootstrap: { B: nriBoot.length,
               nri: m1.nri, nri_ci: nriCI,
               idi: m1.idi, idi_ci: idiCI,
               dAuc: aucB-aucA, dAuc_ci: dAucCI },
  dca: { thresholds: thr, treat_all: nbAll2, model_A: nbA, model_B: nbB }
};
fs.writeFileSync('validation_output/nri_idi_result.json', JSON.stringify(out, null, 2));
console.log(`\n保存: validation_output/nri_idi_result.json`);
console.log(`\n决策曲线净获益@阈值15%: 影像=${nbA[10].toFixed(3)} | 影像+WBC=${nbB[10].toFixed(3)} | 全治=${nbAll2[10].toFixed(3)} | 全不治=0`);
