# ============================================================
# ⚠️  DE-IDENTIFIED VERSION (PHI removed for public release)
# Original 60 patient names replaced with P001–P060.
# See: 60 patients enrolled 2019–2025; 15 Grade≥2 leukopenia events.
# De-identification map is NOT included in this repository.
# ============================================================
"""生成60例四模态完整主表: 临床+剂量+影像+基线(含WBC/Hb/BMI), 整合NRI/IDI结果, 加敏感性分析 sheet"""
import json, openpyxl, math
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

# 加载全部数据
v2 = json.load(open("validation_output/input_data_v2.json"))
fn, pts, X, y = v2["feature_names"], v2["patients"], v2["X"], v2["y"]
dosi = json.load(open("validation_output/dosi/dosi_merged_v2_60.json"))
dose_by = {d["name"]: d for d in dosi["data"]}
bl = json.load(open("validation_output/baseline_60.json"))
nri = json.load(open("validation_output/nri_idi_result.json"))

# 临床表(本批次)
clin = {
"P001":("直肠癌",66,"F",3,0,1,45),"P002":("子宫内膜癌",70,"F",2,0,1,45),
"P003":("子宫内膜癌",74,"F",1,0,1,45),"P004":("宫颈癌",46,"F",1,0,1,45),
"P005":("子宫内膜癌",51,"F",1,0,1,45),"P006":("宫颈癌",59,"F",3,0,1,45),
"P007":("宫颈癌",70,"F",3,0,1,45),"P008":("宫颈癌",48,"F",1,0,1,45),
"P009":("宫颈癌",42,"F",1,0,1,45),"P010":("宫颈癌",59,"F",1,0,1,45),
"P011":("宫颈癌",68,"F",2,0,1,45),"P012":("宫颈癌",66,"F",1,0,1,45),
"P013":("宫颈癌",56,"F",1,0,1,45),"P014":("直肠癌",70,"F",3,0,1,45),
"P015":("直肠癌",60,"M",2,0,1,45),"P016":("直肠癌",68,"F",2,0,1,45),
"P017":("宫颈癌",45,"F",2,0,1,45),"P018":("宫颈癌",54,"F",3,0,1,45),
"P019":("宫颈癌",46,"F",2,0,1,45),"P020":("宫颈癌",54,"F",3,0,1,45),
"P021":("宫颈癌",53,"F",1,0,1,45),"P022":("宫颈癌",65,"F",1,0,1,45),
"P023":("宫颈癌",72,"F",1,1,1,45),"P024":("宫颈癌",57,"F",1,0,1,45),
"P025":("子宫内膜癌",57,"F",3,0,1,45),"P026":("宫颈癌",52,"F",1,0,1,45),
"P027":("宫颈癌",69,"F",1,0,1,45),"P028":("宫颈癌",54,"F",1,0,1,45),
"P029":("宫颈癌",39,"F",2,0,1,45),"P030":("宫颈癌",61,"F",1,0,1,45),
"P031":("宫颈癌",76,"F",3,0,1,45),"P032":("子宫内膜癌",75,"F",1,0,1,45),
"P033":("子宫内膜癌",68,"F",1,0,1,45),"P034":("子宫内膜癌",62,"F",1,0,1,45),
"P035":("宫颈癌",53,"F",1,0,1,45),"P036":("宫颈癌",59,"F",1,0,1,45),
"P037":("直肠癌",60,"M",3,0,1,45),"P038":("宫颈癌",58,"F",1,0,1,45),
"P039":("直肠癌",72,"M",2,0,1,45),"P040":("子宫内膜癌",63,"F",1,0,1,45),
"P041":("子宫内膜癌",61,"F",1,0,1,45),"P042":("直肠癌",66,"M",3,0,1,45),
"P043":("直肠癌",64,"M",3,0,1,45),"P044":("子宫内膜癌",65,"F",1,0,1,45),
"P045":("直肠癌",41,"M",2,0,1,45),"P046":("直肠癌",78,"M",3,0,1,45),
"P047":("子宫内膜癌",45,"F",3,0,1,45),"P048":("直肠癌",83,"M",3,0,1,45),
"P049":("子宫内膜癌",71,"F",1,0,1,45),"P050":("直肠癌",61,"M",3,0,1,45),
"P051":("子宫内膜癌",63,"F",1,0,1,45),"P052":("直肠癌",72,"M",4,0,1,45),
"P053":("直肠癌",54,"F",4,0,1,45),"P054":("子宫内膜癌",74,"F",1,0,1,45),
"P055":("宫颈癌",41,"F",3,0,1,45),"P056":("宫颈癌",63,"F",1,0,1,45),
"P057":("宫颈癌",72,"F",1,0,1,45),"P058":("宫颈癌",48,"F",2,0,1,45),
"P059":("宫颈癌",45,"F",1,0,1,45),"P060":("宫颈癌",72,"F",1,0,1,45)
}
stage_cn = {1:"Ⅰ",2:"Ⅱ",3:"Ⅲ",4:"Ⅳ"}

# 影像预测分数
iDV = fn.index("DependenceVariance"); iMV = fn.index("MeshVolume")
def stdC(arr):
    n=len(arr); mu=sum(arr)/n
    sd=math.sqrt(sum((v-mu)**2 for v in arr)/n)
    return [(v-mu)/(sd or 1) for v in arr]
zDV = stdC([X[i][iDV] for i in range(60)])
zMV = stdC([X[i][iMV] for i in range(60)])
# 影像模型系数(从nri_idi.js的B=12~48子集复制: β=[-1.511, -0.658, -1.270], 但那是子集; 60例全量不同; 重新用全60例拟合)
def sigmoid(z): return 1/(1+math.exp(-z))
n=60; p=3
b = [0.0]*p
for it in range(200):
    pr=[]; W=[]
    for i in range(n):
        z=b[0]+zDV[i]*b[1]+zMV[i]*b[2]
        s=sigmoid(z); pr.append(s); W.append(max(s*(1-s),1e-8))
    A=[[0.0]*p for _ in range(p)]; g=[0.0]*p
    for i in range(n):
        r=y[i]-pr[i]
        for a in range(p):
            g[a]+=r*(1 if a==0 else (zDV[i] if a==1 else zMV[i]))
            for c in range(p):
                v = (1 if a==0 else (zDV[i] if a==1 else zMV[i])) * (1 if c==0 else (zDV[i] if c==1 else zMV[i]))
                A[a][c]+=v*W[i]
    for j in range(p): A[j][j]+=1e-9
    M=[A[i][:]+[g[i]] for i in range(p)]
    for c in range(p):
        piv=c
        for r in range(c+1,p):
            if abs(M[r][c])>abs(M[piv][c]): piv=r
        M[c],M[piv]=M[piv],M[c]
        d=M[c][c] or 1e-12
        for j in range(c,p+1): M[c][j]/=d
        for r in range(p):
            if r!=c and M[r][c]!=0:
                f=M[r][c]
                for j in range(c,p+1): M[r][j]-=f*M[c][j]
    step=[M[i][p] for i in range(p)]
    b=[b[j]+step[j] for j in range(p)]
    if max(abs(s) for s in step)<1e-9: break
print(f"影像模型全60例系数: β0={b[0]:.3f}, βDV={b[1]:.3f}, βMV={b[2]:.3f}")
lp_imaging = [b[0]+zDV[i]*b[1]+zMV[i]*b[2] for i in range(60)]

wb = openpyxl.Workbook()
ws = wb.active; ws.title = "60例四模态主表"
hdr = ["姓名","Group","瘤种","年龄","性别","分期","医嘱(Gy/F)",
       "Volume(cm3)","mean dose(cGy)","V5(%)","V10(%)","V20(%)","V35(%)","V40(%)",
       "DepVar","MeshVolume","影像模型lp",
       "基线WBC","基线Hb","BMI","基线淋巴","采血→放疗|天数|","备注"]
ws.append(hdr)
hf = PatternFill("solid", fgColor="1F4E79")
ev_fill = PatternFill("solid", fgColor="FCE4EC")
new_fill = PatternFill("solid", fgColor="FFF2CC")
thin = Side(style="thin", color="BFBFBF")
border = Border(left=thin,right=thin,top=thin,bottom=thin)
for c in ws[1]:
    c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

for i, pt in enumerate(pts):
    ca, age, sex, st, cc, _, rx = clin[pt]
    d = dose_by[pt]["dosi"]
    b_ = bl.get(pt)
    if b_:
        wbc, hb, bmi, lym, _, _, _, note = b_
        days = None
        try:
            from datetime import datetime
            if b_[4] and b_[5]:
                days = abs((datetime.strptime(b_[5],"%Y-%m-%d") - datetime.strptime(b_[4],"%Y-%m-%d")).days)
        except: pass
        note_full = note or ""
        if days is not None and days > 14: note_full = (note_full + " | 过窗>14天").strip(" |")
    else:
        wbc = hb = bmi = lym = days = None
        note_full = "基线全缺"
    row = [pt, y[i], ca, age, sex, stage_cn.get(st,""), f"{rx}Gy/25F",
           d[0], d[3], d[4], d[5], d[6], d[7], d[8],
           round(X[i][iDV], 4), round(X[i][iMV], 1), round(lp_imaging[i], 3),
           wbc, hb, bmi, lym, days, note_full]
    ws.append(row)
    r = ws.max_row
    for c in ws[r]:
        c.border = border
        c.alignment = Alignment(horizontal="center" if c.column <= 21 else "left")
    if y[i] == 1:
        for c in ws[r]: c.fill = ev_fill
    if pt in ("P059","P060"):
        for c in ws[r]: c.fill = new_fill if y[i]==0 else PatternFill("solid", fgColor="F4B183")
ws.freeze_panes = "A2"
widths = [10,6,11,7,5,6,9,11,12,8,8,8,8,8,11,11,11,9,9,8,9,11,28]
for i,w in enumerate(widths,1): ws.column_dimensions[get_column_letter(i)].width = w

# ============== Sheet 2: 基线组间比较结果 ==============
ws2 = wb.create_sheet("基线变量组间比较")
ws2.append(["基线变量(基线WBC完整者n=48, 事件12/无事件36)","事件组(n=12)","无事件组(n=36)","Wilcoxon p","单变量反向AUC","结论"])
for c in ws2[1]: c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
data2 = [
 ["基线WBC (×10⁹/L) 中位数(IQR)","4.03 (3.70-6.14)","5.45 (4.16-6.54)","0.194","0.627","弱信号,方向一致"],
 ["基线Hb (g/L) 中位数(IQR)","120.5 (94.8-127.5)","120 (114.8-127.3)","0.551","0.559","几乎无差异"],
 ["BMI (kg/m²) 中位数(IQR)","22.04 (21.28-23.54)","22.58 (20.62-24.91)","0.905","0.487","无信号"],
]
for r in data2: ws2.append(r)
ws2.append([])
ws2.append(["时间窗核查(n=48双日期齐全)","","","","",""])
ws2.append(["中位数 |天数|","5天","","","",""])
ws2.append(["≤7天内","31(64.6%)","","","",""])
ws2.append(["8-14天","5(10.4%)","","","",""])
ws2.append([">14天(过窗)","12(25.0%) — 其中3例事件","","","",""])
ws2.append([])
ws2.append(["数据完整度","基线WBC/Hb/BMI/淋巴 48/60 (80.0%) | 放疗日期 56/60 (93.3%) | 手术日期 53/60 (88.3%)","","","",""])
ws2.append([])
ws2.append(["12例过窗患者(均为采血早于放疗,反映术前/新辅助状态)","","","","",""])
ws2.append(["P005(事件)|P018(事件)|P034(事件)|P019|P020|P024|P028|P029|P032|P033|P035|P059","","","","",""])
ws2.append([])
ws2.append(["特殊备注患者","#21P021 尿路感染 | #31P031 根治性放疗 | #43P043、#48P048 新辅助","","","",""])
ws2.column_dimensions["A"].width = 38
for col in "BCDEF": ws2.column_dimensions[col].width = 18

# ============== Sheet 3: 影像模型 ± WBC 敏感性(NRI/IDI/DCA) ==============
ws3 = wb.create_sheet("影像vs影像+WBC敏感性")
ws3.append(["指标","M_imaging(影像 only)","M_imaging+WBC","差值(M+B-A)","95%CI(B=1000分层bootstrap)"])
for c in ws3[1]: c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
ws3.append(["样本量","48(基线WBC完整者)","48","",""])
ws3.append(["事件数","12","12","",""])
ws3.append(["表观AUC",f"{nri['model_A_imaging']['auc']:.4f}",f"{nri['model_B_imaging_WBC']['auc']:.4f}",f"{nri['apparent']['delta_auc']:.4f}",
            f"[{nri['bootstrap']['dAuc_ci'][0]:.4f}, {nri['bootstrap']['dAuc_ci'][1]:.4f}]"])
ws3.append([])
ws3.append(["NRI_continuous","—","—",f"{nri['apparent']['nri']:.3f}",
            f"[{nri['bootstrap']['nri_ci'][0]:.3f}, {nri['bootstrap']['nri_ci'][1]:.3f}]"])
ws3.append(["  事件上移","—","—",f"{nri['apparent']['nri_events_up']}/12",""])
ws3.append(["  非事件下移","—","—",f"{nri['apparent']['nri_nonevents_down']}/36",""])
ws3.append([])
ws3.append(["IDI","—","—",f"{nri['apparent']['idi']:.4f}",
            f"[{nri['bootstrap']['idi_ci'][0]:.4f}, {nri['bootstrap']['idi_ci'][1]:.4f}]"])
ws3.append([])
ws3.append(["决策曲线净获益 @ 阈值概率15%","0.174","0.178","+0.004",""])
ws3.append([])
ws3.append(["结论","基线WBC在影像模型之上提供有限的判别改进(ΔAUC +0.002, NRI 0.167 但95%CI跨0); IDI≈0提示平均预测概率无提升",
            "","",""])
ws3.append(["解读","符合EPV=5约束下的预期: 影像模型已通过骨盆体积捕获骨髓储备信号,与基线WBC存在共线性,WBC不增加新信息",
            "","",""])
ws3.append(["建议","主结果仍以影像2特征模型为主; 基线WBC可写入'基线临床特征'描述性表(已是Table 1候选),不进主模型",
            "","",""])
for col in "ABCDE": ws3.column_dimensions[col].width = 26
ws3.column_dimensions["E"].width = 32

# ============== Sheet 4: 完整四模型对照(主结果) ==============
ws4 = wb.create_sheet("主结果_四模型对照")
ws4.append(["模型","变量","n","表观AUC","DeLong对比","ΔAUC","p"])
for c in ws4[1]: c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
ws4.append(["M1 临床(全60)","年龄、宫颈癌",60,"0.611","vs影像","+0.182","0.089"])
ws4.append(["M2 剂量(全60)","Volume、V5",60,"0.785","vs影像","-0.003","0.876"])
ws4.append(["M3 影像(全60)","DepVar、MeshVolume(LASSO λ.1se)",60,"0.793","(参考)","—","—"])
ws4.append(["M4 联合(全60)","Volume+影像lp",60,"0.788","vs影像","-0.005","0.835"])
ws4.append([])
ws4.append(["B=500全流程Bootstrap乐观校正(影像模型)","60例","0.793","校正后","0.695",""])
ws4.append([])
ws4.append(["特征稳定性(500次重抽样,有效185次)","MeshVolume 48.6% | DepVar 37.8% | DependenceEntropy 32.4%","","","",""])
ws4.append(["整体模型LRT","X²=13.49, df=2, p=0.0012","","","",""])
ws4.append(["Hosmer-Lemeshow","p=0.535(校准良好)","","","",""])
for col in "ABCDEFG": ws4.column_dimensions[col].width = 20

# ============== Sheet 5: 数据治理与队列说明 ==============
ws5 = wb.create_sheet("数据治理与队列说明")
notes = [
["参数口径(决定性)","全部60例统一使用汇总表口径(1409特征,含图像类型列); data.xlsx口径作废"],
["重名特征定位","Contrast→glcm(r=0.61) | GLNU→glrlm(0.98) | GLNUN→glszm(0.77) | GLV→glszm(0.67)"],
["入排","剔除:P061(50Gy/25F,事件)/P062(50Gy/25F);纳入:P059/P060(45Gy/25F);P023保留(术后单次顺铂忽略)"],
["结局窗口","放疗开始至结束后三个月内"],
["处方","统一45Gy/25F(60例); 2例50Gy/25F已剔除"],
["基线血","48/60(80.0%)完整, 12例缺失"],
["过窗处理","12例采血早于放疗>14天(过窗),其中3例事件; 主分析仍用全48例, 严格基线≤7d敏感性:31例/9事件"],
["统计方法","Wilcoxon秩和(tie+连续性校正); AUC基于平均秩; LASSO λ.1se (10-fold CV); 全流程Bootstrap乐观校正; 分层bootstrap NRI/IDI; DeLong配对AUC比较"],
["EPV约束","15事件→主模型≤3特征(实际2); 基线WBC不直接进主模型, 仅作NRI/IDI敏感性"],
]
for r in notes: ws5.append(r)
ws5.column_dimensions["A"].width = 22; ws5.column_dimensions["B"].width = 90

wb.save("../results/tables/60case_4model_master_DEID.xlsx")  # DEID
print("已保存: validation_output/60例四模态主表_含基线.xlsx")
