# ============================================================
# ⚠️  DE-IDENTIFIED VERSION (PHI removed for public release)
# Original 60 patient names replaced with P001–P060.
# See: 60 patients enrolled 2019–2025; 15 Grade≥2 leukopenia events.
# De-identification map is NOT included in this repository.
# ============================================================
# -*- coding: utf-8 -*-
# 新队列(58旧+P059/P060)三模态合并主表 v2
import json, math
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

v2 = json.load(open("validation_output/input_data_v2.json"))
pts, y2, X2, fn = v2["patients"], v2["y"], v2["X"], v2["feature_names"]
dosiOld = json.load(open("validation_output/dosi/dosi_merged_v2_60.json"))  # 60例完整口径(含P059/P060)
doseBy = {d["name"]: d["dosi"] for d in dosiOld["data"]}
oldClin = json.load(open("validation_output/merged_60.json"))
clinBy = {m["name"]: m for m in oldClin}
# 新患者临床（去标识化后 ID）: P059 宫颈癌45岁女IB3 术后辅助45Gy/25F 无事件; P060 宫颈癌72岁女IB2 术后辅助45Gy/25F 事件(WBC 2级)
newClin = {"P059": {"ca": "宫颈癌", "age": 45, "sex": "女", "stage": 1},
           "P060": {"ca": "宫颈癌", "age": 72, "sex": "女", "stage": 1}}
iDV, iMV = fn.index("DependenceVariance"), fn.index("MeshVolume")

# 影像lp (标准化 IRLS, 与 model_comparison_v2.js 一致)
def zc(col):
    mu = sum(col)/len(col); sd = (sum((x-mu)**2 for x in col)/len(col))**0.5
    return [(x-mu)/sd for x in col]
zDV, zMV = zc([r[iDV] for r in X2]), zc([r[iMV] for r in X2])
def sigmoid(z): return 1/(1+math.exp(-z))
n, p = 60, 3; beta = [0.0]*p
Xr = [[1, zDV[i], zMV[i]] for i in range(n)]
for it in range(100):
    W = []; pr = []
    for i in range(n):
        z = sum(Xr[i][j]*beta[j] for j in range(p)); s = sigmoid(z); pr.append(s); W.append(max(s*(1-s), 1e-8))
    A = [[0.0]*p for _ in range(p)]; g = [0.0]*p
    for i in range(n):
        gr = y2[i]-pr[i]
        for a in range(p):
            g[a] += Xr[i][a]*gr
            for b in range(p): A[a][b] += Xr[i][a]*Xr[i][b]*W[i]
    for a in range(p): A[a][a] += 1e-9
    M = [A[i][:]+[g[i]] for i in range(p)]
    for c in range(p):
        piv = c
        for r in range(p):
            if abs(M[r][c]) > abs(M[piv][c]): piv = r
        M[c], M[piv] = M[piv], M[c]
        d = M[c][c] or 1e-12
        for j in range(c, p+1): M[c][j] /= d
        for r in range(p):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                for j in range(c, p+1): M[r][j] -= f*M[c][j]
    step = [M[i][p] for i in range(p)]
    beta = [beta[j]+step[j] for j in range(p)]
    if max(abs(s) for s in step) < 1e-9: break
lp = {pts[i]: sum(Xr[i][j]*beta[j] for j in range(p)) for i in range(n)}

wb = openpyxl.Workbook(); ws = wb.active; ws.title = "新队列60例主表"
hf = PatternFill("solid", fgColor="1F4E79")
hdr = ["姓名", "Group(1=骨髓抑制)", "瘤种", "年龄", "性别", "分期(Ⅰ-Ⅳ)", "放疗医嘱",
       "Volume(cm3)", "min dose(cGy)", "max dose(cGy)", "mean dose(cGy)", "V5(%)", "V10(%)", "V20(%)", "V35(%)", "V40(%)",
       "DependenceVariance", "MeshVolume", "影像模型lp"]
ws.append(hdr)
for c in ws[1]:
    c.font = Font(bold=True, color="FFFFFF"); c.fill = hf; c.alignment = Alignment(horizontal="center")
for i, pt in enumerate(pts):
    c_ = clinBy.get(pt)
    if c_:
        ca = {"cervical": "宫颈癌", "endometrial": "子宫内膜癌", "rectal": "直肠癌"}[
            "cervical" if pt in clinBy and clinBy[pt]["cervical"] == 1 else ("endometrial" if False else "rectal")]
        # 由 merged_60 的 cervical 字段判断
        ca = "宫颈癌" if c_["cervical"] == 1 else ("子宫内膜癌" if c_["stage"] else "")
        age, sex, stage, rx = c_["age"], "女" if c_["sex_F"] else "男", c_["stage"], "45Gy/25F"
    else:
        nc = newClin[pt]
        ca, age, sex, stage, rx = nc["ca"], nc["age"], nc["sex"], nc["stage"], "45Gy/25F"
    d = doseBy.get(pt)
    row = [pt, y2[i], ca, age, sex, stage, rx,
           d[0] if d else None, d[1] if d else None, d[2] if d else None, d[3] if d else None,
           d[4] if d else None, d[5] if d else None, d[6] if d else None, d[7] if d else None, d[8] if d else None,
           X2[i][iDV], X2[i][iMV], round(lp[pt], 4)]
    ws.append(row)
# 瘤种列修正: 用诊断映射
ca_map = {}
for m in oldClin:
    ca_map[m["name"]] = None
# 从临床原始信息重建瘤种
diag = {"P001":"直肠癌","P002":"子宫内膜癌","P003":"子宫内膜癌","P004":"宫颈癌","P005":"子宫内膜癌",
"P006":"宫颈癌","P007":"宫颈癌","P008":"宫颈癌","P009":"宫颈癌","P010":"宫颈癌","P011":"宫颈癌",
"P012":"宫颈癌","P013":"宫颈癌","P014":"直肠癌","P015":"直肠癌","P016":"直肠癌","P017":"宫颈癌",
"P018":"宫颈癌","P019":"宫颈癌","P020":"宫颈癌","P021":"宫颈癌","P022":"宫颈癌","P023":"宫颈癌",
"P024":"宫颈癌","P025":"子宫内膜癌","P026":"宫颈癌","P027":"宫颈癌","P028":"宫颈癌","P029":"宫颈癌",
"P030":"宫颈癌","P031":"宫颈癌","P032":"子宫内膜癌","P033":"子宫内膜癌","P034":"子宫内膜癌",
"P035":"宫颈癌","P036":"宫颈癌","P037":"直肠癌","P038":"宫颈癌","P039":"直肠癌","P040":"子宫内膜癌",
"P041":"子宫内膜癌","P042":"直肠癌","P043":"直肠癌","P044":"子宫内膜癌","P045":"直肠癌","P046":"直肠癌",
"P047":"子宫内膜癌","P048":"直肠癌","P049":"子宫内膜癌","P050":"直肠癌","P051":"子宫内膜癌",
"P052":"直肠癌","P053":"直肠癌","P054":"子宫内膜癌","P055":"宫颈癌","P056":"宫颈癌","P057":"宫颈癌",
"P058":"宫颈癌","P062":"子宫内膜癌","P061":"宫颈癌"}
diag["P059"] = "宫颈癌"; diag["P060"] = "宫颈癌"
for r in ws.iter_rows(min_row=2):
    r[2].value = diag.get(r[0].value, r[2].value)
widths = [10, 14, 13, 7, 6, 10, 11] + [12]*9 + [16, 14, 12]
for i, w in enumerate(widths, 1): ws.column_dimensions[get_column_letter(i)].width = w
# 高亮新患者
new_fill = PatternFill("solid", fgColor="FFF2CC")
for r in ws.iter_rows(min_row=2):
    if r[0].value in ("P059", "P060"):
        for c in r: c.fill = new_fill

ws2 = wb.create_sheet("模型结果")
data2 = [
 ["新队列(60例/15事件, 统一45Gy/25F, 汇总表参数口径)", "", ""],
 ["", "", ""],
 ["模型", "变量", "表观AUC"],
 ["M1 临床 (n=60)", "年龄、宫颈癌", 0.6111],
 ["M2 剂量 (n=60)", "Volume、V5", 0.7852],
 ["M3 影像组学 (n=60)", "DependenceVariance、MeshVolume (LASSO λ.1se)", 0.7926],
 ["M4 联合 (n=60)", "Volume + 影像lp", 0.7881],
 ["", "", ""],
 ["DeLong", "ΔAUC", "p"],
 ["影像 vs 临床 (n=60)", "+0.182", "0.089"],
 ["影像 vs 剂量 (n=60)", "+0.007", "0.876"],
 ["联合 vs 临床 (n=60)", "+0.177", "0.109"],
 ["联合 vs 影像 (n=60)", "-0.004", "0.835"],
 ["", "", ""],
 ["单变量", "MeshVolume(反向) 0.776 | 剂量Volume(反向) 0.784 | V5 0.668 | DependenceVariance 0.271(反向0.729)", ""],
 ["", "", ""],
 ["全样本检验", "", ""],
 ["Wilcoxon筛选候选", "17个 (p<0.05)", ""],
 ["λ.1se", 0.1085, ""],
 ["LRT", "X2=13.49, df=2, p=0.0012", ""],
 ["Hosmer-Lemeshow", "X2=7.01, df=8, p=0.535", ""],
]
for r in data2: ws2.append(r)
for c in ws2[3]: c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
for c in ws2[8]: c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
ws2.column_dimensions["A"].width = 20; ws2.column_dimensions["B"].width = 48; ws2.column_dimensions["C"].width = 12
ws2.append([])
ws2.append(["B=500全流程Bootstrap校正结果见 v2_final 目录(乐观校正后AUC待填)。"])

ws3 = wb.create_sheet("队列变更说明")
notes = [
 ["本次队列调整(2026-10-05)", ""],
 ["剔除", "P061(事件,50Gy/25F)、P062(无事件,50Gy/25F) —— 全队列处方统一为45Gy/25F"],
 ["纳入", "P059(宫颈癌,45岁,女,IB3,术后辅助,45Gy/25F,无事件);P060(宫颈癌,72岁,女,IB2,术后辅助,45Gy/25F,事件:白细胞2级)"],
 ["P023", "保留(术后单次同步顺铂按忽略处理)"],
 ["结局窗口", "放疗开始至结束后三个月内"],
 ["", ""],
 ["参数口径(重要)", "全部101特征值统一取自汇总表提取口径(与新患者文件逐行对齐验证,1409行0错配;新患者对汇总表口径z诊断中位|z|=0.56-0.59,无一特征|z|>3)"],
 ["", "data.xlsx口径与新患者不兼容(46/101与39/101特征|z|>3)——已弃用,论文方法学需按汇总表口径描述特征提取参数"],
 ["重名特征定位", "Contrast→glcm(r=0.61); GrayLevelNonUniformity→glrlm(0.98); GLNUN→glszm(0.77); GrayLevelVariance→glszm(0.67)"],
 ["", ""],
 ["与旧结果差异", "旧模型(Busyness/Contrast/SVR, AUC 0.812)基于data.xlsx口径;新口径下LASSO入选DependenceVariance+MeshVolume,表观AUC 0.793——两套参数的纹理信号不同,新口径是内部一致的可发表版本"],
 ["剂量数据", "P059、P060剂量学数据已并入(2026-10-05用户发送), 剂量分析恢复60例完整口径: Volume、min/max dose、V5、V10 p<0.001, V20/V35不显著——结论与58例口径一致"],
]
for r in notes: ws3.append(r)
ws3.column_dimensions["A"].width = 20; ws3.column_dimensions["B"].width = 110

wb.save("../results/tables/60case_3model_master_v2_DEID.xlsx")  # DEID
print("已保存 validation_output/新队列60例三模态主表.xlsx")
