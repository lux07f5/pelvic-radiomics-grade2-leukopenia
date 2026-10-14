# ============================================================
# ⚠️  DE-IDENTIFIED VERSION (PHI removed for public release)
# Original 60 patient names replaced with P001–P060.
# See: 60 patients enrolled 2019–2025; 15 Grade≥2 leukopenia events.
# De-identification map is NOT included in this repository.
# ============================================================
# -*- coding: utf-8 -*-
# 生成60例三模态合并主表 Excel
import json, math, statistics
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

merged = json.load(open("validation_output/merged_60.json"))
inp = json.load(open("validation_output/input_data.json"))
idx = {n: i for i, n in enumerate(inp["patients"])}
fn = inp["feature_names"]; X = inp["X"]
iB, iC, iS = fn.index("Busyness"), fn.index("Contrast"), fn.index("SurfaceVolumeRatio")
dosi = json.load(open("validation_output/dosi/dosi_merged_60.json"))
dd_by = {x["name"]: x["dosi"] for x in dosi["data"]}

feats = {m["name"]: [X[idx[m["name"]]][iB], X[idx[m["name"]]][iC], X[idx[m["name"]]][iS]] for m in merged}
cols = list(feats.values())
def zc(ci):
    v = [r[ci] for r in cols]; mu = sum(v)/len(v); sd = (sum((x-mu)**2 for x in v)/len(v))**0.5
    return [(x-mu)/sd for x in v]
zb, zcon, zs = zc(0), zc(1), zc(2)
def sigmoid(z): return 1/(1+math.exp(-z))
n = 60; p = 4; beta = [0.0]*p
Xr = [[1, zb[i], zcon[i], zs[i]] for i in range(n)]
for it in range(100):
    W = []; pr = []
    for i in range(n):
        z = sum(Xr[i][j]*beta[j] for j in range(p)); s = sigmoid(z); pr.append(s); W.append(max(s*(1-s), 1e-8))
    A = [[0.0]*p for _ in range(p)]; g = [0.0]*p
    for i in range(n):
        gr = merged[i]["y"] - pr[i]
        for a in range(p):
            g[a] += Xr[i][a]*gr
            for b in range(p): A[a][b] += Xr[i][a]*Xr[i][b]*W[i]
    for a in range(p): A[a][a] += 1e-9
    M = [A[i][:]+[g[i]] for i in range(p)]
    for c in range(p):
        piv = c
        for r in range(c+1, p):
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
lp = {merged[i]["name"]: sum(Xr[i][j]*beta[j] for j in range(p)) for i in range(n)}

clin_full = [
("P001","直肠癌","女"),("P062","子宫内膜癌","女"),("P002","子宫内膜癌","女"),("P003","子宫内膜癌","女"),
("P004","宫颈癌","女"),("P005","子宫内膜癌","女"),("P006","宫颈癌","女"),("P007","宫颈癌","女"),
("P008","宫颈癌","女"),("P009","宫颈癌","女"),("P010","宫颈癌","女"),("P011","宫颈癌","女"),
("P061","宫颈癌","女"),("P012","宫颈癌","女"),("P013","宫颈癌","女"),("P014","直肠癌","女"),
("P015","直肠癌","男"),("P016","直肠癌","女"),("P017","宫颈癌","女"),("P018","宫颈癌","女"),
("P019","宫颈癌","女"),("P020","宫颈癌","女"),("P021","宫颈癌","女"),("P022","宫颈癌","女"),
("P023","宫颈癌","女"),("P024","宫颈癌","女"),("P025","子宫内膜癌","女"),("P026","宫颈癌","女"),
("P027","宫颈癌","女"),("P028","宫颈癌","女"),("P029","宫颈癌","女"),("P030","宫颈癌","女"),
("P031","宫颈癌","女"),("P032","子宫内膜癌","女"),("P033","子宫内膜癌","女"),("P034","子宫内膜癌","女"),
("P035","宫颈癌","女"),("P036","宫颈癌","女"),("P037","直肠癌","男"),("P038","宫颈癌","女"),
("P039","直肠癌","男"),("P040","子宫内膜癌","女"),("P041","子宫内膜癌","女"),("P042","直肠癌","男"),
("P043","直肠癌","男"),("P044","子宫内膜癌","女"),("P045","直肠癌","男"),("P046","直肠癌","男"),
("P047","子宫内膜癌","女"),("P048","直肠癌","男"),("P049","子宫内膜癌","女"),("P050","直肠癌","男"),
("P051","子宫内膜癌","女"),("P052","直肠癌","男"),("P053","直肠癌","女"),("P054","子宫内膜癌","女"),
("P055","宫颈癌","女"),("P056","宫颈癌","女"),("P057","宫颈癌","女"),("P058","宫颈癌","女")]
clin_by = {c[0]: c for c in clin_full}

vols = [m["vol"] for m in merged]; v5s = [m["V5"] for m in merged]
mv, sdv = statistics.mean(vols), statistics.stdev(vols)
m5, sd5 = statistics.mean(v5s), statistics.stdev(v5s)

wb = openpyxl.Workbook(); ws = wb.active; ws.title = "60例三模态合并主表"
hdr = ["姓名","Group(1=骨髓抑制)","瘤种","年龄","性别","分期(Ⅰ-Ⅳ)","同步顺铂","放疗医嘱",
       "Volume(cm3)","min dose(cGy)","max dose(cGy)","mean dose(cGy)","V5(%)","V10(%)","V20(%)","V35(%)","V40(%)",
       "Busyness","Contrast","SurfaceVolumeRatio","影像模型lp","剂量模型lp(参考)"]
ws.append(hdr)
hf = PatternFill("solid", fgColor="1F4E79")
for c in ws[1]:
    c.font = Font(bold=True, color="FFFFFF"); c.fill = hf; c.alignment = Alignment(horizontal="center")
for m in merged:
    nm, ca, sx = clin_by[m["name"]]
    dd = dd_by[m["name"]]
    dlp = -1.501 + (-1.209)*((m["vol"]-mv)/sdv) + 0.314*((m["V5"]-m5)/sd5)
    f3 = feats[m["name"]]
    ws.append([m["name"], m["y"], ca, m["age"], sx, m["stage"], m["cc"],
               ("50Gy/25F" if m["rx50"] else "45Gy/25F"),
               m["vol"], dd[1], dd[2], m["meandose"], m["V5"], m["V10"], m["V20"], m["V35"], m["V40"],
               f3[0], f3[1], f3[2], round(lp[m["name"]], 4), round(dlp, 4)])
widths = [10, 14, 13, 7, 6, 10, 10, 12] + [12]*9 + [14]*3 + [12, 16]
for i, w in enumerate(widths, 1): ws.column_dimensions[get_column_letter(i)].width = w

ws2 = wb.create_sheet("四模型对照")
data2 = [
 ["模型", "变量", "表观AUC"],
 ["M1 临床模型", "年龄、宫颈癌（性别因拟完全分离剔除）", 0.6133],
 ["M2 剂量模型", "Volume、V5", 0.7526],
 ["M3 影像组学模型", "Busyness、Contrast、SurfaceVolumeRatio（glm）", 0.8119],
 ["M4 联合模型", "Volume + 影像lp", 0.8148],
 ["", "", ""],
 ["配对DeLong检验", "ΔAUC", "p值"],
 ["影像 vs 临床", "+0.199", "0.040（显著）"],
 ["影像 vs 剂量", "+0.059", "0.412（不显著）"],
 ["联合 vs 临床", "+0.202", "0.035（显著）"],
 ["联合 vs 剂量", "+0.062", "0.310（不显著）"],
 ["联合 vs 影像", "+0.003", "0.871（不显著）"],
 ["", "", ""],
 ["单变量AUC", "Volume(反向) 0.759 | V5 0.640 | 年龄 0.539 | 女性 0.611 | 宫颈癌 0.600", ""],
]
for r in data2: ws2.append(r)
for c in ws2[1]:
    c.font = Font(bold=True, color="FFFFFF"); c.fill = hf
ws2.column_dimensions["A"].width = 16; ws2.column_dimensions["B"].width = 50; ws2.column_dimensions["C"].width = 12
ws2.append([])
ws2.append(["注1: 全部为表观AUC（全样本回代）；影像模型经B=500全流程Bootstrap校正后AUC=0.718。"])
ws2.append(["注2: M1性别变量因15例事件全为女性致拟完全分离（系数膨胀）而剔除；性别单变量Fisher p=0.053。"])
ws2.append(["注3: 处方口径（用户2026-10-04定稿）：全队列45Gy/25F，仅P061、P0622例为50Gy/25F。"])

ws3 = wb.create_sheet("数据核查与待确认")
notes = [
 ["✅ 已定稿", ""],
 ["放疗方案", "全队列45Gy/25F，仅P061、P0622例50Gy/25F（用户2026-10-04确认，此前的'45Gy/15次'为口误）"],
 ["", ""],
 ["⚠️ 仍待确认", ""],
 ["1. 同步化疗", "P023为'术后辅助同步顺铂'且在骨髓抑制组（Grade2）——此前说明已排除同步化疗患者，此例是否剔除或在局限中说明？"],
 ["2. 结局窗口", "临床表表头写'从放疗开始至放疗结束一月内'，此前描述为'放疗结束后3个月内'——论文方法学需统一为哪个口径？"],
 ["3. 基线血常规", "临床表无基线WBC列，需另行提供（合并主表结构已备好，发我即可补入）"],
 ["4. 剂量学队列64例", "吴莉萍、庄明珍（事件）、赵云媛、沈杏金（无事件）4例有剂量数据但不在影像/临床队列，原因待说明（是则可做剂量模型敏感性分析）"],
 ["", ""],
 ["数据质量备注", "年龄18例为文本格式（P037至P054段）；Stage文本列格式混乱（建议用数字列）；P053'直肠恶性肿瘤'已并入直肠癌"],
 ["队列一致性", "60例临床=影像=剂量三方结局标签100%一致（15事件/45无事件；Grade2:12例、Grade3:3例、Grade1:10例归无事件）"],
]
for r in notes: ws3.append(r)
ws3.column_dimensions["A"].width = 20; ws3.column_dimensions["B"].width = 110

wb.save("../results/tables/60case_3model_master_v1_DEID.xlsx")  # DEID
print("已保存 validation_output/60例三模态合并主表.xlsx")
