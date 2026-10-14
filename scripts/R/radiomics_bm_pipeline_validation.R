# =====================================================================
# 全流程内部验证：单因素筛选 -> LASSO -> 多因素Logistic（Bootstrap 乐观校正）
# ---------------------------------------------------------------------
# 项目   ：定位CT影像组学预测盆腔放疗后 Grade>=2 白细胞减少（骨髓抑制）
# 数据   ：data.xlsx（第1列 Group=结局0/1；其余列为 PyRadiomics 影像特征）
# 用途   ：把"特征筛选 + LASSO选lambda + 定特征 + 拟模型"整条管线装进
#          Bootstrap 重抽样，得到乐观校正后的 AUC —— 替代只对固定3特征
#          公式做校正的 rms::calibrate，回答审稿人"是否全体回代"的质疑
# 依赖包 ：readxl, glmnet, pROC, rms（缺失则 install.packages("包名")）
# 运行   ：Rscript radiomics_bm_pipeline_validation.R
# 说明   ：① 与原分析一致：单因素 p<0.05 筛候选 -> cv.glmnet 取 lambda.1se
#          ② 每次重抽样内【完整重跑上述管线】，乐观度=样本内AUC-原样本AUC
#          ③ 同时输出每次抽样"选中了哪些特征"的频次表，检验3特征是否稳定
# =====================================================================

suppressMessages({
  library(readxl); library(glmnet); library(pROC); library(rms)
})

# ---------------- 0. 参数 ----------------
# Set data_path to point to your DE-IDENTIFIED feature matrix.
# The original analysis used a per-center file; here we expect a generic
# (n × p) feature matrix with a "Group" column (0/1 outcome).
data_path <- "./data/radiomics_features_deid.xlsx"  # <-- ADJUST TO YOUR PATH
out_dir   <- dirname(data_path)                     # Results saved alongside data
B         <- 500                                    # Bootstrap replicates
p_screen  <- 0.05                                   # Univariate screening threshold
lambda_s  <- "lambda.1se"                           # LASSO lambda rule
set.seed(20251207)                                  # Reproducibility

# ---------------- 1. 读数据 ----------------
dat <- as.data.frame(read_excel(data_path, sheet = 1))
dat <- dat[complete.cases(dat), ]
if (!names(dat)[1] %in% c("Group", "group")) {
  warning("第一列不是 Group/group，请检查数据结构！")
}
y   <- as.numeric(dat[[1]])                     # 结局 0/1
X   <- dat[, -1, drop = FALSE]                  # 全部影像特征
n   <- nrow(dat)
stopifnot(all(y %in% c(0, 1)), sum(y) >= 5, ncol(X) >= 10)
cat(sprintf("样本 n=%d | 事件 %d (%.1f%%) | 特征数 %d\n",
            n, sum(y), 100 * mean(y), ncol(X)))

# ---------------- 2. 工具函数 ----------------
# 公式构造：特征名含数字开头（如 10Percentile），必须加反引号
fml_make <- function(resp, vars) {
  as.formula(paste(resp, "~", paste(sprintf("`%s`", vars), collapse = " + ")))
}
# 单因素 p 值：两组均正态 -> t 检验；否则 Wilcoxon（与 CBCgrps 逻辑一致）
pval_one <- function(x, g) {
  if (length(unique(x)) < 3) return(NA_real_)
  x0 <- x[g == 0]; x1 <- x[g == 1]
  ok0 <- tryCatch(shapiro.test(x0)$p.value > 0.05, error = function(e) FALSE)
  ok1 <- tryCatch(shapiro.test(x1)$p.value > 0.05, error = function(e) FALSE)
  if (ok0 && ok1) {
    tryCatch(t.test(x0, x1)$p.value, error = function(e) NA_real_)
  } else {
    tryCatch(wilcox.test(x0, x1)$p.value, error = function(e) NA_real_)
  }
}
# 在给定样本上跑完整管线：筛选 -> LASSO -> 多因素 logistic
fit_pipeline <- function(Xd, g) {
  cand <- names(which(!is.na(vapply(Xd, pval_one, numeric(1), g = g)) &
                        vapply(Xd, pval_one, numeric(1), g = g) < p_screen))
  if (length(cand) < 2) return(NULL)            # 候选太少，无法稳定跑 CV
  Xc <- as.matrix(Xd[, cand, drop = FALSE])
  cvf <- tryCatch(cv.glmnet(Xc, g, family = "binomial", alpha = 1,
                            type.measure = "deviance"),
                  error = function(e) NULL)
  if (is.null(cvf)) return(NULL)
  cf  <- coef(cvf, s = lambda_s)
  sel <- rownames(cf)[as.numeric(cf[, 1]) != 0]
  sel <- setdiff(sel, "(Intercept)")            # 去掉截距行
  if (length(sel) < 1) return(NULL)
  dd  <- data.frame(y = g, Xd[, sel, drop = FALSE])
  fit <- tryCatch(glm(fml_make("y", sel), data = dd, family = binomial),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(fit = fit, sel = sel)
}
# 用已拟合模型预测新样本概率
pred_prob <- function(fit, Xnew, sel) {
  predict(fit, newdata = Xnew[, sel, drop = FALSE], type = "response")
}
# AUC（预测概率需至少 2 个不同取值）
auc_of <- function(pred, truth) {
  if (length(unique(pred)) < 2) return(NA_real_)
  as.numeric(auc(roc(truth, pred, quiet = TRUE,
                     levels = c(0, 1), direction = "<")))
}

# ---------------- 3. 全样本"表观"模型 + LRT + HL ----------------
full     <- fit_pipeline(X, y)
stopifnot(!is.null(full))
auc_app  <- auc_of(pred_prob(full$fit, X, full$sel), y)
cat(sprintf("全样本最终入选特征（%d个）：%s\n", length(full$sel),
            paste(full$sel, collapse = ", ")))

d_full <- data.frame(y = y, X[, full$sel, drop = FALSE])
fit_null <- glm(y ~ 1, data = d_full, family = binomial)
fit_full <- glm(fml_make("y", full$sel), data = d_full, family = binomial)
an       <- anova(fit_null, fit_full, test = "Chisq")
lrt_X2   <- an$Deviance[2]; lrt_p <- an$`Pr(>Chi)`[2]

f_lrm <- lrm(fml_make("y", full$sel), data = d_full, x = TRUE, y = TRUE)
gof   <- residuals(f_lrm, type = "gof")          # c(X2, df, P)
cat(sprintf("整体模型 LRT: X2=%.2f, p=%.4f | Hosmer-Lemeshow: X2=%.2f, df=%g, p=%.4f\n",
            lrt_X2, lrt_p, gof[1], gof[2], gof[3]))

# ---------------- 4. Bootstrap 全流程乐观校正 ----------------
cat(sprintf("开始 Bootstrap 全流程验证（B=%d，每次抽样内完整重跑筛选+LASSO+glm）...\n", B))
O_vec   <- numeric(B); ok <- 0
sel_tab <- list()
for (b in 1:B) {
  idx <- sample(n, n, replace = TRUE)
  pb  <- fit_pipeline(X[idx, , drop = FALSE], y[idx])
  if (is.null(pb)) next
  a_boot <- auc_of(pred_prob(pb$fit, X[idx, , drop = FALSE], pb$sel), y[idx])
  a_orig <- auc_of(pred_prob(pb$fit, X, pb$sel), y)
  if (is.na(a_boot) || is.na(a_orig)) next
  ok <- ok + 1
  O_vec[ok]   <- a_boot - a_orig                  # 单次乐观度
  sel_tab[[ok]] <- pb$sel
  if (b %% 100 == 0) cat(sprintf("  bootstrap %d/%d（成功 %d）\n", b, B, ok))
}
if (ok < 100) warning(sprintf("有效重抽样仅 %d 次，建议增大 B。", ok))
O_vec    <- O_vec[1:ok]
optimism <- mean(O_vec)
auc_corr <- auc_app - optimism
ci_o     <- quantile(O_vec, c(0.025, 0.975))

cat("\n================ 验证结果汇总 ================\n")
cat(sprintf("表观 AUC          : %.3f\n", auc_app))
cat(sprintf("乐观度均值        : %.3f (95%% CI: %.3f-%.3f)\n",
            optimism, ci_o[1], ci_o[2]))
cat(sprintf("乐观校正后 AUC    : %.3f\n", auc_corr))
cat(sprintf("有效 bootstrap 次数: %d / %d\n", ok, B))

# 特征入选频次：检验"3特征"是否稳定，还是每次换一套
freq    <- sort(table(unlist(sel_tab)), decreasing = TRUE)
freq_df <- data.frame(feature = names(freq), 入选次数 = as.integer(freq),
                      入选比例 = round(as.integer(freq) / ok, 3))
cat("\n各特征在 bootstrap 中的入选频次（越高越稳）：\n")
print(head(freq_df, 15))

# ---------------- 5. 保存结果 ----------------
res <- data.frame(
  metric = c("apparent_AUC", "optimism_mean", "optimism_2.5pct", "optimism_97.5pct",
             "optimism_corrected_AUC", "valid_bootstrap_n", "LRT_X2", "LRT_p",
             "HL_X2", "HL_df", "HL_p", "final_features"),
  value  = c(round(auc_app, 4), round(optimism, 4), round(ci_o[1], 4),
             round(ci_o[2], 4), round(auc_corr, 4), ok,
             round(lrt_X2, 4), round(lrt_p, 4), round(gof[1], 4),
             gof[2], round(gof[3], 4), paste(full$sel, collapse = ";"))
)
write.csv(res, file.path(out_dir, "full_pipeline_validation_summary.csv"),
          row.names = FALSE)
write.csv(freq_df, file.path(out_dir, "feature_selection_frequency.csv"),
          row.names = FALSE)

pdf(file.path(out_dir, "optimism_distribution.pdf"), width = 6, height = 5)
hist(O_vec, breaks = 30, col = "steelblue", border = "white",
     main = "Bootstrap optimism (AUC) 分布",
     xlab = "单次重抽样的乐观度 = AUC(样本内) - AUC(原样本)")
abline(v = optimism, col = "red", lwd = 2, lty = 2)
text(optimism, max(hist(O_vec, plot = FALSE)$counts) * 0.9,
     paste("mean =", round(optimism, 3)), col = "red", pos = 4)
dev.off()

cat("\n结果已保存：\n")
cat("  full_pipeline_validation_summary.csv  （表观/校正AUC、LRT、HL）\n")
cat("  feature_selection_frequency.csv       （特征入选频次）\n")
cat("  optimism_distribution.pdf             （乐观度分布）\n")
