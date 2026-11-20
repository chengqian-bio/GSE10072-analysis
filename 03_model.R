# 03_model.R —— LASSO 特征基因筛选 + ROC 曲线
# 评估方式：10 折交叉验证的 out-of-fold（OOF）预测，比单次留出法更稳健、诚实
# 前置：先运行 download_data.R 生成 data/GSE10072_eset.rds

suppressPackageStartupMessages({
  library(Biobase)   # 提供 exprs() / pData() / fData()
  library(glmnet)
  library(pROC)
})

eset <- readRDS("data/GSE10072_eset.rds")
expr <- exprs(eset)
pd   <- pData(eset)

# ---------- 识别分组（与 01 脚本一致）----------
group <- NULL
for (col in colnames(pd)) {
  vals <- tolower(as.character(pd[[col]]))
  has_tumor  <- any(grepl("tumor|adenocarcinoma|carcinoma", vals))
  has_normal <- any(grepl("normal|non.?tumor|adjacent", vals))
  if (has_tumor && has_normal) {
    group <- ifelse(grepl("tumor|adenocarcinoma|carcinoma", vals), "tumor", "normal")
    break
  }
}
if (is.null(group)) stop("未自动识别出 tumor/normal 分组。")

group <- factor(group, levels = c("normal", "tumor"))
y <- as.numeric(group == "tumor")   # 1 = tumor, 0 = normal
x <- t(expr)                        # 行为样本，列为基因
n <- nrow(x)

# ---------- 10 折交叉验证，收集 OOF 预测 ----------
set.seed(42)
folds <- sample(rep(1:10, length.out = n))
pred <- rep(NA, n)

for (k in 1:10) {
  tr <- which(folds != k)
  te <- which(folds == k)
  cvk <- cv.glmnet(x[tr, ], y[tr], family = "binomial", alpha = 1, nfolds = 10)
  pred[te] <- as.numeric(predict(cvk, newx = x[te, ], s = "lambda.min",
                                 type = "response"))
}

roc_obj <- roc(y, pred)

dir.create("results", showWarnings = FALSE)
png("results/roc.png", width = 1600, height = 1600, res = 200)
plot(roc_obj, col = "#d7191c", lwd = 3,
     main = paste0("ROC (10-fold CV, AUC = ", round(auc(roc_obj), 3), ")"))
dev.off()

cat("10 折交叉验证 AUC：", round(auc(roc_obj), 3), "\n")

# ---------- 用全部数据拟合，输出特征基因 ----------
cv_full <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
cf <- as.matrix(coef(cv_full, s = "lambda.min"))
genes <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")

cat("LASSO 特征基因数：", length(genes), "\n")
print(genes)
writeLines(genes, "results/lasso_genes.txt")

cat("已保存 results/roc.png 和 lasso_genes.txt\n")
cat("\n第 3 步完成。\n")
