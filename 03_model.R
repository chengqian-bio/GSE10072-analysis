# 03_model.R —— LASSO 特征基因筛选 + ROC 曲线
# 方法：留出法，70% 训练 / 30% 测试，AUC 为独立测试集结果（诚实评估）
# 前置：先运行 download_data.R 生成 data/GSE10072_eset.rds

suppressPackageStartupMessages({
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

# ---------- 70/30 划分 ----------
set.seed(42)
n <- nrow(x)
idx <- sample(seq_len(n), size = round(0.7 * n))
train_x <- x[idx, ];  train_y <- y[idx]
test_x  <- x[-idx, ]; test_y  <- y[-idx]

cat("训练集样本数：", length(train_y), " 测试集样本数：", length(test_y), "\n")

# ---------- LASSO（10 折交叉验证选 lambda）----------
cv <- cv.glmnet(train_x, train_y, family = "binomial", alpha = 1, nfolds = 10)

cf <- as.matrix(coef(cv, s = "lambda.min"))
genes <- rownames(cf)[cf[, 1] != 0]
genes <- setdiff(genes, "(Intercept)")

cat("\nLASSO 选出的特征基因数：", length(genes), "\n")
print(genes)

dir.create("results", showWarnings = FALSE)
writeLines(genes, "results/lasso_genes.txt")

# ---------- 测试集 ROC ----------
pred <- predict(cv, newx = test_x, s = "lambda.min", type = "response")
roc_obj <- roc(test_y, as.numeric(pred))

png("results/roc.png", width = 1600, height = 1600, res = 200)
plot(roc_obj, col = "#d7191c", lwd = 3,
     main = paste0("ROC（留出法测试集，AUC = ", round(auc(roc_obj), 3), "）"))
dev.off()

cat("\n测试集 AUC：", round(auc(roc_obj), 3), "\n")
cat("已保存 results/lasso_genes.txt 和 roc.png\n")
cat("\n第 3 步完成。\n")
