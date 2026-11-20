# 01_deg_analysis.R —— 差异表达分析（limma）+ 火山图 + 热图
# 前置：先运行 download_data.R 生成 data/GSE10072_eset.rds

suppressPackageStartupMessages({
  library(Biobase)   # 提供 exprs() / pData() / fData()
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

eset <- readRDS("data/GSE10072_eset.rds")
expr <- exprs(eset)
pd   <- pData(eset)

# ---------- 自动识别 肿瘤/正常 分组 ----------
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

if (is.null(group)) {
  stop("未自动识别出 tumor/normal 分组。请运行 download_data.R 并把样本信息发给我，手动指定分组列。")
}

group <- factor(group, levels = c("normal", "tumor"))
cat("分组统计：\n")
print(table(group))

# ---------- limma 差异分析 ----------
design <- model.matrix(~ group)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/DEG_all.csv")

cat("\n差异基因统计：\n")
cat("  P < 0.05 的基因数：", sum(res$P.Value < 0.05), "\n")
cat("  P < 0.05 且 |logFC| > 1 的基因数：",
    sum(res$P.Value < 0.05 & abs(res$logFC) > 1), "\n")

# ---------- 火山图 ----------
res$logP <- -log10(res$P.Value)
res$sig <- "No"
res$sig[res$P.Value < 0.05 & res$logFC > 1]  <- "Up"
res$sig[res$P.Value < 0.05 & res$logFC < -1] <- "Down"

p <- ggplot(res, aes(x = logFC, y = logP, color = sig)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c(Up = "#d7191c", Down = "#2c7bb6", No = "grey75"),
                     name = "差异") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  labs(x = "log2 Fold Change", y = "-log10 P-value",
       title = "GSE10072：肺腺癌 vs 正常组织（limma）") +
  theme_minimal(base_size = 12)

ggsave("results/volcano.png", p, width = 7, height = 6, dpi = 150)
cat("\n已保存 results/volcano.png\n")

# ---------- Top 50 热图 ----------
top50 <- head(res[order(res$P.Value), ], 50)
mat <- expr[rownames(top50), , drop = FALSE]
mat <- t(scale(t(mat)))   # 按基因做 z-score 标准化

pdf("results/heatmap_top50.pdf", width = 10, height = 12)
pheatmap(mat,
         show_rownames = FALSE,
         annotation_col = data.frame(Group = group, row.names = colnames(mat)),
         main = "Top 50 DEGs")
dev.off()
cat("已保存 results/heatmap_top50.pdf\n")

cat("\n第 1 步完成。产物：DEG_all.csv、volcano.png、heatmap_top50.pdf\n")
