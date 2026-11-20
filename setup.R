# setup.R —— 首次运行：安装所有依赖包
# 用法：在 R 控制台运行  source("setup.R")
# 注意：首次安装 Bioconductor 包较慢（约几分钟到十几分钟），属正常现象。

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor 包
bioc_pkgs <- c(
  "GEOquery",       # 下载 GEO 数据
  "limma",          # 差异表达分析
  "clusterProfiler",# 富集分析
  "org.Hs.eg.db"    # 基因 ID 映射（人类）
)

# CRAN 包
cran_pkgs <- c(
  "ggplot2",        # 绘图
  "pheatmap",       # 热图
  "glmnet",         # LASSO 回归
  "pROC"            # ROC 曲线
)

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

cat("\n所有依赖已就绪。接下来可以运行 source('run_all.R')\n")
