# run_all.R —— 一键运行全流程
# 首次使用请先运行 source("setup.R") 安装依赖

cat("========== 开始运行 GSE10072 分析全流程 ==========\n\n")

source("download_data.R")      # 第 0 步：下载数据
source("01_deg_analysis.R")    # 第 1 步：差异表达 + 火山图 + 热图
source("02_enrichment.R")      # 第 2 步：GO/KEGG 富集
source("03_model.R")           # 第 3 步：LASSO + ROC

cat("\n========== 全流程完成，产物都在 results/ 目录 ==========\n")
