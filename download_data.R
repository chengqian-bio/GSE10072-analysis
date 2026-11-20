# download_data.R —— 下载 GSE10072 并保存
# 数据来源：GEO GSE10072（58 例肺腺癌 + 49 例正常肺组织，GPL96）
# 运行后会在 data/ 目录生成 GSE10072_eset.rds

suppressPackageStartupMessages(library(GEOquery))

dir.create("data", showWarnings = FALSE)

cat("开始下载 GSE10072，请耐心等待（受网络影响，可能较慢）...\n")
gse <- getGEO("GSE10072", GSEMatrix = TRUE, destdir = "data")

# GSEMatrix=TRUE 返回一个列表，通常第一个元素就是主表达矩阵
cat("下载到", length(gse), "个表达矩阵对象\n")

eset <- gse[[1]]
saveRDS(eset, "data/GSE10072_eset.rds")

# 打印样本信息，方便核对肿瘤/正常分组
cat("\n==== 样本信息列名 ====\n")
print(colnames(pData(eset)))

cat("\n==== 前 6 个样本的信息 ====\n")
print(head(pData(eset), 6))

cat("\n==== 表达矩阵维度（基因 x 样本）====\n")
print(dim(exprs(eset)))

cat("\n已保存到 data/GSE10072_eset.rds\n")
cat("请把上面打印的样本信息发给我，我来确认分组字段。\n")
