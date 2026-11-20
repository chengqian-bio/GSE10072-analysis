# 02_enrichment.R —— GO / KEGG 功能富集分析
# 前置：先运行 01_deg_analysis.R 生成 results/DEG_all.csv

suppressPackageStartupMessages({
  library(Biobase)   # 提供 exprs() / pData() / fData()
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

eset <- readRDS("data/GSE10072_eset.rds")
res  <- read.csv("results/DEG_all.csv", row.names = 1)
fd   <- fData(eset)

# ---------- 探针 ID 映射到基因符号 ----------
sym_col <- grep("gene.?symbol", colnames(fd), ignore.case = TRUE, value = TRUE)

if (length(sym_col) == 0) {
  # fData 里没有 Gene Symbol 列，改用注释包
  if (!requireNamespace("hgu133a.db", quietly = TRUE)) {
    stop("未找到 Gene Symbol 列。请先运行：BiocManager::install('hgu133a.db')，再重跑本脚本。")
  }
  sym <- AnnotationDbi::mapIds(hgu133a.db::hgu133a.db, keys = rownames(fd),
                               column = "SYMBOL", keytype = "PROBEID")
  map <- data.frame(probe = rownames(fd), symbol = sym, stringsAsFactors = FALSE)
} else {
  map <- data.frame(probe = rownames(fd),
                    symbol = fd[[sym_col[1]]],
                    stringsAsFactors = FALSE)
}

# 去掉空值和无效符号
map <- map[!is.na(map$symbol) & map$symbol != "" & map$symbol != "---", ]

# ---------- 取差异基因 ----------
deg <- res[res$P.Value < 0.05 & abs(res$logFC) > 1, , drop = FALSE]
probes <- rownames(deg)
sym <- map$symbol[match(probes, map$probe)]
sym <- unique(sym[!is.na(sym)])
cat("用于富集的差异基因数：", length(sym), "\n")

ids <- bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("成功映射到 ENTREZID 的基因数：", nrow(ids), "\n")

dir.create("results", showWarnings = FALSE)

# ---------- GO 富集（生物过程 BP）----------
ego <- enrichGO(gene          = ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  png("results/GO_dotplot.png", width = 2000, height = 1600, res = 200)
  print(dotplot(ego, showCategory = 15))
  dev.off()
  write.csv(as.data.frame(ego), "results/GO_BP.csv")
  cat("已保存 results/GO_dotplot.png 和 GO_BP.csv\n")
} else {
  cat("GO 富集没有显著结果，请检查差异基因筛选条件。\n")
}

# ---------- KEGG 通路富集 ----------
ekegg <- tryCatch(
  enrichKEGG(gene = ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05),
  error = function(e) NULL
)

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  png("results/KEGG_dotplot.png", width = 2000, height = 1600, res = 200)
  print(dotplot(ekegg, showCategory = 15))
  dev.off()
  write.csv(as.data.frame(ekegg), "results/KEGG.csv")
  cat("已保存 results/KEGG_dotplot.png 和 KEGG.csv\n")
} else {
  cat("KEGG 富集无结果或网络请求失败（KEGG 需要联网），可稍后单独重试。\n")
}

cat("\n第 2 步完成。\n")
