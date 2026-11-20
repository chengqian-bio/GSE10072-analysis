# GSE10072-analysis

肺腺癌（LUAD）与正常肺组织的基因表达差异分析及诊断模型构建。

## 数据来源

- GEO GSE10072：58 例肺腺癌 + 49 例正常肺组织，Affymetrix HG-U133A（GPL96）
- 原始数据**不随仓库分发**，由 `download_data.R` 自动下载到 `data/` 目录

## 分析内容

1. **差异表达分析**（limma）：筛选差异基因，输出火山图、Top50 热图
2. **功能富集**（clusterProfiler）：GO 生物过程 + KEGG 通路富集
3. **诊断模型**（LASSO）：筛选特征基因，留出法测试集 ROC 评估

## 环境依赖

- R ≥ 4.4
- 依赖包：GEOquery、limma、clusterProfiler、org.Hs.eg.db、ggplot2、pheatmap、glmnet、pROC

## 运行方式

首次运行先安装依赖：

```r
source("setup.R")
```

一键跑全流程：

```r
source("run_all.R")
```

也可以分步运行，便于排查问题：

```r
source("download_data.R")      # 下载数据
source("01_deg_analysis.R")    # 差异表达 + 火山图 + 热图
source("02_enrichment.R")      # GO/KEGG 富集
source("03_model.R")           # LASSO + ROC
```

## 结果产物（results/ 目录）

| 文件 | 说明 |
| --- | --- |
| DEG_all.csv | 全部基因的差异分析结果 |
| volcano.png | 火山图 |
| heatmap_top50.pdf | Top50 差异基因热图 |
| GO_BP.csv / GO_dotplot.png | GO 富集结果 |
| KEGG.csv / KEGG_dotplot.png | KEGG 富集结果 |
| lasso_genes.txt | LASSO 筛选出的特征基因 |
| roc.png | 测试集 ROC 曲线 |
