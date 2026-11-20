# setup.R —— 首次运行：安装所有依赖包
# 策略：CRAN 包全部强制用预编译 binary（快），BiocManager 只装 Bioconductor 自己的包，
#       避免 stringi / RSQLite / igraph 等从源码编译（极慢）。

options(install.packages.check.source = "no")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 自动检测本机 7890 端口是否有代理在运行；有则让 R 走代理，无代理环境自动跳过
proxy_ok <- tryCatch({
  con <- socketConnection(host = "127.0.0.1", port = 7890, timeout = 1, blocking = FALSE)
  close(con)
  TRUE
}, error = function(e) FALSE, warning = function(w) FALSE)
if (proxy_ok) {
  Sys.setenv(http_proxy = "http://127.0.0.1:7890")
  Sys.setenv(https_proxy = "http://127.0.0.1:7890")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", type = "win.binary")
}

# ---------- 第 1 步：CRAN 依赖（全部预编译）----------
cran_pkgs <- c(
  "ggplot2", "pheatmap", "glmnet", "pROC",
  "stringi", "stringr", "httr", "curl", "openssl", "xml2",
  "Rcpp", "data.table", "igraph", "systemfonts", "RSQLite",
  "fastmap", "formatR", "base64enc", "cachem", "digest",
  "lambda.r", "futile.options", "gridGraphics", "tweenr", "polyclip",
  "lazyeval", "blob", "memoise", "png", "futile.logger", "snow", "BH",
  "cowplot", "fastmatch", "ggplotify", "patchwork", "ggforce", "ape",
  "DBI", "reshape2", "aplot", "ggfun", "ggnewscale", "ggrepel",
  "ggtangle", "scatterpie", "fs", "rappdirs", "gson", "plyr", "yulab.utils"
)

missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cat("安装 CRAN 依赖（binary）：", paste(missing, collapse = ", "), "\n")
  install.packages(missing, type = "win.binary")
} else {
  cat("CRAN 依赖已全部就绪。\n")
}

# ---------- 第 2 步：Bioconductor 包 ----------
bioc_pkgs <- c(
  "GEOquery",        # 下载 GEO 数据
  "limma",           # 差异表达分析
  "clusterProfiler", # 富集分析
  "org.Hs.eg.db"     # 基因 ID 映射（人类）
)

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("安装 Bioconductor 包：", p, "\n")
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

cat("\n所有依赖已就绪。接下来可以运行 source('run_all.R')\n")
