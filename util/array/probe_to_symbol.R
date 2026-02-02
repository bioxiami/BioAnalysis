gpl_number <- "GPL23126"

# 方法一：GEO 官网下载 GPL 注释表
get_gpl_table <- function(gpl_number, destdir = getwd(), download = FALSE) {
  if (!str_starts(gpl_number, "GPL|gpl")) {
    stop("Invalid GPL Number")
  }
  gpl_number <- str_to_upper(gpl_number)
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
    gpl_number, "&targ=self&form=text&view=data"
  )
  message(url)
  if (download) {
    utils::download.file(url, destfile = paste0(destdir, gpl_number, ".txt"))
    message("gpl table download successfully")
  }
}

options(timeout = 600000)
# get_gpl_table(gpl_number, destdir = "~/bio-project/BioAnalysis/dataset/", download = T)
get_gpl_table(gpl_number, download = T)

read_gpl_data <- function(gpl_file_path) {
  # 读取文件，跳过注释行，找到数据起始位置
  lines <- readLines(gpl_file_path)

  # 找到数据表开始的位置（!platform_table_begin之后）
  table_start <- grep("!platform_table_begin", lines) + 1

  # 从table_start位置开始读取数据
  read.table(gpl_file_path,
    header = TRUE,
    sep = "\t",
    skip = table_start - 1,
    quote = "",
    fill = TRUE,
    comment.char = ""
  )
}

gpl_data <- read_gpl_data("~/bio-project/BioAnalysis/dataset/GPL23126.txt")

gpl_table <- gpl_data %>%
  select(probeset_id, gene_assignment) %>%
  separate(gene_assignment,
    into = c("EnsemblId", "GeneSymbol"),
    sep = " // ", fill = "right", extra = "drop"
  ) %>%
  rename(ProbeKey = probeset_id) %>%
  select(ProbeKey, EnsemblId, GeneSymbol)


# 方法二：用 GEOquery 自动获取
# 如果报错parsing failed--expected only one '!series_data_table_begin'，说明可能是芯片太新了，或者下载失败！
# 如果下载下来的文件是3.8K，内容是HTML文件，说明是下载失败
# 会在指定目录下载GPL570.soft.gz文件
gpl_table <- getGEO(gpl_number, destdir = "~/bio-project/BioAnalysis/dataset")
gpl_table <- Table(gpl_table)
names(gpl_table)
probe_symbol <- gpl_table[, c("ID", "Gene Symbol")]

# 方法三：使用 Bioconductor 官方注释包
gpl_anno_pkg <- tinyarray::pkg_all
anno_pkg_name <- gpl_anno_pkg[gpl_anno_pkg$gpl == gpl_number, 2]

if (!require(hgu133plus2.db)) {
  BiocManager::install("hgu133plus2.db", ask = F, update = F)
}
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
probe_symbol <- toTable(hgu133plus2SYMBOL)

# 方法四：使用geneExpressionFromGEO注释包
# install.packages("geneExpressionFromGEO")
# https://cran.r-project.org/web/packages/geneExpressionFromGEO/index.html
# GPL11532, GPL23126, GPL6244, GPL8300, GPL80, GPL96, GPL570, GPL571, GPL20115, GPL1293, GPL6102, GPL6104, GPL6883, GPL6884, GPL13497, GPL14550, GPL17077, GPL6480
gse_number <- "GSE315114"
dataset <- geneExpressionFromGEO::getDatasetFeaturesFromGEOcode(gse_number, verbose = FALSE)
