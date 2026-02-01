library(GEOquery)
gse_number <- "GSE315114"
# 加载数据，如果本地有文件，则优先使用本地的数据
# getGPL控制是否下载对应的 GPL（平台）注释文件
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315114
geo_source_dataset <- getGEO("GSE315114", destdir = "~/bio-project/BioAnalysis/dataset/", getGPL = F)
class(geo_source_dataset)

dataset <- geo_source_dataset[[1]]
class(dataset)

express_matrix <- exprs(dataset)

# 查看表达矩阵行数
dim(express_matrix)

# 查看数据范围
range(express_matrix)

# 查看是否有异常样本
# las 参数的作用是控制坐标轴刻度标签（即坐标轴上的数字或文字）的方向
# NA、Inf、负值均为异常样本（可以有少量的负数）
# 正常的情况下，各个样本里面箱线图高度应该相等。
boxplot(express_matrix, las = 2)

# 依据上一步的值判断是否需要取log（+1是为了避免出现负值）
if (max(express_matrix) > 20) {
  express_matrix <- log2(express_matrix + 1)
}

# 提取临床信息（用于构建分组）
pheno_data <- Biobase::pData(dataset)

# 确保临床信息与表达矩阵一一对应
# identical函数精确地判断两个对象是否完全相等
if (!identical(rownames(pheno_data), colnames(express_matrix))) {
  common_sample <- intersect(rownames(pheno_data), colnames(express_matrix))
  express_matrix <- express_matrix[, common_sample]
  pheno_data <- pheno_data[common_sample, ]
}

# 构建分组信息
control_group_key <- "Y"
control_group_name <- "Young"
treatment_group_name <- "Old"


group <- ifelse(str_detect(pheno_data$title, control_group_key), control_group_name, treatment_group_name)

# 将分组转换为因子，对照组在前，处理组在后
group <- factor(group, levels = c(control_group_name, treatment_group_name))

# 提取探针平台编号
gpl_number <- dataset@annotation
