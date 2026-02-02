library(GEOquery)
gse_number <- "GSE315114"
# 加载数据，如果本地有文件，则优先使用本地的数据
# getGPL控制是否下载对应的 GPL（平台）注释文件
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315114

geo_source_dataset <- getGEO("GSE315114", destdir = "~/bio-project/BioAnalysis/dataset/", getGPL = F)

dataset <- geo_source_dataset[[1]]

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
  message("logs precess ...")
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

gpl_data <- read_gpl_data("~/bio-project/BioAnalysis/dataset/GPL23126.txt")

gpl_table <- gpl_data %>%
  select(probeset_id, gene_assignment) %>%
  separate(gene_assignment,
    into = c("EnsemblId", "GeneSymbol"),
    sep = " // ", fill = "right", extra = "drop"
  ) %>%
  rename(ProbeKey = probeset_id) %>%
  select(ProbeKey, EnsemblId, GeneSymbol)

get_deg_by_limma <- function(express_matrix, group) {
  # 差异分析
  library(limma)

  # 构建设计矩阵（Design Matrix）。这是告诉计算机你的实验设计是什么样的
  design_matrix <- model.matrix(~group)
  # 使用最小二乘法拟合线性模型，根据提供的设计矩阵，计算出每个基因在不同组间的表达差异（即回归系数）
  fit <- lmFit(express_matrix, design_matrix)
  # 对拟合的模型进行经验贝叶斯（Empirical Bayes）调整
  fit <- eBayes(fit)
  # 提取差异表达基因的统计结果，deg 是 Differentially Expressed Genes 的缩写，中文意思是“差异表达基因”
  deg <- topTable(fit, coef = 2, number = Inf)
}


deg <- get_deg_by_limma(express_matrix, group)

# 为deg数据框添加几列
# 加probeKey列，把行名变成一列
library(dplyr)
deg <- mutate(deg, ProbeKey = rownames(deg))
# 加上探针注释（有可能多个探针对应一个基因，进行了去重）
gpl_table <- distinct(gpl_table, GeneSymbol, .keep_all = T)
deg <- inner_join(deg, gpl_table, by = "ProbeKey")
# 如果行数为0说明找的探针注释是错的
nrow(deg)


# 标记上下调基因
logFC_threshold <- 1
p_threshold <- 0.05
condition1 <- (deg$P.Value < p_threshold) & (deg$logFC < -logFC_threshold)
condition2 <- (deg$P.Value < p_threshold) & (deg$logFC > logFC_threshold)
deg <- mutate(deg, change = ifelse(condition1, "down", ifelse(condition2, "up", "stable")))
table(deg$change)


up_genes <- deg[deg$change == "up", ]
down_genes <- deg[deg$change == "down", ]

# 按padj排序
up_genes <- up_genes[order(up_genes$P.Value), ]
down_genes <- down_genes[order(down_genes$P.Value), ]

# 火山图
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
  library(ggplot2)
}

# 准备火山图数据
volcano_data <- deg
volcano_data$neg_log10_padj <- -log10(volcano_data$P.Value)

# 为了更好的可视化，设置padj的上限
max_y <- max(volcano_data$neg_log10_padj[!is.infinite(volcano_data$neg_log10_padj)], na.rm = TRUE)
volcano_data$neg_log10_padj[is.infinite(volcano_data$neg_log10_padj)] <- max_y + 10

# 设置颜色
volcano_data$color <- factor(volcano_data$change, levels = c("down", "stable", "up"))

# 标记top基因（按padj排序的前10个上调和下调基因）
top_up_genes <- head(up_genes, 10)
top_down_genes <- head(down_genes, 10)
volcano_data$label <- NA
volcano_data$label[match(top_up_genes$ProbeKey, volcano_data$ProbeKey)] <- top_up_genes$GeneSymbol
volcano_data$label[match(top_down_genes$ProbeKey, volcano_data$ProbeKey)] <- top_down_genes$GeneSymbol

# 绘制火山图
p <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_padj, color = color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("down" = "#2E86C1", "stable" = "grey", "up" = "#E74C3C"),
    name = "Differential Expression"
  ) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Volcano Plot of Differential Gene Expression",
    subtitle = paste0(
      "Up: ", sum(volcano_data$color == "up"),
      " | Down: ", sum(volcano_data$color == "down")
    ),
    x = "log2 Fold Change",
    y = "-log10 (adjusted p-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 添加基因标签（如果安装了ggrepel）
if (require("ggrepel", quietly = TRUE)) {
  p <- p + geom_text_repel(aes(label = label),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  )
} else {
  # 如果没有ggrepel，使用基础的geom_text
  p <- p + geom_text(aes(label = label),
    size = 3,
    hjust = -0.1,
    vjust = 0.5,
    show.legend = FALSE,
    na.rm = TRUE
  )
}
print(p)

# 差异基因热图
plot_deg_heatmap <- function(diff_express_matrix, group) {
  library(pheatmap)
  annotation_col <- data.frame(group = group)
  rownames(annotation_col) <- colnames(diff_express_matrix)
  pheatmap(diff_express_matrix,
    show_colnames = F,
    show_rownames = F,
    scale = "row",
    annotation_col = annotation_col,
    breaks = seq(-3, 3, length.out = 100)
  )
}

express_matrix <- express_matrix[deg$ProbeKey, ]
rownames(express_matrix) <- deg$GeneSymbol
diff_gene <- deg$GeneSymbol[deg$change != "stable"]
diff_express_matrix <- express_matrix[rownames(express_matrix) %in% diff_gene, ]

plot_deg_heatmap(diff_express_matrix, group)
