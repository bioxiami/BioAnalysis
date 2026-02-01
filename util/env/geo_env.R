if (!require("BiocManager")) {
  install.packages("BiocManager", update = F, ask = F)
}

cran_packages <- c(
  "tidyr",
  "tibble",
  "dplyr",
  "stringr",
  "ggplot2",
  "ggpubr",
  "factoextra",
  "FactoMineR",
  "devtools",
  "cowplot",
  "patchwork",
  "basetheme",
  "paletteer",
  "AnnoProbe",
  "tinyarray",
  "ggthemes",
  "VennDiagram",
  "survminer"
)

biocductor_packages <- c(
  "GEOquery",
  "GO.db",
  "hgu133plus2.db",
  "ggnewscale",
  "limma",
  "impute",
  "GSEABase",
  "GSVA",
  "clusterProfiler",
  "org.Hs.eg.db",
  "preprocessCore",
  "enrichplot"
)

for (pkg in cran_packages) {
  if (!require(pkg, character.only = T, quietly = T)) {
    install.packages(pkg, ask = F, update = F)
    require(pkg, character.only = T)
  }
}


for (pkg in biocductor_packages) {
  if (!require(pkg, character.only = T, quietly = T)) {
    BiocManager::install(pkg, ask = F, update = F)
    require(pkg, character.only = T)
  }
}

# BiocManager::install("ggtree") 报错object 'check_linewidth' not found
# remotes::install_github("YuLab-SMU/ggtree")

# 检查R包是否安装成功
for (pkg in c(biocductor_packages, cran_packages)) {
  require(pkg, character.only = T)
}

