iseu <- "result_out/01_process_DR_Cluster/celltype_ker/subcell.rds"

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplotify)
library(ggplot2)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 48)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS(iseu)

# 定义抽样比例（例如 10%）
prop <- 0.1

# 分层抽样
set.seed(123)
sampled_cells <- seurat_obj@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(`SCT_snn_res.0.2`) %>%
  sample_frac(size = prop, replace = FALSE) %>%
  pull(cell)

# 创建子集
seurat_subset <- subset(seurat_obj, cells = sampled_cells)

# 验证比例
table(seurat_subset$celltype)

DimPlot(seurat_subset, group.by='SCT_snn_res.0.2', label=TRUE) + ggtitle('Zhou et al Control Cortex') + NoLegend()

seurat_obj <- SetupForWGCNA(
  seurat_subset,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("SCT_snn_res.0.2"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 15, # nearest-neighbors parameter
  max_shared = 8, # maximum number of shared cells between two metacells
  ident.group = "SCT_snn_res.0.2" # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = '0', # the name of the group of interest in the group.by column
  group.by='SCT_snn_res.0.2', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'SCT', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_outdir = "result_out/02_hdWGCNA/celltype_ker/TOM.rda",
  tom_name = 'Ker' # name of the topoligical overlap matrix written to disk
)

p <- PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

