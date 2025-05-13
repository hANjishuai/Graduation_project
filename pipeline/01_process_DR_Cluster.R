#' 01_process_DR_Cluster
#' Single-cell RNA-seq Processing Pipeline
#' 
#' @description 该脚本用于处理Scanpy转换后的Seurat对象，包含质控、过滤、标准化、降维和聚类。
#' @param input_rds 输入文件路径 (Scanpy转换的Seurat对象)
#' @param output_rds 输出文件路径 (处理后的Seurat对象)
#' @param config 包含以下参数的列表：
#'               - species: 物种 ("human"/"mouse")
#'               - dims: PCA降维维度 (默认20)
#'               - subsets_param: 过滤参数 (nFeature_RNA_min, nFeature_RNA_max, percent_mt_max)
#'               - output_figures: 图表输出路径列表 (vln_raw, scatter_raw, vln_qc, scatter_qc, elbow)
#' @return 处理后的Seurat对象并保存结果
#' @examples
#' config <- list(
#'   species = "human",
#'   dims = 20,
#'   subsets_param = list(nFeature_RNA_min=50, nFeature_RNA_max=7500, percent_mt_max=5),
#'   output_figures = list(
#'     vln_raw = "figures/vln_raw.png",
#'     scatter_raw = "figures/scatter_raw.png",
#'     vln_qc = "figures/vln_qc.png",
#'     scatter_qc = "figures/scatter_qc.png",
#'     elbow = "figures/elbow.pdf"
#'   )
#' )
#' process_seurat_pipeline("input.rds", "output.rds", config)

# 加载依赖包
.libPaths("~/R/4.4.1/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(glmGamPoi)
  library(harmony)
  library(dplyr)
})

#------------------- 核心函数 -------------------#
load_seurat_data <- function(input_rds, species) {
 message("\n[1/7] Loading data from: ", input_rds)
  subcells <- readRDS(input_rds)
  metaData <- subcells@meta.data %>% select(c("sample", "cell_type_lvl1", "condition"))
  counts <- LayerData(subcells, assay = "RNA", layer = "counts")
  
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metaData,
    project = unique(metaData$cell_type_lvl1)
  )
  seurat_obj$orig.ident <- seurat_obj$sample
  message("✓ Data loaded. Cells: ", ncol(seurat_obj), ", Features: ", nrow(seurat_obj))
  return(seurat_obj)
}

calculate_qc_metrics <- function(seurat_obj, species) {
  message("\n[2/7] Calculating QC metrics for species: ", species)
  
  mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  rb_pattern <- ifelse(species == "human", "^RP[SL]", "^Rp[sl]")
  rb_genes <- grep(rb_pattern, rownames(seurat_obj), value = TRUE)
  
  counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
  percent.ribo <- Matrix::colSums(counts[rb_genes, ]) / Matrix::colSums(counts) * 100
  seurat_obj <- AddMetaData(seurat_obj, percent.ribo, col.name = "percent.ribo")
  
  message("✓ QC metrics added: percent.mt & percent.ribo")
  return(seurat_obj)
}

generate_qc_plots <- function(seurat_obj, output_paths) {
  message("\n[3/7] Generating QC plots")
  
  p1 <- VlnPlot(seurat_obj, features = c("percent.ribo", "percent.mt", "nFeature_RNA", "nCount_RNA"),
                group.by = "sample", ncol = 4, pt.size = 0.01, alpha = 0.5)
  ggsave(output_paths[[1]], p1, width = 10, height = 6)
  
  plot1 <- FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA")
  ggsave(output_paths[[2]], plot1 + plot2, width = 10, height = 5)
  
  message("✓ Plots saved to: ", paste(output_paths, collapse = ", "))
}

filter_cells <- function(seurat_obj, subsets_param, species) {
  message("\n[4/7] Filtering cells with params: ", toString(subsets_param))
  seurat_sub <- subset(seurat_obj,
    subset = nFeature_RNA > subsets_param$nFeature_RNA_min &
             nFeature_RNA < subsets_param$nFeature_RNA_max &
             percent.mt < subsets_param$percent_mt_max)
  message("✓ Cells after filtering: ", ncol(seurat_sub), " (Removed ", ncol(seurat_obj) - ncol(seurat_sub), ")")
  seurat_sub <- calculate_qc_metrics(seurat_sub, species)  # 重新计算
  return(seurat_sub)
}

process_seurat <- function(seurat_obj, dims, species, elbow_path) {
  message("\n[5/7] Normalizing data using SCTransform")
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  
  message("[6/7] Running PCA and Harmony integration")
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  p_elbow <- ElbowPlot(seurat_obj)
  ggsave(elbow_path, p_elbow, width = 8, height = 6)
  seurat_obj <- RunHarmony(seurat_obj, "orig.ident")
  
  message("[7/7] Clustering with dims=1:", dims)
  seurat_obj <- RunHarmony(seurat_obj, "orig.ident")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:dims, seed.use = 10)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = seq(0, 1, 0.1), random.seed = 10)
  
  return(seurat_obj)
}

#------------------- 执行主函数 -------------------#
process_seurat_pipeline <- function(input_rds, output_rds, config) {
  options(future.globals.maxSize = 5000000 * 1024^2)
  
  # 数据加载与质控
  seurat_data <- load_seurat_data(input_rds, config$species) %>%
    calculate_qc_metrics(config$species) 
  
  # 生成质控图
  generate_qc_plots(seurat_data, config$output_figures[1:2])
  
  # 细胞过滤
  seurat_filtered <- filter_cells(seurat_data, config$subsets_param,config$species)
  generate_qc_plots(seurat_filtered, config$output_figures[3:4])
  
  # 数据处理与保存
  seurat_processed <- process_seurat(seurat_filtered, config$dims, config$species, config$output_figures)
  saveRDS(seurat_processed, output_rds)
  message("\n✓ Pipeline completed. Output saved to: ", output_rds)
}

# Snakemake集成入口
if (exists("snakemake")) {
  params <- snakemake@params
  
  # 显式转换并匹配参数名称
  config <- list(
    species = params$species,
    dims = as.integer(params$dims),
    subsets_param = list(
      nFeature_RNA_min = as.integer(params$subsets$nFeature_RNA_min),
      nFeature_RNA_max = as.integer(params$subsets$nFeature_RNA_max),
      percent_mt_max = as.numeric(params$subsets$percent_mt_max)
    ),
    output_figures = snakemake@config$output_figures
  )
  
  # 确保目录存在
  dir.create(dirname(snakemake@output$rds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(config$output_figures$vln_raw), recursive = TRUE, showWarnings = FALSE)
  
  # 调用主函数
  process_seurat_pipeline(
    input_rds = snakemake@input$rds,
    output_rds = snakemake@output$rds,
    config = config
  )
}
