#!/usr/bin/env Rscript
# 单细胞数据转换工具 - Scanpy转Seurat
# 版本: 1.0.0
# 命令行运行
#Rscript scanpy2seurat.R \
#  -i output/s07_scanpy2seurat/celltype_ker \
#  -o output/s08_seurat/celltype_ker.rds
#单细胞数据转换流程启动
#🔹 [2023-10-02 14:30:00] 开始时间: 2023-10-02 14:30:00
#
#✅ [2023-10-02 14:30:01] 读取表达矩阵: 32421基因 x 5000细胞
#✅ [2023-10-02 14:30:02] 读取细胞元数据: 8列
#✅ [2023-10-02 14:30:03] 读取降维坐标: pca (50D)
#✅ [2023-10-02 14:30:03] 读取降维坐标: umap (2D)
#⚠️ [2023-10-02 14:30:04] 未找到降维文件: tsne.csv
#
#🔹 [2023-10-02 14:30:05] 创建Seurat对象...
#✅ [2023-10-02 14:30:06] 初始化对象: 5000细胞
#✅ [2023-10-02 14:30:07] 添加降维坐标: pca
#✅ [2023-10-02 14:30:08] 添加降维坐标: umap
#
#🔹 [2023-10-02 14:30:09] 保存结果文件...
#✅ [2023-10-02 14:30:10] 成功保存Seurat对象: output/s08_seurat/celltype_ker.rds
#
#An object of class Seurat 
#5000 features across 5000 samples 
#Active assay: RNA 
#3 dimensional reductions calculated: pca, umap

.libPaths("~/R/4.4.1/library/")

library(Seurat)
library(Matrix)
library(argparse)
library(cli)
library(glue)

# --------------------------
# 日志系统
# --------------------------
log_step <- function(message, symbol = "🔹") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_info("{.timestamp [{ts}]} {symbol} {message}")
}

log_success <- function(message) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_success("{.timestamp [{ts}]} ✅ {message}")
}

log_warning <- function(message) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_warning("{.timestamp [{ts}]} ⚠️ {message}")
}

# --------------------------
# 核心功能函数
# --------------------------
read_sc_data <- function(input_dir) {
  log_step("加载输入文件...")
  
  # 读取表达矩阵并转置
  mat_path <- file.path(input_dir, "matrix.mtx")
  if(!file.exists(mat_path)) stop("未找到matrix.mtx文件")
  raw_counts <- readMM(mat_path)
  counts <- t(raw_counts)  # 关键修复：转置为 genes x cells
  log_success(glue("读取表达矩阵: {nrow(counts)}基因 x {ncol(counts)}细胞"))
  
  # 读取基因名称
  gene_path <- file.path(input_dir, "genes.tsv")
  if(file.exists(gene_path)){
    gene_names <- read.delim(gene_path, header=FALSE)$V1
    rownames(counts) <- gene_names
    log_success(glue("加载基因名称: {length(gene_names)}个"))
  } else {
    rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))
    log_warning("未找到genes.tsv，使用自动生成基因名")
  }
  
  # 读取细胞元数据
  meta_path <- file.path(input_dir, "metadata.csv")
  if(!file.exists(meta_path)) stop("未找到metadata.csv文件")
  metadata <- read.csv(meta_path, row.names = 1)
  colnames(counts) <- rownames(metadata)  # 设置细胞名称
  log_success(glue("读取细胞元数据: {ncol(metadata)}列"))
  
  # 读取降维数据
  read_reduction <- function(name) {
    path <- file.path(input_dir, glue("{name}.csv"))
    if(file.exists(path)) {
      df <- read.csv(path, row.names = 1)
      # 验证细胞名称一致性
      if(!all(rownames(df) %in% rownames(metadata))) {
        stop(glue("{name}降维数据包含未知细胞"))
      }
      log_success(glue("读取降维坐标: {name} ({ncol(df)}D)"))
      return(as.matrix(df))
    }
    log_warning(glue("未找到降维文件: {name}.csv"))
    return(NULL)
  }
  
  reductions <- list(
    pca = read_reduction("pca"),
    umap = read_reduction("umap"),
    tsne = read_reduction("tsne")
  )
  
  return(list(
    counts = counts,
    metadata = metadata,
    reductions = reductions
  ))
}

create_seurat <- function(sc_data) {
  log_step("创建Seurat对象...")
  
  # 创建基础对象
  seu <- CreateSeuratObject(
    counts = sc_data$counts,
    meta.data = sc_data$metadata,
    assay = "RNA"
  )
  log_success(glue("初始化对象: {ncol(seu)}细胞 x {nrow(seu)}基因"))
  
  # 添加降维数据
  for (reduc_name in names(sc_data$reductions)) {
    reduc <- sc_data$reductions[[reduc_name]]
    if(!is.null(reduc)) {
      colnames(reduc) <- paste0(reduc_name, "_", 1:ncol(reduc))
      seu[[reduc_name]] <- CreateDimReducObject(
        embeddings = reduc[colnames(seu), ], # 确保顺序一致
        key = glue("{toupper(reduc_name)}_")
      )
      log_success(glue("添加降维坐标: {reduc_name}"))
    }
  }
  
  return(seu)
}

# --------------------------
# 主函数
# --------------------------
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-i", "--input_dir", required=TRUE, 
                     help="输入目录路径（需包含matrix.mtx, metadata.csv等）")
  parser$add_argument("-o", "--output", required=TRUE,
                     help="输出RDS文件路径（*.rds）")
  args <- parser$parse_args()
  
  input_dir <- args$input_dir
  output <- args$output

  cli_h1("单细胞数据转换流程启动")
  log_step(glue("开始时间: {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
  
  tryCatch({
    sc_data <- read_sc_data(input_dir)
    seu <- create_seurat(sc_data)
    
    log_step("保存结果文件...")
    saveRDS(seu, output)
    log_success(glue("成功保存Seurat对象: {output}"))
    
    cli_h1("转换完成")
    print(seu)
  }, error = function(e) {
    cli_alert_danger("流程异常终止: {e$message}")
    quit(status = 1)
  })
}

if (!interactive()) {
  main()
}
