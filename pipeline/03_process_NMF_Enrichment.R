.libPaths("~/R/4.4.1/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(NMF)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(ggprism)
  library(glue)
  library(cli)
})

#------------------- 日志系统 -------------------#
log_info <- function(message, symbol = "🔹") {
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

log_error <- function(message) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_danger("{.timestamp [{ts}]} ❌ {message}")
}

#------------------- NMF 核心函数 (支持跳过计算) -------------------#
run_nmf_analysis <- function(seu, config) {
  # 创建NMF专用对象
  sub_sce <- CreateSeuratObject(
    counts = seu@assays$RNA$counts,
    meta.data = seu@meta.data,
    assay = "RNA"
  )
  
  # 数据预处理
  sub_sce <- NormalizeData(sub_sce) %>% 
    FindVariableFeatures() %>% 
    ScaleData(do.center = FALSE)
  
  vm <- sub_sce@assays$RNA$scale.data
  
  # 检查是否提供了已有的NMF结果
  use_existing_nmf <- !is.null(config$existing_nmf_model) && file.exists(config$existing_nmf_model) &&
    !is.null(config$existing_final_nmf) && file.exists(config$existing_final_nmf)
  
  if (use_existing_nmf) {
    # 使用已有NMF结果
    log_info("Loading existing NMF models")
    res <- readRDS(config$existing_nmf_model)
    estimate <- readRDS(config$existing_final_nmf)
    
    # 确定最佳rank
    coph_diffs <- diff(res$measures$cophenetic)
    rank <- which.max(abs(coph_diffs)) + (config$rank_start - 1)
    
    log_success(glue("Loaded existing NMF models with rank: {rank}"))
    
    return(list(
      nmf_result = res,
      final_model = estimate,
      seurat_obj = sub_sce,
      optimal_rank = rank
    ))
  }
  
  # 从头计算NMF
  ranks <- seq(config$rank_start, config$rank_stop, config$rank_step)
  
  # 初始NMF分析
  log_info(glue("Running NMF with ranks: {min(ranks)}-{max(ranks)}"))
  res <- nmf(
    vm,
    rank = ranks,
    nrun = config$nrun,
    .options = "v+m+p10"
  )
  
  # 保存模型和评估图
  saveRDS(res, file = config$nmf_model_path)
  ggsave(
    config$rank_estimate_plot,
    plot(res),
    width = 12,
    height = 8
  )
  
  # 确定最佳rank
  coph_diffs <- diff(res$measures$cophenetic)
  rank <- which.max(abs(coph_diffs)) + (config$rank_start - 1)
  log_success(glue("Selected optimal rank: {rank}"))
  
  # 最终NMF模型
  log_info("Running final NMF estimation")
  estimate <- nmf(
    vm,
    rank = rank,
    nrun = config$nrun,
    .options = "v+m+p10",
    stop = TRUE
  )
  saveRDS(estimate, file = config$final_nmf_path)
  
  return(list(
    nmf_result = res,
    final_model = estimate,
    seurat_obj = sub_sce,
    optimal_rank = rank
  ))
}

integrate_nmf_results <- function(nmf_results, resolution) {
  sub_sce <- nmf_results$seurat_obj
  rank <- nmf_results$optimal_rank
  estimate <- nmf_results$final_model
  
  # 创建NMF降维对象
  sub_sce <- RunPCA(sub_sce,verbose = F)
  sub_sce@reductions$nmf <- sub_sce@reductions$pca
  # 修复：正确提取系数矩阵
  coef_matrix <- as.matrix(coef(estimate))
  basis_matrix <- as.matrix(basis(estimate))

  sub_sce@reductions$nmf@cell.embeddings <- t(coef_matrix)
  sub_sce@reductions$nmf@feature.loadings <- basis_matrix 
  
  # 基于NMF的UMAP和聚类
  sub_sce <- RunUMAP(sub_sce, 
                     reduction = "nmf",
                     dims = 1:rank,
                     reduction.name = "umap_nmf")
  
  sub_sce <- FindNeighbors(sub_sce, 
                           reduction = "nmf", 
                           dims = 1:rank) %>% 
    FindClusters(resolution = resolution)
  
  # 添加NMF聚类信息
  cluster_assignments <- apply(coef_matrix[1:rank, ], 2, which.max)
  sub_sce$nmf_cluster <- cluster_assignments
  sub_sce$NMF_Cluster <- paste0("NMF_Cluster_", cluster_assignments) %>% 
    as.factor()
  
  Idents(sub_sce) <- 'NMF_Cluster'
  log_success(glue("Identified {length(levels(sub_sce))} NMF clusters"))
  
  return(sub_sce)
}

extract_top_genes <- function(nmf_results, topN) {
  estimate <- nmf_results$final_model
  sig.order <- extractFeatures(estimate, topN)
  sig.order <- lapply(sig.order, function(x) rownames(estimate)[x])
  names(sig.order) <- paste0("NMF_Cluster_", seq_along(sig.order))
  
  log_success(glue("Extracted top {topN} genes per cluster"))
  return(sig.order)
}

#------------------- 富集分析函数 -------------------#
run_cluster_enrichment <- function(gene_list, cluster_name, results_dir, species = "human") {
  # 确保目录存在
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  # 基因ID转换
  convert_gene_symbols <- function(symbols) {
    gene_list <- clusterProfiler::bitr(
      symbols, 
      fromType = "SYMBOL", 
      toType = "ENTREZID", 
      OrgDb = if (species == "human") org.Hs.eg.db else org.Mm.eg.db
    )
    readr::write_csv(gene_list, file.path(results_dir, glue("GeneConversion.csv")))
    return(gene_list)
  }
  
  # KEGG富集
  run_kegg_enrichment <- function(entrez_ids) {
    enrich_res <- enrichKEGG(
      gene = entrez_ids,
      organism = if (species == "human") 'hsa' else 'mmu',
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
    
    if (!is.null(enrich_res) && nrow(enrich_res) > 0) {
      result_df <- enrich_res@result %>% 
        dplyr::mutate(DB = "KEGG")
      readr::write_csv(result_df, file.path(results_dir, glue("KEGG_Enrichment.csv")))
    }
    return(if (!is.null(enrich_res)) enrich_res@result else NULL)
  }
  
  # GO富集
  run_go_enrichment <- function(entrez_ids) {
    enrich_res <- enrichGO(
      gene = entrez_ids,
      OrgDb = if (species == "human") org.Hs.eg.db else org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH"
    )
    
    if (!is.null(enrich_res) && nrow(enrich_res) > 0) {
      result_df <- enrich_res@result %>% 
        dplyr::mutate(DB = "GO")
      readr::write_csv(result_df, file.path(results_dir, glue("GO_Enrichment.csv")))
    }
    return(if (!is.null(enrich_res)) enrich_res@result else NULL)
  }
  
  # 处理富集结果数据
  process_enrichment_data <- function(kegg_res, go_res) {
    combined <- bind_rows(
      if (!is.null(kegg_res)) kegg_res %>% dplyr::mutate(DB = "KEGG") else tibble(),
      if (!is.null(go_res)) go_res %>% dplyr::mutate(DB = "GO") else tibble()
    )
    
    if (nrow(combined) == 0) return(combined)
    
    combined %>%
      dplyr::mutate(
        # 计算富集分数 = GeneCount / BgCount
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
        # 缩短通路名称
        Description_short = stringr::str_trunc(Description, width = 50),
        # 计算-log10(p.adjust)用于颜色映射
        neg_log10_padj = -log10(p.adjust)
      ) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::filter(!duplicated(Description_short)) %>%
      dplyr::mutate(
        # 确保通路名称顺序正确
        Description_short = factor(Description_short, levels = Description_short[order(p.adjust, decreasing = TRUE)])
      )
  }
  
  # 绘制转置条形图
  plot_enrichment_bars <- function(plot_data, color_palette = NULL) {
    # 只取前10个显著通路
    top10 <- plot_data %>% 
      dplyr::slice_min(p.adjust, n = 10) %>%
      dplyr::mutate(Description_short = forcats::fct_reorder(Description_short, GeneRatio_num))
    
    # 默认配色方案
    if (is.null(color_palette)) {
      p <- ggplot(top10, aes(x = GeneRatio_num, y = Description_short, fill = neg_log10_padj)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(option = "D", name = "-log10(adj.p)") +
        labs(title = glue("Top Enriched Pathways: {cluster_name}"),
             x = "Enrichment Score (Gene Ratio)",
             y = "") +
        theme_prism(base_size = 12) +
        theme(axis.text.y = element_text(size = 10))
    } else {
      # 自定义配色
      p <- ggplot(top10, aes(x = GeneRatio_num, y = Description_short, fill = neg_log10_padj)) +
        geom_bar(stat = "identity") +
        scale_fill_gradientn(colors = color_palette, name = "-log10(adj.p)") +
        labs(title = glue("Top Enriched Pathways: {cluster_name}"),
             x = "Enrichment Score (Gene Ratio)",
             y = "") +
        theme_prism(base_size = 12) +
        theme(axis.text.y = element_text(size = 10))
    }
    
    return(p)
  }
  
  # 保存配色方案
  save_color_palette <- function(plot_data, output_path) {
    # 创建颜色映射元数据
    color_meta <- list(
      min_value = min(plot_data$neg_log10_padj, na.rm = TRUE),
      max_value = max(plot_data$neg_log10_padj, na.rm = TRUE),
      default_palette = viridisLite::viridis(10, option = "D")
    )
    
    saveRDS(color_meta, output_path)
  }

  # 执行富集流程
  log_info(glue("Processing cluster: {cluster_name}"))
  converted_genes <- convert_gene_symbols(gene_list)
  
  kegg_res <- run_kegg_enrichment(converted_genes$ENTREZID)
  go_res <- run_go_enrichment(converted_genes$ENTREZID)
  
  combined_res <- process_enrichment_data(kegg_res, go_res)

  if (nrow(combined_res) == 0) {
    log_warning(glue("No enrichment results for cluster {cluster_name}"))
    return()
  }
  
  # 保存完整结果
  readr::write_csv(combined_res, file.path(results_dir, "Enrich_summary.csv"))
  
  # 提取前10通路数据并保存
  top10_data <- combined_res %>% 
    dplyr::slice_min(p.adjust, n = 10) %>%
    dplyr::select(DB, ID, Description, GeneRatio, GeneRatio_num, p.adjust, neg_log10_padj)
  
  readr::write_csv(top10_data, file.path(results_dir, "top10_pathways.csv"))
  
  # 绘制并保存条形图
  p <- plot_enrichment_bars(combined_res)
  ggsave(
    file.path(results_dir, "top10_pathways.pdf"),
    plot = p,
    width = 10,
    height = 6
  )
  
  ggsave(
    file.path(results_dir, "top10_pathways.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # 保存默认配色方案
  save_color_palette(
    combined_res,
    file.path(results_dir, "color_palette.rds")
  )
  
  log_success(glue("Completed enrichment for {cluster_name}"))
}

#------------------- 配色调整函数 -------------------#
replot_with_custom_colors <- function(results_dir, custom_colors = NULL) {
  # 加载数据
  plot_data <- readr::read_csv(file.path(results_dir, "top10_pathways.csv"))
  color_meta <- readRDS(file.path(results_dir, "color_palette.rds"))
  
  # 获取簇名
  cluster_name <- basename(results_dir)
  
  # 创建新图
  if (is.null(custom_colors)) {
    p <- ggplot(plot_data, aes(x = GeneRatio_num, y = Description, fill = neg_log10_padj)) +
      geom_bar(stat = "identity") +
      scale_fill_gradientn(
        colors = color_meta$default_palette,
        limits = c(color_meta$min_value, color_meta$max_value),
        name = "-log10(adj.p)"
      ) +
      labs(title = glue("Top Enriched Pathways: {cluster_name}"),
           x = "Enrichment Score (Gene Ratio)",
           y = "") +
      theme_prism(base_size = 12)
  } else {
    p <- ggplot(plot_data, aes(x = GeneRatio_num, y = Description, fill = neg_log10_padj)) +
      geom_bar(stat = "identity") +
      scale_fill_gradientn(
        colors = custom_colors,
        limits = c(color_meta$min_value, color_meta$max_value),
        name = "-log10(adj.p)"
      ) +
      labs(title = glue("Top Enriched Pathways: {cluster_name}"),
           x = "Enrichment Score (Gene Ratio)",
           y = "") +
      theme_prism(base_size = 12)
  }
  
  # 保存新图
  ggsave(
    file.path(results_dir, "custom_color_pathways.pdf"),
    plot = p,
    width = 10,
    height = 6
  )
  
  ggsave(
    file.path(results_dir, "custom_color_pathways.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # 返回绘图对象用于进一步调整
  return(p)
}

#------------------- 主流程函数 (支持跳过NMF计算) -------------------#
process_nmf_pipeline <- function(input_rds, config) {
  # 加载数据
  log_info(glue("Loading Seurat object: {input_rds}"))
  seu <- readRDS(input_rds)
  
  # 执行NMF分析（支持跳过）
  nmf_results <- run_nmf_analysis(seu, config)
  
  # 整合结果到Seurat
  sub_sce <- integrate_nmf_results(nmf_results, config$resolution)
  
  # 提取top基因
  top_genes <- extract_top_genes(nmf_results, config$topN)
  
  # 富集分析
  log_info("Starting enrichment analysis")
  for (cluster_name in names(top_genes)) {
    cluster_dir <- file.path(config$enrich_dir, cluster_name)
    run_cluster_enrichment(
      gene_list = top_genes[[cluster_name]],
      cluster_name = cluster_name,
      results_dir = cluster_dir,
      species = config$species
    )
  }
  
  # 保存最终对象
  saveRDS(sub_sce, config$output_rds)
  log_success(glue("Pipeline completed. Results saved to {config$output_rds}"))
  return(sub_sce)
}

#------------------- Snakemake集成入口 (支持跳过NMF计算) -------------------#
if (exists("snakemake")) {
  # 构建配置参数
  config <- list(
    rank_start = snakemake@params$rank_start,
    rank_stop = snakemake@params$rank_stop,
    rank_step = snakemake@params$rank_step,
    nrun = snakemake@params$nrun,
    topN = snakemake@params$topN,
    resolution = snakemake@params$resolution,
    species = snakemake@params$species,
    nmf_model_path = snakemake@output$nmf_model,
    final_nmf_path = snakemake@output$final_model,
    rank_estimate_plot = snakemake@output$rank_plot,
    enrich_dir = snakemake@params$enrich_dir,
    output_rds = snakemake@output$rds,
    # 添加可选的已有NMF结果路径
    existing_nmf_model = if ("existing_nmf_model" %in% names(snakemake@params)) 
      snakemake@params$existing_nmf_model else NULL,
    existing_final_nmf = if ("existing_final_nmf" %in% names(snakemake@params)) 
      snakemake@params$existing_final_nmf else NULL
  )
  
  # 运行主流程
  process_nmf_pipeline(
    input_rds = snakemake@input[[1]],
    config = config
  )
} else {
  # 直接运行示例配置
  config <- list(
    rank_start = 5,
    rank_stop = 8,
    rank_step = 1,
    nrun = 10,
    topN = 50,
    resolution = 0.2,
    species = "human",
    nmf_model_path = "nmf_models/initial_model.rds",
    final_nmf_path = "nmf_models/final_model.rds",
    rank_estimate_plot = "results/rank_estimate.pdf",
    enrich_dir = "results/enrichment",
    output_rds = "results/nmf_seurat.rds",
    # 可选：提供已有NMF结果路径
    existing_nmf_model = NULL,  # "existing_models/nmf_initial.rds"
    existing_final_nmf = NULL   # "existing_models/nmf_final.rds"
  )
  
  process_nmf_pipeline("input.rds", config)
}
