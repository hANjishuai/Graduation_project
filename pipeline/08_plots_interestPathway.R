# 加载必要的包
.libPaths("~/R/4.4.1/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(GSEABase)
  library(GSVA)
  library(ggplot2)
  library(ComplexHeatmap)
  library(GO.db)
  library(circlize)
  library(glue)
  library(cli)
  library(purrr)
  library(readr)
  library(BiocParallel)
  library(dplyr)
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  library(openxlsx)
  library(Seurat)
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

#------------------- 核心函数 -------------------#
parse_gmt <- function(gmt_file) {
  log_info(glue("Parsing GMT file: {gmt_file}"))
  
  tryCatch({
    lines <- readLines(gmt_file)
    pathway_info <- list()
    
    for (line in lines) {
      parts <- strsplit(line, "\t")[[1]]
      pathway_id <- parts[1]
      
      # 判断通路类型并获取名称
      if (grepl("^\\d{5}$", pathway_id)) {
        kegg_id <- paste0("hsa", pathway_id)
        pathway_name <- tryCatch({
          kegg_name <- KEGGPATHID2NAME[[kegg_id]]
          if (is.null(kegg_name)) parts[2] else kegg_name
        }, error = function(e) parts[2])
      } else if (grepl("^GO:", pathway_id)) {
        pathway_name <- tryCatch({
          go_term <- Term(GOTERM[[pathway_id]])
          if (is.na(go_term)) parts[2] else go_term
        }, error = function(e) parts[2])
      } else {
        pathway_name <- pathway_id
      }
      
      genes <- parts[-(1:2)] %>% 
        keep(~ .x != "")  # 移除空基因
      
      pathway_info[[pathway_id]] <- list(
        name = pathway_name,
        genes = genes
      )
    }
    
    log_success(glue("Parsed {length(pathway_info)} pathways from GMT file"))
    return(pathway_info)
  }, error = function(e) {
    log_error(glue("Failed to parse GMT file: {e$message}"))
    stop(e)
  })
}

save_gsva_results <- function(gsva_results, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    # 保存为CSV
    csv_file <- file.path(output_dir, "gsva_scores.csv")
    write.csv(gsva_results, csv_file)
    log_info(glue("GSVA scores saved to {csv_file}"))
    
    # 保存为RDS
    rds_file <- file.path(output_dir, "gsva_results.rds")
    saveRDS(gsva_results, rds_file)
    log_info(glue("GSVA results saved to {rds_file}"))
    
    return(TRUE)
  }, error = function(e) {
    log_error(glue("Failed to save GSVA results: {e$message}"))
    return(FALSE)
  })
}

str_sim <- function(a, b) {
  1 - adist(a, b) / max(nchar(a), nchar(b))   # 0 完全不同，1 完全相同
}

#------------------- 新增辅助函数 -------------------#
add_gsva_scores_to_seurat <- function(gsva_results, seurat_obj, pathway_info, output_dir) {
  log_info("Adding GSVA scores to Seurat object metadata")
  
  # 创建通路ID-名称-清理名称映射
  pathway_map <- tibble(
    pathway_id = names(pathway_info),
    pathway_name = map_chr(pathway_info, ~ .x$name)
  )
  
  # 保存映射表
  write_csv(pathway_map, file.path(output_dir, "pathway_id_name_map.csv"))
  
  # 转置GSVA结果并添加清理后的列名
  gsva_df <- as.data.frame(gsva_results)
  gsva_t <- as.data.frame(t(gsva_df))

  # 清理通路名称（转换为有效的R变量名）
  clean_pathway_names <- make.names(rownames(gsva_df))
  colnames(gsva_t) <- clean_pathway_names

  # 添加细胞ID列
  gsva_t$Cell_ID <- rownames(gsva_t)

  # 将GSVA评分添加到metadata
  # 确保细胞ID匹配（Seurat默认使用细胞barcode作为行名）
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = gsva_t
  )
  
  # 保存带注释的通路评分
  extracted_data <- FetchData(seurat_obj, vars = c("NMF_Cluster", clean_pathway_names))
  write.csv(extracted_data, file.path(output_dir, "pathway_scores_with_annotations.csv"), row.names = TRUE)
  saveRDS(seurat_obj, file.path(output_dir, "pathway_scores_with_annotations_seu.rds"))
  
  log_success("GSVA scores added to Seurat object")
  return(list(seurat_obj = seurat_obj, pathway_map = pathway_map))
}

analyze_single_pathway <- function(seurat_obj, pathway_map, target_pathway, target_group= 'NMF_Cluster_5',output_dir, output_plot_dir) {
  log_info(glue("Analyzing pathway: {target_pathway}"))
  
  # 准备分析数据
  gsva_data <- FetchData(seurat_obj, vars = c("NMF_Cluster", target_pathway)) %>%
    rownames_to_column("Cell_ID") %>%
    rename(Group = NMF_Cluster, Score = !!target_pathway) %>%
    select(Cell_ID, Group, Score)
  
  # 执行统计检验
  stat_test <- gsva_data %>%
    dunn_test(Score ~ Group, p.adjust.method = "BH") %>%
    filter(group1 == target_group | group2 == target_group) %>%
    mutate(comparison = ifelse(group1 == target_group, 
                               paste(group2, "vs", target_group),
                               paste(group1, "vs", target_group)))
  
  # 准备Prism数据
  prepare_prism_data <- function(gsva_data) {
    gsva_data %>%
      group_by(Group) %>%
      mutate(row_id = row_number()) %>%
      pivot_wider(
        id_cols = row_id,
        names_from = Group,
        values_from = Score
      ) %>%
      select(-row_id) %>%
      select(sort(colnames(.)))
  }
  
  prism_bar <- prepare_prism_data(gsva_data)
  prism_violin <- prism_bar  # 相同数据结构
  
  # 导出Prism数据
  wb <- createWorkbook()
  addWorksheet(wb, "BarPlot_Data")
  writeData(wb, "BarPlot_Data", prism_bar)
  addWorksheet(wb, "ViolinPlot_Data")
  writeData(wb, "ViolinPlot_Data", prism_violin)
  addWorksheet(wb, "Statistics")
  writeData(wb, "Statistics", stat_test)
  
  output_file <- file.path(output_dir, glue("Prism_Plotting_Data_{make.names(target_pathway)}.xlsx"))
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # 获取通路信息
  pathway_id <- sapply(pathway_map$pathway_name,function(x) str_sim(x,target_pathway)) %>% sort(decreasing=T) %>% names() %>% head(n=1)
  pathway_desc <- pathway_map$pathway_name[pathway_map$pathway_id==pathway_id][[1]]
  
  # 生成图形
  generate_pathway_plot <- function(data, plot_type = "bar") {
    p <- ggplot(data, aes(x = Group, y = Score, fill = Group)) +
      scale_fill_manual(values = scales::hue_pal()(n_distinct(data$Group))) +
      labs(title = glue("Pathway Activity: {pathway_id}"),
           subtitle = pathway_desc,
           x = "B-cell Subsets", y = "GSVA Score") +
      theme_pubr() +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (plot_type == "bar") {
      p + 
        geom_bar(stat = "summary", fun = "mean", alpha = 0.7) +
        geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, shape = 16) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
    } else {
      p +
        geom_violin(alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 1, alpha = 0.5, shape = 16)
    }
  }
  
  # 保存图形
  bar_plot <- generate_pathway_plot(gsva_data, "bar")
  violin_plot <- generate_pathway_plot(gsva_data, "violin")
  
  plot_prefix <- make.names(target_pathway)
  ggsave(file.path(output_plot_dir, glue("BarPlot_{plot_prefix}.png")), bar_plot, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_plot_dir, glue("ViolinPlot_{plot_prefix}.png")), violin_plot, width = 10, height = 6, dpi = 300)
  
  log_success(glue("Analysis completed for pathway: {target_pathway}"))
  return(stat_test)
}

#------------------- 主流程 -------------------#
run_gsva_pipeline <- function(gmt_file, 
                              seurat_rds, 
                              output_dir, 
                              output_plot_dir, 
                              target_group = "NMF_Cluster_5",
                              target_pathway = NULL) {
  # 初始化
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. 解析GMT
  pathway_info <- parse_gmt(gmt_file)
  
  # 2. 加载Seurat对象
  seurat_obj <- readRDS(seurat_rds)
  
  # 3. 创建基因集集合
  gene_sets <- map(pathway_info, ~ {
    GeneSet(setName = .x$name, geneIds = .x$genes, geneIdType = SymbolIdentifier())
  })
  gsc <- GeneSetCollection(gene_sets)
  
  # 4. 运行GSVA
  expr_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  param <- gsvaParam(
    exprData = as.matrix(expr_data),
    geneSets = gsc,
    minSize = 5,
    maxSize = 500,
    kcdf = "Gaussian"
  )
  
  parallel.sz <- 144
  
  if (packageVersion("GSVA") >= "1.40.0") {
    param <- gsvaParam(
      exprData = as.matrix(expr_data),
      geneSets = gsc,
      minSize = 5,
      maxSize = 500,
      kcdf = "Gaussian"
    )
    
    # 设置并行
    if (parallel.sz > 1) {
      if (.Platform$OS.type == "windows") {
        BiocParallel::register(SnowParam(workers = parallel.sz))
      } else {
        BiocParallel::register(MulticoreParam(workers = parallel.sz))
      }
    }
  }
  gsva_results <- gsva(param, verbose = TRUE)
  
  # 5. 保存原始结果
  save_gsva_results(gsva_results, output_dir)
  log_success("GSVA analysis completed successfully")

  # 6. 添加评分到Seurat对象
  result <- add_gsva_scores_to_seurat(gsva_results, seurat_obj, pathway_info, output_dir)
  seurat_obj <- result$seurat_obj
  pathway_map <- result$pathway_map
  
  # 7. 分析特定通路
  if(is.null(target_pathway)){
    pathwayForcircle <- intersect(pathway_map$pathway_name %>% make.names(),colnames(seurat_obj@meta.data))
    for(i in pathwayForcircle){
      log_info(glue("select {i}"))
      target_pathway = i
      analyze_single_pathway(
          seurat_obj = seurat_obj,
          pathway_map = pathway_map,
          target_pathway = target_pathway,
          target_group = target_group ,
          output_dir = output_dir,
          output_plot_dir= output_plot_dir
      )
    }   
    log_success("GSVA pipeline completed successfully")
  } else {
    analyze_single_pathway(
        seurat_obj = seurat_obj,
        pathway_map = pathway_map,
        target_pathway = target_pathway,
        target_group = target_group ,
        output_dir = output_dir,
        output_plot_dir= output_plot_dir
    )
  }
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  run_gsva_pipeline(
    gmt_file = snakemake@input$gmt_file,
    seurat_rds = snakemake@input$seurat_rds,
    output_dir = snakemake@params$output_dir,
    output_plot_dir = snakemake@params$output_plot_dir,
    target_group = snakemake@params$target_group,
    target_pathway = snakemake@params$target_pathway  # 从Snakemake传递
  )
  file.create(snakemake@output[[1]])
}
