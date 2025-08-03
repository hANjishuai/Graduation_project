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
  log_info("Parsing GMT file: {gmt_file}")
  
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
    
    log_success("Parsed {length(pathway_info)} pathways from GMT file")
    return(pathway_info)
  }, error = function(e) {
    log_error("Failed to parse GMT file: {e$message}")
    stop(e)
  })
}

create_pseudobulk_matrix <- function(seurat_obj, cluster_col) {
  log_info("Creating pseudobulk matrix by cell clusters")
  
  tryCatch({
    cluster_avg <- AverageExpression(
      seurat_obj,
      assays = "RNA",
      features = rownames(seurat_obj),
      group.by = cluster_col,
      slot = "data"
    )$RNA
    
    log_success(glue("Created pseudobulk matrix: {ncol(cluster_avg)} clusters x {nrow(cluster_avg)} genes"))
    return(cluster_avg)
  }, error = function(e) {
    log_error("Failed to create pseudobulk matrix: {e$message}")
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
    log_error("Failed to save GSVA results: {e$message}")
    return(FALSE)
  })
}

generate_pathway_heatmap <- function(gsva_results, pathway_map,output_plot_dir) {
  log_info("Generating pathway activity heatmap")
  
  tryCatch({
    # 修改行名为对应的ID
    rownames(gsva_results) <- pathway_map$pathway_id[
      match(rownames(gsva_results), pathway_map$pathway_name)
    ]
    
    # PNG版本
    heatmap_png <- file.path(output_plot_dir, "pathway_activity_heatmap.png")
    png(heatmap_png, width = 400, height = 650)
    
    ht <- Heatmap(gsva_results,
                  cluster_rows = FALSE,
                  rect_gp = gpar(width = 3,length = 3),
                  show_heatmap_legend = TRUE,
                  name = "GSVA Score",
                  col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                  cluster_columns = FALSE,
                  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                  heatmap_legend_param = list( 
                    title = "GSVA Score",
                    title_gp = gpar(fontsize = 8, fontface = "bold"),
                    labels_gp = gpar(fontsize = 6),
                    direction = "vertical",
                    legend_width = unit(3, "cm"),
                    at = c(-1, 0, 1)
                  ))
    draw(ht, heatmap_legend_side = "right")
    dev.off()
    
    # PDF版本
    heatmap_pdf <- file.path(output_plot_dir, "pathway_activity_heatmap.pdf")
    pdf(heatmap_pdf, width = 10, height = 16)
    ht <- Heatmap(gsva_results,
              cluster_rows = FALSE,
              rect_gp = gpar(width = 3,length = 3),
              show_heatmap_legend = TRUE,
              name = "GSVA Score",
              col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 15, fontface = "bold"),
              column_names_gp = gpar(fontsize = 15, fontface = "bold"),
              heatmap_legend_param = list( 
                title = "GSVA Score",
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 9),
                direction = "vertical",
                legend_width = unit(3, "cm"),
                at = c(-1, 0, 1)
              ))
    draw(ht, heatmap_legend_side = "right")
    dev.off()
    
    log_success(glue("Pathway heatmaps saved to {heatmap_png} and {heatmap_pdf}"))
    return(TRUE)
  }, error = function(e) {
    log_error(glue("Failed to generate pathway heatmap: {e$message}"))
    return(FALSE)
  })
}

generate_gene_heatmaps <- function(pathway_info, cluster_avg, output_plot_dir) {
  log_info(glue("Generating gene expression heatmaps for {length(pathway_info)} pathways"))
  
  n_top_genes <- 50
  success_count <- 0
  
  walk(names(pathway_info), function(pathway_id) {
    genes <- pathway_info[[pathway_id]]$genes
    pathway_name <- pathway_info[[pathway_id]]$name
    
    tryCatch({
      valid_genes <- intersect(genes, rownames(cluster_avg))
      
      if (length(valid_genes) > 1) {
        # 计算基因表达的变异系数
        expr_matrix <- cluster_avg[valid_genes, , drop = FALSE]
        gene_variability <- apply(expr_matrix, 1, function(x) {
          if (mean(x) > 0) sd(x) / mean(x) else 0
        })
        
        # 选择高变基因
        n_show <- min(n_top_genes, length(valid_genes))
        selected_genes <- names(sort(gene_variability, decreasing = TRUE))[1:n_show]
        gene_expr <- cluster_avg[selected_genes, , drop = FALSE]
        
        # 处理NA值
        mat <- as.matrix(gene_expr)
        mat[!is.finite(mat)] <- 0
        mat <- mat[rowSums(mat != 0) >= 2, , drop = FALSE]
        
        # 格式化标题
        if(grepl("\\d{5}", pathway_id)){
          main_title <- paste(
          strwrap(pathway_name, width = 70) %>% paste(collapse = "\n"),
          "\nPathway ID:", paste0("hsa",pathway_id)
        )
        }else {
          main_title <- paste(
            strwrap(pathway_name, width = 70) %>% paste(collapse = "\n"),
            "\nPathway ID:", pathway_id
          )
        }


        p <- pheatmap(
          mat,
          main = main_title,
          scale = "row",
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          clustering_method = "ward.D2",
          fontsize_row = 10,
          show_rownames = TRUE,
          fontface = "bold",
          heatmap_legend_param= list(title="")
        )

        
        # 绘制热图
        output_file <- file.path(output_plot_dir, paste0("gene_expression_", pathway_id, ".pdf"))
        pdf(output_file, width = 7, height = max(6, n_show * 0.35))
        draw(p)
        dev.off()
        
        success_count <<- success_count + 1
      }
    }, error = function(e) {
      log_warning(glue("Failed to generate heatmap for {pathway_id}: {e$message}"))
    })
  })
  
  log_success(glue("Generated {success_count} gene expression heatmaps"))
  return(success_count)
}

#------------------- 主流程 -------------------#
run_gsva_pipeline <- function(gmt_file, seurat_rds, output_dir, output_plot_dir) {
  # 设置日志开始
  log_info("Starting GSVA analysis pipeline")
  log_info(glue("Output directory: {output_dir}"))
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 步骤1: 解析GMT文件
  pathway_info <- parse_gmt(gmt_file)
  
  # 步骤2: 加载Seurat对象
  log_info(glue("Loading Seurat object: {seurat_rds}"))
  seurat_obj <- readRDS(seurat_rds)
  
  # 步骤3: 创建通路ID-名称映射
  pathway_map <- tibble(
    pathway_id = names(pathway_info),
    pathway_name = map_chr(pathway_info, ~ .x$name)
  )
  write_csv(pathway_map, file.path(output_dir, "pathway_id_name_map.csv"))
  
  # 步骤4: 创建基因集集合
  gene_sets <- map(pathway_info, ~ {
    GeneSet(setName = .x$name, geneIds = .x$genes, geneIdType = SymbolIdentifier())
  })
  gsc <- GeneSetCollection(gene_sets)
  
  # 步骤5: 创建伪bulk矩阵
  cluster_avg <- create_pseudobulk_matrix(seurat_obj, cluster_col = "NMF_Cluster")
  
  # 步骤6: 运行GSVA分析
  log_info("Running GSVA analysis with parallel processing")
  parallel.sz <- 96
  
  if (packageVersion("GSVA") >= "1.40.0") {
    param <- gsvaParam(
      exprData = as.matrix(cluster_avg),
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
  log_success("GSVA analysis completed successfully")
  
  # 步骤7: 保存结果
  save_gsva_results(gsva_results, output_dir)
  
  # 步骤8: 生成通路活性热图
  generate_pathway_heatmap(gsva_results, pathway_map, output_plot_dir)
  
  # 步骤9: 生成基因表达热图
  generate_gene_heatmaps(pathway_info, cluster_avg, output_plot_dir)
  
  log_success("GSVA pipeline completed successfully")
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")){
  # 解析配置参数
  config <- list(
    gmt_file = snakemake@input$gmt_file,
    seurat_rds = snakemake@input$seurat_rds,
    output_dir = snakemake@params$output_dir,
    output_plot_dir = snakemake@params$output_plot_dir
  )

  # 运行分析
  run_gsva_pipeline(
    gmt_file = config$gmt_file,
    seurat_rds = config$seurat_rds,
    output_dir = config$output_dir,
    output_plot_dir = config$output_plot_dir
  )
  
  # 创建完成标记文件
  file.create(snakemake@output[[1]])
  log_success("Gene set collection built successfully")
}


