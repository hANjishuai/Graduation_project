# 整合脚本：整合UMAP绘图、GSVA分析和热图生成
.libPaths("~/R/4.4.1/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(GSVA)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(GO.db)
  library(tidyverse)
  library(ggprism)
  library(glue)
  library(cli)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
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

#--------------------来源其他脚本的函数-------------------#
get_gene_set <- function(pathway_id, species = "human") {
  # 确定物种
  species_code <- switch(
    species,
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
    "rat" = "Rattus norvegicus",
    "Homo sapiens"
  )

  get_kegg_species_code <- function(species) {
    switch(
      species,
      "human" = "hsa",
      "mouse" = "mmu",
      "rat" = "rno",
      "hsa"  # 默认是人类
    )
  }
  
  # 处理不同格式的ID
  clean_id <- pathway_id
  
  # 处理GO ID格式
  if (grepl("^GO:", clean_id) || grepl("^GO\\d+$", clean_id)) {
    # 提取数字部分
    go_id <- gsub("^GO:|^GO", "", clean_id)
    go_id <- paste0("GO:", go_id)
    
    # 使用AnnotationDbi获取GO基因集
    if (species == "human") {
      gene_ids <- get(go_id, org.Hs.egGO2ALLEGS)
    } else if (species == "mouse") {
      gene_ids <- get(go_id, org.Mm.egGO2ALLEGS)
    } else if (species == "rat") {
      gene_ids <- get(go_id, org.Rn.egGO2ALLEGS)
    }
    
    if (!is.null(gene_ids) && length(gene_ids) > 0) {
      # 将Entrez ID转换为基因符号
      if (species == "human") {
        symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID")
      } else if (species == "mouse") {
        symbols <- mapIds(org.Mm.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID")
      } else if (species == "rat") {
        symbols <- mapIds(org.Rn.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID")
      }
      
      genes <- unname(symbols[!is.na(symbols)])
      
      if (length(genes) > 0) {
        # 获取GO描述
        term <- tryCatch({
          GOTERM[[go_id]]@Term
        }, error = function(e) go_id)
        
        return(list(
          genes = unique(genes),
          db = "GO",
          id = go_id,
          description = term
        ))
      }
    }
    
    # 如果AnnotationDbi失败，尝试msigdbr（但只获取整个GO数据库）
    log_warning(glue("Could not retrieve GO gene set {go_id} through AnnotationDbi. Trying msigdbr..."))
    go_df <- msigdbr(species = species_code, category = "C5")
    
    if (!is.null(go_df) && nrow(go_df) > 0) {
      go_df <- go_df %>% 
        filter(gs_exact_source == go_id)
      
      if (nrow(go_df) > 0) {
        genes <- unique(go_df$gene_symbol)
        description <- unique(go_df$gs_name)
        
        return(list(
          genes = genes,
          db = "GO",
          id = go_id,
          description = ifelse(length(description) > 0, description[1], go_id)
        ))
      }
    }
    
    # 最后尝试GO API
    log_warning(glue("Could not retrieve GO gene set {go_id} through standard methods. Trying API..."))
    go_url <- glue("https://api.geneontology.org/api/bioentity/function/{go_id}/genes?rows=1000")
    go_response <- tryCatch({
      jsonlite::fromJSON(go_url)
    }, error = function(e) NULL)
    
    if (!is.null(go_response) && "associations" %in% names(go_response)) {
      genes <- unique(go_response$associations$gene$symbol)
      
      if (length(genes) > 0) {
        return(list(
          genes = genes,
          db = "GO",
          id = go_id,
          description = go_response$associations$subject$label[1]
        ))
      }
    }
  }
  
  # 处理KEGG ID格式
  if (grepl("^hsa", clean_id) || grepl("^\\d+$", clean_id)) {
    # 提取数字部分
    kegg_id <- gsub("^hsa", "", clean_id)

    # 获取KEGG物种代码
    kegg_species_code <- get_kegg_species_code(species)
    full_kegg_id <- paste0(kegg_species_code, kegg_id)

    # 尝试msigdbr (最可靠的方法)
    log_info(glue("Retrieving KEGG gene set for {full_kegg_id} through msigdbr"))
    kegg_df <- msigdbr(species = species_code, category = "C2", subcategory = "CP:KEGG")

    if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
      # 查找匹配的KEGG通路
      kegg_match <- kegg_df %>% 
        filter(gs_exact_source == full_kegg_id)

      if (nrow(kegg_match) == 0) {
        # 宽松匹配
        kegg_match <- kegg_df %>% 
          filter(str_detect(gs_exact_source, fixed(full_kegg_id, ignore_case = TRUE)))
      }

      if (nrow(kegg_match) > 0) {
        genes <- unique(kegg_match$gene_symbol)
        description <- unique(kegg_match$gs_name)

        return(list(
          genes = genes,
          db = "KEGG",
          id = kegg_id,
          description = ifelse(length(description) > 0, description[1], full_kegg_id)
        ))
      }
    }

    # 使用KEGG API (更可靠)
    log_warning(glue("Could not retrieve KEGG gene set {full_kegg_id} through msigdbr. Trying KEGG API..."))

    # 第一步：获取通路信息
    kegg_info_url <- glue("https://rest.kegg.jp/get/{full_kegg_id}")
    kegg_info <- tryCatch({
      readLines(kegg_info_url)
    }, error = function(e) NULL)

    description <- full_kegg_id
    if (!is.null(kegg_info)) {
      # 提取通路名称
      name_line <- grep("^NAME", kegg_info, value = TRUE)
      if (length(name_line) > 0) {
        description <- gsub("^NAME\\s+", "", name_line[1]) %>% 
          gsub("\\s+\\(.*$", "", .)  # 移除括号中的额外信息
      }
    }

    # 第二步：获取通路基因
    kegg_genes_url <- glue("https://rest.kegg.jp/link/{kegg_species_code}/{full_kegg_id}")
    kegg_genes_response <- tryCatch({
      readLines(kegg_genes_url)
    }, error = function(e) NULL)

    if (!is.null(kegg_genes_response)) {
      # 解析基因ID (格式: path:hsa04920\thsa:10000)
      gene_ids <- str_extract(kegg_genes_response, "hsa:\\d+$") %>%
        gsub("hsa:", "", .) %>%
        unique()

      if (length(gene_ids) > 0) {
        # 转换Entrez ID为基因符号
        if (species == "human") {
          entrez_to_symbol <- suppressMessages(
            bitr(gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
          )
        } else if (species == "mouse") {
          entrez_to_symbol <- suppressMessages(
            bitr(gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
          )
        } else if (species == "rat") {
          entrez_to_symbol <- suppressMessages(
            bitr(gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Rn.eg.db)
          )
        }

        if (nrow(entrez_to_symbol) > 0) {
          return(list(
            genes = unique(entrez_to_symbol$SYMBOL),
            db = "KEGG",
            id = kegg_id,
            description = description
          ))
        }
      }
    }

    # 最终尝试：使用clusterProfiler
    log_warning(glue("Could not retrieve KEGG gene set through API. Trying clusterProfiler..."))

    tryCatch({
      kegg_pathways <- download_KEGG(kegg_species_code, keggType = "KEGG")

      if (!is.null(kegg_pathways)) {
        # 查找特定通路
        pathway_entry <- kegg_pathways$KEGGPATHID2NAME %>% 
          filter(from == full_kegg_id)

        if (nrow(pathway_entry) > 0) {
          kegg_genes <- kegg_pathways$KEGGPATHID2EXTID %>% 
            filter(from == full_kegg_id) %>% 
            pull(to)

          if (length(kegg_genes) > 0) {
            # 转换Entrez ID为基因符号
            if (species == "human") {
              entrez_to_symbol <- bitr(kegg_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
            } else if (species == "mouse") {
              entrez_to_symbol <- bitr(kegg_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
            } else if (species == "rat") {
              entrez_to_symbol <- bitr(kegg_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Rn.eg.db)
            }

            if (nrow(entrez_to_symbol) > 0) {
              return(list(
                genes = unique(entrez_to_symbol$SYMBOL),
                db = "KEGG",
                id = kegg_id,
                description = pathway_entry$Description[1]
              ))
            }
          }
        }
      }
    }, error = function(e) {
      log_warning(glue("clusterProfiler failed: {e$message}"))
    })
  }


  log_error(glue("Failed to retrieve gene set for pathway: {pathway_id}"))
  return(NULL)
}

run_gsva_analysis <- function(seurat_obj,
                              pathway_ids, 
                              species = "human", 
                              output_dir = "results/gsva", 
                              pathway_df = NULL) {
  # 创建输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 获取Seurat对象中的所有基因
  all_genes <- rownames(seurat_obj)
  
  # 构建基因集列表
  log_info(glue("Retrieving gene sets for {length(pathway_ids)} pathways"))
  
  genesets <- list()
  pathway_info <- tibble()
  valid_pathways <- 0
  
  for (id in pathway_ids) {
    gene_set <- get_gene_set(id, species)
    
    if (!is.null(gene_set)) {
      # 过滤掉Seurat对象中不存在的基因
      valid_genes <- intersect(gene_set$genes, all_genes)
      
      if (length(valid_genes) >= 5) {
        geneset_name <- glue("{gene_set$description}")
        genesets[[geneset_name]] <- valid_genes
        
        pathway_info <- bind_rows(
          pathway_info,
          tibble(
            Pathway = geneset_name,
            DB = gene_set$db,
            ID = gene_set$id,
            Description = gene_set$description,
            Genes = paste(valid_genes, collapse = ","),
            GeneCount = length(valid_genes)
          )
        )
        
        valid_pathways <- valid_pathways + 1
      } else {
        log_warning(glue("Pathway {id} has only {length(valid_genes)} valid genes (minimum 5 required)"))
      }
    }
  }
  
  if (length(genesets) == 0) {
    log_error("No valid gene sets found for GSVA analysis")
    stop("Gene sets not found")
  }
  
  log_success(glue("Retrieved {valid_pathways} valid gene sets"))
  
  # 保存基因集信息
  write_csv(pathway_info, file.path(output_dir, "pathway_info.csv"))
  
  # 准备伪bulk表达矩阵
  log_info("Creating pseudobulk expression matrix")
  
  Idents(seurat_obj) <- "NMF_Cluster"
  clusters <- levels(seurat_obj)
  
  # 获取所有基因集涉及的所有基因
  all_gs_genes <- unique(unlist(genesets))
  
  pb_matrix <- matrix(
    nrow = length(all_gs_genes),
    ncol = length(clusters),
    dimnames = list(all_gs_genes, clusters)
  )
  
  for (cluster in clusters) {
    cluster_cells <- WhichCells(seurat_obj, idents = cluster)
    if (length(cluster_cells) > 0) {
      expr <- GetAssayData(seurat_obj, slot = "data")[all_gs_genes, cluster_cells, drop = FALSE]
      pb_matrix[, cluster] <- Matrix::rowMeans(expr)
    } else {
      pb_matrix[, cluster] <- 0
    }
  }
  
  # 过滤掉在少于5%细胞中表达的基因
  log_info("Filtering low-expression genes")
  expr_frac <- rowSums(pb_matrix > 0) / ncol(pb_matrix)
  keep_genes <- expr_frac > 0.05
  pb_matrix <- pb_matrix[keep_genes, ]
  
  # 过滤掉太小的基因集
  filtered_genesets <- list()
  for (name in names(genesets)) {
    valid_genes <- genesets[[name]][genesets[[name]] %in% rownames(pb_matrix)]
    if (length(valid_genes) >= 5) {
      filtered_genesets[[name]] <- valid_genes
    }
  }
  
  log_info(glue("Using {length(filtered_genesets)} gene sets after filtering"))
  
  # 运行GSVA
  # 设置并行参数
  parallel.sz <- 4 

  log_info("Running GSVA analysis")
  if (packageVersion("GSVA") >= "1.40.0") {
    param <- GSVA::gsvaParam(
      exprData = as.matrix(pb_matrix),
      geneSets = filtered_genesets,
      minSize = 5,
      maxSize = 500,
      kcdf = "Gaussian"
    )

    # 设置并行
    if (parallel.sz > 1) {
      if (.Platform$OS.type == "windows") {
        BiocParallel::register(BiocParallel::SnowParam(workers = parallel.sz))
      } else {
        BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel.sz))
      }
    }

    gsva_scores <- GSVA::gsva(param, verbose = TRUE)
  }
  
  # 准备结果数据
  log_info("Preparing results")
  
  # GSVA打分矩阵
  gsva_df <- as.data.frame(gsva_scores) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(
      cols = -Pathway,
      names_to = "Cluster",
      values_to = "GSVA_Score"
    ) %>%
    left_join(
      pathway_info %>% 
        dplyr::select(Pathway, DB, ID, Description),
      by = "Pathway"
    ) %>%
    dplyr::select(DB, ID, Description, Cluster, GSVA_Score)
  
  # 添加cluster信息到结果中
  if (!is.null(pathway_df)) {
    gsva_df <- gsva_df %>%
      left_join(
        pathway_df %>% 
          dplyr::select(ID, Origin_Cluster = Cluster),
        by = "ID"
      ) %>%
      dplyr::relocate(Origin_Cluster, .after = Description)
  } else {
    gsva_df <- gsva_df %>%
      mutate(Origin_Cluster = "All") %>%
      dplyr::relocate(Origin_Cluster, .after = Description)
  }

  # 保存结果
  write_csv(gsva_df, file.path(output_dir, "gsva_scores.csv"))
  
  # 创建默认配色
  default_colors <- scales::hue_pal()(length(clusters))
  names(default_colors) <- clusters
  saveRDS(default_colors, file.path(dirname(output_dir), "gsva_color_palette.rds"))
  
  log_success("GSVA analysis completed")
  
  return(list(
    scores = gsva_df,
    color_palette = default_colors
  ))
}

run_gsva_pipeline <- function(seurat_rds,
                              pathway_df,
                              species = "human",
                              top_n = 10,
                              color_palette = NULL,
                              output_dir = "result_out/GSVA",
                              output_plot_dir = "figure_out/GSVA/plot") {
  # 加载Seurat对象
  log_info(glue("Loading Seurat object: {seurat_rds}"))
  seu <- readRDS(seurat_rds)
  
  for(c in unique(pathway_df$Cluster)){
    sub_pathway_df <- subset(pathway_df, subset=Cluster==c)
    pathway_ids <- unique(sub_pathway_df$ID)

    cluster_dir <- glue("{output_dir}/{c}")
    if(!dir.exists(cluster_dir)) dir.create(cluster_dir,recursive = T)
    
    # 运行GSVA分析
    results <- run_gsva_analysis(seurat_obj = seu, 
                                pathway_ids = pathway_ids, 
                                species = species, 
                                output_dir = cluster_dir,
                                pathway_df = sub_pathway_df  # 传递cluster信息
                                )

    # 绘制条形图
    cluster_plot_dir <- glue("{output_plot_dir}/{c}")
    if(missing(color_palette)){
      plot_gsva_bars_prism(results$scores, results$color_palette, top_n, cluster_plot_dir)
    } else {
      plot_gsva_bars_prism(results$scores, color_palette, top_n, cluster_plot_dir)
    }

  }

  log_success(glue("GSVA pipeline completed. Results saved to {output_dir} and {output_plot_dir}"))
  
  return(results)
}

#------------------- 绘图函数 -------------------#
plot_gsva_bars_prism <- function(gsva_scores, color_palette, top_n = 10, output_plot_dir = "results/gsva") {
  
  # 加载必要包
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(dplyr)) install.packages("dplyr")
  if (!require(forcats)) install.packages("forcats")
  if (!require(ggprism)) install.packages("ggprism")
  
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(ggprism)
  
  # 获取所有唯一的目标聚类
  all_clusters <- unique(gsva_scores$Cluster)
  origin_cluster <- unique(gsva_scores$Origin_Cluster %>% na.omit())
  gsva_scores$Origin_Cluster <- origin_cluster
  
  # 验证配色方案
  if (missing(color_palette)) {
    stop("必须提供color_palette参数")
  }
  
  if (!all(origin_cluster %in% names(color_palette))) {
    missing_clusters <- setdiff(origin_cluster, names(color_palette))
    stop(paste("配色方案中缺少以下Origin_Cluster:", paste(missing_clusters, collapse = ", ")))
  }
  
  # 数据预处理：选择每个来源聚类的top通路
  plot_data <- gsva_scores %>%
    group_by(Origin_Cluster, Description) %>%
    mutate(Max_Score = max(abs(gsva_scores$GSVA_Score))) %>%
    ungroup() %>%
    group_by(gsva_scores$Origin_Cluster) %>%
    arrange(desc(Max_Score)) %>%
    distinct(Description, .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    dplyr::select(-Max_Score) %>%
    inner_join(gsva_scores, by = c("Origin_Cluster", "Description")) %>%
    mutate(
      # 标记是否为Origin Cluster
      is_origin = gsva_scores$Cluster == origin_cluster,
      # 按最大GSVA得分排序通路
      Description = fct_reorder(gsva_scores$Description, gsva_scores$GSVA_Score, .fun = max, .desc = TRUE),
      Cluster.y = factor(Cluster.y, levels = rev(all_clusters))
    ) %>%
    # 添加颜色信息
    mutate(
      bar_color = ifelse(is_origin, 
                         color_palette[as.character(Origin_Cluster)], 
                         "grey80")
    )
  
  plot_data <- plot_data %>% 
    mutate(Description_short = stringr::str_trunc(as.character(Description), width = 75, side ="left"))
  
  log_info(glue("{plot_data$Description_short}"))

  manual_colors <- set_names(plot_data$bar_color,plot_data$Cluster.y)
  if(manual_colors[[origin_cluster]]!=color_palette[[origin_cluster]]){
    log_error(glue("未设定正确的颜色{color_palette[[origin_cluster]]},现在是{manual_colors[[origin_cluster]]}"))
  }
    

  # 创建条形图 - Prism风格
  p <- ggplot(plot_data, aes(x = Description_short, y = GSVA_Score.y, fill = Cluster.y)) +
    geom_bar(
      aes(color = Cluster.y, ),
      stat = "identity", 
      position = position_dodge(width = 0.8), 
      width = 0.7,
      size = 0.3
    ) +
    # 使用指定的颜色方案
    scale_fill_manual(values = manual_colors) +
    scale_color_manual(values = manual_colors) +
    coord_flip() +
    labs(
      title = paste("GSVA Scores -", origin_cluster),
      x = "Pathway",
      y = "GSVA Score"
    ) +
    theme_prism(
      base_size = 12,
      base_fontface = "plain",
      base_family = "sans"
    ) +
    theme(
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none",  # 移除图例
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.background = element_rect(fill = "white")
    ) +
    # 添加零线
    geom_hline(yintercept = 0, color = "black", size = 0.5)
  
  # 添加分面（如果有多个来源聚类）
  if (length(unique(plot_data$Origin_Cluster)) > 1) {
    p <- p + 
      facet_wrap(~ Origin_Cluster, scales = "free_y", ncol = 1) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", color = "black")
      )
  }
  
  # 创建输出目录
  if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)
  
  # 保存图片
  ggsave(
    file.path(output_plot_dir, "gsva_barplot_prism_style.pdf"),
    plot = p,
    width = 10,
    height = max(6, top_n * 0.4),
    device = cairo_pdf
  )
  
  ggsave(
    file.path(output_plot_dir, "gsva_barplot_prism_style.png"),
    plot = p,
    width = 10,
    height = max(6, top_n * 0.4),
    dpi = 300,
    bg = "white"
  )
  
  message("Prism-style GSVA barplot created successfully!")
  
  return(p)
}


#------------------- 1. UMAP绘图函数 --------------------#
generate_umap_plot <- function(seu, 
                               color_palette = NULL, 
                               output_dir,
                               output_plot_dir) {
  # 如果没有提供配色方案，创建默认配色
  if (is.null(color_palette)) {
    clusters <- sort(levels(Idents(seu)))
    message = paste0(clusters,collapse=",")
    log_info(glue("使用的seurat对象的Idents是{message}"))
    default_colors <- scales::hue_pal()(length(clusters))
    color_palette <- setNames(default_colors, clusters)
    
    # 保存默认配色供后续修改
    color_df <- data.frame(cluster = names(color_palette), colors = color_palette, row.names = NULL)
    write.csv(color_df, file.path(output_dir, "default_colors.csv"), row.names = FALSE)
    message("\n")
    log_info(glue("Default color palette saved to:  {(file.path(output_dir, 'default_colors.csv'))}"))
    message("\n")
    log_info("Edit this file and re-run with --color_palette argument to use custom colors.csv")
  } else {
    # 加载自定义配色
    if (file.exists(color_palette)) {
      color_df <- read.csv(color_palette)
      color_palette <- setNames(color_df$colors, color_df$cluster)
    } else {
      warning("Color palette file not found: ", color_palette, ". Using default colors.")
      clusters <- sort(levels(Idents(seu)))
      color_palette <- setNames(scales::hue_pal()(length(clusters)), clusters)
    }
  }
  
  # 保存最终使用的配色方案
  saveRDS(color_palette, file.path(output_dir, "umap_color_palette.rds"))
  
  # 生成UMAP图
  p <- DimPlot(seu,
               label = FALSE,
               label.color = "black",
               alpha = 0.8,
               pt.size = 0.8,
               repel = TRUE,
               label.size = 2.5) + 
    scale_color_manual(values = color_palette) +
    theme(legend.justification = c(0, 1),
          legend.position = c(0.05, 0.99),
          legend.byrow = FALSE,
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.1, "pt"))
  
  # 保存图片
  umap_pdf <- file.path(output_dir, "umap_plot.pdf")
  ggsave(umap_pdf, p, width = 4, height = 4)
  log_info(glue("UMAP plot saved to: {umap_pdf}"))
  
  return(color_palette)
}

#------------------- 2. GSVA分析函数 --------------------#
run_gsva_analysis_re <- function(seurat_rds,
                              pathway_df,
                              species,
                              top_n,
                              color_palette,
                              output_dir,
                              output_plot_dir
                              ) {
  # 确保目录存在
  gsva_output_dir <- file.path(output_dir, "gsva_results")
  dir.create(gsva_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  gsva_plot_output_dir <- file.path(output_plot_dir, "gsva_results")
  dir.create(gsva_plot_output_dir, showWarnings = FALSE, recursive = TRUE)

  # 加载Seurat对象和通路数据
  seu <- readRDS(seurat_rds)
  pathway_df <- read.csv(pathway_df)
  
  # 运行GSVA分析（这里简化了实际分析过程）
  message("Running GSVA analysis...")
  run_gsva_pipeline(seurat_rds = seurat_rds,
                    pathway_df = pathway_df,
                    species = species,
                    top_n = top_n,
                    color_palette = color_palette,
                    output_dir = gsva_output_dir,
                    output_plot_dir = gsva_plot_output_dir
                  )
}

#------------------- 3. 热图生成函数 --------------------#
generate_integrated_heatmap <- function(base_dir, 
                                        color_palette, 
                                        output_dir,
                                        output_plot_dir) {
  # 获取所有NMF_Cluster目录路径
  cluster_dirs <- list.dirs(base_dir, recursive = FALSE) %>% 
    str_subset("NMF_Cluster_\\d+")
  
  # 读取所有gsva_scores.csv文件并合并
  all_data <- map_dfr(cluster_dirs, ~ {
    file_path <- file.path(.x, "gsva_scores.csv")
    if (file.exists(file_path)) {
      df <- read_csv(file_path, show_col_types = FALSE) 
      Origin_Cluster <- unique(na.omit(df$Origin_Cluster))
      df$Origin_Cluster <- Origin_Cluster
      df
    }
  })
  
  # 创建表达矩阵 (Description × Cluster)
  expression_matrix <- all_data %>%
    dplyr::select(Description, Cluster, GSVA_Score) %>%
    group_by(Description, Cluster) %>%
    summarise(GSVA_Score = mean(GSVA_Score), .groups = "drop") %>%  # 处理重复值
    pivot_wider(
      names_from = Cluster,
      values_from = GSVA_Score,
      values_fill = 0  # 填充缺失值为0
    ) %>%
    column_to_rownames("Description") %>%
    as.matrix()
  
  # 创建通路-Origin_Cluster映射表
  pathway_origin_map <- all_data %>%
    dplyr::select(Description, Origin_Cluster) %>%
    distinct(Description, .keep_all = TRUE) %>%
    arrange(Origin_Cluster, Description)
  
  # 保存结果
  expr_matrix_file <- file.path(output_dir, "gsva_expression_matrix.csv")
  map_file <- file.path(output_dir, "pathway_origin_mapping.csv")
  write.csv(expression_matrix, expr_matrix_file)
  write.csv(pathway_origin_map, map_file)
  message("Expression matrix saved to: ", expr_matrix_file)
  message("Pathway origin mapping saved to: ", map_file)
  
  # 创建行注释向量
  row_split_vec <- pathway_origin_map$Origin_Cluster[
    match(rownames(expression_matrix), pathway_origin_map$Description)
  ]
  
  # 确保颜色映射包含所有Origin_Cluster
  all_clusters <- unique(row_split_vec)
  missing_clusters <- setdiff(all_clusters, names(color_palette))
  if (length(missing_clusters) > 0) {
    warning("Color palette missing clusters: ", paste(missing_clusters, collapse = ", "),
            ". Assigning default colors.")
    default_colors <- scales::hue_pal()(length(missing_clusters))
    color_palette[missing_clusters] <- default_colors
  }
  
  # 创建行注释
  ha_row <- rowAnnotation(
    Origin = row_split_vec,
    col = list(Origin = color_palette),
    annotation_legend_param = list(
      title = "Origin Cluster",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold"),
      grid_height = unit(0.7, "cm"),
      grid_width = unit(0.7, "cm"),
      direction = "horizontal",
      nrow = 2),
    show_annotation_name = FALSE
  )
  
  # 绘制热图
  heatmap_png <- file.path(output_plot_dir, "gsva_heatmap.png")
  heatmap_pdf <- file.path(output_plot_dir, "gsva_heatmap.pdf")
  
  # PNG版本
  png(heatmap_png, width = 400, height = 1050)
  ht <- Heatmap(expression_matrix,
                cluster_rows = FALSE,
                rect_gp = gpar(width = 3),
                show_heatmap_legend = FALSE,
                name = "GSVA Score",
                col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                right_annotation = ha_row,
                row_split = row_split_vec,
                cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                heatmap_legend_param = list( 
                  title = "GSVA Score",
                  title_gp = gpar(fontsize = 8, fontface = "bold"),
                  labels_gp = gpar(fontsize = 6),
                  direction = "horizontal",
                  legend_width = unit(1.5, "cm"),
                  at = c(-1, 0, 1)
                ),
                row_names_max_width = max_text_width(rownames(expression_matrix), gp = gpar(fontsize = 4))
  )
  draw(ht, 
       annotation_legend_side = "top",
       heatmap_legend_side = "top",
       merge_legend = TRUE)
  dev.off()
  
  # PDF版本
  pdf(heatmap_pdf, width = 6.5, height = 21)
  draw(ht, 
       annotation_legend_side = "top",
       heatmap_legend_side = "top",
       merge_legend = TRUE)
  dev.off()
  
  message("Heatmap saved to: ", heatmap_png, " and ", heatmap_pdf)
}

#------------------- 主执行流程 --------------------#
main <- function(seurat_rds,
                 color_palette,
                 species = "human",
                 pathway_df,
                 top_n,
                 base_dir,
                 output_dir,
                 output_plot_dir
                 ) {
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_plot_dir,showWarnings = FALSE, recursive = TRUE)
  # 1. 读取Seurat对象
  message("\nReading Seurat object from: ", seurat_rds)
  seu <- readRDS(seurat_rds)

  # 2. 生成UMAP图并获取/创建配色方案
  message("\nGenerating UMAP plot...")
  color_palette <- generate_umap_plot(seu, 
                                      color_palette = color_palette, 
                                      output_dir = output_dir,
                                      output_plot_dir = output_plot_dir)

  # 3. 运行GSVA分析
  message("\nRunning GSVA analysis...")
  run_gsva_analysis_re(seurat_rds = seurat_rds,
                    pathway_df = pathway_df,
                    species = species,
                    top_n = top_n,
                    color_palette = color_palette,
                    output_dir = output_dir,
                    output_plot_dir = output_plot_dir
                    )

  # 4. 生成整合热图
  message("\nGenerating integrated heatmap...")
  generate_integrated_heatmap(base_dir = base_dir,
                              color_palette = color_palette,
                              output_dir = output_dir,
                              output_plot_dir = output_plot_dir)

  message("\nAnalysis completed successfully! All results saved to: ", output_dir)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")){
  # 构建配置参数
  config <- list(
    species = snakemake@params$species,
    top_n = snakemake@params$top_n,
    base_dir = snakemake@params$base_dir,
    output_dir = snakemake@params$output_dir,
    output_plot_dir =  snakemake@params$output_plot_dir 
  )

  # 运行主流程
  main(
    seurat_rds = snakemake@input$seurat,
    pathway_df = snakemake@input$pathway_df,
    color_palette = snakemake@input$color_palette,
    species = config$species,
    top_n = config$top_n,
    base_dir = config$base_dir,
    output_dir = config$output_dir,
    output_plot_dir = config$output_plot_dir
  )
}