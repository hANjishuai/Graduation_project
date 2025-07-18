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
  library(tidyverse)
  library(ggprism)
  library(glue)
  library(cli)
  library(purrr)
  library(readr)
  library(stringr)
  library(jsonlite)
  library(httr)
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

#------------------- 基因集管理核心函数 -------------------#
get_kegg_species_code <- function(species) {
  switch(
    species,
    "human" = "hsa",
    "mouse" = "mmu",
    "rat" = "rno",
    "hsa"  # 默认是人类
  )
}

get_go_gene_set <- function(go_id, species = "human") {
  species_code <- switch(
    species,
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
    "rat" = "Rattus norvegicus",
    "Homo sapiens"
  )
  
  # 1. 尝试通过AnnotationDbi获取
  if (species == "human") {
    gene_ids <- tryCatch({
      get(go_id, org.Hs.egGO2ALLEGS)
    }, error = function(e) NULL)
  } else if (species == "mouse") {
    gene_ids <- tryCatch({
      get(go_id, org.Mm.egGO2ALLEGS)
    }, error = function(e) NULL)
  } else if (species == "rat") {
    gene_ids <- tryCatch({
      get(go_id, org.Rn.egGO2ALLEGS)
    }, error = function(e) NULL)
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
  
  # 2. 尝试通过msigdbr获取
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
  
  # 3. 尝试通过GO API获取
  log_warning(glue("Could not retrieve GO gene set {go_id} through msigdbr. Trying API..."))
  go_url <- glue("https://api.geneontology.org/api/bioentity/function/{go_id}/genes?rows=1000")
  
  tryCatch({
    go_response <- fromJSON(go_url)
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
  }, error = function(e) {
    log_warning(glue("GO API failed: {e$message}"))
  })
  
  return(NULL)
}

get_kegg_gene_set <- function(kegg_id, species = "human") {
  species_code <- switch(
    species,
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
    "rat" = "Rattus norvegicus",
    "Homo sapiens"
  )
  
  kegg_species_code <- get_kegg_species_code(species)
  full_kegg_id <- paste0(kegg_species_code, kegg_id)
  
  # 1. 尝试通过msigdbr获取
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
  
  # 2. 尝试通过KEGG API获取
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
    gene_lines <- str_trim(kegg_genes_response)
    gene_ids <- sapply(strsplit(gene_lines, "\t"), function(x) x[2])
    gene_ids <- unique(gsub(paste0(kegg_species_code, ":"), "", gene_ids))
    
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
  
  # 3. 尝试通过clusterProfiler获取
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
  
  return(NULL)
}

get_gene_set <- function(pathway_id, species = "human") {
  # 处理不同格式的ID
  clean_id <- pathway_id
  
  # 处理GO ID格式
  if (grepl("^GO:", clean_id) || grepl("^GO\\d+$", clean_id)) {
    # 提取数字部分
    go_id <- gsub("^GO:|^GO", "", clean_id)
    go_id <- paste0("GO:", go_id)
    return(get_go_gene_set(go_id, species))
  }
  
  # 处理KEGG ID格式
  if (grepl("^hsa", clean_id) || grepl("^\\d+$", clean_id)) {
    # 提取数字部分
    kegg_id <- gsub("^hsa", "", clean_id)
    return(get_kegg_gene_set(kegg_id, species))
  }
  
  log_error(glue("Unsupported pathway ID format: {pathway_id}"))
  return(NULL)
}

#------------------- GMT文件处理函数 -------------------#
load_custom_gene_sets <- function(custom_tsv) {
  log_info(glue("Loading custom gene sets from {custom_tsv}"))
  
  custom_df <- read.table(custom_tsv,header = T)
  
  # 检查列名
  if (!all(c("Pathway", "Gene") %in% colnames(custom_df))) {
    stop("Custom TSV must contain 'Pathway' and 'Gene' columns")
  }
  
  # 按通路分组
  gene_sets <- custom_df %>%
    group_by(Pathway) %>%
    summarise(genes = list(unique(Gene))) %>%
    mutate(
      db = "Custom",
      id = Pathway,
      description = Pathway
    ) %>%
    dplyr::select(db, id, description, genes)
  
  log_success(glue("Loaded {nrow(gene_sets)} custom gene sets"))
  return(gene_sets)
}

convert_to_gmt <- function(gene_sets, output_file) {
  log_info(glue("Converting gene sets to GMT format: {output_file}"))
  
  # 创建输出目录
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  # 写入GMT文件
  con <- file(output_file, "w")
  for (i in 1:nrow(gene_sets)) {
    genes <- paste(gene_sets$genes[[i]], collapse = "\t")
    line <- paste(gene_sets$id[i], gene_sets$description[i], genes, sep = "\t")
    writeLines(line, con)
  }
  close(con)
  
  log_success(glue("Saved GMT file: {output_file}"))
}

build_gene_set_collection <- function(kegg_ids = NULL, 
                                      go_ids = NULL, 
                                      custom_tsv = NULL, 
                                      species = "human",
                                      output_dir = "gene_sets") {
  log_info("Building gene set collection")
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 存储所有基因集的元数据和基因列表
  all_gene_sets <- tibble()
  all_genes_list <- list()
  
  # 1. 处理KEGG基因集
  if (!is.null(kegg_ids)) {
    log_info(glue("Processing {length(kegg_ids)} KEGG pathways"))
    
    for (id in kegg_ids) {
      gene_set <- get_kegg_gene_set(id, species)
      if (!is.null(gene_set)) {
        all_gene_sets <- bind_rows(
          all_gene_sets,
          tibble(
            db = gene_set$db,
            id = gene_set$id,
            description = gene_set$description,
            genes = list(gene_set$genes)
          )
        )
        all_genes_list[[gene_set$description]] <- gene_set$genes
      }
    }
  }
  
  # 2. 处理GO基因集
  if (!is.null(go_ids)) {
    log_info(glue("Processing {length(go_ids)} GO terms"))
    
    for (id in go_ids) {
      gene_set <- get_go_gene_set(id, species)
      if (!is.null(gene_set)) {
        all_gene_sets <- bind_rows(
          all_gene_sets,
          tibble(
            db = gene_set$db,
            id = gene_set$id,
            description = gene_set$description,
            genes = list(gene_set$genes)
          )
        )
        all_genes_list[[gene_set$description]] <- gene_set$genes
      }
    }
  }
  
  # 3. 处理自定义基因集
  if (!is.null(custom_tsv)) {
    custom_sets <- load_custom_gene_sets(custom_tsv)
    all_gene_sets <- bind_rows(all_gene_sets, custom_sets)
    for (i in 1:nrow(custom_sets)) {
      all_genes_list[[custom_sets$description[i]]] <- custom_sets$genes[[i]]
    }
  }
  
  # 4. 保存元数据
  metadata_file <- file.path(output_dir, "gene_set_metadata.csv")
  write_csv(all_gene_sets, metadata_file)
  log_success(glue("Saved gene set metadata: {metadata_file}"))
  
  # 5. 保存GMT文件
  # 5.1 总GMT文件
  total_gmt_file <- file.path(output_dir, "all_gene_sets.gmt")
  convert_to_gmt(all_gene_sets, total_gmt_file)
  
  # 5.2 分类型GMT文件
  # 按数据库类型分组
  db_types <- unique(all_gene_sets$db)
  for (db_type in db_types) {
    db_sets <- all_gene_sets %>% filter(db == db_type)
    db_gmt_file <- file.path(output_dir, glue("{db_type}_gene_sets.gmt"))
    convert_to_gmt(db_sets, db_gmt_file)
  }
  
  log_success(glue("Gene set collection built in {output_dir}"))
  
  return(list(
    metadata = all_gene_sets,
    gene_sets = all_genes_list
  ))
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 解析配置参数
  config <- list(
    custom_tsv = snakemake@input$custom_tsv,
    kegg_ids = snakemake@params$kegg_ids,
    go_ids = snakemake@params$go_ids,
    species = snakemake@params$species,
    output_dir = snakemake@params$output_dir
  )
  
  # 运行基因集构建
  gene_set_collection <- build_gene_set_collection(
    kegg_ids = config$kegg_ids,
    go_ids = config$go_ids,
    custom_tsv = config$custom_tsv,
    species = config$species,
    output_dir = config$output_dir
  )
  
  # 创建完成标记文件
  file.create(snakemake@output[[1]])
  log_success("Gene set collection built successfully")
}
