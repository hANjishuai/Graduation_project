# 02_process_scRepertoire.R

# 加载依赖包
source("~/R_development/ImmuneRepertoireVis/R/Vis_scRepertoire.R")

#------------------- 核心函数 -------------------#
workflow <- function(config) {
    flow_item <- create_visual_workflow(config) %>%
        add_visual_task(
          func_name = "clonal_quant",
          func = clonal_quant,
          params = list(
            group.by = c("sample", "ID"),
            cloneCall = c("gene", "nt", "aa", "strict"),
            chain = c("IGH", "IGL", "both"),
            scale = c(TRUE, FALSE)
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_abundance",
          func = clonal_abundance,
          params = list(
            group.by = c("sample", "ID"),
            cloneCall = c("gene", "nt", "aa", "strict"),
            chain = c("IGH", "IGL", "both"),
            scale = c(TRUE, FALSE)
          )
        ) %>% 
        add_visual_task(
          func_name = "clonal_compare",
          func = clonal_compare,
          params = list(
            cloneCall = c("nt", "aa"),
            top.clones = seq(4, 10, 2),
            graph = c("alluvial", "area")
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_homeostasis",
          func = clonal_homeostasis,
          params = list(
            cloneCall = c("nt", "aa", "gene", "strict"),
            group.by = c("ID", "sample")
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_proportion",
          func = clonal_proportion,
          params = list(
            cloneCall = c("nt", "aa", "gene", "strict"),
            group.by = c("ID", "sample")
          )
        ) %>% 
        add_visual_task(
          func_name = "percent_aa",
          func = percent_aa,
          params = list(
            group.by = c("ID", "sample"),
            aa.length = seq(20, 30, 2),
            chain = c("IGH", "IGL", "both")
          )
        ) %>%
        add_visual_task(
          func_name = "positional_entropy",
          func = positional_entropy,
          params = list(
            group.by = c("ID", "sample"),
            method = c("shannon", "inv.simpson", "norm.entropy"),
            aa.length = 20,
            chain = c("IGH", "IGL", "both")
          )
        ) %>%
        add_visual_task(
          func_name = "positional_property",
          func = positional_property,
          params = list(
            group.by = c("ID", "sample"),
            method = c("Atchley", "Kidera", "stScales", "tScales", "VHSE"),
            aa.length = 20,
            chain = c("IGH", "IGL", "both")
          )
        ) %>%
        add_visual_task(
          func_name = "viz_genes",
          func = viz_genes,
          params = list(
            x.axis = c("IGHV", "IGHJ"),
            y.axis = "IGLV",
            field = c("DLE", "SLE"),
            plot = "heatmap",
            scale = c(TRUE, FALSE)
          )
        ) %>%
        add_visual_task(
          func_name = "percent_genes",
          func = percent_genes,
          params = list(
            group.by = c("ID", "sample"),
            gene = c("Vgene", "Jgene"),
            chain = c("IGH", "IGL")
          )
        ) %>%
        add_visual_task(
          func_name = "pca_percent_genes",
          func = pca_percent_genes,
          params = list(
            group.by = c("ID", "sample"),
            gene = c("Vgene", "Jgene"),
            chain = c("IGH", "IGL")
          )
        ) %>%
        add_visual_task(
          func_name = "percent_vj",
          func = percent_vj,
          params = list(
            group.by = c("ID", "sample"),
            chain = c("IGH", "IGL")
          )
        ) %>%
        add_visual_task(
          func_name = "pca_percent_vj",
          func = pca_percent_vj,
          params = list(
            group.by = c("ID", "sample"),
            chain = c("IGH", "IGL")
          )
        ) %>%
        add_visual_task(
          func_name = "percent_kmer",
          func = percent_kmer,
          params = list(
            group.by = c("ID", "sample"),
            cloneCall = c("nt", "aa"),
            chain = c("IGH", "IGL", "both"),
            motif.length = seq(3, 9, 3),
            top.motifs = 25
          )
        ) %>%
        # Comparing Clonal Diversity and Overlap
        add_visual_task(
          func_name = "clonal_diversity",
          func = clonal_diversity,
          params = list(
            x.axis = "sample",
            cloneCall = c("nt", "aa"),
            chain = c("IGH", "IGL", "both"),
            metrics = c("shannon", "inv.simpson", "ACE", "norm.entropy", "gini.simpson", "chao1")
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_cluster",
          func = clonal_cluster,
          params = list(
            sequence = c("aa", "nt"),
            chain = c("IGH", "IGL"),
            group.by = c("ID", "sample")
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_rarefaction",
          func = clonal_rarefaction,
          params = list(
            cloneCall = c("nt", "aa"),
            chain = c("IGH", "IGL", "both"),
            group.by = c("ID", "sample"),
            plot.type = 1:3,
            hill.numbers = 0:2
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_size_distribution",
          func = clonal_size_distribution,
          params = list(
            cloneCall = c("aa"),
            chain = c("IGH", "IGL"),
            method = c("ward.D2"),
            group.by = c("ID")
          )
        ) %>%
        add_visual_task(
          func_name = "clonal_overlap",
          func = clonal_overlap,
          params = list(
            cloneCall = c("aa", "nt", "strict"),
            chain = c("IGH", "IGL", "both"),
            method = c("overlap", "morisita", "jaccard", "cosine", "raw"),
            group.by = c("ID")
          )
        )
    return(flow_item)
}
#------------------- 执行主函数 -------------------#
process_scRepertoire_pipeline <- function(
    rawdatadir = "data/Skin_IR",
    output_dir = "results/",
    visual_output_dir = "results/visual",
    parallel = TRUE) {
    
    # 环境初始化
    initiate()
    
    # 数据处理阶段
    processed_data <- process_immune_repertoire(
        rawdatadir = rawdatadir,
        output_dir = output_dir
        )
    message("\n data processed. Output saved to: ", output_dir)
    
    # 初始化配置
    config <- VisualizationConfig$new()

    # 创建可视化流程
    workflow <- workflow(config)

    # 执行流程
    results <- run_visual_pipeline(
        processed_data$data,
        workflow,
        output_dir = visual_output_dir,
        parallel = parallel
        )
    
    message("\n✓ Pipeline completed. Visual output saved to: ", visual_output_dir)
} 

# Snakemake集成入口
if (exists("snakemake")) {
    params <- snakemake@params

    # 显式转换并匹配参数名称
    config <- list(
        parallel = as.logical(params$parallel)
    )

    # 确保目录存在
    dir.create(dirname(snakemake@output$output_dir), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(snakemake@output$visual_output_dir), recursive = TRUE, showWarnings = FALSE)
    
    # 调用主函数
    process_scRepertoire_pipeline(
      rawdatadir = snakemake@input$rawdatadir,
      output_dir = snakemake@output$output_dir,
      visual_output_dir = snakemake@output$visual_output_dir,
      parallel = config$parallel[[1]]
    ) 
}
