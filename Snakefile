configfile: "config.yaml"

# seurat亚群指控，分群聚类
rule all_process_DR_Cluster:
    input:
        rds = config["output"]["result"].format(celltype=config["celltype"]),
        vln_raw = config["output_figures"]["vln_raw"].format(celltype=config["celltype"]),
        scatter_raw = config["output_figures"]["scatter_raw"].format(celltype=config["celltype"]),
        vln_qc = config["output_figures"]["vln_qc"].format(celltype=config["celltype"]),
        scatter_qc = config["output_figures"]["scatter_qc"].format(celltype=config["celltype"]),
        elbow = config["output_figures"]["elbow"].format(celltype=config["celltype"])

rule process_DR_Cluster:
    input:
        rds = config["input"]["seurat"].format(celltype=config["celltype"])
    output:
        rds = config["output"]["result"].format(celltype=config["celltype"]),
        vln_raw = config["output_figures"]["vln_raw"].format(celltype=config["celltype"]),
        scatter_raw = config["output_figures"]["scatter_raw"].format(celltype=config["celltype"]),
        vln_qc = config["output_figures"]["vln_qc"].format(celltype=config["celltype"]),
        scatter_qc = config["output_figures"]["scatter_qc"].format(celltype=config["celltype"]),
        elbow = config["output_figures"]["elbow"].format(celltype=config["celltype"])
    params:
        species = config["params"]["species"],
        dims = config["params"]["dims"],
        subsets = config["params"]["subsets"]
    script:
        "pipeline/01_process_DR_Cluster.R"
    
# 单细胞免疫组库可视化流程
rule process_scRepertoire:
    input: 
        rawdatadir = config["input02"]["contig"].format(group=config["group"])
    output:
        output_dir = directory(config["output02"]["result"].format(group=config["group"])),
        visual_output_dir = directory(config["output02"]["figures"].format(group=config["group"]))
    params:
        parallel = config["params02"]["parallel"]
    script:
        "pipeline/02_process_scRepertoire.R"

# B细胞亚群NMF分析流程
# NMF分析流程
rule all_NMF_analysis:
    input:
        rds = config["output03"]["result"].format(celltype=config["celltype"]),
        nmf_model = config["output03"]["nmf_model"].format(celltype=config["celltype"]),
        final_model = config["output03"]["final_model"].format(celltype=config["celltype"]),
        rank_plot = config["output03"]["rank_plot"].format(celltype=config["celltype"])

rule NMF_analys is:
    input:
        rds = config["input03"]["seurat"].format(celltype=config["celltype"])
    output:
        rds = config["output03"]["result"].format(celltype=config["celltype"]),
        nmf_model = config["output03"]["nmf_model"].format(celltype=config["celltype"]),
        final_model = config["output03"]["final_model"].format(celltype=config["celltype"]),
        rank_plot = config["output03"]["rank_plot"].format(celltype=config["celltype"])
    params:
        species = config["params03"]["species"],
        rank_start = config["params03"]["rank_start"],
        rank_stop = config["params03"]["rank_stop"],
        rank_step = config["params03"]["rank_step"],
        nrun = config["params03"]["nrun"],
        topN = config["params03"]["topN"],
        resolution = config["params03"]["resolution"],
        enrich_dir = config["output03"]["enrich_dir"].format(celltype=config["celltype"]),
        # 可选参数：已有的NMF模型
        existing_nmf_model = config["params03"]["existing_nmf_model"],
        existing_final_nmf = config["params03"]["existing_final_nmf"]
    script:
        "pipeline/03_process_NMF_Enrichment.R"

# NMF亚群GSVA分析
rule all_GSVA_analysis:
    input:
        expand("result_out/04_process_runGSVA/{celltype}/GSVA_completed.txt",
                celltype=config["celltype"])

rule run_GSVA_analysis:
    input:
        seurat = config["output03"]["result"].format(celltype="{celltype}"),
        # 依赖所有cluster的富集结果
        enrich_dirs = directory(config["output03"]["enrich_dir"].format(celltype="{celltype}"))
    output:
        touch("result_out/04_process_runGSVA/{celltype}/GSVA_completed.txt")
    params:
        species = config["params04"]["species"],
        top_n = config["params04"]["top_n"],
        output_dir = config["output04"]["gsva_result_dir"],
        output_plot_dir = config["output04"]["gsva_plot_dir"]
    script:
        "pipeline/04_process_runGSVA.R"

# 对以上NMF亚群进行作图
rule all_integrate_analysis:
    input:
        expand("result_out/05_plots_nmfAnno_enrichment/{celltype}/plot_completed.txt",
                celltype=config["celltype"])

rule integrate_analysis:
    input:
        seurat = config["input05"]["seurat"],
        pathway_df = config["input05"]["pathway_df"],
        color_palette = config["input05"]["color_palette"]
    output:
        touch("result_out/05_plots_nmfAnno_enrichment/{celltype}/plot_completed.txt")
    params:
        species = config["params05"]["species"],
        top_n = config["params05"]["top_n"],
        base_dir = config["params05"]["base_dir"],
        output_dir = config["output05"]["output_dir"],
        output_plot_dir = config["output05"]["output_plot_dir"]
    script:
        "pipeline/05_plots_nmfAnno_enrichment.R"



# NMF亚群筛选感兴趣的通路
rule all_build_gene_sets:
    input:
        expand("result_out/06_process_Build_Geneset/{celltype}/gene_sets_completed.txt", 
               celltype=config["celltype"])

rule build_gene_sets:
    input:
        custom_tsv = config["input06"]["custom_tsv"]
    output:
        touch("result_out/06_process_Build_Geneset/{celltype}/gene_sets_completed.txt")
    params:
        kegg_ids = config["params06"]["kegg_ids"],
        go_ids = config["params06"]["go_ids"],
        species = config["params06"]["species"],
        output_dir = config["output06"]["gmt_files_dir"]
    script:
        "pipeline/06_process_Build_Geneset.R"


