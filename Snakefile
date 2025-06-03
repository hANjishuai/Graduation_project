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
