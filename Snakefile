configfile: "config.yaml"

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
