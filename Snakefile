configfile: "config.yaml"

rule all_process_DR_Cluster:
    input:
        config["output"]["result"],
        expand(
            "figure_out/01_process_DR_Cluster/celltype_ker/{plot}.png",
            plot=["vln_raw", "scatter_raw", "vln_qc", "scatter_qc", "elbow"]
        )

rule process_DR_Cluster:
    input:
        rds = config["input"]["seurat"]
    output:
        rds = config["output"]["result"],
        figures = expand(
            "figure_out/01_process_DR_Cluster/celltype_ker/{plot}.png",
            plot=["vln_raw", "scatter_raw", "vln_qc", "scatter_qc", "elbow"]
        )
    params:
        species = config["params"]["species"],
        dims = config["params"]["dims"],
        subsets = config["params"]["subsets"]
    script:
        "pipeline/01_process_DR_Cluster.R"
