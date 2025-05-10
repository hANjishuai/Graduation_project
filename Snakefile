# Snakefile

# 加载配置文件
configfile:"config.yaml"

# 定义通配符
#SUFFIXES = ['dle', 'sle', 'hc']#疾病状态
#SUFFIXES2 = ['C00_Bn_IL4R+','C01_Bn_IGKV3_20+','C04_Bm_AIM2+','C05_Bm_ZEB2+',
#            'C07_B1b_CD5-','C09_Plasma_JCHAIN+','C02_Bn_IGHV3_23+','C06_B1a_CD5+']

# 定义规则
rule rawdataqc:
    input:
        rawdatadir = config["rawdatadir"]
    output:
        # 输出文件
        "figure_out/01.QC/integrated.pdf",
        "result_out/01.QC/2.raw.Rdata"
    script:
        # R 脚本路径
        "pipeline/00_rawdatamerge.R"

#rule scRNA_QC:
#    input:
#        "result_out/01.QC/2.raw.Rdata"
#    output:
#        #输出文件
#        "figure_out/01.QC/0.raw.vlnplot.png",
#        "result_out/01.QC/3.QC.Rdata"
#    script:
#        #R 脚本路径
#        "scripts/01_scRNA.QC.R" 
#
#rule cluster:
#    input:
#        "result_out/01.QC/3.QC.Rdata"
#    output:
#        #输出文件
#        "result_out/02.cluster/cluster.markers.top10.csv",
#        "result_out/01.QC/5.Dle_b1.cluster.Rdata",
#        "figure_out/02.cluster/04.Visualization.cluster.png"
#    script:
#        #R 脚本路径
#        "scripts/02_0_cluster.R"
#
#rule cellmarkerplot:
#    input:
#        "result_out/01.QC/5.Dle_b1.cluster.Rdata"
#    output:
#        #输出文件
#        "figure_out/03.celltype/0.manual_all_marker_heatmap.pdf",
#        "figure_out/03.celltype/1.manual_single_marker_heatmap.pdf"
#    script:
#        #R脚本路径
#        "scripts/03_0_cellmarkerplot.R"
#
#rule cellannobysingleR:
#    input:
#        "result_out/01.QC/5.Dle_b1.cluster.Rdata"
#    output:
#        #输出文件
#        "figure_out/03.celltype/2.cellplot.pdf",
#        "result_out/01.QC/6.Dle_b1.celltype.Rdata"
#    script:
#        #R脚本路径
#        "scripts/04_0_cellanno.R"
#
#rule cellannobyman:
#    input:
#        "result_out/01.QC/6.Dle_b1.celltype.Rdata",
#        "result_out/03.celltype/0.singleR_cellplot.csv"
#    output:
#        #输出文件
#        "figure_out/03.celltype/3.cellplot.manual_L1.pdf",
#        "result_out/01.QC/7.Dle_b1.celltype.check.Rdata"
#    script:
#        #R脚本路径
#        "scripts/04_1_cellanno_check.R"
#
#rule subsetcell:
#    input:
#        "result_out/01.QC/7.Dle_b1.celltype.check.Rdata"
#    output:
#        #输出文件
#        "result_out/01.QC/8.subset_Myeloid_cells.Rdata"
#    script:
#        #R脚本路径
#        "scripts/05_0_subsetcell.R"
#
#rule subset:
#    input:
#        "result_out/01.QC/8.subset_Myeloid_cells.Rdata"
#    output:
#        #输出文件
#        "figure_out/04.subcell/01.QC/4.Elbowplot.pdf",
#        "result_out/04.subcell/02.cluster/1.Myeloid_cells.cluster.Rdata"
#    script:
#        #R脚本路径
#        "scripts/06_0_subset.R"
#
#rule anno_marker:
#    input:
#        imarker = config["anno_level2"],
#        iseu = "result_out/04.subcell/02.cluster/1.Myeloid_cells.cluster.Rdata"
#    output:
#        omarker_dotplot_whole = "figure_out/04.subcell/02.celltype/marker_dotplot_whole.pdf",
#        omarker_dotplot_bypage = "figure_out/04.subcell/02.celltype/marker_dotplot_bypage.pdf"
#    params:
#        ispecies = config["species"],
#        pnsheet = config["anno_marker_table"]
#    shell:
#        """
#        Rscript scripts/08_0_marker_anno_find.R \
#            --ispecies {params.ispecies} \
#            --imarker {input.imarker} \
#            --iseu {input.iseu} \
#            --pnsheet {params.pnsheet} \
#            --omarker_dotplot_whole {output.omarker_dotplot_whole} \
#            --omarker_dotplot_bypage {output.omarker_dotplot_bypage}
#        """
#
#rule subcellanno:
#    input:
#        "result_out/04.subcell/02.cluster/1.Myeloid_cells.cluster.Rdata"
#    output:
#        "figure_out/04.subcell/04.annotation/1.cellplot.manual_L2.pdf",
#        "result_out/04.subcell/02.cluster/2.Myeloid_cells.cluster.anno.Rdata"
#    script:
#        #R脚本路径
#        "scripts/07_0_subcellanno.R"
#
#rule plot:
#    input:
#        iseu = "result_out/01.QC/7.Dle_b1.celltype.check.Rdata",
#        isubseu = "result_out/04.subcell/02.cluster/2.Myeloid_cells.cluster.anno.Rdata"
#    output:
#        ocolorbar = "result_out/04.subcell/04.annotation/colorbar.RData"
#    params:
#        pfig_path =config['pfig_path']
#    shell:
#        """
#        Rscript scripts/09_0_plot_all.R \
#        --iseu {input.iseu} \
#        --isubseu {input.isubseu} \
#        --pfig_path {params.pfig_path} \
#        --ocolorbar {output.ocolorbar} \
#        """

