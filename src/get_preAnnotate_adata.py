"""
单细胞注释流程（模块化工程版）
支持多种注释方法和灵活的参数配置
参数	类型	默认值	说明
-i/--input	路径	必填	输入h5ad文件路径
-o/--output	路径	必填	输出h5ad文件路径
-g/--groupby	字符串	leiden_harmony_0.2	聚类分组字段
-s/--species	选项	human	物种信息 (human/mouse)
-c/--color-features	列表	Epithelial cells	可视化特征（可多个）
-m/--marker-file	路径	可选	标记基因文件（CSV/TSV）
-d/--diff-genes-file	路径	可选	差异基因文件
--model	字符串	Immune_All_Low.pkl	CellTypist模型名称
-p/--plot-dir	路径	annotation_plots	绘图输出目录

# 输入文件格式
标记基因文件 (CSV/TSV)：
cell_type,genes
T cells,CD3D;CD4
B cells,CD19;MS4A1

差异基因文件：
LYZ
ACTB
S100A6
S100A4
CST3

# 典型工作流
# 基本运行
python annotate.py -i input.h5ad -o annotated.h5ad

# 完整参数运行
python annotate.py -i input.h5ad -o output.h5ad \
  -g leiden_0.5 -s human \
  -m markers.csv -d diff_genes.txt \
  -p results/plots --model Immune_All_High.pkl

# 注意事项
    GPU加速：
        需要NVIDIA显卡和匹配的CUDA版本
        自动检测GPU可用性，无需手动配置

    内存管理：
        处理大数据时建议分配至少32GB内存
        可使用--subset参数进行测试运行

    模型下载：
        首次使用CellTypist模型时会自动下载
        模型缓存位置：~/.celltypist/

"""

import click
import time
import scanpy as sc
import rapids_singlecell as rsc
import celltypist as ct
import decoupler as dc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
from typing import Dict, List, Optional

# 初始化配置
warnings.filterwarnings("ignore")
#sc.settings.autosave = True
#sc.settings.verbosity = 2
plt.rcParams['figure.max_open_warning'] = 0  # 关闭绘图警告

def init_gpu():
    """初始化GPU加速配置"""
    try:
        import cupy as cp
        import rapids_singlecell as rsc
        from rmm.allocators.cupy import rmm_cupy_allocator
        import rmm
        
        rmm.reinitialize(
            managed_memory=False,
            pool_allocator=False,
            devices=0,
        )
        cp.cuda.set_allocator(rmm_cupy_allocator)
        print("✅ GPU加速已启用")
        return True
    except ImportError:
        print("⚠️ 未检测到GPU支持，将使用CPU运行")
        return False

def load_markers(marker_file: str) -> Dict[str, List[str]]:
    """从文件加载标记基因字典"""
    df = pd.read_csv(marker_file, sep='\t' if marker_file.endswith('.tsv') else ',')
    return {
        row['cell_type']: row['genes'].split(';')
        for _, row in df.iterrows()
    }

def manual_annotation(
    adata: sc.AnnData,
    marker_file: str,
    groupby: str,
    plot_dir: Path
) -> sc.AnnData:
    """基于标记基因的手动注释"""
    markers = load_markers(marker_file)
    
    # 生成点图
    plt.figure(figsize=(12, 6))
    sc.pl.dotplot(adata, markers, groupby=groupby, show=False)
    plt.savefig(str(f"{plot_dir}/manual_dotplot.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def auto_annotation(
    adata: sc.AnnData,
    model_name: str,
    groupby: str,
    plot_dir: Path
) -> sc.AnnData:
    """CellTypist自动注释"""
    if not Path(model_name).exists():
        ct.models.download_models(model=[model_name])
    
    model = ct.models.Model.load(model_name)
    predictions = ct.annotate(
        adata, 
        model=model,
        majority_voting=True,
        over_clustering=groupby
    )
    
    adata = predictions.to_adata()
    plt.figure(figsize=(12, 12))
    sc.pl.tsne(adata, color="majority_voting", ncols=2, show=False,legend_loc="on data")
    plt.savefig(plot_dir/"auto_annotation.pdf", dpi=300)
    plt.close()
    
    return adata

def enrichment_annotation(
    adata: sc.AnnData,
    species: str,
    groupby: str,
    color_features: List[str],
    plot_dir: Path
) -> sc.AnnData:
    """富集分析注释"""
    markers = dc.get_resource("PanglaoDB", organism=species)
    markers = markers[markers["canonical_marker"]].drop_duplicates(["cell_type", "genesymbol"])
    #建议用户在运行run_mlm之前，将数据显式转换为CPU格式，并确保数据类型正确。这可能涉及到在函数内部添加数据转换步骤，或者在调用run_mlm之前手动转换数据，以避免类型不匹配的问题
    rsc.get.anndata_to_CPU(adata)
    dc.run_mlm(adata, net=markers, source="cell_type", weight=None,target="genesymbol", use_raw=False)
    acts = dc.get_acts(adata, obsm_key="mlm_estimate")
    
    # 可视化特征分布
    sc.pl.tsne(acts, color=color_features, wspace=0.5, ncols=2, show=False)
    plt.savefig(plot_dir/"enrichment.png", dpi=300)
    plt.close()
    
    # 分配注释
    mean_enr = dc.summarize_acts(acts, groupby=groupby)
    annotation_dict = dc.assign_groups(mean_enr)
    adata.obs["dc_anno"] = [annotation_dict[clust] for clust in adata.obs[groupby]]
    
    sc.pl.tsne(adata, color=["dc_anno"], wspace=0.5, ncols=2, show=False,legend_loc="on data")
    plt.savefig(plot_dir/"enrichment_dc_anno.png", dpi=300)
    plt.close()
    return adata

def diff_gene_analysis(
    adata: sc.AnnData,
    gene_file: str,
    groupby: str,
    plot_dir: Path
) -> sc.AnnData:
    """差异基因分析"""
    with open(gene_file) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    sc.tl.rank_genes_groups(adata, groupby=groupby)
    sc.tl.filter_rank_genes_groups(adata, min_fold_change=1.5)
    
    # 差异基因点图
    sc.pl.rank_genes_groups_dotplot(adata, groupby=groupby, n_genes=5, show=False)
    plt.savefig(plot_dir/"diff_genes.png", dpi=300)
    plt.close()
    
    # 基因表达可视化
    sc.pl.tsne(adata, color=genes[:5]+[groupby], legend_loc="on data", ncols=3, show=False)
    plt.savefig(plot_dir/"gene_expression.png", dpi=300)
    plt.close()
    
    return adata

@click.command()
@click.option('--input', '-i', required=True, help='输入h5ad文件路径')
@click.option('--output', '-o', required=True, help='输出h5ad文件路径')
@click.option('--groupby', '-g', default='leiden_harmony_0.2', show_default=True, 
             help='聚类分组字段')
@click.option('--species', '-s', default='human', type=click.Choice(['human', 'mouse']),
             help='物种信息')
@click.option('--color-features', '-c', multiple=True, 
             default=['Epithelial cells', 'Endothelial cells'],
             help='可视化特征列表')
@click.option('--marker-file', '-m', type=click.Path(exists=True),
             help='标记基因文件路径（CSV/TSV）')
@click.option('--diff-genes-file', '-d', type=click.Path(exists=True),
             help='差异基因文件路径（每行一个基因）')
@click.option('--model', default='Immune_All_Low.pkl', 
             help='CellTypist预训练模型名称')
@click.option('--plot-dir', '-p', default='annotation_plots',
             help='绘图输出目录')
def main(
    input: str,
    output: str,
    groupby: str,
    species: str,
    color_features: List[str],
    marker_file: Optional[str],
    diff_genes_file: Optional[str],
    model: str,
    plot_dir: str
):
    """
    🧬 单细胞注释流程

    📌 典型运行时间:
    - 10k细胞: 10-15分钟
    - 50k细胞: 20-30分钟 
    - 100k细胞: 40-60分钟（建议GPU）

    🚀 示例命令:
    python annotate.py -i input.h5ad -o output.h5ad \
    -g leiden_0.5 -s human -m markers.csv -d genes.txt \
    -p results/plots --model Immune_All_High.pkl
    """
    start_time = time.time()
    plot_path = Path(plot_dir)
    plot_path.mkdir(parents=True, exist_ok=True)
    
    # GPU初始化
    if init_gpu():
        import rapids_singlecell as rsc
    
    try:
        click.echo(f"\n📥 加载数据: {Path(input).resolve()}")
        adata = sc.read(input)
        
        if 'anndata_to_GPU' in dir(rsc.get):
            rsc.get.anndata_to_GPU(adata)
        click.echo(f"✅ 加载完成 | 细胞数: {adata.n_obs:,}")

        # 执行注释流程
        if marker_file:
            click.echo("\n🔍 执行手动注释...")
            adata = manual_annotation(adata, marker_file, groupby, plot_path)
            
        click.echo("\n🤖 执行自动注释...")
        adata = auto_annotation(adata, model, groupby, plot_path)
        
        click.echo("\n📊 执行富集分析...")
        adata = enrichment_annotation(adata, species, groupby, color_features, plot_path)
        
        click.echo(f"\n💾 保存结果到: {output}")
        rsc.get.anndata_to_CPU(adata)
        adata.write_h5ad(output)       

        if diff_genes_file:
            click.echo("\n🧬 执行差异分析...")
            adata = diff_gene_analysis(adata, diff_genes_file, groupby, plot_path)
        
    except Exception as e:
        click.echo(f"\n❌ 流程执行失败: {str(e)}", err=True)
        raise
    
    click.echo(f"\n✅ 流程完成! 总耗时: {(time.time()-start_time)/60:.1f}分钟")

if __name__ == "__main__":
    main()
