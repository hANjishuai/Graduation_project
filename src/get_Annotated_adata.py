"""
单细胞数据预处理流程（工程化改造版）

功能概述：
1. 模块化设计，支持命令行参数和Snakemake集成
2. GPU加速支持，自动内存管理
3. 可视化结果自动保存
4. 支持外部映射文件和标记基因文件
5. 细胞亚群智能提取

输入文件格式：
mapping_file.csv:
cluster,cell_type（这一行不要，只是在这说明列内容）
0,Epithelial cells
1,Epithelial cells
...

marker_file.csv:
cell_type,genes
Epithelial cells,KRT5;KRT14
Immune cells,PTPRC;CD45

典型运行示例：
python preprocess.py \
  -i skin_epi.h5ad \
  -o subdata_Epi.h5ad \
  -g leiden_harmony_0.2 \
  -c cell_type_lvl1 \
  -t "Epithelial cells" \
  -m cluster_mapping.csv \
  --marker-file marker_genes.csv \
  -p results/plots

运行监控输出样例：
[12:34:56] 📥 加载数据: skin_epi.h5ad (cells: 50,000)
[12:35:01] ✅ 数据预处理完成 (耗时: 5.2s)
[12:35:15] 🎨 生成3种可视化图表 (保存至: results/plots)
[12:35:20] ✂️ 提取Epithelial细胞亚群 (cells: 35,000)
[12:35:25] 💾 结果保存至: subdata_Epi.h5ad
[12:35:25] 🎉 总耗时: 29.3s
"""

import click
import time
import numpy as np
import scanpy as sc
import rapids_singlecell as rsc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
from typing import Dict, List

# 初始化配置
warnings.filterwarnings("ignore")
plt.rcParams['figure.max_open_warning'] = 0  # 关闭Matplotlib警告

def init_gpu() -> bool:
    """初始化GPU加速环境"""
    try:
        import rmm
        import cupy as cp
        from rmm.allocators.cupy import rmm_cupy_allocator
        
        rmm.reinitialize(
            managed_memory=True,
            pool_allocator=True,
            devices=0
        )
        cp.cuda.set_allocator(rmm_cupy_allocator)
        click.echo("✅ GPU加速已启用 (RAPIDS内存池初始化)")
        return True
    except ImportError:
        click.echo("⚠️ 未检测到GPU支持，使用CPU模式运行")
        return False

def load_cluster_mapping(mapping_file: str) -> Dict[str, str]:
    """加载聚类-ID到细胞类型的映射表"""
    try:
        df = pd.read_csv(mapping_file, sep='\t' if mapping_file.endswith('.tsv') else ',')
        assert df.shape[1] >= 2, "映射文件需要至少包含两列"
        return {str(k): v for k, v in zip(df.iloc[:,0], df.iloc[:,1])}
    except Exception as e:
        raise ValueError(f"映射文件解析失败: {str(e)}")

def load_marker_genes(marker_file: str) -> Dict[str, List[str]]:
    """加载标记基因字典"""
    try:
        df = pd.read_csv(marker_file, sep='\t' if marker_file.endswith('.tsv') else ',')
        return {
            row['cell_type']: row['genes'].split(';') 
            for _, row in df.iterrows()
        }
    except Exception as e:
        raise ValueError(f"标记基因文件解析失败: {str(e)}")

def generate_qc_plots(adata, markers: Dict,cell_type: str, plot_dir: Path):
    """生成质量控制可视化图表"""
    with plt.rc_context({"figure.figsize": (12, 8)}):
        # 点图展示标记基因表达
        marker_list=list(markers.values())
        flattened_list = [gene for sublist in marker_list for gene in sublist]
        sc.pl.dotplot(
            adata, 
            flattened_list, 
            groupby=cell_type,
            standard_scale='var',
            show=False,
            title='Marker Gene Expression'
        )
        plt.savefig(plot_dir/'marker_dotplot.pdf', bbox_inches='tight')
        plt.close()
        
        # t-SNE可视化
        sc.pl.tsne(
            adata,
            color=[cell_type],
            legend_loc='on data',
            frameon=False,
            show=False,
            add_outline=False,
            title='t-SNE Projection'
        )
        plt.savefig(plot_dir/'tsne_projection.pdf', bbox_inches='tight')
        plt.close()
        
        # UMAP可视化
        adata.obsm['X_umap'] = adata.obsm['X_umap_harmony']
        sc.pl.umap(
            adata,
            color=[cell_type],
            legend_loc='on data',
            frameon=False,
            show=False,
            title='UMAP Projection'
        )
        plt.savefig(plot_dir/'umap_projection.pdf', bbox_inches='tight')
        plt.close()

def process_dataset(
    input_path: str,
    output_path: str,
    groupby: str,
    cell_type_col: str,
    target_cell: str,
    mapping_file: str,
    marker_file: str,
    plot_dir: str
) -> None:
    """主处理流程"""
    timer = time.time()
    
    # 阶段1: 数据加载
    click.echo(f"\n📥 加载数据: {Path(input_path).name}")
    adata = sc.read(input_path)
    click.echo(f"✅ 初始数据: {adata.n_obs:,} 个细胞 | {adata.n_vars:,} 个基因")
    
    # 阶段2: GPU加速
    if init_gpu():
        rsc.get.anndata_to_GPU(adata)
        click.echo(f"🔧 数据已转至GPU格式")
    
    # 阶段3: 细胞类型注释
    click.echo(f"\n🗺️ 应用细胞类型映射...")
    cluster_map = load_cluster_mapping(mapping_file)
    adata.obs[cell_type_col] = adata.obs[groupby].map(cluster_map)
    click.echo(f"🔍 细胞类型分布:\n{adata.obs[cell_type_col].value_counts()}")
    
    # 阶段4: 可视化
    click.echo(f"\n🎨 生成质量控制图表...")
    plot_path = Path(plot_dir)
    plot_path.mkdir(parents=True, exist_ok=True)
    markers = load_marker_genes(marker_file)
    generate_qc_plots(adata, markers, cell_type_col, plot_path)
    
    # 阶段5: 亚群提取
    click.echo(f"\n✂️ 提取目标细胞类型: {target_cell}")
    subdata = adata[adata.obs[cell_type_col] == target_cell].copy()
    click.echo(f"✅ 亚群数据: {subdata.n_obs:,} 个细胞")
    
    # 阶段6: 数据保存
    if hasattr(rsc.get, 'anndata_to_CPU'):
        rsc.get.anndata_to_CPU(subdata)
    subdata.obs['condition_logi'] = subdata.obs['sample'].str.contains(r'SLE', regex=True, na=False)
    subdata.obs['condition'] = np.where(subdata.obs['condition_logi'], 'Ctr', 'Exp')
    
    # 数据校验步骤
    if not np.isfinite(subdata.X.data).all():
        click.echo("发现非数值数据，正在清理...")
        subdata.X = np.nan_to_num(subdata.X, nan=0, posinf=0, neginf=0)

    # 确保数据为非负数
    if (subdata.X < 0).sum() > 0:
        click.echo("检测到负值，进行非负处理...")
        subdata.X = subdata.X.clip(min=0)

    # 处理零值问题（避免log计算错误）
    click.echo("处理零值（替换为1e-9）...")
    subdata.X = subdata.X + 1e-9

    subdata.write_h5ad(output_path)
    click.echo(f"\n💾 结果已保存至: {output_path}")
    
    # 性能监控
    total_time = time.time() - timer
    click.echo(f"\n🕒 总运行时间: {total_time:.1f}秒")

@click.command()
@click.option('-i', '--input-path', required=True, help='输入h5ad文件路径')
@click.option('-o', '--output-path', required=True, help='输出h5ad文件路径')
@click.option('-g', '--groupby', default='leiden_harmony_0.2', show_default=True,
             help='原始聚类结果列名')
@click.option('-c', '--cell-type-col', default='cell_type_lvl1', show_default=True,
             help='新注释的细胞类型列名')
@click.option('-t', '--target-cell', default='Epithelial cells', show_default=True,
             help='需要提取的目标细胞类型')
@click.option('-m', '--mapping-file', required=True, 
             help='聚类ID到细胞类型的映射文件 (CSV/TSV)')
@click.option('--marker-file', required=True,
             help='标记基因定义文件 (CSV/TSV)')
@click.option('-p', '--plot-dir', default='qc_plots', show_default=True,
             help='可视化结果输出目录')
def main(input_path, output_path, groupby, cell_type_col, target_cell, mapping_file, marker_file, plot_dir):
    """
    🧬 单细胞数据预处理流程
    
    核心功能：
    1. 基于外部映射文件进行细胞类型注释
    2. 根据标记基因进行质量验证
    3. 自动生成诊断图表
    4. 提取指定细胞亚群(若要提取多个细胞亚群，请在映射文件中统一命名，如“KC_pDC”)
    
    技术亮点：
    • 支持GPU加速，提升大数据处理效率
    • 自动化内存管理，优化资源使用
    • 模块化设计，方便扩展和维护
    • 生产级错误处理，保障流程稳定性
    """
    try:
        process_dataset(
            input_path,
            output_path,
            groupby,
            cell_type_col,
            target_cell,
            mapping_file,
            marker_file,
            plot_dir
        )
    except Exception as e:
        click.echo(f"\n❌ 流程执行失败: {str(e)}", err=True)
        raise

if __name__ == "__main__":
    main()