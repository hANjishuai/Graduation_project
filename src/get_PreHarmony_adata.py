"""
单细胞分析流程（GPU加速版）
████████████████████████████████████████████████
python get_PreHarmony_adata --help

# 输出示例：
Options:
  -i, --input TEXT         🖥️ 输入文件路径（h5ad格式）  [required]
  -o, --output TEXT        💾 输出h5ad文件路径  [required]
  -p, --params TEXT        ⚙️ 参数文件路径（JSON格式）  [required]
  --plot_dir TEXT          📊 图表输出目录（默认：./plots）
  --n_pcs INTEGER          🔧 PCA成分数（默认：100）
  --n_top_genes INTEGER    🔧 高变基因数量（默认：5000）
  --copy_params            🔗 是否复制参数文件到输出目录

████████████████████████████████████████████████
功能描述：
1. 数据加载与GPU内存初始化
2. 基于QC参数的细胞/基因过滤
3. 双联体检测(Scrublet)
4. 数据标准化与高变基因选择
5. PCA分析及方差解释率可视化

输入输出说明：
████████████████████████████████████████████████
输入要求：
- 输入文件：标准h5ad格式单细胞数据
- 必需字段：
  adata.obs中需包含'sample'分组信息
  adata.var中需有MT/RIBO基因标记
- 参数文件：JSON格式，包含以下键：
  "n_genes_max", "mt_pct_max", "min_gene_count"

输出内容：
- 分析结果：标准h5ad文件（包含全部分析结果）
- 参数文件：与输入参数相同（自动复制）
- 可视化图表：PCA方差解释率图（PDF格式）

典型处理时间（NVIDIA A100）：
数据集规模       | 预估耗时
---------------------------------
1-5万细胞      | 8-15分钟
5-10万细胞     | 15-30分钟
10万+细胞      | 30-60分钟
████████████████████████████████████████████████

运行示例：
# 基础用法
python sc_analysis.py \
  -i raw_data.h5ad \
  -o processed.h5ad \
  -p params.json \
  --plot_dir ./plots

# 自定义分析参数
python sc_analysis.py \
  -i ./data/raw.h5ad \
  -o ./results/processed.h5ad \
  -p ./params/analysis.json \
  --plot_dir ./qc_plots \
  --n_pcs 50 \
  --n_top_genes 3000

# Snakemake集成示例
rule single_cell:
    input: 
        data="raw.h5ad",
        params="params.json"
    output: 
        result="processed.h5ad",
        plot="plots/pca_variance.pdf"
    params:
        n_pcs=100,
        n_top_genes=5000
    threads: 8
    resources:
        mem_mb=64000,
        gpu=1
    shell:
        "python sc_analysis.py --input {input.data} "
        "--output {output.result} --params {input.params} "
        "--plot_dir {output.plot} --n_pcs {params.n_pcs} "
        "--n_top_genes {params.n_top_genes}"
"""

import click
import time
import json
import scanpy as sc
import cupy as cp
import matplotlib.pyplot as plt
import rapids_singlecell as rsc
import warnings
from pathlib import Path
from shutil import copyfile

# GPU内存管理
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

warnings.filterwarnings("ignore")

def init_gpu():
    """初始化GPU配置"""
    rmm.reinitialize(
        managed_memory=False,
        pool_allocator=False,
        devices=0,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

@click.command()
@click.option('--input', '-i', required=True,
              help='🖥️ 输入文件路径（h5ad格式）', type=click.Path(exists=True))
@click.option('--output', '-o', required=True,
              help='💾 输出h5ad文件路径', type=click.Path())
@click.option('--params', '-p', required=True,
              help='⚙️ 参数文件路径（JSON格式）', type=click.Path(exists=True))
@click.option('--plot_dir', default="./plots",
              help='📊 图表输出目录（默认：./plots）', type=click.Path())
@click.option('--n_pcs', default=100, 
              help='🔧 PCA成分数（默认：100）', type=int)
@click.option('--n_top_genes', default=5000,
              help='🔧 高变基因数量（默认：5000）', type=int)
@click.option('--copy_params', is_flag=True,
              help='🔗 是否复制参数文件到输出目录')
def main(input, output, params, plot_dir, n_pcs, n_top_genes, copy_params):
    """单细胞分析流程命令行接口"""
    # 初始化GPU
    init_gpu()
    
    try:
        # 初始化流程
        start_time = time.time()
        Path(plot_dir).mkdir(parents=True, exist_ok=True)
        
        # === 参数处理 ===
        with open(params) as f:
            qc_params = json.load(f)
        
        click.echo("\n⚙️ 当前分析参数：")
        click.echo(f"|-- 质量控制参数 --")
        click.echo(f"| 最大基因数: {qc_params['n_genes_max']}")
        click.echo(f"| 线粒体阈值: {qc_params['mt_pct_max']}%")
        click.echo(f"| 最小基因计数: {qc_params['min_gene_count']}")
        click.echo(f"|-- 分析参数 --")
        click.echo(f"| PCA成分数: {n_pcs}")
        click.echo(f"| 高变基因数: {n_top_genes}")

        # === 数据加载 ===
        click.echo("\n📥 正在加载数据...")
        adata = sc.read(input)
        rsc.get.anndata_to_GPU(adata)
        click.echo(f"✅ 初始数据维度：{adata.shape}")

        # === 质量控制 ===
        click.echo("\n🔍 执行质量控制...")
        adata = adata[adata.obs["n_genes_by_counts"] < qc_params["n_genes_max"]]
        adata = adata[adata.obs["pct_counts_MT"] < qc_params["mt_pct_max"]]
        rsc.pp.filter_genes(adata, min_count=qc_params["min_gene_count"])
        click.echo(f"✅ 过滤后维度：{adata.shape}")

        # === 双联体检测 ===
        click.echo("\n🕵️ 执行双联体检测...")
        rsc.pp.scrublet(adata, batch_key="sample")

        # === 数据预处理 ===
        click.echo("\n⚡ 执行预处理...")
        adata.layers["counts"] = adata.X.copy()
        rsc.pp.normalize_total(adata, target_sum=1e4)
        rsc.pp.log1p(adata)

        # 高变基因选择
        adata.raw = adata
        rsc.pp.highly_variable_genes(
            adata, 
            n_top_genes=n_top_genes,
            flavor="cell_ranger",
            batch_key="sample"
        )
#        adata.raw = adata
        rsc.pp.filter_highly_variable(adata)

        # 回归校正
        rsc.pp.regress_out(adata, keys=["total_counts", "pct_counts_MT"])
        rsc.pp.scale(adata, max_value=10)

        # === PCA分析 ===
        click.echo(f"\n📊 执行PCA分析（n_pcs={n_pcs}）...")
        rsc.tl.pca(adata, n_comps=n_pcs)
        
        # 保存方差解释率图
        plot_path = Path(plot_dir)
        plot_path.mkdir(parents=True, exist_ok=True)
        rsc.get.anndata_to_CPU(adata)
        sc.pl.pca_variance_ratio(
            adata, 
            log=True, 
            n_pcs=n_pcs,
            show=False
        )
        plt.savefig(str(plot_path/"pca_variance.pdf"))
        # === 结果保存 ===
        adata.write_h5ad(output)
        if copy_params:
            copyfile(params, Path(output).parent / "used_params.json")
        
        # 时间统计
        total_time = time.time() - start_time
        click.echo(f"\n✅ 分析完成！总耗时：{total_time//60:.0f}分{total_time%60:.1f}秒")

    except Exception as e:
        click.echo(f"\n❌ 错误发生：{str(e)}", err=True)
        raise click.Abort()

if __name__ == "__main__":
    main()
