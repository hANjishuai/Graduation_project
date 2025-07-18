"""
单细胞数据预处理流水线（GPU加速版）
████████████████████████████████████████████████
# 查看完整帮助
python get_filteredParam_for_qc.py --help

# 输出示例：

Options:
  -i, --input TEXT        🖥️ 输入文件路径（h5ad格式）  [required]
  -o, --output TEXT       💾 输出h5ad文件路径  [required]
  -p, --params TEXT       ⚙️ 过滤参数保存路径（JSON格式）  [required]
  --plot_dir TEXT         📊 图表输出目录（默认：./qc_plots）
  --n_genes_max INTEGER   🔧 最大基因数阈值（默认：8000）
  --mt_pct_max INTEGER    🔧 最大线粒体基因百分比（默认：25）
  --min_gene_count INTEGER 🔧 基因最小表达量（默认：3）

████████████████████████████████████████████████
功能描述：
1. GPU内存初始化配置
2. 数据加载与转换到GPU显存
3. 计算质量控制指标（线粒体/核糖体基因）
4. 生成质控可视化图表
5. 保存过滤参数供下游分析使用

输入输出说明：
████████████████████████████████████████████████
输入要求：
- 输入文件：标准h5ad格式单细胞数据
- 必需字段：adata.obs中需包含'sample'分组信息

输出内容：
- 预处理数据：标准h5ad文件（保留原始数据）
- 过滤参数：JSON格式阈值文件
- 质控图表：PNG格式（3散点图 + 3小提琴图）

参数文件示例：
{
    "n_genes_max": 8000,
    "mt_pct_max": 25,
    "min_gene_count": 3
}

典型处理时间：
数据集规模       | 预估耗时
---------------------------------
1-5万细胞      | 2-5分钟
5-10万细胞     | 5-10分钟
10万+细胞      | 10-20分钟
████████████████████████████████████████████████

运行示例：
# 基础用法
python preprocess.py \
  --input raw_data.h5ad \
  --output processed.h5ad \
  --params filter_params.json

# 自定义参数
python preprocess.py \
  --input ./data/raw.h5ad \
  --output ./results/processed.h5ad \
  --params ./params/filters.json \
  --plot_dir ./qc_plots \
  --n_genes_max 7500 \
  --mt_pct_max 20 \
  --min_gene_count 5

# Snakemake集成示例
rule preprocess:
    input: "raw.h5ad"
    output: 
        data="processed.h5ad",
        params="params.json"
    params:
        plot_dir="qc_plots",
        mt_pct=25
    shell:
        "python preprocess.py --input {input} --output {output.data} "
        "--params {output.params} --plot_dir {params.plot_dir} --mt_pct_max {params.mt_pct}"
"""

import click
import time
import json
import scanpy as sc
import cupy as cp
import rapids_singlecell as rsc
import warnings
from pathlib import Path
import matplotlib.pyplot as plt
# GPU内存管理
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

warnings.filterwarnings("ignore")

def init_gpu():
    """初始化GPU内存配置"""
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
              help='⚙️ 过滤参数保存路径（JSON格式）', type=click.Path())
@click.option('--plot_dir', default="./qc_plots",
              help='📊 图表输出目录（默认：./qc_plots）', type=click.Path())
@click.option('--n_genes_max', default=8000, 
              help='🔧 最大基因数阈值（默认：8000）', type=int)
@click.option('--mt_pct_max', default=25,
              help='🔧 最大线粒体基因百分比（默认：25）', type=int)
@click.option('--min_gene_count', default=3,
              help='🔧 基因最小表达量（默认：3）', type=int)
def main(input, output, params, plot_dir, n_genes_max, mt_pct_max, min_gene_count):
    """单细胞数据预处理工具（GPU加速版）"""
    
    # 初始化GPU
    init_gpu()
    
    try:
        total_start = time.time()
        click.echo(f"\n🚀 启动预处理流程 [{time.ctime()}]")
        click.echo(f"📥 输入文件：{Path(input).resolve()}")
        click.echo(f"📤 输出文件：{Path(output).resolve()}")
        click.echo(f"⚙️ 参数文件：{Path(params).resolve()}")

        # === 数据加载阶段 ===
        click.echo("\n🔌 阶段1/3：数据加载...")
        load_start = time.time()
        
        adata = sc.read(input)
        rsc.get.anndata_to_GPU(adata)
        
        load_time = time.time() - load_start
        click.echo(f"✅ 加载完成 | 细胞数：{adata.n_obs:,} | 基因数：{adata.n_vars:,}")
        click.echo(f"⏱️ 加载耗时：{load_time:.1f}秒")

        # === 质控计算阶段 ===
        click.echo("\n🧪 阶段2/3：质控指标计算...")
        qc_start = time.time()
        
        rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="MT")
        rsc.pp.flag_gene_family(adata, gene_family_name="RIBO", gene_family_prefix="RPS")
        rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT", "RIBO"])
        
        qc_time = time.time() - qc_start
        click.echo(f"✅ 质控计算完成 | 新增指标：MT%, RIBO%")
        click.echo(f"⏱️ 计算耗时：{qc_time:.1f}秒")

        # === 可视化阶段 ===
        click.echo(f"\n🎨 阶段3/3：生成质控图表到目录：{plot_dir}")
        viz_start = time.time()
        
        plot_path = Path(plot_dir)
        plot_path.mkdir(parents=True, exist_ok=True)
        
        # 散点图
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_MT", 
                     show=False)
        plt.savefig(str(plot_path/"total_vs_MT.png"))
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", 
                     show=False)
        plt.savefig(str(plot_path/"total_vs_genes.png"))
        
        # 小提琴图
        for metric in ["n_genes_by_counts", "total_counts", "pct_counts_MT"]:
            sc.pl.violin(adata, metric, jitter=0.4, groupby="sample", #这里的sample是开始的样本名列
                        show=False)
            plt.savefig(str(plot_path/f"{metric}_violin.png"))        
        viz_time = time.time() - viz_start
        click.echo(f"✅ 生成图表：{len(list(plot_path.glob('*.png')))}个文件")
        click.echo(f"⏱️ 可视化耗时：{viz_time:.1f}秒")

        # === 保存结果 ===
        click.echo("\n💾 正在保存结果...")
        save_start = time.time()
        
        # 保存参数
        params = Path(params)
        params.mkdir(parents=True, exist_ok=True)
        # 指定要保存的 JSON 文件的路径
        params_file = params/"parameters.json"  # 在目录中创建一个名为 

        with open(params_file, 'w') as f:
            param_dict = {
            "n_genes_max": n_genes_max,
            "mt_pct_max": mt_pct_max,
            "min_gene_count": min_gene_count
            }
            json.dump(param_dict, f, indent=4)
        
        click.echo("✅ 参数已保存，包含以下内容：")
        click.echo(json.dumps(param_dict, indent=4))
        # 保存数据
        adata.write_h5ad(output)
        
        save_time = time.time() - save_start
        click.echo(f"✅ 数据保存完成：{Path(output).stat().st_size/1024/1024:.1f}MB")
        click.echo(f"⏱️ 保存耗时：{save_time:.1f}秒")

        # === 总耗时统计 ===
        total_time = time.time() - total_start
        click.echo(f"\n🎉 处理完成！总耗时：{total_time//60:.0f}分{total_time%60:.1f}秒")

    except Exception as e:
        click.echo(f"\n❌ 错误发生：{str(e)}", err=True)
        raise click.Abort()

if __name__ == "__main__":
    main()
