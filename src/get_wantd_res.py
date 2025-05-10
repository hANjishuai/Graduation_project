"""
单细胞分析流程（参数化+工程化版本）
支持通过起始/结束/步长生成分辨率参数
Usage: sc_analysis.py [OPTIONS]

  单细胞分析主流程

  典型运行时间:
  - 10k细胞: 15-20分钟
  - 50k细胞: 30-45分钟
  - 100k细胞: 1-2小时

  示例命令:
  python sc_analysis.py -i input.h5ad -o output.h5ad \
  --res-start 0.0 --res-end 1.0 --res-step 0.2 \
  -n 40 -g leiden_harmony_0.2 -p results/plots

Options:
  -i, --input TEXT      输入h5ad文件路径  [required]
  -o, --output TEXT     输出h5ad文件路径  [required]
  --res-start FLOAT     Leiden分辨率起始值  [default: 0.0]
  --res-end FLOAT       Leiden分辨率结束值（包含）  [default: 1.0]
  --res-step FLOAT      Leiden分辨率步长  [default: 0.2]
  -n, --n-pcs INTEGER   t-SNE使用的PC数量  [default: 40]
  -g, --groupby TEXT    差异表达分析的分组字段  [default: leiden_harmony_0.2]
  -p, --plot-dir TEXT   绘图输出目录  [default: plots]
  --help                Show this message and exit.

内存管理：
# 处理完成后释放内存
adata = None
import gc
gc.collect()

并行计算
# 设置Scanpy并行参数
sc.settings.n_jobs = 4  # 使用4个CPU核心

自定义绘图样式
# 在脚本开头设置全局绘图参数
sc.settings.set_figure_params(
    dpi=300,
    facecolor='white',
    frameon=False,
    fontsize=8
)
"""

import time
from pathlib import Path
import click
import numpy as np
import cupy as cp
import scanpy as sc
import rapids_singlecell as rsc
import matplotlib
matplotlib.use('Agg')  # 无头模式
import matplotlib.pyplot as plt
# GPU配置
def setup_gpu(pool_allocator=False):
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    
    rmm.reinitialize(
        managed_memory=False,
        pool_allocator=pool_allocator,
        devices=0,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

def generate_resolutions(start, end, step):
    """生成分辨率列表（包含结束值）"""
    if start > end:
        raise ValueError("起始值不能大于结束值")
    if step <= 0:
        raise ValueError("步长必须大于0")
    
    # 添加微小增量确保包含终点
    resolutions = np.arange(
        start=start,
        stop=end + step*0.5,  # 防止浮点精度问题
        step=step
    )
    resolutions = np.round(resolutions, 2)
    return sorted(np.unique(resolutions[resolutions <= end]))

@click.command()
@click.option('--input', '-i', required=True, help='输入h5ad文件路径')
@click.option('--output', '-o', required=True, help='输出h5ad文件路径')
@click.option('--res-start', type=float, default=0.0, show_default=True,
             help='Leiden分辨率起始值')
@click.option('--res-end', type=float, default=1.0, show_default=True,
             help='Leiden分辨率结束值（包含）')
@click.option('--res-step', type=float, default=0.2, show_default=True,
             help='Leiden分辨率步长')
@click.option('--n-pcs', '-n', default=40, show_default=True,
             help='t-SNE使用的PC数量')
@click.option('--groupby', '-g', default='leiden_harmony_0.2',
             show_default=True, help='差异表达分析的分组字段')
@click.option('--plot-dir', '-p', default='plots',
             show_default=True, help='绘图输出目录')
def pipeline(input, output, res_start, res_end, res_step, n_pcs, groupby, plot_dir):
    """
    单细胞分析主流程
    
    典型运行时间:
    - 10k细胞: 15-20分钟
    - 50k细胞: 30-45分钟
    - 100k细胞: 1-2小时
    
    示例命令:
    python sc_analysis.py -i input.h5ad -o output.h5ad \
    --res-start 0.0 --res-end 1.0 --res-step 0.2 \
    -n 40 -g leiden_harmony_0.2 -p results/plots
    """
    # 初始化环境
    start_time = time.time()
    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    setup_gpu()
    
    try:
        # 数据加载
        adata = sc.read(input)
        print(f"加载数据: {adata.n_obs} 细胞 | {adata.n_vars} 基因")
        
        # 生成分辨率列表
        resolutions = generate_resolutions(res_start, res_end, res_step)
        print(f"生成分辨率参数: {resolutions}")
        
        # Leiden聚类
        for res in resolutions:
            print(f"正在处理分辨率 {res:.1f}...")
            rsc.tl.leiden(
                adata,
                resolution=res,
                neighbors_key="harmony",
                key_added=f"leiden_harmony_{res:.1f}",
                random_state=42
            )
            
            # UMAP可视化
            sc.pl.embedding(
                adata,
                basis="X_umap_harmony",
                color=[f"leiden_harmony_{res:.1f}"],
                legend_loc="on data"
            )
            plt.savefig(f"{plot_dir}/umap_res{res:.1f}.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        # 样本分布可视化
        sc.pl.embedding(adata, basis="X_umap_harmony", color=["sample"], ncols=1)
        plt.savefig(f"{plot_dir}/umap_samples.png", dpi=300)
        plt.close()
        
        # t-SNE分析
        print(f"运行t-SNE (n_pcs={n_pcs})")
        rsc.tl.tsne(
            adata,
            n_pcs=n_pcs,
            perplexity=30,
            early_exaggeration=12,
            learning_rate=200,
            use_rep="X_pca_harmony"
        )
        
        # t-SNE可视化
        for res in resolutions:
            sc.pl.tsne(
                adata,
                color=[f'leiden_harmony_{res:.1f}'],
                legend_loc="on data"
            )
            plt.savefig(f"{plot_dir}/tsne_res{res:.1f}.png", dpi=300)
            plt.close()
        
        # 样本t-SNE可视化
        sc.pl.tsne(adata, color=["sample"], ncols=1)
        plt.savefig(f"{plot_dir}/tsne_samples.png", dpi=300)
        plt.close()
        
        # 差异表达分析
        print(f"差异表达分析 (决定分辨度的分组: {groupby})")
        rsc.tl.rank_genes_groups_logreg(adata, groupby=groupby, use_raw=False)
        sc.pl.rank_genes_groups(adata)
        plt.savefig(f"{plot_dir}/de_{groupby}.png", dpi=300)
        plt.close()
        
        # 保存结果
        adata.write_h5ad(output)
        print(f"结果已保存至 {output}")
        
    except Exception as e:
        print(f"流程执行失败: {str(e)}")
        raise
    
    # 计时
    total_time = time.time() - start_time
    print(f"流程完成，总耗时: {total_time/60:.2f} 分钟")

if __name__ == "__main__":
    pipeline()
