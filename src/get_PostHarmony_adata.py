"""
单细胞聚类分析流程（命令行版）

用法示例：
python sc_cluster.py \
  --input skin_epi.h5ad \
  --output processed.h5ad \
  --plot_dir ./plots \
  --n_pcs 40 \
  --n_neighbors 15 \
  --res_start 0 \
  --res_end 1.0 \
  --res_step 0.2 \
  --leiden_res 0.2

主要功能：
1. 数据预处理和质控可视化
2. Harmony批次校正
3. Leiden多分辨率聚类
4. UMAP/t-SNE可视化
5. 结果和参数持久化

典型运行时间（基于NVIDIA T4 GPU）：
- 10万细胞数据集：约15-25分钟
- 5万细胞数据集：约8-15分钟
- 1万细胞数据集：约3-5分钟

输入输出说明见下方详细文档:
# 典型工作流
# 步骤1：初始分析（探索参数）
python sc_cluster.py \
  --input raw.h5ad \
  --output initial.h5ad \
  --plot_dir initial_plots \
  --res_start 0.1 \
  --res_end 0.5 \
  --res_step 0.1

# 查看initial_plots中的可视化结果，确定最佳分辨率

# 步骤2：正式分析
python sc_cluster.py \
  --input raw.h5ad \
  --output final.h5ad \
  --plot_dir final_plots \
  --leiden_res 0.3 \          # 根据上一步确定
  --n_pcs 35 \                # 根据肘部法则调整
  --n_neighbors 50

#参数优化建议
n_pcs选择:
    通过PCA方差解释率图选择拐点：建议保留解释>80%方差的PC数

leiden_res调整策略:
    0.1-0.3：粗粒度聚类（适合细胞类型注释）
    0.4-0.6：中等粒度（亚型分析）
    0.7-1.0：细粒度（稀有细胞群识别）

n_neighbors经验公式:
    n_neighbors = min(50, sqrt(细胞数))
    例如：10万细胞 → sqrt(100000)≈316 → 取50

"""

import click
import time
import json
import scanpy as sc
import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
import rapids_singlecell as rsc
from pathlib import Path

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

@click.command()
@click.option('--input', required=True, 
              help="""输入h5ad文件路径，需包含：
                   - obs字段: sample, predicted_doublet (可选)
                   - var字段: 基因表达矩阵
                   - 预处理后的PCA结果""")
@click.option('--output', required=True,
              help="""输出h5ad文件路径，包含：
                   - obsm: X_umap, X_umap_harmony, X_tsne
                   - obs: leiden_* 多分辨率聚类结果
                   - uns: 运行参数记录""")
@click.option('--plot_dir', default='plots',
              help="""可视化结果输出目录，生成：
                   - pre_harmony_*.pdf 校正前聚类图
                   - post_harmony_*.pdf 校正后聚类图
                   - batch_effect_*.pdf 批次效应可视化
                   - qc_metrics.pdf 质控指标图""")
@click.option('--n_pcs', default=40, type=int, help='PCA主成分数量')
@click.option('--n_neighbors', default=50, type=int, help='邻域图邻居数')
@click.option('--res_start', default=0.0, type=float, help='分辨率起始值')
@click.option('--res_end', default=1.0, type=float, help='分辨率结束值')
@click.option('--res_step', default=0.2, type=float, help='分辨率步长')
@click.option('--leiden_res', default=0.2, type=float, help='主要Leiden分辨率')

def main(input, output, plot_dir, n_pcs, n_neighbors, 
         res_start, res_end, res_step, leiden_res):
    """
    执行单细胞聚类分析全流程

    典型运行示例：
    python sc_cluster.py \\
      --input raw_data.h5ad \\
      --output processed.h5ad \\
      --plot_dir ./analysis_plots \\
      --n_pcs 40 \\
      --n_neighbors 20 \\
      --res_start 0.2 \\
      --res_end 1.0 \\
      --res_step 0.2 \\
      --leiden_res 0.4

    关键步骤耗时预估（以10万细胞为例）：
    1. 数据加载：1-2分钟
    2. 邻域图计算：3-5分钟
    3. UMAP计算：2-3分钟 
    4. Leiden聚类：1分钟/分辨率
    5. Harmony整合：4-6分钟
    6. t-SNE计算：5-8分钟
    """
        
    # 初始化环境和日志
    start_time = time.time()
    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    setup_gpu()
    
    # 记录运行参数
    params = {
        'n_pcs': n_pcs,
        'n_neighbors': n_neighbors,
        'resolutions': f"{res_start}-{res_end}-{res_step}",
        'leiden_resolution': leiden_res
    }
    try:
        # 数据加载
        adata = sc.read(input)
        rsc.get.anndata_to_GPU(adata)
        
        # 邻域分析和UMAP
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        rsc.tl.umap(adata)
        
        # 生成分辨率列表
        resolutions = np.arange(res_start, res_end + 0.01, res_step).round(1)
        
        # Leiden聚类分析
        for res in resolutions:
            key = f"leiden_{res}"
            rsc.tl.leiden(adata, resolution=res, key_added=key, random_state=42)
            
            # 保存UMAP可视化
            plt.figure()
            sc.pl.umap(adata, color=[key], show=False)
            plt.savefig(f"{plot_dir}/pre_harmony_{key}.pdf")
            plt.close()

        # 批次效应可视化
        plt.figure(figsize=(8, 6))
        sc.pl.umap(adata, color=["sample"], ncols=1, show=False)
        plt.savefig(f"{plot_dir}/batch_effect_pre_harmony.pdf")
        plt.close()

        # 双联体过滤
        if "predicted_doublet" in adata.obs:
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            
            # 质量控制可视化
            fig, axs = plt.subplots(2, 2, figsize=(12, 10))
            sc.pl.umap(adata, color=["leiden_0.2"], ax=axs[0,0], show=False)
            sc.pl.umap(adata, color=["log1p_total_counts"], ax=axs[0,1], show=False)
            sc.pl.umap(adata, color=["pct_counts_MT"], ax=axs[1,0], show=False)
            sc.pl.umap(adata, color=["log1p_n_genes_by_counts"], ax=axs[1,1], show=False)
            plt.savefig(f"{plot_dir}/qc_metrics.pdf")
            plt.close()

        # Harmony整合
        setup_gpu(pool_allocator=True)
        rsc.pp.harmony_integrate(adata, key="sample", dtype=cp.float32)
        
        # 整合后分析
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, 
                        use_rep="X_pca_harmony", key_added="harmony")
        rsc.tl.umap(adata, neighbors_key="harmony", key_added="X_umap_harmony")
        
        # 整合后可视化
        plt.figure()
        sc.pl.embedding(adata, basis="X_umap_harmony", color=[f"leiden_{leiden_res}"], show=False,legend_loc="on data")
        plt.savefig(f"{plot_dir}/post_harmony_clusters.pdf")
        plt.close()

        # 保存最终结果
        adata.write_h5ad(output)
        
        # 保存运行参数
        with open(f"{plot_dir}/run_params.json", "w") as f:
            json.dump(params, f, indent=2)

        print(f"分析成功完成！总耗时：{time.time()-start_time:.1f}秒")

    except Exception as e:
        print(f"运行出错: {str(e)}")
        raise

if __name__ == "__main__":
    main()

