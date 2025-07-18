"""
🌟 [2023-10-02 09:30:00] 单细胞数据转换流程启动
🕒 开始时间: 2023-10-02 09:30:00

🌟 [2023-10-02 09:30:01] 步骤1/4: 加载数据
📥 输入文件: pbmc3k_processed.h5ad

🌟 [2023-10-02 09:30:02] 输入数据解剖报告

📦 对象概览
  细胞数: 2638 | 基因数: 1838
  X矩阵类型: csr_matrix (稀疏)

🔹 核心数据结构
  obs元数据 (2638, 3):
    ⋯ 列: n_genes, percent_mito, louvain
  var元数据 (1838, 2):
    ⋯ 列: gene_ids, highly_variable

🔹 降维坐标
  X_pca: (2638, 50) → 50D
  X_umap: (2638, 2) → 2D

🔹 原始数据备份
  基因数: 1838 (当前使用: 1838)
  var示例: MIR1302-10, FAM138A, OR4F5...

🔹 附加元数据
  存储的键: pca, umap, louvain

🌟 [2023-10-02 09:30:03] 步骤2/4: 数据源选择
✅ 原始基因数: 1838 | 当前基因数: 1838
是否使用原始完整基因数据？ [y/N]: y

🌟 [2023-10-02 09:30:04] 步骤3/4: 准备输出
✅ 创建输出目录: ./output

🌟 [2023-10-02 09:30:05] 步骤4/4: 导出数据
✅ 保存表达矩阵 (1838 基因)
  路径: ./output/matrix.mtx
  示例基因: MIR1302-10, FAM138A, OR4F5...
✅ 保存细胞元数据 (3 列)
✅ 已保存 X_pca (50D)
✅ 已保存 X_umap (2D)

🌟 [2023-10-02 09:30:06] 转换完成
🎉 成功生成以下文件:
./output/
├── matrix.mtx
├── metadata.csv
├── pca.csv
└── umap.csv

"""

import click
import scanpy as sc
import pandas as pd
import scipy.sparse
from scipy.io import mmwrite
import os
from pathlib import Path
from datetime import datetime

# --------------------------
# 日志系统
# --------------------------
def log_header(message):
    """步骤开始标题"""
    click.echo(click.style(f"\n[{datetime.now()}] 🌟 {message}", fg="magenta", bold=True))

def log_subheader(message):
    """子步骤标题"""
    click.echo(click.style(f"[{datetime.now()}] 🔹 {message}", fg="blue"))

def log_success(message):
    """成功提示"""
    click.echo(click.style(f"[{datetime.now()}] ✅ {message}", fg="green"))

def log_warning(message):
    """警告提示"""
    click.echo(click.style(f"[{datetime.now()}] ⚠️ {message}", fg="yellow"))

def log_error(message):
    """错误提示"""
    click.echo(click.style(f"[{datetime.now()}] ❌ {message}", fg="red"))

# --------------------------
# 数据验证与展示
# --------------------------
def inspect_anndata(adata):
    """详细展示AnnData对象结构"""
    log_header("输入数据解剖报告")
    
    # 基础信息
    click.echo(f"📦 对象概览")
    click.echo(f"  细胞数: {adata.n_obs} | 基因数: {adata.n_vars}")
    click.echo(f"  X矩阵类型: {type(adata.X).__name__} ({'稀疏' if scipy.sparse.issparse(adata.X) else '稠密'})")
    
    # 核心数据
    log_subheader("核心数据结构")
    click.echo(f"  obs元数据 ({adata.obs.shape}):")
    click.echo(f"    ⋯ 列: {', '.join(adata.obs.columns[:])}" + ("..." if len(adata.obs.columns)>3 else ""))
    click.echo(f"  var元数据 ({adata.var.shape}):")
    click.echo(f"    ⋯ 列: {', '.join(adata.var.columns[:])}" + ("..." if len(adata.var.columns)>3 else ""))
    
    # 降维信息
    log_subheader("降维坐标")
    if len(adata.obsm) == 0:
        click.echo("  ∅ 未检测到任何降维坐标")
    else:
        for key in adata.obsm:
            dims = adata.obsm[key].shape[1]
            click.echo(f"  {key}: {adata.obsm[key].shape} → {dims}D")
    
    # 原始数据状态
    log_subheader("原始数据备份")
    if adata.raw is None:
        click.echo("  ∅ 无原始数据 (adata.raw is None)")
    else:
        click.echo(f"  基因数: {adata.raw.n_vars} (当前使用: {adata.n_vars})")
        click.echo(f"  var示例: {', '.join(adata.raw.var.index[:3])}...")
    
    # 非结构化数据
    log_subheader("附加元数据")
    if len(adata.uns) == 0:
        click.echo("  ∅ 无附加数据")
    else:
        keys = list(adata.uns.keys())[:3]
        click.echo(f"  存储的键: {', '.join(keys)}" + ("..." if len(adata.uns)>3 else ""))

# --------------------------
# 核心功能函数
# --------------------------
# 修改数据源选择逻辑
def check_raw_data(adata, use_raw):
    if adata.raw is None:
        return adata.X, adata.var
    else:
        if use_raw:  # 自动确认逻辑
            click.echo("✅ 自动选择原始基因数据")
            return adata.raw.X, adata.raw.var
        else:
            if click.confirm("是否使用原始完整基因数据？", default=True):
                return adata.raw.X, adata.raw.var
            else:
                return adata.X, adata.var

def save_sparse_matrix(matrix, genes, output_dir):
    """保存表达矩阵"""
    # 保存基因名
    gene_path = Path(output_dir) / "genes.tsv"
    pd.Series(genes).to_csv(gene_path, index=False, header=False)

    # 保存表达矩阵
    output_path = Path(output_dir) / "matrix.mtx"
    try:
        mmwrite(output_path, matrix)
        log_success(f"保存表达矩阵 ({matrix.shape[1]} 基因)")
        click.echo(f"  路径: {output_path}")
        click.echo(f"  示例基因: {', '.join(genes[:3])}...")
    except Exception as e:
        log_error(f"保存失败: {str(e)}")
        raise

def save_reductions(adata, output_dir, reductions):
    """处理降维结果"""
    log_subheader("正在导出降维坐标")
    for key in reductions.split(","):
        if key in adata.obsm:
            df = pd.DataFrame(adata.obsm[key], index=adata.obs.index)
            output_path = Path(output_dir) / f"{key.replace('X_', '')}.csv"
            df.to_csv(output_path)
            log_success(f"已保存 {key} ({df.shape[1]}D)")
        else:
            log_warning(f"跳过未找到的降维坐标: {key}")

# --------------------------
# 主函数
# --------------------------
@click.command()
@click.option("--input_h5ad", required=True, help="输入AnnData文件路径 (.h5ad)")
@click.option("--output_dir", required=True, help="输出目录路径")
@click.option("--reductions", default="X_pca,X_umap", help="逗号分隔的降维坐标名称")
@click.option("--use_raw", is_flag=True, help="自动使用原始基因数据")  # 新增参数
def main(input_h5ad, output_dir, reductions, use_raw):
    # 初始化
    log_header("单细胞数据转换流程启动")
    click.echo(f"🕒 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 加载数据
    log_header("步骤1/4: 加载数据")
    click.echo(f"📥 输入文件: {click.format_filename(input_h5ad)}")
    adata = sc.read_h5ad(input_h5ad)
    inspect_anndata(adata)  # 关键新增功能
    
    # 数据源确认
    log_header("步骤2/4: 数据源选择")
    matrix, gene_meta = check_raw_data(adata, use_raw)
    
    # 创建输出目录
    log_header("步骤3/4: 准备输出")
    os.makedirs(output_dir, exist_ok=True)
    log_success(f"创建输出目录: {output_dir}")
    
    # 导出数据
    log_header("步骤4/4: 导出数据")
    save_sparse_matrix(matrix, gene_meta.index.tolist(), output_dir)
    
    # 保存元数据
    pd.DataFrame(adata.obs).to_csv(Path(output_dir)/"metadata.csv")
    log_success(f"保存细胞元数据 ({adata.obs.shape[1]} 列)")
    
    # 处理降维
    save_reductions(adata, output_dir, reductions)
    
    # 完成
    log_header("转换完成")
    click.echo(click.style("🎉 成功生成以下文件:", fg="bright_green"))
    os.system(f"tree {output_dir}")

if __name__ == "__main__":
    main()
