"""
单细胞RNA-seq数据预处理流水线
████████████████████████████████████████████
# 查看帮助文档
python get_anndata_from_10xMtx.py --help

# 运行示例
python get_anndata_from_10xMtx.py \
  --samples_table samples.tsv \
  --output combined.h5ad \
  --delimiter $'\t'
████████████████████████████████████████████

功能描述：
1. 从表格文件批量读取10x Genomics单细胞数据
2. 合并多个样本数据为统一Anndata对象
3. 添加样本来源标签
4. 输出标准h5ad格式文件

输入输出说明：
████████████████████████████████████████████
输入要求：
- 样本表格式：2列CSV/TSV，无表头
  示例：
    HC1,./data/HC2-1118
    SLE4,./data/SLE4
    （支持#开头的注释行）

- 数据目录结构：
  每个样本路径应包含：
  - matrix.mtx.gz
  - features.tsv.gz
  - barcodes.tsv.gz

输出文件：
  H5AD格式标准单细胞数据集，包含：
  - 所有样本的合并数据
  - obs中包含样本来源标签
████████████████████████████████████████████

典型分析时间：
  样本数 | 预计耗时
  ----------------
  5个样本 | 2-5分钟
  10个样本 | 5-10分钟
  50个样本 | 20-30分钟
"""

import click
import scanpy as sc
import anndata as ad
from pathlib import Path
import csv
from datetime import datetime

sc.settings.set_figure_params(dpi=50, facecolor="white")

def get_absolute_path(file_path):
    """将相对路径转换为绝对路径"""
    return str(Path(file_path).resolve())

def read_samples_table(table_path, delimiter=','):
    """
    从表格文件读取样本信息
    
    参数：
    - table_path: 样本表路径
    - delimiter: 分隔符（默认逗号）
    
    输入格式要求：
    每行格式：<样本ID><分隔符><数据路径>
    示例：
        HC1,/data/HC1
        SLE4,/data/SLE4
    """
    samples = {}
    with open(table_path, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row_num, row in enumerate(reader, 1):
            if not row or row[0].startswith("#"):
                continue
            try:
                sample_id = row[0].strip()
                path = row[1].strip()
                samples[sample_id] = path
            except IndexError:
                raise ValueError(f"第{row_num}行格式错误，应为2列，实际{len(row)}列")
    if not samples:
        raise ValueError("样本表中未找到有效样本")
    return samples

@click.command(help="单细胞数据合并工具 | 版本：1.1 | 作者：jifanghan@lab")
@click.option('--samples_table', required=True, 
              help='样本表路径（CSV/TSV）| 示例：samples.csv')
@click.option('--output', required=True, 
              help='输出文件路径 | 示例：combined_data.h5ad')
@click.option('--delimiter', default=',', 
              help='样本表分隔符 | 默认：逗号 | 可选：制表符=> $"\t"')
def main(samples_table, output, delimiter):
    """主处理流程（记录运行时间）"""
    start_time = datetime.now()
    click.echo(f"⏰ 开始时间：{start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 读取样本表
    try:
        samples = read_samples_table(samples_table, delimiter)
        click.echo(f"📋 成功读取 {len(samples)} 个样本")
    except Exception as e:
        click.echo(f"❌ 读取样本表失败: {e}", err=True)
        raise click.Abort()

    # 处理样本
    adatas = {}
    for sample_id, filename in samples.items():
        try:
            abs_path = get_absolute_path(filename)
            click.echo(f"\n🔄 [{datetime.now().strftime('%H:%M:%S')}] 处理 {sample_id}...")
            
            sample_adata = sc.read_10x_mtx(abs_path, cache=False)
            sample_adata.var_names_make_unique()
            adatas[sample_id] = sample_adata
            
            click.echo(f"   ✅ 细胞数：{sample_adata.n_obs:,} | 基因数：{sample_adata.n_vars:,}")
        except Exception as e:
            click.echo(f"❌ 处理样本 {sample_id} 失败: {e}", err=True)
            raise click.Abort()

    # 合并保存
    try:
        click.echo(f"\n🔗 合并{len(adatas)}个样本...")
        adata = ad.concat(adatas, label="sample")
        adata.obs_names_make_unique()
        
        click.echo(f"📊 合并后维度：{adata.n_obs:,}细胞 × {adata.n_vars:,}基因")
        click.echo(adata.obs["sample"].value_counts().to_string())
        
        adata.write_h5ad(output)
        click.echo(f"\n💾 保存成功：{get_absolute_path(output)}")
    except Exception as e:
        click.echo(f"❌ 合并保存失败: {e}", err=True)
        raise click.Abort()
    
    # 计算耗时
    time_cost = datetime.now() - start_time
    click.echo(f"\n⏱️ 总耗时：{time_cost}")

if __name__ == "__main__":
    main()
