import pandas as pd
import click
from igfold import IgFoldRunner
from igfold.refine.pyrosetta_ref import init_pyrosetta
from Bio import SeqIO
import os
import re
from tqdm import tqdm
import time

# ========================= CLI 命令组 =========================
@click.group()
def cli():
    """单细胞BCR分析工具: 准备数据、配对抗体、结构预测"""
    pass

# ========================= 数据准备命令 =========================
@cli.command()
@click.option('--input', '-i', required=True, help='输入CSV文件路径')
@click.option('--output', '-o', default='prepared_sequences.csv', help='输出CSV文件路径')
def prepare(input, output):
    """准备数据: 提取生产性序列并拼接V区"""
    click.echo(f"准备数据: 从 {input} 读取数据...")
    df = pd.read_csv(input, sep=",")
    
    # 统一productive列处理
    if df["productive"].dtype == bool:
        df = df[df["productive"]]
    else:
        df = df[df["productive"].str.lower() == "true"]
    
    # 拼接V区序列
    df["V_nt"] = df["fwr1_nt"] + df["cdr1_nt"] + df["fwr2_nt"] + df["cdr2_nt"] + df["fwr3_nt"] + df["cdr3_nt"] + df["fwr4_nt"]
    df["V_aa"] = df["fwr1"] + df["cdr1"] + df["fwr2"] + df["cdr2"] + df["fwr3"] + df["cdr3"] + df["fwr4"]
    
    # 确保输出目录存在
    os.makedirs(os.path.dirname(output), exist_ok=True)
    df.to_csv(output, index=False)
    click.echo(f"✓ 数据准备完成! 保存到: {output}")
    click.echo(f"生产性序列数量: {len(df)}")

# ========================= 抗体配对命令 =========================
@cli.command()
@click.option('--input', '-i', required=True, help='输入CSV文件路径')
@click.option('--output', '-o', default='paired_antibodies.fasta', help='输出FASTA文件路径')
def pair(input, output):
    """配对抗体: 生成重链和轻链配对的FASTA文件"""
    click.echo(f"配对抗体: 从 {input} 读取数据...")
    df = pd.read_csv(input)
    
    # 统一productive列处理
    if df["productive"].dtype == bool:
        df = df[df["productive"]]
    else:
        df = df[df["productive"].str.lower() == "true"]
    
    # 确保存在V_aa列
    if "V_aa" not in df.columns:
        click.echo("⚠ 警告: 未找到V_aa列，将自动拼接V区...")
        df["V_aa"] = df["fwr1"] + df["cdr1"] + df["fwr2"] + df["cdr2"] + df["fwr3"] + df["cdr3"] + df["fwr4"]
    
    # 配对处理
    paired_antibodies = []
    for (barcode, clonotype), group in df.groupby(["barcode", "raw_clonotype_id"]):
        chains = group["chain"].tolist()
        if "IGH" in chains and ("IGK" in chains or "IGL" in chains):
            heavy = group[group["chain"] == "IGH"].iloc[0]
            light = group[group["chain"].isin(["IGK", "IGL"])].iloc[0]
            paired_antibodies.append({
                "barcode": barcode,
                "clonotype": clonotype,
                "heavy_V_aa": heavy["V_aa"],
                "light_V_aa": light["V_aa"]
            })
    
    # 生成FASTA
    fasta_content = ""
    for ab in paired_antibodies:
        fasta_content += f">{ab['clonotype']}_{ab['barcode']}_heavy\n{ab['heavy_V_aa']}\n"
        fasta_content += f">{ab['clonotype']}_{ab['barcode']}_light\n{ab['light_V_aa']}\n"
    
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as f:
        f.write(fasta_content)
    
    click.echo(f"✓ 配对完成! 保存到: {output}")
    click.echo(f"配对抗体数量: {len(paired_antibodies)}")

# ========================= 结构预测命令 =========================
@cli.command()
@click.option('--fasta', '-f', required=True, help='输入FASTA文件路径')
@click.option('--output_dir', '-o', default='igfold_models', help='输出目录')
@click.option('--refine', is_flag=True, help='启用PyRosetta优化')
@click.option('--renum', is_flag=True, help='启用IMGT重编号')
@click.option('--log', '-l', default='igfold.log', help='日志文件路径')
def fold(fasta, output_dir, refine, renum, log):
    """使用IgFold预测抗体结构"""
    # 初始化PyRosetta
    init_pyrosetta()
    os.makedirs(output_dir, exist_ok=True)
    
    # 日志记录
    log_file = open(log, 'w')
    log_file.write(f"IgFold预测开始: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    log_file.write(f"参数: refine={refine}, renum={renum}\n\n")
    
    # 克隆型分组
    clonotypes = {}
    for record in SeqIO.parse(fasta, "fasta"):
        match = re.search(r'^(clonotype\d+)', record.id)
        if not match:
            continue
        clonotype_id = match.group(1)
        
        if clonotype_id not in clonotypes:
            clonotypes[clonotype_id] = {"H": None, "L": None}
        
        if "_heavy" in record.id:
            clonotypes[clonotype_id]["H"] = str(record.seq)
        elif "_light" in record.id:
            clonotypes[clonotype_id]["L"] = str(record.seq)
    
    # 进度条处理
    total = len(clonotypes)
    processed = 0
    skipped = 0
    progress_bar = tqdm(total=total, desc="预测进度", unit="clonotype")
    
    # 处理每个克隆型
    for clonotype_id, chains in clonotypes.items():
        progress_bar.set_description(f"处理 {clonotype_id}")
        
        # 检查链完整性
        if not chains["H"] or not chains["L"]:
            msg = f"跳过 {clonotype_id} - 缺少链数据"
            log_file.write(f"[跳过] {msg}\n")
            skipped += 1
            progress_bar.update(1)
            continue
        
        # 创建克隆型目录
        clonotype_dir = os.path.join(output_dir, clonotype_id)
        os.makedirs(clonotype_dir, exist_ok=True)
        output_pdb = os.path.join(clonotype_dir, f"{clonotype_id}.pdb")
        
        # 结构预测
        try:
            start_time = time.time()
            igfold = IgFoldRunner()
            igfold.fold(
                output_pdb,
                sequences={"H": chains["H"], "L": chains["L"]},
                do_refine=refine,
                do_renum=renum
            )
            elapsed = time.time() - start_time
            log_file.write(f"[成功] {clonotype_id} 预测完成 ({elapsed:.1f}s)\n")
            processed += 1
        except Exception as e:
            log_file.write(f"[错误] {clonotype_id} 失败: {str(e)}\n")
            skipped += 1
        
        progress_bar.update(1)
    
    # 清理和报告
    progress_bar.close()
    log_file.write(f"\n总计: {total}, 成功: {processed}, 跳过: {skipped}\n")
    log_file.close()
    
    click.echo(f"\n预测完成! 成功: {processed}/{total} (跳过: {skipped})")
    click.echo(f"模型保存在: {output_dir}")
    click.echo(f"日志文件: {log}")

if __name__ == '__main__':
    cli()
