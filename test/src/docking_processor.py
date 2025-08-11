import click
from pathlib import Path
import pandas as pd
import csv
import sys
import re
import datetime
import json
from typing import Dict, List, Tuple, Set

@click.group()
def cli():
    """对接结果处理工具"""
    pass

def create_run_info(output_dir: Path, params: dict, processed_count: int):
    """创建运行信息文件"""
    run_info = {
        "timestamp": datetime.datetime.now().isoformat(),
        "parameters": params,
        "processed_records": processed_count,
        "output_file": str(output_dir.resolve())
    }
    
    info_file = output_dir.parent / "run_info.json"
    with info_file.open("w") as f:
        json.dump(run_info, f, indent=2)
    
    return info_file

def print_run_summary(output_path: Path, processed_count: int, success_count: int):
    """打印运行摘要"""
    click.echo("\n" + "=" * 50)
    click.echo("对接结果处理摘要:")
    click.echo(f"• 处理的成功对接记录: {success_count}")
    click.echo(f"• 总处理记录: {processed_count}")
    click.echo(f"• 输出文件: {output_path}")
    click.echo(f"• 运行记录: {output_path.parent / 'run_info.json'}")
    click.echo("=" * 50)

def extract_metadata(output_dir: str) -> dict:
    """从output_dir路径提取元数据"""
    # 标准化路径并分割
    parts = Path(output_dir).parts
    if len(parts) < 2:
        return {}
    
    antibody_full = parts[-2]
    antigen_full = parts[-1]
    
    # 拆分抗体名
    antibody_parts = antibody_full.split('_')
    antibody_part1 = antibody_parts[0] if len(antibody_parts) > 0 else ""
    antibody_part2 = antibody_parts[1] if len(antibody_parts) > 1 else ""
    
    # 拆分抗原名
    antigen_parts = antigen_full.split('_')
    # 抗原类型（除最后两部分外的所有部分）
    antigen_type = '_'.join(antigen_parts[:-2]) if len(antigen_parts) >= 3 else ""
    # 完整PDB ID和链
    pdb_chain = '_'.join(antigen_parts[-2:]) if len(antigen_parts) >= 2 else ""
    # 仅PDB ID（不带链）
    pdb_id = antigen_parts[-2] if len(antigen_parts) >= 2 else ""
    
    return {
        'antibody_full': antibody_full,
        'antigen_full': antigen_full,
        'antibody_part1': antibody_part1,
        'antibody_part2': antibody_part2,
        'antigen_type': antigen_type,
        'antigen_pdb_chain': pdb_chain,
        'antigen_pdb_id': pdb_id
    }

def parse_capri_tsv(capri_path: Path) -> dict:
    """解析capri_ss.tsv文件并返回rank=1的结果"""
    if not capri_path.exists():
        return {}
    
    try:
        with capri_path.open('r') as f:
            # 读取文件头
            header = f.readline().strip().split('\t')
            # 查找caprieval_rank列位置
            rank_idx = header.index('caprieval_rank') if 'caprieval_rank' in header else -1
            
            if rank_idx == -1:
                return {}
            
            # 查找rank=1的行
            for line in f:
                values = line.strip().split('\t')
                if len(values) <= rank_idx:
                    continue
                
                try:
                    if int(values[rank_idx]) == 1:
                        return dict(zip(header, values))
                except ValueError:
                    continue
    except Exception as e:
        click.echo(f"错误解析文件 {capri_path}: {str(e)}")
    
    return {}

@cli.command()
@click.argument("summary_file", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--output", "-o", default="docking_results_with_metadata.csv", 
              type=click.Path(path_type=Path),
              help="输出文件路径 (默认: docking_results_with_metadata.csv)")
@click.option("--dry-run", is_flag=True, help="模拟运行而不实际输出文件")
@click.option("--verbose", "-v", is_flag=True, help="显示详细处理信息")
def process_docking_results(summary_file: Path, output: Path, dry_run: bool, verbose: bool):
    """
    处理对接总结文件并生成元数据和对接结果的总文件
    
    SUMMARY_FILE: 对接总结文件路径 (如: test/registry/docking_logs/docking_summary_master.csv)
    """
    # 收集运行参数
    params = {
        "summary_file": str(summary_file.resolve()),
        "output_file": str(output.resolve()),
        "command": " ".join(sys.argv)
    }
    
    if verbose:
        click.echo(f"开始处理对接总结文件: {summary_file}")
    
    # 读取总结文件
    try:
        with summary_file.open('r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
    except Exception as e:
        raise click.ClickException(f"读取总结文件失败: {str(e)}")
    
    if verbose:
        click.echo(f"• 读取到 {len(rows)} 条记录")
    
    # 筛选成功状态
    success_rows = [row for row in rows if row.get('status', '') == 'success']
    success_count = len(success_rows)
    
    if verbose:
        click.echo(f"• 找到 {success_count} 条成功对接记录")
    
    if not success_rows:
        click.echo("警告: 没有找到成功的对接记录")
        return
    
    # 处理每条成功记录
    processed_count = 0
    results = []
    
    for row in success_rows:
        output_dir = row.get('output_dir', '')
        if not output_dir:
            continue
        
        # 提取元数据
        metadata = extract_metadata(output_dir)
        if not metadata:
            if verbose:
                click.echo(f"  警告: 无法从路径提取元数据: {output_dir}")
            continue
        
        # 构建capri_ss.tsv路径
        capri_path = Path(output_dir) / "analysis" / "7_caprieval_analysis" / "capri_ss.tsv"
        
        if verbose:
            click.echo(f"  处理: {metadata['antibody_full']} vs {metadata['antigen_full']}")
            click.echo(f"    capri路径: {capri_path}")
        
        # 解析capri结果
        capri_data = parse_capri_tsv(capri_path) if capri_path.exists() else {}
        
        if not capri_data and verbose:
            click.echo(f"    警告: 未找到rank=1的对接结果或文件不存在")
        
        # 合并数据
        merged_row = {**row, **metadata, **capri_data}
        results.append(merged_row)
        processed_count += 1
        
        if verbose:
            if capri_data:
                score = capri_data.get('score', 'N/A')
                irmsd = capri_data.get('irmsd', 'N/A')
                click.echo(f"    成功提取对接结果: score={score}, irmsd={irmsd}")
    
    if verbose:
        click.echo(f"• 成功处理 {processed_count} 条记录")
    
    # 准备输出
    if not results:
        click.echo("错误: 没有有效结果可输出")
        return
    
    # 收集所有可能的字段
    all_fields = set()
    for row in results:
        all_fields.update(row.keys())
    
    # 创建输出目录
    if not dry_run:
        output.parent.mkdir(parents=True, exist_ok=True)
    
    # 写入结果
    if dry_run:
        click.echo("\n模拟运行详情:")
        click.echo(f"• 将创建输出文件: {output}")
        click.echo(f"• 将写入 {processed_count} 条记录")
        click.echo(f"• 字段: {', '.join(sorted(all_fields))}")
    else:
        with output.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=sorted(all_fields))
            writer.writeheader()
            writer.writerows(results)
        
        # 创建运行记录
        info_file = create_run_info(output, params, processed_count)
        
        if verbose:
            click.echo(f"• 结果已保存至: {output}")
            click.echo(f"• 运行记录已创建: {info_file}")
    
    # 打印摘要
    print_run_summary(output, len(rows), processed_count)

if __name__ == "__main__":
    cli()


#脚本说明：
#
#    命令结构：
#        主命令组：process_docking_results
#        参数：对接总结文件路径（必需）
#        选项：
#            --output：指定输出文件路径（默认：docking_results_with_metadata.csv）
#            --dry-run：模拟运行模式
#            --verbose：显示详细处理信息
#
#    核心功能：
#        读取并解析对接总结CSV文件
#        筛选status为"success"的记录
#        从output_dir路径提取元数据（抗体名、抗原名等）
#        查找并解析capri_ss.tsv文件，提取rank=1的结果
#        合并所有数据并输出到CSV文件
#
#    辅助功能：
#        创建详细的运行记录（run_info.json）
#        提供模拟运行模式（dry-run）
#        详细的进度和错误报告
#        处理结束时的摘要统计
#
#    元数据提取逻辑：
#
#    python
#
#    复制代码
#    # 示例路径: ./dock/DLE_DLE1_clonotype131/Skin_celltype_3CW1_X
#    antibody_full = "DLE_DLE1_clonotype131"  # 完整抗体名
#    antigen_full = "Skin_celltype_3CW1_X"    # 完整抗原名
#
#    # 拆分结果：
#    antibody_part1 = "DLE"           # 抗体第一部分
#    antibody_part2 = "DLE1"          # 抗体第二部分
#    antigen_type = "Skin_celltype"   # 抗原类型
#    antigen_pdb_chain = "3CW1_X"     # PDB ID和链
#    antigen_pdb_id = "3CW1"          # PDB ID（不含链）

#使用示例：
#
#    基本用法：
#
#    bash
#
#复制代码
#python docking_processor.py process-docking-results test/registry/docking_logs/docking_summary_master.csv
#
#指定输出文件：
#
#bash
#
#复制代码
#python docking_processor.py process-docking-results summary.csv --output results/final_results.csv
#
#详细模式+模拟运行：
#
#bash
#
#复制代码
#python docking_processor.py process-docking-results summary.csv --verbose --dry-run
#
#查看帮助：
#
#bash
#
#    复制代码
#    python docking_processor.py --help
#    python docking_processor.py process-docking-results --help
#
#输出文件内容：
#
#输出CSV文件将包含以下列：
#
#    原始总结文件的所有列（cfg_path, status, time, output_dir, timestamp）
#    提取的元数据列（antibody_full, antigen_full, antibody_part1, antibody_part2, antigen_type, antigen_pdb_chain, antigen_pdb_id）
#    capri_ss.tsv文件中rank=1行的所有列（caprieval_rank, score, irmsd, fnat等）
