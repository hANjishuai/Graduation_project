import click
from pathlib import Path
import csv
import shutil
from typing import Set
import datetime
import json
import sys
import re

@click.group()
def cli():
    """抗原数据处理工具"""
    pass

def create_run_info(output_dir: Path, params: dict, antigens: Set[str]):
    """创建运行信息文件"""
    run_info = {
        "timestamp": datetime.datetime.now().isoformat(),
        "parameters": params,
        "antigens": sorted(list(antigens)),
        "antigen_count": len(antigens)
    }
    
    info_file = output_dir / "run_info.json"
    with info_file.open("w") as f:
        json.dump(run_info, f, indent=2)
    
    return info_file

def print_run_summary(output_dir: Path, processed_files: int, antigen_count: int):
    """打印运行摘要"""
    click.echo("\n" + "=" * 50)
    click.echo("运行摘要:")
    click.echo(f"• 处理的抗原数量: {antigen_count}")
    click.echo(f"• 处理的文件数量: {processed_files}")
    click.echo(f"• 输出目录: {output_dir}")
    click.echo(f"• 运行记录: {output_dir / 'run_info.json'}")
    click.echo("=" * 50)

def safe_rmtree(path: Path):
    """安全删除目录"""
    if path.exists() and path.is_dir():
        shutil.rmtree(path)

@cli.command()
@click.argument("input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("--threshold", default=0.15, type=float, 
              help="Ratio阈值 (默认: 0.15)")
@click.option("--move", is_flag=True, 
              help="移动文件而不是复制")
@click.option("--dry-run", is_flag=True, 
              help="模拟运行而不实际执行操作")
@click.option("--force", is_flag=True,
              help="强制覆盖现有输出目录")
def filter_antigens(input_dir: Path, threshold: float, move: bool, 
                    dry_run: bool, force: bool):
    """
    根据epitope_ratios.csv筛选抗原并处理文件
    
    INPUT_DIR: 包含epitope_ratios.csv和valid_prot_ids的目录
    """
    # 验证必要文件存在
    csv_file = input_dir / "epitope_ratios.csv"
    valid_dir = input_dir / "valid_prot_ids"
    output_dir = input_dir / "fitted_antigens"
    
    if not csv_file.exists():
        raise click.ClickException(f"文件不存在: {csv_file}")
    if not valid_dir.exists():
        raise click.ClickException(f"目录不存在: {valid_dir}")
    
    # 检查输出目录
    if output_dir.exists():
        if dry_run:
            click.echo(f"! 输出目录已存在: {output_dir} (模拟运行不会删除)")
        elif force:
            click.echo(f"! 删除现有输出目录: {output_dir}")
            safe_rmtree(output_dir)
        else:
            click.confirm(
                f"输出目录 {output_dir} 已存在。是否删除并继续？", 
                abort=True
            )
            safe_rmtree(output_dir)
    
    # 读取并筛选抗原
    selected_antigens = read_and_filter_antigens(csv_file, threshold)
    
    # 收集运行参数
    params = {
        "threshold": threshold,
        "move_files": move,
        "input_dir": str(input_dir.resolve()),
        "command": " ".join(sys.argv)
    }
    
    if dry_run:
        click.echo("\n模拟运行详情:")
        click.echo(f"• 阈值: {threshold}")
        click.echo(f"• 操作模式: {'移动' if move else '复制'}")
        click.echo(f"• 选择的抗原数量: {len(selected_antigens)}")
        click.echo(f"• 输出目录: {output_dir}")
        click.echo(f"• 将创建运行记录文件: {output_dir / 'run_info.json'}")
        return
    
    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 处理文件
    action = "移动" if move else "复制"
    click.echo(f"\n开始处理: {input_dir.name}")
    click.echo(f"• 阈值: {threshold}")
    click.echo(f"• 操作模式: {action}文件")
    click.echo(f"• 选择的抗原数量: {len(selected_antigens)}")
    
    processed_files = process_files(selected_antigens, valid_dir, output_dir, move)
    
    # 创建运行记录
    info_file = create_run_info(output_dir, params, selected_antigens)
    click.echo(f"• 创建运行记录: {info_file}")
    
    # 打印摘要
    print_run_summary(output_dir, processed_files, len(selected_antigens))

def read_and_filter_antigens(csv_file: Path, threshold: float) -> Set[str]:
    """读取CSV文件并返回符合条件的抗原集合"""
    selected = set()
    
    with csv_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=",")
        for i, row in enumerate(reader, 1):
            try:
                antigen = row["Antigen"]
                ratio = float(row["Ratio"])
                
                if ratio > threshold:
                    selected.add(antigen)
                    click.echo(f"  ✓ 选择: {antigen} (比率: {ratio:.4f})")
            except (KeyError, ValueError):
                continue
    
    if not selected:
        click.echo("警告: 没有找到符合条件的抗原")
    
    return selected

def find_matching_files(antigen: str, source_dir: Path) -> list:
    """查找与抗原匹配的文件（忽略大小写）"""
    # 创建忽略大小写的正则表达式
    pattern = re.compile(rf"^{re.escape(antigen)}_discotope3\.(csv|pdb)$", re.IGNORECASE)
    
    matching_files = []
    for file in source_dir.iterdir():
        if file.is_file() and pattern.match(file.name):
            matching_files.append(file)
    
    return matching_files

def process_files(antigens: Set[str], source_dir: Path, dest_dir: Path, move: bool) -> int:
    """处理匹配的文件并返回处理的数量（忽略大小写）"""
    processed = 0
    
    for antigen in antigens:
        # 创建目标文件夹
        target_folder = dest_dir / antigen
        target_folder.mkdir(exist_ok=True)
        
        # 查找匹配的文件
        matching_files = find_matching_files(antigen, source_dir)
        
        if not matching_files:
            click.echo(f"  警告: 未找到 {antigen} 的文件")
            continue
            
        # 处理找到的文件
        for source_file in matching_files:
            target_file = target_folder / source_file.name
            
            if move:
                source_file.rename(target_file)
                action = "移动"
            else:
                shutil.copy2(source_file, target_file)
                action = "复制"
            
            click.echo(f"    {action}: {source_file.name}")
            processed += 1
    
    return processed

if __name__ == "__main__":
    cli()
