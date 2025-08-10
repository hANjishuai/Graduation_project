import click
from pathlib import Path
import pandas as pd
import re
import os
from collections import OrderedDict

def clean_sequence(seq):
    """清理序列中的非常规字符"""
    return re.sub(r'[^A-Z]', '', seq.upper()) if seq else ""

def find_cdr_position(sequence, cdr_seq):
    """在序列中查找CDR位置，处理可能的插入符号"""
    cdr_clean = clean_sequence(cdr_seq)
    if not cdr_clean:
        return -1, -1
    
    # 尝试直接匹配
    start = sequence.find(cdr_clean)
    if start != -1:
        return start, start + len(cdr_clean) - 1
    
    # 处理可能的插入符号（如haddock中的格式）
    cdr_pattern = re.compile(re.escape(cdr_clean).replace('X', '.'))
    match = cdr_pattern.search(sequence)
    if match:
        return match.start(), match.end() - 1
    
    return -1, -1

def generate_paratope_formats(cdr_positions):
    """生成HADDOCK和PyMOL格式的paratope定义"""
    # 收集所有位置
    all_positions = []
    for (start, end) in cdr_positions.values():
        all_positions.extend(range(start, end + 1))
    
    # 排序并去重
    all_positions = sorted(set(all_positions))
    
    # 生成HADDOCK格式 (空格分隔)
    haddock_format = " ".join(map(str, all_positions))
    
    # 生成PyMOL格式 (加号分隔)
    pymol_format = "+".join(map(str, all_positions))
    
    return haddock_format, pymol_format

def write_haddock_actpass_file(positions, output_path):
    """写入HADDOCK active-passive文件"""
    with output_path.open('w') as f:
        # 第一行: active residues
        f.write(positions + "\n")
        # 第二行: passive residues (留空)
        f.write("\n")
    
    return output_path

@click.group()
def cli():
    """抗体序列处理工具集 - FASTA/CDR分析"""
    pass

@cli.command()
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option('--heavy', '-H', default='heavy.fasta', 
              help='重链输出文件路径', show_default=True, type=Path)
@click.option('--light', '-L', default='light.fasta', 
              help='轻链输出文件路径', show_default=True, type=Path)
def split_chains(input_file: Path, heavy: Path, light: Path):
    """
    将FASTA文件中的重链(>H)和轻链(>L)分离到不同文件
    
    \b
    示例:
      split_chains input.fasta
      split_chains input.fasta -H heavy.fa -L light.fa
      split_chains input.fasta --heavy output/H.fasta --light output/L.fasta
    """
    # 确保输出目录存在
    heavy.parent.mkdir(parents=True, exist_ok=True)
    light.parent.mkdir(parents=True, exist_ok=True)
    
    # 初始化变量
    current_chain = None
    heavy_data = []
    light_data = []
    
    # 读取并处理FASTA文件
    with input_file.open('r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:  # 跳过空行
                continue
                
            if stripped.startswith('>H'):
                current_chain = 'heavy'
                heavy_data.append(line)  # 保留原始行（包含换行符）
            elif stripped.startswith('>L'):
                current_chain = 'light'
                light_data.append(line)
            elif current_chain == 'heavy':
                heavy_data.append(line)
            elif current_chain == 'light':
                light_data.append(line)
    
    # 写入输出文件
    heavy.write_text(''.join(heavy_data), encoding='utf-8')
    light.write_text(''.join(light_data), encoding='utf-8')
    
    # 结果反馈
    click.secho("✓ FASTA文件分离成功!", fg='green', bold=True)
    click.echo(f"• 重链序列数: {len([l for l in heavy_data if l.startswith('>')])}")
    click.echo(f"• 轻链序列数: {len([l for l in light_data if l.startswith('>')])}")
    click.echo(f"• 重链输出: {heavy.absolute()}")
    click.echo(f"• 轻链输出: {light.absolute()}")

@cli.command()
@click.argument('fasta_file', type=click.Path(exists=True, path_type=Path))
@click.option('--csv', '-c', default='paired_antibody_sequences.csv', 
              help='抗体序列CSV文件路径', show_default=True, type=Path)
@click.option('--clonotype', '-t', default=None,
              help='指定clonotype ID (默认为FASTA文件名)')
@click.option('--output', '-o', default=None,
              help='输出文件路径', type=Path)
@click.option('--pdb', '-p', default=None,
              help='PDB文件路径(用于PyMOL命令生成)', type=Path)
@click.option('--actpass', '-a', default='antibody-paratope.act-pass',
              help='HADDOCK active-passive文件路径', show_default=True, type=Path)
def locate_cdr(fasta_file: Path, csv: Path, clonotype: str, output: Path, pdb: Path, actpass: Path):
    """
    定位FASTA文件中CDR区域的位置并生成paratope定义
    
    \b
    功能:
      1. 读取FASTA文件并合并重链和轻链
      2. 从CSV文件中提取指定clonotype的CDR序列
      3. 在合并序列中定位CDR位置
      4. 生成HADDOCK active-passive文件
      5. 生成PyMOL可视化命令
    
    \b
    示例:
      locate_cdr clonotype54.fasta
      locate_cdr clonotype54.fasta -c my_sequences.csv -o positions.txt
      locate_cdr ab.fasta --clonotype clonotype123 --pdb antibody.pdb
    """
    # 确定clonotype ID
    if not clonotype:
        clonotype = fasta_file.stem
        click.echo(f"• 使用FASTA文件名作为clonotype ID: {clonotype}")
    
    # 1. 读取并合并FASTA序列
    chains = {}
    with fasta_file.open('r') as fasta:
        current_chain = None
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                current_chain = line[1:]
                chains[current_chain] = ""
            elif current_chain:
                chains[current_chain] += line
    
    # 识别链类型（H为重链，L为轻链）
    heavy_key = next((k for k in chains if k.startswith("H") or "heavy" in k.lower()), None)
    light_key = next((k for k in chains if k.startswith("L") or "light" in k.lower()), None)
    
    if not heavy_key or not light_key:
        click.secho("错误: FASTA文件中必须包含重链(H)和轻链(L)", fg='red', bold=True)
        click.echo(f"检测到的链: {list(chains.keys())}")
        return
    
    heavy_seq = clean_sequence(chains[heavy_key])
    light_seq = clean_sequence(chains[light_key])
    combined_seq = heavy_seq + light_seq
    heavy_length = len(heavy_seq)
    
    click.echo(f"\n• 重链序列 ({heavy_key}): {heavy_seq[:30]}... (长度: {heavy_length})")
    click.echo(f"• 轻链序列 ({light_key}): {light_seq[:30]}... (长度: {len(light_seq)})")
    click.echo(f"• 合并序列: {combined_seq[:30]}... (总长度: {len(combined_seq)})\n")
    
    # 2. 读取CSV文件并处理数据
    try:
        df = pd.read_csv(csv)
    except FileNotFoundError:
        click.secho(f"错误: 找不到CSV文件 '{csv}'", fg='red', bold=True)
        return
    
    if "raw_clonotype_id" not in df.columns:
        click.secho("错误: CSV文件中缺少 'raw_clonotype_id' 列", fg='red', bold=True)
        click.echo(f"可用列: {list(df.columns)}")
        return
    
    clonotype_data = df[df["raw_clonotype_id"] == clonotype]
    
    if clonotype_data.empty:
        click.secho(f"警告: 在CSV文件中未找到 '{clonotype}' 的条目", fg='yellow', bold=True)
        available_clonotypes = df["raw_clonotype_id"].unique()[:5]  # 只显示前5个
        click.echo(f"可用的clonotype_id: {', '.join(available_clonotypes)}" + 
                   ("..." if len(df["raw_clonotype_id"].unique()) > 5 else ""))
        return
    
    # 3. 提取CDR序列并定位位置
    cdr_positions = OrderedDict()
    chain_mapping = {"IGH": "H", "IGK": "L", "IGL": "L"}
    
    for _, row in clonotype_data.iterrows():
        chain_type = row["chain"]
        simple_chain = chain_mapping.get(chain_type)
        
        if not simple_chain:
            click.secho(f"警告: 未知链类型 '{chain_type}'，跳过", fg='yellow')
            continue
        
        # 获取CDR序列
        cdr1 = str(row.get("cdr1", ""))
        cdr2 = str(row.get("cdr2", ""))
        cdr3 = str(row.get("cdr3", ""))
        
        # 确定搜索序列和偏移量
        search_seq = heavy_seq if simple_chain == "H" else light_seq
        offset = 0 if simple_chain == "H" else heavy_length
        
        # 查找CDR位置
        for i, cdr in enumerate([cdr1, cdr2, cdr3], 1):
            if not cdr or cdr.lower() == "nan":
                click.secho(f"警告: {simple_chain}_CDR{i} 序列缺失", fg='yellow')
                continue
                
            start, end = find_cdr_position(search_seq, cdr)
            if start == -1:
                click.secho(f"错误: 在{simple_chain}链中未找到CDR{i}序列 '{cdr}'", fg='red')
                continue
                
            # 转换为1-based索引并添加偏移量
            start_1based = start + 1 + offset
            end_1based = end + 1 + offset
            cdr_name = f"{simple_chain}_CDR{i}"
            
            cdr_positions[cdr_name] = (start_1based, end_1based)
            
            # 验证找到的序列
            found_seq = combined_seq[start+offset:end+offset+1]
            click.secho(f"✓ {cdr_name} 定位成功: 位置 {start_1based}-{end_1based}", fg='green')
            click.echo(f"  预期序列: {cdr}")
            click.echo(f"  实际序列: {found_seq}\n")
    
    # 4. 输出结果
    if not cdr_positions:
        click.secho("错误: 未找到任何CDR位置", fg='red', bold=True)
        return
    
    click.secho("\nCDR位置汇总 (1-based索引):", bold=True)
    for cdr, (start, end) in cdr_positions.items():
        click.echo(f"{cdr}: {start}-{end}")
    
    # 5. 生成paratope格式
    haddock_paratope, pymol_paratope = generate_paratope_formats(cdr_positions)
    
    click.secho("\nHADDOCK paratope格式 (直接用于配置文件):", bold=True)
    click.echo(haddock_paratope)
    
    click.secho("\nPyMOL paratope格式 (用于可视化):", bold=True)
    click.echo(pymol_paratope)
    
    # 6. 创建HADDOCK active-passive文件
    actpass_path = write_haddock_actpass_file(haddock_paratope, actpass)
    click.secho(f"\n✓ HADDOCK active-passive文件已创建: {actpass_path.absolute()}", fg='green', bold=True)
    
    # 7. 生成PyMOL命令
    pdb_name = pdb.name if pdb else "antibody.pdb"
    pymol_cmds = f"""
# PyMOL可视化命令
load {pdb_name}
select paratope, resi {pymol_paratope}
show surface, {pdb_name}
show cartoon, {pdb_name}
color gray, {pdb_name}
color red, paratope
set surface_color, red, paratope
set transparency, 0.5
zoom paratope
"""
    
    click.secho("\nPyMOL可视化命令:", bold=True)
    click.echo(pymol_cmds)
    
    # 8. 保存位置信息到文件
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('w') as pos_file:
            pos_file.write(f"# CDR位置报告 - {clonotype}\n")
            pos_file.write(f"FASTA文件: {fasta_file}\n")
            pos_file.write(f"CSV文件: {csv}\n")
            pos_file.write(f"HADDOCK active-passive文件: {actpass_path.absolute()}\n\n")
            pos_file.write(f"重链序列 ({heavy_key}): {heavy_seq}\n")
            pos_file.write(f"轻链序列 ({light_key}): {light_seq}\n")
            pos_file.write(f"合并序列: {combined_seq}\n\n")
            pos_file.write(f"总长度: {len(combined_seq)}\n")
            pos_file.write(f"重链长度: {heavy_length}\n")
            pos_file.write(f"轻链长度: {len(light_seq)}\n\n")
            pos_file.write("CDR位置 (1-based):\n")
            for cdr, (start, end) in cdr_positions.items():
                pos_file.write(f"{cdr}: {start}-{end}\n")
            
            pos_file.write("\nHADDOCK active-passive文件内容:\n")
            pos_file.write(f"{haddock_paratope}\n\n")
            
            pos_file.write("PyMOL paratope格式:\n")
            pos_file.write(f"{pymol_paratope}\n\n")
            
            pos_file.write("PyMOL可视化命令:\n")
            pos_file.write(pymol_cmds)
        
        click.secho(f"\n✓ 所有信息已保存到 '{output}'", fg='green', bold=True)

if __name__ == '__main__':
    cli()
