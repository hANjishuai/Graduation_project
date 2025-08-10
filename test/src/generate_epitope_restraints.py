import click
from pathlib import Path
import pandas as pd
import subprocess
import shutil
import logging
import sys
from datetime import datetime
import tomlkit
from tomlkit import TOMLDocument, table, array, comment, nl

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)
logger = logging.getLogger(__name__)

@click.group()
def epitope():
    """DiscoTope3表位预测结果处理工具"""
    pass

def run_passive_from_active(pdb_file, active_residues):
    """运行haddock3-restraints passive_from_active命令"""
    residues_str = ",".join(map(str, active_residues))
    cmd = f"haddock3-restraints passive_from_active {pdb_file} {residues_str}"
    
    try:
        logger.info(f"执行命令: {cmd}")
        result = subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True,
            timeout=30
        )
        
        output_lines = result.stdout.splitlines()
        if not output_lines:
            logger.error("命令未返回任何输出")
            return None
        
        passive_line = output_lines[0].strip()
        
        if ":" in passive_line:
            passive_line = passive_line.split(":", 1)[1].strip()
        
        return list(map(int, passive_line.split()))
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败: {e.stderr}")
        return None
    except subprocess.TimeoutExpired:
        logger.error("命令执行超时")
        return None

def write_actpass_file(active_res, passive_res, output_path):
    """写入HADDOCK active-passive文件"""
    with output_path.open('w') as f:
        f.write(" ".join(map(str, active_res)) + "\n")
        if passive_res:
            f.write(" ".join(map(str, passive_res)) + "\n")
        else:
            f.write("\n")
    return output_path

def update_actpass_file(active_res, passive_res, output_path):
    """更新HADDOCK active-passive文件，保留原始内容"""
    tmp_file = output_path.with_suffix('.tmp')
    
    with tmp_file.open('w') as f:
        f.write(" ".join(map(str, active_res)) + "\n")
        if passive_res:
            f.write(" ".join(map(str, passive_res)) + "\n")
        else:
            f.write("\n")
    
    shutil.move(tmp_file, output_path)
    return output_path

def generate_passive_residues(pdb_file, active_residues):
    """生成被动残基的封装函数"""
    logger.info("尝试生成被动残基...")
    return run_passive_from_active(pdb_file, active_residues)

@epitope.command()
@click.argument('input_csv', type=click.Path(exists=True, path_type=Path))
@click.option('--output', '-o', default='antigen-epitope.act-pass', 
              help='输出文件路径', show_default=True, type=Path)
@click.option('--pdb', '-p', required=True, 
              help='对应的PDB文件路径', type=Path)
@click.option('--min_score', '-s', type=float, default=None,
              help='最低DiscoTope分数阈值')
@click.option('--min_rsa', '-r', type=float, default=None,
              help='最低相对可及性(RSA)阈值')
@click.option('--min_plddt', '-l', type=float, default=None,
              help='最低pLDDT置信度阈值')
@click.option('--generate_passive', '-g', is_flag=True, default=False,
              help='使用passive_from_active生成passive残基')
@click.option('--verbose', '-v', is_flag=True, default=False,
              help='显示详细日志')
def extract_epitope(input_csv: Path, output: Path, pdb: Path, 
                    min_score, min_rsa, min_plddt, generate_passive, verbose):
    """
    从DiscoTope3预测结果中提取表位残基并生成HADDOCK active-passive文件
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    if not pdb.exists():
        logger.error(f"PDB文件 '{pdb}' 不存在")
        return
    
    try:
        logger.info(f"读取CSV文件: {input_csv}")
        df = pd.read_csv(input_csv)
        
        required_columns = ['res_id', 'epitope', 'DiscoTope-3.0_score', 'rsa', 'pLDDTs']
        missing_cols = [col for col in required_columns if col not in df.columns]
        
        if missing_cols:
            logger.error(f"CSV文件中缺少必要的列: {', '.join(missing_cols)}")
            return
        
        epitope_df = df[df['epitope'] == True]
        
        if epitope_df.empty:
            logger.warning("未找到任何表位残基 (epitope=True)")
            if 'DiscoTope-3.0_score' in df.columns:
                logger.info("尝试使用DiscoTope分数最高的残基...")
                epitope_df = df.sort_values('DiscoTope-3.0_score', ascending=False).head(10)
                epitope_df['epitope'] = True
            else:
                logger.error("无法确定表位残基")
                return
        
        filtered_df = epitope_df.copy()
        initial_count = len(filtered_df)
        
        if min_score is not None:
            filtered_df = filtered_df[filtered_df['DiscoTope-3.0_score'] >= min_score]
            logger.info(f"应用DiscoTope分数阈值 (≥{min_score}): 剩余 {len(filtered_df)}/{initial_count} 个残基")
        
        if min_rsa is not None:
            filtered_df = filtered_df[filtered_df['rsa'] >= min_rsa]
            logger.info(f"应用RSA阈值 (≥{min_rsa}): 剩余 {len(filtered_df)}/{initial_count} 个残基")
        
        if min_plddt is not None:
            filtered_df = filtered_df[filtered_df['pLDDTs'] >= min_plddt]
            logger.info(f"应用pLDDT阈值 (≥{min_plddt}): 剩余 {len(filtered_df)}/{initial_count} 个残基")
        
        active_residues = sorted(filtered_df['res_id'].unique().tolist())
        
        if not active_residues:
            logger.error("筛选后未找到任何表位残基")
            return
        
        logger.info("\n表位残基汇总:")
        logger.info(f"• 残基数量: {len(active_residues)}")
        logger.info(f"• 残基列表: {active_residues}")
        
        passive_residues = []
        if generate_passive:
            logger.info("\n尝试生成passive残基...")
            passive_residues = generate_passive_residues(pdb, active_residues)
            
            if passive_residues:
                logger.info("passive残基生成成功!")
                logger.info(f"• Passive残基数量: {len(passive_residues)}")
                logger.info(f"• Passive残基列表: {passive_residues}")
            else:
                logger.warning("被动残基生成失败，将使用原始活性残基作为被动残基")
                passive_residues = active_residues.copy()
        
        output_path = write_actpass_file(active_residues, passive_residues, output)
        logger.info(f"\n✓ HADDOCK active-passive文件已创建: {output_path.absolute()}")
        
        with output_path.open('r') as f:
            logger.info("\n文件内容:")
            logger.info(f.read())
        
        residues_str = ",".join(map(str, active_residues))
        passive_cmd = f"haddock3-restraints passive_from_active {pdb} {residues_str}"
        logger.info("\npassive_from_active命令参考:")
        logger.info(passive_cmd)
        
    except Exception as e:
        logger.exception(f"处理过程中发生错误: {str(e)}")

@epitope.command()
@click.argument('pdb_file', type=click.Path(exists=True, path_type=Path))
@click.argument('act_pass_file', type=click.Path(exists=True, path_type=Path))
@click.option('--backup/--no-backup', default=True,
              help='是否创建备份文件', show_default=True)
@click.option('--verbose', '-v', is_flag=True, default=False,
              help='显示详细日志')
def generate_passive(pdb_file, act_pass_file, backup, verbose):
    """
    根据已有的active-passive文件生成被动残基并更新文件
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    if not pdb_file.exists():
        logger.error(f"PDB文件 '{pdb_file}' 不存在")
        return
    
    if not act_pass_file.exists():
        logger.error(f"act-pass文件 '{act_pass_file}' 不存在")
        return
    
    try:
        if backup:
            backup_file = act_pass_file.with_suffix(act_pass_file.suffix + '.bak')
            shutil.copyfile(act_pass_file, backup_file)
            logger.info(f"已创建备份文件: {backup_file}")
        
        with open(act_pass_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            logger.error("act-pass文件为空")
            return
        
        active_line = lines[0].strip()
        if not active_line:
            logger.error("第一行未找到活性残基")
            return
        
        try:
            active_residues = list(map(int, active_line.split()))
        except ValueError as e:
            logger.error(f"解析活性残基时出错: {str(e)}")
            return
        
        logger.info(f"从文件读取的活性残基: {active_residues}")
        passive_residues = generate_passive_residues(pdb_file, active_residues)
        
        if not passive_residues:
            logger.warning("被动残基生成失败，将使用原始活性残基作为被动残基")
            passive_residues = active_residues.copy()
        
        update_actpass_file(active_residues, passive_residues, act_pass_file)
        logger.info(f"\n✓ 文件已更新: {act_pass_file.absolute()}")
        
        logger.info("文件内容:")
        with open(act_pass_file, 'r') as f:
            logger.info(f.read())
        
        logger.info("\n摘要:")
        logger.info(f"• 活性残基数量: {len(active_residues)}")
        logger.info(f"• 被动残基数量: {len(passive_residues)}")
        
    except Exception as e:
        logger.exception(f"处理过程中发生错误: {str(e)}")
        if backup and 'backup_file' in locals() and backup_file.exists():
            logger.info(f"已恢复备份文件: {backup_file}")
            shutil.move(backup_file, act_pass_file)

@epitope.command()
@click.argument('config_file', type=click.Path(exists=True, path_type=Path))
@click.option('--run_dir', default=None, 
              help='运行目录名称 (默认为run_<timestamp>)')
@click.option('--antibody', required=True, help='抗体PDB文件路径', type=Path)
@click.option('--antigen', required=True, help='抗原PDB文件路径', type=Path)
@click.option('--ambig_fname', required=True, 
              help='ambig-paratope-NMR-epitope.tbl文件路径', type=Path)
@click.option('--unambig_fname', required=True, 
              help='antibody-unambig.tbl文件路径', type=Path)
@click.option('--reference_str',help='参考复合物结构cluster_1_model_1.pdb.gz文件路径', type=Path)
@click.option('--sampling', type=int, default=100, 
              help='rigidbody采样数', show_default=True)
@click.option('--select', type=int, default=40, 
              help='seletop选择模型数', show_default=True)
@click.option('--top_models', type=int, default=4, 
              help='每个簇选择模型数', show_default=True)
@click.option('--output', '-o', default=None, 
              help='输出配置文件路径 (默认为输入文件加"_updated"后缀)', type=Path)
@click.option('--backup/--no-backup', default=True,
              help='是否创建备份文件', show_default=True)
@click.option('--verbose', '-v', is_flag=True, default=False,
              help='显示详细日志')
def configure_docking(config_file, run_dir, antibody, antigen, ambig_fname, 
                     unambig_fname, sampling, select, top_models, output, reference_str,
                     backup, verbose):
    """
    配置HADDOCK分子对接参数文件
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    # 验证输入文件存在
    required_files = [antibody, antigen, ambig_fname, unambig_fname]
    for file in required_files:
        if not file.exists():
            logger.error(f"文件不存在: {file}")
            return
    
    # 设置默认值
    if run_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = f"run_{timestamp}"
    
    if output is None:
        output = config_file.with_name(f"{config_file.stem}_updated.cfg")
    
    # 创建备份
    if backup:
        backup_file = config_file.with_name(f"{config_file.stem}.bak")
        shutil.copyfile(config_file, backup_file)
        logger.info(f"已创建备份文件: {backup_file}")
    
    try:
        # 读取原始配置文件
        logger.info(f"读取配置文件: {config_file}")
        with open(config_file, 'r') as f:
            config_content = f.read()
        
        # 将配置文件内容转换为字典
        config_dict = {}
        current_section = None
        
        for line in config_content.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
                
            # 处理节标题
            if line.startswith("[") and line.endswith("]"):
                current_section = line[1:-1]
                config_dict[current_section] = {}
                continue
                
            # 处理键值对
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                
                # 处理数组值
                if value.startswith("[") and value.endswith("]"):
                    # 提取数组元素
                    items = value[1:-1].split(",")
                    value = [item.strip().strip('"').strip("'") for item in items]
                else:
                    # 处理引号包裹的字符串值
                    if (value.startswith('"') and value.endswith('"')) or \
                       (value.startswith("'") and value.endswith("'")):
                        value = value[1:-1]
                
                if current_section:
                    config_dict[current_section][key] = value
                else:
                    config_dict[key] = value
        
        # 更新顶层参数
        logger.info("更新顶层参数...")
        config_dict["run_dir"] = run_dir
        config_dict["molecules"] = [str(antibody), str(antigen)]
        
        # 确保顶层参数有值
        config_dict.setdefault("mode", "local")
        config_dict.setdefault("ncores", 50)
        config_dict.setdefault("postprocess", True)
        config_dict.setdefault("clean", True)
        
        # 更新rigidbody部分
        logger.info("更新rigidbody部分...")
        config_dict.setdefault("rigidbody", {})
        config_dict["rigidbody"]["ambig_fname"] = str(ambig_fname)
        config_dict["rigidbody"]["unambig_fname"] = str(unambig_fname)
        config_dict["rigidbody"]["sampling"] = sampling
        
        # 更新seletop部分
        logger.info("更新seletop部分...")
        config_dict.setdefault("seletop", {})
        config_dict["seletop"]["select"] = select
        
        # 更新flexref部分
        logger.info("更新flexref部分...")
        config_dict.setdefault("flexref", {})
        config_dict["flexref"]["ambig_fname"] = str(ambig_fname)
        config_dict["flexref"]["unambig_fname"] = str(unambig_fname)
        
        # 更新emref部分
        logger.info("更新emref部分...")
        config_dict.setdefault("emref", {})
        config_dict["emref"]["ambig_fname"] = str(ambig_fname)
        config_dict["emref"]["unambig_fname"] = str(unambig_fname)
        
        # 更新seletopclusts部分
        logger.info("更新seletopclusts部分...")
        config_dict.setdefault("seletopclusts", {})
        config_dict["seletopclusts"]["top_models"] = top_models
        
        # 更新caprieval部分
        logger.info("更新caprieval部分...")
        config_dict.setdefault("caprieval", {})
        if reference_str is not None:
            config_dict["caprieval"]["reference_fname"] = str(reference_str)

        # 确保所有HADDOCK模块都存在
        for section in ["topoaa", "rigidbody", "seletop", "flexref", 
                        "emref", "clustfcc", "seletopclusts", "caprieval", "contactmap"]:
            config_dict.setdefault(section, {})
        
        # 生成新的配置文件内容
        new_content = f"""# ====================================================================
# HADDOCK 分子对接配置文件
# 生成时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
# ====================================================================

# directory in which the scoring will be done
run_dir = "{run_dir}"

# compute mode
mode = "{config_dict.get('mode', 'local')}"
ncores = {config_dict.get('ncores', 50)}

# Post-processing to generate statistics and plots
postprocess = {str(config_dict.get('postprocess', True)).lower()}
clean = {str(config_dict.get('clean', True)).lower()}

molecules = [
    "{str(antibody)}",
    "{str(antigen)}"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
"""

        # 添加所有模块 - 确保所有模块都出现，即使没有参数
        sections = [
            "topoaa", "rigidbody", "seletop", "flexref", 
            "emref", "clustfcc", "seletopclusts", "caprieval", "contactmap"
        ]
        
        for section in sections:
            new_content += f"\n[{section}]\n"
            if section in config_dict and config_dict[section]:
                for key, value in config_dict[section].items():
                    if isinstance(value, list):
                        # 格式化数组
                        array_str = ", ".join([f'"{v}"' for v in value])
                        new_content += f"{key} = [{array_str}]\n"
                    else:
                        # 确保字符串值有引号
                        if isinstance(value, str):
                            # 路径值加引号
                            if "/" in value or "\\" in value:
                                new_content += f'{key} = "{value}"\n'
                            # 布尔值不加引号
                            elif value.lower() in ["true", "false"]:
                                new_content += f"{key} = {value.lower()}\n"
                            # 数字值不加引号
                            elif value.isdigit():
                                new_content += f"{key} = {value}\n"
                            else:
                                new_content += f'{key} = "{value}"\n'
                        else:
                            new_content += f"{key} = {value}\n"
            # 即使没有参数，也保留模块标题
            new_content += "\n"
        
        # 写入更新后的配置文件
        logger.info(f"写入更新后的配置文件: {output}")
        with open(output, 'w') as f:
            f.write(new_content)
        
        logger.info(f"\n✓ 配置文件已更新: {output.absolute()}")
        logger.info("\n摘要:")
        logger.info(f"• 运行目录: {run_dir}")
        logger.info(f"• 抗体文件: {antibody}")
        logger.info(f"• 抗原文件: {antigen}")
        logger.info(f"• Ambig约束文件: {ambig_fname}")
        logger.info(f"• Unambig约束文件: {unambig_fname}")
        logger.info(f"• 采样数: {sampling}")
        logger.info(f"• 选择模型数: {select}")
        logger.info(f"• 每簇模型数: {top_models}")
        logger.info(f"• 参考结构: {reference_str}")

        
        # 显示文件内容预览
        logger.info("\n配置文件预览:")
        with open(output, 'r') as f:
            for i, line in enumerate(f):
                if i < 25:  # 显示前25行
                    logger.info(line.strip())
                else:
                    logger.info("...")
                    break
        
    except Exception as e:
        logger.exception(f"配置更新失败: {str(e)}")
        if backup and 'backup_file' in locals() and backup_file.exists():
            logger.info(f"已恢复备份文件: {backup_file}")
            shutil.move(backup_file, config_file)

if __name__ == '__main__':
    epitope()