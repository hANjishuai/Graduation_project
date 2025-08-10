#!/bin/bash

# 设置基础目录
ANTIBODY_DIR="./antibody"

# 创建处理单个克隆型的函数
process_clonotype() {
    local clonotype_dir="$1"
    
    # 提取路径组件
    local sample_dir=$(dirname "$(dirname "$clonotype_dir")")
    local sample_name=$(basename "$sample_dir")
    local group_dir=$(dirname "$sample_dir")
    local group_name=$(basename "$group_dir")
    local clonotype_name=$(basename "$clonotype_dir")
    
    # 设置输入输出路径
    local fasta_file="$clonotype_dir/$clonotype_name.fasta"
    local csv_file="$sample_dir/prepared_sequences.csv"
    local output_dir="$clonotype_dir/paratope_meta"
    local output_file="$output_dir/cdr_positions.txt"
    local act_pass_file="$output_dir/antibody-paratope.act-pass"
    local prefix="$output_dir/clonetype_renumbered"
    
    # 确保输入文件存在
    if [ ! -f "$fasta_file" ]; then
        echo "⚠ 错误: FASTA文件不存在 - $fasta_file"
        return 1
    fi
    
    if [ ! -f "$csv_file" ]; then
        echo "⚠ 错误: CSV文件不存在 - $csv_file"
        return 1
    fi
    
    # 创建输出目录
    mkdir -p "$output_dir"
    
    # 执行CDR定位命令
    echo "处理: $ANTIBODY_DIR/$group_name/$sample_name/all_igfold_models/$clonotype_name"
    python ./src/fasta_tools.py locate-cdr \
        "$fasta_file" \
        -t "$clonotype_name" \
        -c "$csv_file" \
        -o "$output_file" \
        -a "$act_pass_file" \
        -p "$prefix"
    
    # 检查执行状态
    if [ $? -eq 0 ]; then
        echo "✓ 成功: $clonotype_name"
    else
        echo "⚠ 失败: $clonotype_name"
    fi
}

# 导出函数以便在parallel中使用
export -f process_clonotype

# 查找所有克隆型目录
find "$ANTIBODY_DIR" -type d -path "*/all_igfold_models/clonotype*" -prune | \
parallel -j 64 --progress --bar --joblog cdr_joblog.txt \
    process_clonotype {}

echo "所有CDR定位处理完成!"
