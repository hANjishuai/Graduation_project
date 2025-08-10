#!/bin/bash

# 抗原处理脚本
# 用于并行处理 Skin_celltype 和 Skin_tissue 目录下的所有抗原

# 设置基础目录
BASE_DIR="./antigen"
SCRIPT_PATH="./src/generate_epitope_restraints.py"  # 确保这是正确的脚本路径

# 要处理的目录列表
DIRECTORIES=("Skin_celltype" "Skin_tissue")

# 设置并行任务数（根据CPU核心数调整）
PARALLEL_JOBS=64

# 创建处理单个抗原的函数
process_antigen() {
    local antigen_dir="$1"
    local antigen_name=$(basename "$antigen_dir")
    
    # 查找CSV文件（不区分大小写）
    local csv_file=$(find "$antigen_dir" -maxdepth 1 -type f -iname "${antigen_name}_discotope3.csv" | head -1)
    
    # 查找PDB文件（不区分大小写）
    local pdb_file=$(find "$antigen_dir" -maxdepth 1 -type f -iname "${antigen_name}_clean.pdb" | head -1)
    
    # 检查文件是否存在
    if [ -z "$csv_file" ]; then
        echo "⚠ 未找到CSV文件: ${antigen_name}_discotope3.csv (在 $antigen_dir)"
        return 1
    fi
    
    if [ -z "$pdb_file" ]; then
        echo "⚠ 未找到PDB文件: ${antigen_name}_clean.pdb (在 $antigen_dir)"
        return 1
    fi
    
    # 设置输出文件路径
    local output_act_pass="${antigen_dir}/antigen-epitope.act-pass"
    local output_positions="${antigen_dir}/antigen-epitope-positions.txt"
    
    # 运行命令
    echo "处理抗原: $antigen_name"
    echo "• 使用CSV: $(basename "$csv_file")"
    echo "• 使用PDB: $(basename "$pdb_file")"
    
    python "$SCRIPT_PATH" extract-epitope \
        "$csv_file" \
        -p "$pdb_file" \
        -o "$output_act_pass" \
        -g --generate_passive > "$output_positions"
    
    # 检查执行状态
    if [ $? -eq 0 ]; then
        echo "✓ 成功处理: $antigen_name"
        echo "• 输出文件: antigen-epitope.act-pass"
        echo "• 输出文件: antigen-epitope-positions.txt"
        return 0
    else
        echo "⚠ 处理失败: $antigen_name"
        return 1
    fi
}

# 导出函数以便在parallel中使用
export -f process_antigen
export BASE_DIR SCRIPT_PATH

# 主处理循环
for dir_name in "${DIRECTORIES[@]}"; do
    dir_path="${BASE_DIR}/${dir_name}/fitted_antigens"
    
    echo "开始处理目录: $dir_name"
    echo "=========================================="
    
    # 检查目录是否存在
    if [ ! -d "$dir_path" ]; then
        echo "⚠ 目录不存在: $dir_path"
        continue
    fi
    
    # 查找所有抗原目录
    antigen_dirs=($(find "$dir_path" -mindepth 1 -maxdepth 1 -type d))
    antigen_count=${#antigen_dirs[@]}
    
    if [ $antigen_count -eq 0 ]; then
        echo "⚠ 未找到抗原目录: $dir_path"
        continue
    fi
    
    echo "找到 $antigen_count 个抗原"
    echo "使用 $PARALLEL_JOBS 个并行任务处理..."
    
    # 使用parallel并行处理所有抗原
    printf "%s\n" "${antigen_dirs[@]}" | \
        parallel -j $PARALLEL_JOBS --progress --eta --bar --joblog "${dir_path}/epitope_joblog_${dir_name}.txt" \
            "process_antigen {}"
    
    # 检查处理结果
    success_count=$(grep -c "成功处理" "${dir_path}/epitope_joblog_${dir_name}.txt" || true)
    failed_count=$(grep -c "处理失败" "${dir_path}/epitope_joblog_${dir_name}.txt" || true)
    
    echo "处理完成:"
    echo "• 成功: $success_count"
    echo "• 失败: $failed_count"
    echo "• 作业日志: ${dir_path}/epitope_joblog_${dir_name}.txt"
    
    echo "=========================================="
    echo ""
done

echo "所有抗原处理完成!"
