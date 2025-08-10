#!/bin/bash

# 抗原PDB标准化处理脚本
# 用于处理 Skin_celltype 和 Skin_tissue 目录下的所有抗原PDB文件

# 设置基础目录
BASE_DIR="./antigen"

# 要处理的目录列表
DIRECTORIES=("Skin_celltype" "Skin_tissue")

# 创建处理单个抗原PDB的函数
process_antigen() {
    local antigen_dir="$1"
    local antigen_name=$(basename "$antigen_dir")
    
    # 查找匹配的PDB文件（不区分大小写）
    local input_pdb=$(find "$antigen_dir" -maxdepth 1 -type f -iname "${antigen_name}_discotope3.pdb" | head -1)
    
    # 检查输入文件是否存在
    if [ -z "$input_pdb" ]; then
        echo "⚠ 未找到PDB文件: ${antigen_name}_discotope3.pdb (尝试查找不同大小写)"
        return 1
    fi
    
    # 设置输出文件路径
    local output_pdb="${antigen_dir}/${antigen_name}_clean.pdb"
    
    # 执行PDB处理管道
    cat "$input_pdb" | \
      pdb_tidy -strict | \
      pdb_delhetatm | \
      pdb_selaltloc | \
      pdb_keepcoord | \
      pdb_chain -B | \
      pdb_chainxseg | \
      pdb_tidy -strict > "$output_pdb"
    
    # 检查是否成功
    if [ $? -eq 0 ]; then
        # 检查输出文件是否非空
        if [ -s "$output_pdb" ]; then
            echo "✓ 成功处理: $antigen_name (使用: $(basename "$input_pdb"))"
            return 0
        else
            echo "⚠ 输出文件为空: $output_pdb"
            return 1
        fi
    else
        echo "⚠ 处理失败: $antigen_name"
        return 1
    fi
}

# 导出函数以便在parallel中使用
export -f process_antigen

# 主处理函数
process_directory() {
    local dir_name="$1"
    local dir_path="${BASE_DIR}/${dir_name}/fitted_antigens"
    
    echo "开始处理目录: $dir_name"
    echo "=========================================="
    
    # 检查目录是否存在
    if [ ! -d "$dir_path" ]; then
        echo "⚠ 目录不存在: $dir_path"
        return 1
    fi
    
    # 查找所有抗原目录
    local antigen_dirs=($(find "$dir_path" -mindepth 1 -maxdepth 1 -type d))
    local antigen_count=${#antigen_dirs[@]}
    
    if [ $antigen_count -eq 0 ]; then
        echo "⚠ 未找到抗原目录: $dir_path"
        return 1
    fi
    
    echo "找到 $antigen_count 个抗原"
    
    # 使用parallel并行处理所有抗原
    printf "%s\n" "${antigen_dirs[@]}" | \
    parallel -j 64 --progress --eta --bar --joblog "${dir_path}/pdb_joblog_${dir_name}.txt" \
        process_antigen {}
    
    echo "=========================================="
    echo "完成处理目录: $dir_name"
    echo "作业日志保存在: ${dir_path}/pdb_joblog_${dir_name}.txt"
    echo ""
}

# 导出函数和变量以便在parallel中使用
export -f process_directory
export BASE_DIR

# 使用parallel并行处理所有目录
echo "开始处理 ${#DIRECTORIES[@]} 个抗原目录..."
echo "并行处理中..."

# 为每个目录创建处理任务
for dir in "${DIRECTORIES[@]}"; do
    process_directory "$dir" &
done

# 等待所有后台任务完成
wait

echo "所有抗原PDB标准化处理完成!"
