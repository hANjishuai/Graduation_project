#!/bin/bash

# 抗原处理脚本
# 用于处理 Skin_celltype 和 Skin_tissue 目录

# 设置基础目录
BASE_DIR="./antigen"
SCRIPT_PATH="./src/Get_epitope.py"  # 确保这是正确的脚本路径

# 要处理的目录列表
DIRECTORIES=("Skin_celltype" "Skin_tissue")

# 阈值设置
THRESHOLD=0.45

# 创建处理单个目录的函数
process_directory() {
    local dir_name="$1"
    local dir_path="${BASE_DIR}/${dir_name}"
    
    echo "开始处理目录: $dir_name"
    echo "=========================================="
       
    # 执行抗原筛选命令
    python "${SCRIPT_PATH}" filter-antigens \
        "${dir_path}" \
        --threshold $THRESHOLD
    
    # 检查执行状态
    if [ $? -eq 0 ]; then
        echo "✓ 成功处理: $dir_name"
    else
        echo "⚠ 处理失败: $dir_name"
    fi
    
    echo "=========================================="
    echo ""
}

# 导出函数和变量以便在parallel中使用
export -f process_directory
export BASE_DIR SCRIPT_PATH THRESHOLD

# 使用parallel并行处理所有目录
# 根据CPU核心数调整并行度
echo "开始处理 ${#DIRECTORIES[@]} 个目录..."
echo "阈值: $THRESHOLD"
echo "并行处理..."

parallel -j 2 --progress --eta --bar --joblog antigen_joblog.txt \
    process_directory ::: "${DIRECTORIES[@]}"

echo "所有抗原目录处理完成!"
echo "作业日志保存在: antigen_joblog.txt"
