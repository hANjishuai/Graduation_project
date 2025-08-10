#!/bin/bash

# 抗体非歧义性约束生成脚本
# 为所有克隆型的重编号PDB文件生成约束

# 设置基础目录
ANTIBODY_DIR="./antibody"

# 创建处理单个克隆型的函数
generate_restraints() {
    local renumber_pdb="$1"
    local output_tbl="${renumber_pdb%/*}/antibody-unambig.tbl"
    
    # 提取克隆型名称用于日志
    local clonotype_dir=$(dirname "$(dirname "$renumber_pdb")")
    local clonotype_name=$(basename "$clonotype_dir")
    
    echo "处理克隆型: $clonotype_name"
    echo "• 输入文件: $(basename "$renumber_pdb")"
    
    # 运行约束生成命令
    haddock3-restraints restrain_bodies "$renumber_pdb" > "$output_tbl"
    
    # 检查执行状态
    if [ $? -eq 0 ]; then
        # 检查输出文件是否非空
        if [ -s "$output_tbl" ]; then
            echo "✓ 成功生成约束: $clonotype_name"
            echo "• 输出文件: antibody-unambig.tbl ($(wc -l < "$output_tbl") 行约束)"
            return 0
        else
            echo "⚠ 输出文件为空: $output_tbl"
            return 1
        fi
    else
        echo "⚠ 约束生成失败: $clonotype_name"
        return 1
    fi
}

# 导出函数以便在parallel中使用
export -f generate_restraints

# 查找所有重编号PDB文件
echo "搜索重编号PDB文件..."
renumber_files=($(find "$ANTIBODY_DIR" -type f -path "*/paratope_meta/clonetype_renumber.pdb"))
file_count=${#renumber_files[@]}

if [ $file_count -eq 0 ]; then
    echo "⚠ 未找到任何 clonetype_renumber.pdb 文件"
    exit 1
fi

echo "找到 $file_count 个重编号PDB文件"
echo "使用 64 个并行任务处理..."

# 使用parallel并行处理所有文件
printf "%s\n" "${renumber_files[@]}" | \
parallel -j 64 --progress --eta --bar --joblog restraints_joblog.txt \
    "generate_restraints {}"

# 检查处理结果
success_count=$(grep -c "成功生成约束" restraints_joblog.txt || true)
failed_count=$(grep -c "约束生成失败" restraints_joblog.txt || true)
empty_count=$(grep -c "输出文件为空" restraints_joblog.txt || true)

echo ""
echo "处理完成:"
echo "• 成功: $success_count"
echo "• 失败: $failed_count"
echo "• 空输出: $empty_count"
echo "• 作业日志: restraints_joblog.txt"

echo "所有抗体约束生成完成!"
