#!/bin/bash

# 设置基础目录
ANTIBODY_DIR="./antibody"

# 创建处理单个克隆型的函数
process_clonotype() {
    local clonotype_dir="$1"
    local clonotype_name=$(basename "$clonotype_dir")
    local pdb_file="${clonotype_dir}/${clonotype_name}.pdb"
    local output_dir="${clonotype_dir}/paratope_meta"
    
    # 创建输出目录
    mkdir -p "$output_dir"
    
    # 执行PDB处理管道
    pdb_reres -1 "$pdb_file" | \
    pdb_chain -A | \
    pdb_chainxseg | \
    pdb_tidy > "${output_dir}/clonetype_renumber.pdb"
    
    # 检查是否成功
    if [ $? -eq 0 ]; then
        echo "✓ 成功处理: $clonotype_dir"
    else
        echo "⚠ 处理失败: $clonotype_dir"
    fi
}

# 导出函数以便在parallel中使用
export -f process_clonotype

# 查找所有需要处理的克隆型目录
find "$ANTIBODY_DIR" -type d -path "*/all_igfold_models/clonotype*" | \
parallel -j 64 --progress --bar --joblog pdb_joblog.txt \
    process_clonotype {}

echo "所有PDB处理完成!"
