#!/bin/bash

# 注册表生成脚本（并行优化版）
# 使用 GNU Parallel 加速抗体和抗原的文件路径注册表创建

# 设置基础目录
ANTIBODY_BASE="./antibody"
ANTIGEN_BASE="./antigen"
REGISTRY_DIR="./registry"

# 创建注册表目录
mkdir -p "$REGISTRY_DIR"

# 初始化注册表文件
ANTIBODY_REGISTRY="$REGISTRY_DIR/antibody_registry.csv"
ANTIGEN_REGISTRY="$REGISTRY_DIR/antigen_registry.csv"

echo "clonotype_id,renumber_pdb,act_pass,unambig_tbl" > "$ANTIBODY_REGISTRY"
echo "antigen_id,clean_pdb,act_pass" > "$ANTIGEN_REGISTRY"

# 获取绝对路径函数
get_abs_path() {
    echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}

# 导出函数以便在子进程中使用
export -f get_abs_path

# 并行处理抗体文件
process_antibodies() {
    # 查找所有重编号PDB文件
    find "$ANTIBODY_BASE" -type f -name "clonetype_renumber.pdb" | \
    parallel --progress --bar -j $(nproc) '
        renumber_pdb={}
        meta_dir=$(dirname "$renumber_pdb")
        clonotype_dir=$(dirname "$meta_dir")
        clonotype_id=$(basename "$clonotype_dir")
        patient_id=$(basename "$(dirname "$(dirname "$clonotype_dir")")")
        disease_type=$(basename "$(dirname "$(dirname "$(dirname "$clonotype_dir")")")")
        
        full_clonotype_id="${disease_type}_${patient_id}_${clonotype_id}"
        
        act_pass="${meta_dir}/antibody-paratope.act-pass"
        unambig_tbl="${meta_dir}/antibody-unambig.tbl"
        
        if [[ ! -f "$act_pass" ]]; then
            echo "⚠ 缺少抗体paratope文件: $act_pass" >&2
            exit 0
        fi
        
        if [[ ! -f "$unambig_tbl" ]]; then
            echo "⚠ 缺少抗体非歧义约束文件: $unambig_tbl" >&2
            exit 0
        fi
        
        abs_renumber=$(get_abs_path "$renumber_pdb")
        abs_act_pass=$(get_abs_path "$act_pass")
        abs_unambig_tbl=$(get_abs_path "$unambig_tbl")
        
        echo "${full_clonotype_id},${abs_renumber},${abs_act_pass},${abs_unambig_tbl}"
    ' >> "$ANTIBODY_REGISTRY"
    
    local count=$(wc -l < "$ANTIBODY_REGISTRY")
    echo "抗体注册表完成: 添加 $((count-1)) 个克隆型"  # 减1排除标题行
}

# 并行处理抗原文件
process_antigens() {
    # 查找所有抗原act-pass文件
    find "$ANTIGEN_BASE" -type f -name "antigen-epitope.act-pass" | \
    parallel --progress --bar -j $(nproc) '
        act_pass={}
        antigen_dir=$(dirname "$act_pass")
        antigen_id=$(basename "$antigen_dir")
        antigen_type=$(basename "$(dirname "$(dirname "$antigen_dir")")")
        
        full_antigen_id="${antigen_type}_${antigen_id}"
        
        clean_pdb=$(find "$antigen_dir" -maxdepth 1 -type f \( -iname "*clean.pdb" -o -iname "${antigen_id}_clean.pdb" \) | head -1)
        
        if [[ -z "$clean_pdb" ]]; then
            echo "⚠ 缺少抗原clean.pdb文件: $antigen_dir" >&2
            exit 0
        fi
        
        abs_clean_pdb=$(get_abs_path "$clean_pdb")
        abs_act_pass=$(get_abs_path "$act_pass")
        
        echo "${full_antigen_id},${abs_clean_pdb},${abs_act_pass}"
    ' >> "$ANTIGEN_REGISTRY"
    
    local count=$(wc -l < "$ANTIGEN_REGISTRY")
    echo "抗原注册表完成: 添加 $((count-1)) 个抗原"  # 减1排除标题行
}

# 主执行流程
echo "开始并行生成抗体注册表..."
process_antibodies

echo -e "\n开始并行生成抗原注册表..."
process_antigens

echo -e "\n注册表生成完成!"
echo "抗体注册表: $ANTIBODY_REGISTRY ($(wc -l < "$ANTIBODY_REGISTRY") 行)"
echo "抗原注册表: $ANTIGEN_REGISTRY ($(wc -l < "$ANTIGEN_REGISTRY") 行)"
