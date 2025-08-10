#!/bin/bash

# 设置输入和输出目录
RAW_DIR="./raw"
ANTIBODY_DIR="./antibody"

# 创建样本列表
SAMPLES=("DLE1" "DLE2" "DLE3" "DLE4" "SLE1" "SLE2" "SLE3")

# 创建处理单个样本的函数
process_sample() {
    local sample=$1
    local group=${sample:0:3}  # 提取前3个字符作为组名 (DLE/SLE)
    
    echo "开始处理样本: $sample (组: $group)"
    
    # 步骤1: 准备数据
    python ./src/GetBcr_seq_information.py prepare \
        -i "$RAW_DIR/$sample/matrix.csv" \
        -o "$ANTIBODY_DIR/$group/$sample/prepared_sequences.csv"
    
    # 步骤2: 配对抗体
    python ./src/GetBcr_seq_information.py pair \
        -i "$ANTIBODY_DIR/$group/$sample/prepared_sequences.csv" \
        -o "$ANTIBODY_DIR/$group/$sample/paired_antibodies.fasta"
    
    # 步骤3: 结构预测
    python ./src/GetBcr_seq_information.py fold \
        -f "$ANTIBODY_DIR/$group/$sample/paired_antibodies.fasta" \
        -o "$ANTIBODY_DIR/$group/$sample/all_igfold_models" \
        --refine \
        --renum \
        -l "$ANTIBODY_DIR/$group/$sample/fold.log"
    
    echo "完成处理样本: $sample"
}

# 导出函数以便在parallel中使用
export -f process_sample
export RAW_DIR ANTIBODY_DIR

# 使用parallel并行处理所有样本
# 根据CPU核心数调整并行度 (192核系统可设置较高并行度)
parallel -j 24 --eta --progress --bar --joblog joblog.txt \
    process_sample ::: "${SAMPLES[@]}"

echo "所有样本处理完成!"
