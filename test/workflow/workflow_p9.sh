#!/bin/bash

# 设置基础目录
REGISTRY_DIR="./registry"
AMBIG_DIR="./registry/ambig_restraints"
SAMPLE_SIZE=${1:-1500}  # 可从命令行指定抽样数量，默认1500
SEED=${2:-42}           # 可从命令行指定随机种子，默认42

# 输入注册表文件
ANTIBODY_REGISTRY="$REGISTRY_DIR/antibody_registry.csv"
ANTIGEN_REGISTRY="$REGISTRY_DIR/antigen_registry.csv"

# 创建输出目录
mkdir -p "$AMBIG_DIR"

# 检查注册表文件是否存在
if [[ ! -f "$ANTIBODY_REGISTRY" ]]; then
    echo "错误: 抗体注册表文件不存在: $ANTIBODY_REGISTRY" >&2
    exit 1
fi

if [[ ! -f "$ANTIGEN_REGISTRY" ]]; then
    echo "错误: 抗原注册表文件不存在: $ANTIGEN_REGISTRY" >&2
    exit 1
fi

# 导出变量用于子进程
export AMBIG_DIR

# 定义处理函数
process_pair() {
    local antibody_line="$1"
    local antigen_line="$2"
    
    # 确保抗体行非空
    if [[ -z "$antibody_line" ]]; then
        echo "错误: 抗体行为空" >&2
        return 1
    fi
    
    # 使用 awk 安全解析 CSV 行（处理带逗号的字段）
    local clonotype_id=$(echo "$antibody_line" | awk -F, '{print $1}')
    local antibody_act_pass=$(echo "$antibody_line" | awk -F, '{print $3}')
    local antigen_id=$(echo "$antigen_line" | awk -F, '{print $1}')
    local antigen_act_pass=$(echo "$antigen_line" | awk -F, '{print $3}')
    
    # 验证文件路径
    if [[ ! -f "$antibody_act_pass" ]]; then
        echo "错误: 抗体 active_passive 文件不存在: $antibody_act_pass" >&2
        return 1
    fi
    
    if [[ ! -f "$antigen_act_pass" ]]; then
        echo "错误: 抗原 active_passive 文件不存在: $antigen_act_pass" >&2
        return 1
    fi
    
    # 创建输出目录
    local output_dir="$AMBIG_DIR/$clonotype_id/$antigen_id"
    mkdir -p "$output_dir"
    
    # 设置输出文件路径
    local output_tbl="$output_dir/ambig-paratope-NMR-epitope.tbl"
    
    # 生成 ambig 约束文件（带错误处理）
    if ! haddock3-restraints active_passive_to_ambig \
        "$antibody_act_pass" \
        "$antigen_act_pass" \
        --segid-one A \
        --segid-two B > "$output_tbl" 2> "$output_dir/haddock-error.log"; then
        
        echo "生成 ambig 约束文件失败: $clonotype_id/$antigen_id" >&2
        echo "详见: $output_dir/haddock-error.log" >&2
        return 1
    fi
    
    # 验证约束文件
    if ! haddock3-restraints validate_tbl "$output_tbl" --silent 2> "$output_dir/validate-error.log"; then
        echo "验证约束文件失败: $output_tbl" >&2
        echo "详见: $output_dir/validate-error.log" >&2
        return 1
    fi
    
    # 检查 TBL 文件是否非空
    if [[ ! -s "$output_tbl" ]]; then
        echo "错误: 生成的 TBL 文件为空: $output_tbl" >&2
        return 1
    fi
    
    return 0
}

# 导出处理函数
export -f process_pair

# 分层抽样函数（带随机种子）
sample_clonotypes() {
    local sample_size=$1
    local seed=$2
    
    # 临时文件
    local patient_temp=$(mktemp)
    local sampled_temp=$(mktemp)
    
    # 提取患者信息（假设第一列为患者ID格式为DLE1_clonotype1）
    echo "提取患者分布信息..." >&2
    # 跳过标题行，只处理数据行
    tail -n +2 "$ANTIBODY_REGISTRY" | awk -F, '{split($1, a, "_"); patient=a[1]; print patient}' | \
        sort | uniq -c > "$patient_temp"
    
    # 计算总克隆数
    total_clonotypes=$(awk '{sum += $1} END {print sum}' "$patient_temp")
    echo "总克隆型数量: $total_clonotypes" >&2
    echo "随机种子: $seed" >&2
    
    # 计算每个患者应抽取的数量
    awk -v total="$total_clonotypes" -v sample="$sample_size" '
    {
        count = $1
        patient = $2
        proportion = count / total
        sample_count = int(sample * proportion + 0.5)  # 四舍五入
        
        # 确保至少抽1个
        if (sample_count < 1) sample_count = 1
        
        # 保存分配结果
        printf "%s\t%d\n", patient, sample_count
    }' "$patient_temp" > "${patient_temp}_alloc"
    
    # 验证抽样总数
    allocated_total=$(awk '{sum += $2} END {print sum}' "${patient_temp}_alloc")
    if [ "$allocated_total" -ne "$sample_size" ]; then
        echo "调整抽样数量: $allocated_total → $sample_size" >&2
        # 调整最后一个患者的数量以匹配总数
        last_patient=$(tail -1 "${patient_temp}_alloc" | cut -f1)
        last_count=$(tail -1 "${patient_temp}_alloc" | cut -f2)
        adjusted=$((last_count + sample_size - allocated_total))
        sed -i "\$s/\t[0-9]*$/\t$adjusted/" "${patient_temp}_alloc"
    fi
    
    # 从每个患者中随机抽样（使用固定种子）
    echo "执行分层抽样（使用种子 $seed）..." >&2
    while read -r patient count; do
        # 提取该患者的所有克隆型（跳过标题行）
        # 使用精确匹配：确保患者ID在第一个字段中
        awk -F, -v patient="$patient" 'NR>1 && $1 ~ "^" patient "_" {print $0}' "$ANTIBODY_REGISTRY" | \
            awk -v seed="$seed" -v patient="$patient" '
            BEGIN {
                # 生成基于患者和主种子的唯一子种子
                srand(seed + length(patient) + index("ABCDEFGHIJKLMNOPQRSTUVWXYZ", substr(patient,1,1)))
            }
            { 
                # 为每行生成随机数
                print rand() "\t" $0 
            }' | \
            sort -k1,1n | \
            head -n "$count" | \
            cut -f2- >> "$sampled_temp"
    done < "${patient_temp}_alloc"
    
    # 返回抽样结果（只包含有效数据行）
    cat "$sampled_temp"
    
    # 清理临时文件
    rm -f "$patient_temp" "${patient_temp}_alloc" "$sampled_temp"
}

# 主处理函数
generate_ambig_restraints() {
    # 创建临时文件
    local antibody_temp=$(mktemp)
    local antigen_temp=$(mktemp)
    
    # 执行分层抽样
    echo "开始分层抽样，目标抽样数: $SAMPLE_SIZE" >&2
    sample_clonotypes "$SAMPLE_SIZE" "$SEED" > "$antibody_temp"
    
    # 检查抽样结果是否为空
    if [[ ! -s "$antibody_temp" ]]; then
        echo "错误: 抽样结果为空，请检查抗体注册表文件" >&2
        exit 1
    fi
    
    # 保存抽样列表（用于复现）
    local sample_list="$AMBIG_DIR/sampled_clonotypes_${SAMPLE_SIZE}_seed${SEED}.txt"
    cp "$antibody_temp" "$sample_list"
    echo "抽样克隆型列表已保存: $sample_list" >&2
    
    # 获取抗原列表（跳过标题行）
    tail -n +2 "$ANTIGEN_REGISTRY" > "$antigen_temp"
    
    # 检查抗原列表是否为空
    if [[ ! -s "$antigen_temp" ]]; then
        echo "错误: 抗原列表为空，请检查抗原注册表文件" >&2
        exit 1
    fi
    
    # 计算任务数
    local antibody_count=$(wc -l < "$antibody_temp")
    local antigen_count=$(wc -l < "$antigen_temp")
    local total_tasks=$((antibody_count * antigen_count))
    
    echo "抽样克隆型数量: $antibody_count" >&2
    echo "抗原数量: $antigen_count" >&2
    echo "总组合数: $total_tasks" >&2
    echo "并行处理中..." >&2
    
    # 使用 parallel 并行处理所有组合
    parallel --progress --bar --eta -j $(nproc) \
        --joblog "$AMBIG_DIR/parallel_joblog.txt" \
        --retries 2 \
        process_pair {1} {2} :::: "$antibody_temp" :::: "$antigen_temp"
    
    # 检查并行处理结果
    local exit_status=$?
    
    # 生成错误报告
    if [ $exit_status -ne 0 ]; then
        echo -e "\n⚠ 部分任务处理失败，生成错误报告..." >&2
        awk -F '\t' '$7 != 0 {print "Job " $1 " failed with exit code " $7 ": " $9}' \
            "$AMBIG_DIR/parallel_joblog.txt" > "$AMBIG_DIR/failed_jobs.txt"
        echo "失败任务列表: $AMBIG_DIR/failed_jobs.txt" >&2
    fi
    
    # 清理临时文件
    rm -f "$antibody_temp" "$antigen_temp"
    
    if [ $exit_status -eq 0 ]; then
        echo -e "\n所有 ambig 约束文件生成完成!" >&2
        echo "输出目录: $AMBIG_DIR" >&2
        echo "随机种子: $SEED" >&2
        echo "抽样列表: $sample_list" >&2
    else
        echo -e "\n⚠ 部分任务处理失败，请检查错误报告" >&2
        return $exit_status
    fi
}

# 执行主函数
generate_ambig_restraints
