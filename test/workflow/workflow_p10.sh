#!/bin/bash

# 设置基础路径
REGISTRY_DIR="./registry"
AMBIG_RESTRAINTS_DIR="$REGISTRY_DIR/ambig_restraints"
ANTIBODY_REGISTRY="$REGISTRY_DIR/antibody_registry.csv"
ANTIGEN_REGISTRY="$REGISTRY_DIR/antigen_registry.csv"
WORKFLOW_TEMPLATE="./db/workflows/docking-antibody-antigen.cfg"
SCRIPT_PATH="./src/generate_epitope_restraints.py"
DOCK_OUTPUT_DIR="./dock"  # 新增：对接结果输出目录

# 检查必要的文件和目录
if [[ ! -f "$ANTIBODY_REGISTRY" ]]; then
    echo "错误: 抗体注册表文件不存在: $ANTIBODY_REGISTRY" >&2
    exit 1
fi

if [[ ! -f "$ANTIGEN_REGISTRY" ]]; then
    echo "错误: 抗原注册表文件不存在: $ANTIGEN_REGISTRY" >&2
    exit 1
fi

if [[ ! -f "$WORKFLOW_TEMPLATE" ]]; then
    echo "错误: 工作流模板文件不存在: $WORKFLOW_TEMPLATE" >&2
    exit 1
fi

if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "错误: Python脚本不存在: $SCRIPT_PATH" >&2
    exit 1
fi

# 创建日志目录
LOG_DIR="$REGISTRY_DIR/configure_logs"
mkdir -p "$LOG_DIR"
JOB_LOG="$LOG_DIR/configure_joblog_$(date +%Y%m%d_%H%M%S).txt"

# 确保对接输出目录存在
mkdir -p "$DOCK_OUTPUT_DIR"

# 导出变量用于子进程
export WORKFLOW_TEMPLATE
export SCRIPT_PATH
export ANTIBODY_REGISTRY
export ANTIGEN_REGISTRY
export DOCK_OUTPUT_DIR

# 定义处理函数
process_ambig() {
    local ambig_file="$1"
    
    # 从路径中提取 clonotype_id 和 antigen_id
    local base_dir=$(dirname "$ambig_file")
    local antigen_id=$(basename "$base_dir")
    local clonotype_dir=$(dirname "$base_dir")
    local clonotype_id=$(basename "$clonotype_dir")
    
    # 记录处理信息
    echo "处理: $clonotype_id / $antigen_id"
    
    # 从抗体注册表中提取信息
    local antibody_info=$(awk -F, -v id="$clonotype_id" '$1 == id' "$ANTIBODY_REGISTRY")
    
    if [[ -z "$antibody_info" ]]; then
        echo "错误: 在抗体注册表中找不到 $clonotype_id" >&2
        return 1
    fi
    
    # 解析抗体信息
    IFS=, read -r clonotype_id renumber_pdb active_passive unambig_tbl <<< "$antibody_info"
    
    # 检查文件是否存在
    if [[ ! -f "$renumber_pdb" ]]; then
        echo "错误: 抗体 PDB 文件不存在: $renumber_pdb" >&2
        return 1
    fi
    
    if [[ ! -f "$unambig_tbl" ]]; then
        echo "错误: 抗体约束文件不存在: $unambig_tbl" >&2
        return 1
    fi
    
    # 从抗原注册表中提取信息
    local antigen_info=$(awk -F, -v id="$antigen_id" '$1 == id' "$ANTIGEN_REGISTRY")
    
    if [[ -z "$antigen_info" ]]; then
        echo "错误: 在抗原注册表中找不到 $antigen_id" >&2
        return 1
    fi
    
    # 解析抗原信息
    IFS=, read -r antigen_id clean_pdb antigen_active_passive <<< "$antigen_info"
    
    # 检查文件是否存在
    if [[ ! -f "$clean_pdb" ]]; then
        echo "错误: 抗原 PDB 文件不存在: $clean_pdb" >&2
        return 1
    fi
    
    # 设置运行目录 - 使用安全路径格式
    local run_dir="${DOCK_OUTPUT_DIR}/${clonotype_id}/${antigen_id}"
    mkdir -p "$run_dir"
    
    # 设置输出配置文件路径 - 使用绝对路径
    local output_cfg="${base_dir}/docking-antibody-antigen_updated.cfg"
    
    # 创建输出目录（如果不存在）
    mkdir -p "$base_dir"
    
    ############################################################
    # 关键修复：使用数组形式传递参数，正确处理路径中的特殊字符
    ############################################################
    local cmd_args=(
        "configure-docking"
        "$WORKFLOW_TEMPLATE"
        "--run_dir" "$run_dir"
        "--output" "$output_cfg"
        "--antibody" "$renumber_pdb"
        "--antigen" "$clean_pdb"
        "--ambig_fname" "$ambig_file"
        "--unambig_fname" "$unambig_tbl"
        "--sampling" "100"
        "--select" "20"
        "--top_models" "3"
    )
    
    # 打印完整命令以便调试
    echo "执行命令: python \"$SCRIPT_PATH\" ${cmd_args[@]}"
    
    # 执行命令
    python "$SCRIPT_PATH" "${cmd_args[@]}"
    
    # 检查命令执行状态
    local status=$?
    if [ $status -ne 0 ]; then
        echo "错误: 配置命令执行失败: $clonotype_id / $antigen_id (状态码: $status)" >&2
        
        # 尝试不带参数的运行以检查脚本是否可执行
        echo "测试脚本可用性: python \"$SCRIPT_PATH\" --version"
        python "$SCRIPT_PATH" --version
        test_status=$?
        echo "脚本测试状态码: $test_status"
        
        return $status
    fi
    
    # 检查输出文件
    if [[ ! -f "$output_cfg" ]]; then
        echo "错误: 配置文件未生成: $output_cfg" >&2
        return 1
    fi
    
    # 检查运行目录是否创建成功
    if [[ ! -d "$run_dir" ]]; then
        echo "错误: 运行目录未创建: $run_dir" >&2
        return 1
    fi
    
    # 记录成功
    echo "✓ 成功生成配置文件: $output_cfg"
    echo "✓ 运行目录创建于: $run_dir"
    
    return 0
}

# 导出处理函数
export -f process_ambig

# 主处理函数
main() {
    # 查找所有 ambig 约束文件
    local ambig_files=($(find "$AMBIG_RESTRAINTS_DIR" -name "ambig-paratope-NMR-epitope.tbl"))
    
    if [[ ${#ambig_files[@]} -eq 0 ]]; then
        echo "警告: 未找到任何 ambig 约束文件" >&2
        return 0
    fi
    
    # 统计文件数量
    local total_files=${#ambig_files[@]}
    echo "找到 $total_files 个 ambig 约束文件"
    echo "开始处理..."
    echo "对接结果将输出到: $DOCK_OUTPUT_DIR"
    
    # 使用 parallel 并行处理
    parallel --progress --bar --eta -j $(nproc) \
        --joblog "$JOB_LOG" \
        --retries 2 \
        process_ambig ::: "${ambig_files[@]}"
    
    # 检查处理结果
    local exit_status=$?
    
    # 生成摘要报告
    echo ""
    echo "===== 处理摘要 ====="
    echo "总约束文件数: $total_files"
    echo "对接结果目录: $DOCK_OUTPUT_DIR"
    
    if [ $exit_status -eq 0 ]; then
        echo "状态: ✅ 所有文件处理成功"
    else
        echo "状态: ⚠ 部分文件处理失败 (退出码: $exit_status)"
        
        # 生成失败任务报告
        local failed_log="${LOG_DIR}/failed_jobs_$(date +%Y%m%d_%H%M%S).txt"
        awk -F '\t' '$7 != 0 {print "任务 " $1 " 失败 (退出码: " $7 ") - " $9}' "$JOB_LOG" > "$failed_log"
        echo "失败任务详情: $failed_log"
    fi
    
    echo "作业日志: $JOB_LOG"
    echo "====================="
    
    return $exit_status
}

# 执行主函数
main
