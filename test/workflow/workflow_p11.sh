#!/bin/bash

# 并行对接脚本
# 用法: ./run_docking.sh [并发任务数] [清理模式]

# 设置基础路径
REGISTRY_DIR="./registry"
AMBIG_RESTRAINTS_DIR="$REGISTRY_DIR/ambig_restraints"
LOG_DIR="$REGISTRY_DIR/docking_logs"
CFG_PATTERN="docking-antibody-antigen_updated.cfg"
MASTER_SUMMARY="$LOG_DIR/docking_summary_master.csv"  # 主摘要文件

# 设置默认并发任务数
DEFAULT_JOBS=$(( $(nproc) / 4 ))  # 降低并发度为CPU核心数的1/4
MAX_JOBS=${1:-$DEFAULT_JOBS}
CLEAN_MODE=${2:-safe}  # 默认安全模式: safe | force

# 创建日志目录
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
JOB_LOG="$LOG_DIR/docking_joblog_$TIMESTAMP.txt"
SESSION_SUMMARY="$LOG_DIR/docking_summary_$TIMESTAMP.csv"  # 本次会话摘要

# 初始化主摘要文件（如果不存在）
if [[ ! -f "$MASTER_SUMMARY" ]]; then
    echo "cfg_path,status,time,output_dir,timestamp" > "$MASTER_SUMMARY"
fi

# 查找所有配置文件
echo "正在搜索配置文件..."
cfg_files=($(find "$AMBIG_RESTRAINTS_DIR" -type f -name "$CFG_PATTERN"))
total_jobs=${#cfg_files[@]}

if [[ $total_jobs -eq 0 ]]; then
    echo "错误: 未找到任何配置文件 ($CFG_PATTERN)" >&2
    exit 1
fi

echo "找到 $total_jobs 个对接任务"
echo "最大并发任务数: $MAX_JOBS"
echo "清理模式: $CLEAN_MODE"

# 从主摘要文件中加载已完成任务
declare -A completed_jobs
if [[ -s "$MASTER_SUMMARY" ]]; then
    echo "从主摘要文件中加载已完成任务: $MASTER_SUMMARY"
    while IFS=',' read -r cfg_path status _ run_dir _; do
        if [[ "$cfg_path" == "cfg_path" ]]; then continue; fi  # 跳过标题行
        completed_jobs["$cfg_path"]=1
        echo "  已记录: $cfg_path ($status)"
    done < "$MASTER_SUMMARY"
fi

# 过滤出待处理任务
pending_jobs=()
for cfg_file in "${cfg_files[@]}"; do
    if [[ -z "${completed_jobs[$cfg_file]}" ]]; then
        pending_jobs+=("$cfg_file")
    else
        echo "跳过已记录任务: $cfg_file"
    fi
done

pending_count=${#pending_jobs[@]}
completed_count=$((total_jobs - pending_count))

echo "待处理任务: $pending_count"
echo "已记录任务: $completed_count"

if [[ $pending_count -eq 0 ]]; then
    echo "所有任务已记录，无需执行"
    exit 0
fi

# 定义带锁机制的对接函数
run_docking() {
    local cfg_file="$1"
    local log_id=$(basename $(dirname $(dirname "$cfg_file")))
    local log_file="$LOG_DIR/${log_id}_${TIMESTAMP}.log"
    
#    # 提取运行目录
#    local run_dir=$(awk -F'=' '/^run_dir/ {gsub(/[[:space:]]+/, "", $2); print $2}' "$cfg_file")

     # 提取运行目录（去掉两侧双引号）
    local run_dir
    run_dir=$(sed -n 's/^[[:space:]]*run_dir[[:space:]]*=[[:space:]]*\(.*\)[[:space:]]*$/\1/p' "$cfg_file" | sed 's/^["'\'']//;s/["'\'']$//')

    local lock_file="$run_dir/.haddock_lock"
    
    # 检查主摘要文件中是否已有记录（双重检查）
    if grep -q "^$cfg_file," "$MASTER_SUMMARY"; then
        echo "✅ 跳过: $cfg_file (已在主摘要中记录)" | tee -a "$log_file"
        return 0
    fi
    
    # 清理残留目录（如果需要）
    if [ -d "$run_dir" ]; then
        if [ "$CLEAN_MODE" == "force" ]; then
            echo "⚠️ 强制清理目录: $run_dir" | tee -a "$log_file"
            rm -rf "$run_dir"
        else
            echo "⚠️ 安全清理目录: $run_dir" | tee -a "$log_file"
            find "$run_dir" -maxdepth 1 -type f \( -name "*.inp" -o -name "*.out" -o -name "*.err" -o -name "*.log" -o -name ".haddock_lock" \) -delete
            find "$run_dir" -type d -empty -delete
        fi
    fi
    
    # 如果清理后目录仍然存在，则可能是非空目录，报错
    if [ -d "$run_dir" ]; then
        echo "❌ 错误: 运行目录非空且无法清理: $run_dir" | tee -a "$log_file"
        # 记录到本次会话摘要
        echo "$cfg_file,failed-clean,0,$run_dir,$TIMESTAMP" >> "$SESSION_SUMMARY"
        return 1
    fi
    
    # 创建运行目录的父目录
    mkdir -p "$(dirname "$run_dir")"
    
    # 目录锁检查
    if [ -f "$lock_file" ]; then
        echo "⚠️ 跳过: $cfg_file (已在运行)" | tee -a "$log_file"
        return 1
    fi
    touch "$lock_file"
    
    # 资源检查 - 避免系统过载
    local load=$(awk '{print $1}' /proc/loadavg | cut -d. -f1)
    local mem_avail=$(free -m | awk '/Mem:/ {print $7}')
    if [ $load -gt $(nproc) ] || [ $mem_avail -lt 2048 ]; then
        echo "⚠️ 资源紧张: 负载 $load | 内存 ${mem_avail}MB - 随机延迟" | tee -a "$log_file"
        sleep $(( RANDOM % 30 + 10 )) # 随机延迟10-40秒
    fi
    
    echo "启动对接: $cfg_file" | tee -a "$log_file"
    echo "运行目录: $run_dir" | tee -a "$log_file"
    
    local start_time=$(date +%s)
    
    # 运行对接（带超时防止卡住）
    timeout 12h haddock3 "$cfg_file" > "$log_file" 2>&1
    local status=$?
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # 释放锁
    rm -f "$lock_file"
    
    # 记录结果到本次会话摘要
    if [[ $status -eq 0 ]]; then
        echo "✅ 成功: $cfg_file (用时: ${duration}秒)" | tee -a "$log_file"
        echo "$cfg_file,success,$duration,$run_dir,$TIMESTAMP" >> "$SESSION_SUMMARY"
    else
        echo "❌ 失败[代码:$status]: $cfg_file (用时: ${duration}秒)" | tee -a "$log_file"
        echo "$cfg_file,failed-$status,$duration,$run_dir,$TIMESTAMP" >> "$SESSION_SUMMARY"
        
        # 提取错误信息
        local error_msg=$(grep -i -m1 -e "error" -e "fail" -e "exception" "$log_file" | head -1)
        echo "   错误信息: ${error_msg:0:120}..." | tee -a "$log_file"
    fi
    
    return $status
}

# 导出函数用于parallel
export -f run_docking
export TIMESTAMP
export LOG_DIR
export SESSION_SUMMARY
export MASTER_SUMMARY
export CLEAN_MODE

echo "开始并行对接..."

# 禁用GNU Parallel的引用提示
export PARALLEL_HOME="$LOG_DIR"
mkdir -p "$PARALLEL_HOME"
echo '--will-cite' > "$PARALLEL_HOME/config"

# 使用parallel并行运行
parallel --link --progress --bar --eta -j "$MAX_JOBS" \
    --joblog "$JOB_LOG" \
    --retries 1 \
    --resume-failed \
    --halt soon,fail=10% \
    --delay 0.1 \
    run_docking ::: "${pending_jobs[@]}" </dev/null

# 将本次会话结果合并到主摘要文件
echo "合并结果到主摘要文件: $MASTER_SUMMARY"
grep -v "cfg_path" "$SESSION_SUMMARY" >> "$MASTER_SUMMARY"  # 追加新结果，跳过标题行

# 生成摘要报告
session_success=$(grep -c "success" "$SESSION_SUMMARY")
session_failed=$(grep -c "failed" "$SESSION_SUMMARY")
session_skipped=$((pending_count - session_success - session_failed))

master_total=$(($(wc -l < "$MASTER_SUMMARY") - 1))  # 减去标题行
master_success=$(grep -c "success" "$MASTER_SUMMARY")
master_failed=$((master_total - master_success))

echo ""
echo "===== 本次会话摘要 ====="
echo "待处理任务: $pending_count"
echo "本次成功: $session_success"
echo "本次失败: $session_failed"
echo "本次跳过: $session_skipped"
echo "本次用时: $(awk -F',' 'NR>1 {sum+=$3} END {print sum "秒 (" int(sum/60) "分)"}' "$SESSION_SUMMARY")"

echo ""
echo "===== 总体进度摘要 ====="
echo "总任务数: $total_jobs"
echo "已记录任务: $master_total"
echo "总体成功: $master_success"
echo "总体失败: $master_failed"
echo "剩余任务: $((total_jobs - master_total))"
echo ""
echo "详细日志: $JOB_LOG"
echo "本次会话摘要: $SESSION_SUMMARY"
echo "主摘要文件: $MASTER_SUMMARY"

if [[ $session_failed -gt 0 ]]; then
    echo ""
    echo "失败任务列表:"
    grep "failed" "$SESSION_SUMMARY" | cut -d',' -f1
    echo ""
    echo "请检查失败任务的日志文件获取详细信息"
    
    # 如果有失败任务，建议使用强制清理模式
    if [ "$CLEAN_MODE" != "force" ]; then
        echo ""
        echo "提示: 部分失败可能是由于残留文件导致，可以尝试使用强制清理模式重新运行:"
        echo "      $0 $MAX_JOBS force"
    fi
fi

exit 0



#使用方式：
#
#    首次运行：
#
#    bash
#
#复制代码
#./run_docking.sh 8
#
#    创建主摘要文件 docking_summary_master.csv
#    记录所有运行过的任务状态
#
#后续运行：
#
#bash
#
#复制代码
## 再次运行，自动跳过主摘要文件中已记录的任务
#./run_docking.sh 8
#
#    只运行新任务或之前失败的任务
#    结果自动追加到主摘要文件
#
#强制重试失败任务：
#
#bash
#
#    复制代码
#    # 首先从主摘要文件中删除失败记录
#    grep -v "failed" $LOG_DIR/docking_summary_master.csv > $LOG_DIR/temp.csv
#    mv $LOG_DIR/temp.csv $LOG_DIR/docking_summary_master.csv
#
#    # 然后使用force模式重新运行
#    ./run_docking.sh 8 force
#
#文件说明：
#
#    主摘要文件 (docking_summary_master.csv)：
#        持久化存储所有任务的状态
#        格式：cfg_path,status,time,output_dir,timestamp
#        每次运行后追加新结果
#
#    会话摘要文件 (docking_summary_<timestamp>.csv)：
#        记录本次运行的任务结果
#        运行结束后合并到主摘要文件
#
#    作业日志 (docking_joblog_<timestamp>.txt)：
#        记录并行作业的详细日志
#        用于调试和监控
#
#这个系统确保：
#
#    每次运行只处理新任务
#    任务状态持久化存储
#    避免重复运行导致目录冲突
#    提供清晰的进度报告
#    支持中断后继续运行
