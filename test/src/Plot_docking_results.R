#!/usr/bin/env Rscript
# 对接结果数据分析工具
# 版本: 1.1.5
# 更新: 调整抗体顺序(SLE在前DLE在后) + 行聚类(两类)
# 命令行运行
# Rscript analyze_docking_results.R \
#   -i docking_results_with_metadata.csv \
#   -o results

.libPaths("~/R/4.4.1/library/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(openxlsx)
  library(cli)
  library(glue)
  library(argparse)
})

#--------------------------
# 日志系统
#--------------------------
log_step <- function(message, symbol = "🔹") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_info("{.timestamp [{ts}]} {symbol} {message}")
}

log_success <- function(message) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_success("{.timestamp [{ts}]} ✅ {message}")
}

log_warning <- function(message) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli_alert_warning("{.timestamp [{ts}]} ⚠️ {message}")
}

#--------------------------
# 核心功能函数
#--------------------------
read_docking_data <- function(input_file) {
  log_step(glue("读取数据: {input_file}"))
  
  data <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # 检查必要列是否存在
  required_cols <- c("antibody_part2", "antigen_pdb_chain", "desolv", "elec", "total", "vdw", "air")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(glue("输入文件缺少必要列: {paste(missing_cols, collapse=', ')}"))
  }
  
  log_success(glue("读取数据: {input_file} ({nrow(data)} 行)"))
  
  # 转换数值列
  processed_data <- data %>%
    select(all_of(required_cols)) %>%
    mutate(across(c(desolv, elec, total, vdw, air), as.numeric))
  
  return(processed_data)
}

compute_average_energy <- function(data) {
  log_step("计算平均能量值...")
  
  avg_table <- data %>%
    group_by(antigen_pdb_chain, antibody_part2) %>%
    summarise(
      desolv = mean(desolv, na.rm = TRUE),
      elec = mean(elec, na.rm = TRUE),
      total = mean(total, na.rm = TRUE),
      vdw = mean(vdw, na.rm = TRUE),
      air = mean(air, na.rm = TRUE),
      .groups = 'drop'
    )
  
  log_success(glue("创建表: {length(unique(avg_table$antigen_pdb_chain))} 抗原 x {length(unique(avg_table$antibody_part2))} 抗体"))
  
  return(avg_table)
}

generate_normalized_matrices <- function(avg_table) {
  log_step("生成标准化能量矩阵...")
  
  energy_types <- c("desolv", "elec", "total", "vdw", "air")
  matrices <- list()
  
  # 获取抗体类型并排序：SLE在前，DLE在后
  antibody_types <- unique(avg_table$antibody_part2)
  sle_antibodies <- sort(grep("^SLE", antibody_types, value = TRUE))
  dle_antibodies <- sort(grep("^DLE", antibody_types, value = TRUE))
  antibody_order <- c(sle_antibodies, dle_antibodies)  # SLE在前，DLE在后
  
  # 抗原按字母顺序排序（热图中将根据聚类重新排序）
  antigen_order <- sort(unique(avg_table$antigen_pdb_chain))
  
  for (energy in energy_types) {
    # 创建原始矩阵（按指定顺序）
    matrix_df <- avg_table %>%
      select(antigen_pdb_chain, antibody_part2, !!sym(energy)) %>%
      mutate(
        antigen_pdb_chain = factor(antigen_pdb_chain, levels = antigen_order),
        antibody_part2 = factor(antibody_part2, levels = antibody_order)
      ) %>%
      pivot_wider(
        names_from = antibody_part2,
        values_from = !!sym(energy),
        values_fill = NA
      ) %>%
      arrange(antigen_pdb_chain) %>%
      column_to_rownames("antigen_pdb_chain")
    
    # 对矩阵进行行标准化
    normalized_matrix <- t(scale(t(matrix_df)))
    
    matrices[[energy]] <- list(
      raw = matrix_df,
      normalized = normalized_matrix
    )
  }
  
  log_success(glue("生成{length(energy_types)}个标准化能量矩阵"))
  
  return(matrices)
}

plot_heatmaps <- function(matrices, output_dir) {
  log_step("绘制热图（行标准化，行聚类两类）...")
  
  for (energy in names(matrices)) {
    heatmap_data <- matrices[[energy]]$normalized
    
    # 跳过全为NA的矩阵
    if (all(is.na(heatmap_data))) {
      log_warning(glue("能量项 '{energy}' 的数据全为NA，跳过绘图"))
      next
    }
    
    # 设置颜色范围（使用分位数避免异常值）
    color_range <- quantile(heatmap_data, probs = c(0.05, 0.95), na.rm = TRUE)
    
    # 创建热图（启用行聚类并分为两类）
    pheatmap(
      heatmap_data,
      main = glue("{toupper(energy)} Energy (Z-score)"),
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(color_range[1], color_range[2], length.out = 101),
      cluster_rows = TRUE,   # 启用行聚类
      cluster_cols = FALSE,  # 禁用列聚类
      cutree_rows = 2,       # 将行聚类分为两类
      fontsize_row = 8,
      fontsize_col = 8,
      filename = file.path(output_dir, glue("heatmap_{energy}.png")),
      show_rownames = TRUE,
      show_colnames = TRUE
    )
  }
  
  log_success(glue("保存{length(matrices)}个热图"))
}

generate_prism_data <- function(avg_table, output_dir) {
  log_step("准备Prism数据...")
  
  prism_dir <- file.path(output_dir, "prism_data")
  if (!dir.exists(prism_dir)) {
    dir.create(prism_dir, recursive = TRUE)
  }
  
  # 获取抗体类型并排序：SLE在前，DLE在后
  antibody_types <- unique(avg_table$antibody_part2)
  sle_antibodies <- sort(grep("^SLE", antibody_types, value = TRUE))
  dle_antibodies <- sort(grep("^DLE", antibody_types, value = TRUE))
  antibody_order <- c(sle_antibodies, dle_antibodies)  # SLE在前，DLE在后
  
  # 获取所有抗原
  antigens <- unique(avg_table$antigen_pdb_chain)
  
  for (ag in antigens) {
    # 清理抗原名称，确保文件名安全
    clean_ag <- gsub("[^[:alnum:]]", "", ag)
    
    ag_data <- avg_table %>%
      filter(antigen_pdb_chain == ag) %>%
      # 按照抗体顺序排序
      mutate(antibody_part2 = factor(antibody_part2, levels = antibody_order)) %>%
      arrange(antibody_part2) %>%
      select(antibody_part2, desolv, elec, total, vdw, air)
    
    # 转置为Prism格式（行：能量类型，列：抗体）
    transposed_data <- ag_data %>%
      pivot_longer(cols = -antibody_part2) %>%
      pivot_wider(names_from = antibody_part2, values_from = value) %>%
      rename(Energy_Type = name)
    
    write.csv(transposed_data, 
              file.path(prism_dir, glue("prism_{clean_ag}.csv")), 
              row.names = FALSE)
  }
  
  log_success(glue("保存{length(antigens)}个Prism数据文件"))
}

save_results <- function(avg_table, matrices, output_dir) {
  log_step("保存汇总Excel文件")
  
  wb <- createWorkbook()
  
  # 添加平均总表
  addWorksheet(wb, "Average_Energy")
  writeData(wb, "Average_Energy", avg_table)
  
  # 添加能量矩阵（原始值和标准化值）
  for (energy in names(matrices)) {
    # 原始矩阵
    addWorksheet(wb, glue("Matrix_{energy}_raw"))
    writeData(wb, glue("Matrix_{energy}_raw"), matrices[[energy]]$raw, rowNames = TRUE)
    
    # 标准化矩阵
    addWorksheet(wb, glue("Matrix_{energy}_zscore"))
    writeData(wb, glue("Matrix_{energy}_zscore"), matrices[[energy]]$normalized, rowNames = TRUE)
  }
  
  # 添加Prism数据（每个抗原一个工作表）
  antigens <- unique(avg_table$antigen_pdb_chain)
  
  # 获取抗体类型并排序：SLE在前，DLE在后
  antibody_types <- unique(avg_table$antibody_part2)
  sle_antibodies <- sort(grep("^SLE", antibody_types, value = TRUE))
  dle_antibodies <- sort(grep("^DLE", antibody_types, value = TRUE))
  antibody_order <- c(sle_antibodies, dle_antibodies)  # SLE在前，DLE在后
  
  for (i in seq_along(antigens)) {
    ag <- antigens[i]
    
    # 创建唯一的工作表名称（使用索引避免冲突）
    sheet_name <- glue("Prism_{i}")
    
    # 确保工作表名称长度不超过31个字符（Excel限制）
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
    }
    
    ag_data <- avg_table %>%
      filter(antigen_pdb_chain == ag) %>%
      # 按照抗体顺序排序
      mutate(antibody_part2 = factor(antibody_part2, levels = antibody_order)) %>%
      arrange(antibody_part2) %>%
      select(antibody_part2, desolv, elec, total, vdw, air)
    
    transposed_data <- ag_data %>%
      pivot_longer(cols = -antibody_part2) %>%
      pivot_wider(names_from = antibody_part2, values_from = value) %>%
      rename(Energy_Type = name)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, transposed_data)
  }
  
  excel_file <- file.path(output_dir, "energy_results.xlsx")
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  # 保存平均总表为CSV
  write.csv(avg_table, file.path(output_dir, "average_energy.csv"), row.names = FALSE)
  
  log_success(glue("成功保存所有结果到: {output_dir}"))
}

#--------------------------
# 主函数
#--------------------------
main <- function() {
  parser <- ArgumentParser(description='分析对接结果数据并生成热图和Prism数据')
  parser$add_argument('-i', '--input', required=TRUE, help='输入CSV文件路径')
  parser$add_argument('-o', '--output_dir', required=TRUE, help='输出目录')
  args <- parser$parse_args()
  
  input_file <- args$input
  output_dir <- args$output_dir
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cli_h1("对接结果数据分析流程启动")
  log_step(glue("开始时间: {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
  
  tryCatch({
    # 1. 读取数据
    data <- read_docking_data(input_file)
    
    # 2. 计算平均能量值
    avg_table <- compute_average_energy(data)
    
    # 3. 生成标准化能量矩阵
    matrices <- generate_normalized_matrices(avg_table)
    
    # 4. 绘制热图（行标准化，行聚类两类）
    plot_heatmaps(matrices, output_dir)
    
    # 5. 生成Prism数据
    generate_prism_data(avg_table, output_dir)
    
    # 6. 保存结果
    save_results(avg_table, matrices, output_dir)
    
    cli_h1("分析完成")
    cli_alert_success(glue("所有结果已保存至: {output_dir}"))
    
  }, error = function(e) {
    cli_alert_danger(glue("流程异常终止: {e$message}"))
    quit(status = 1)
  })
}

if (!interactive()) {
  main()
}


#主要更新：
#
#    行标准化处理：
#        在generate_normalized_matrices函数中，对每个能量矩阵进行行标准化（Z-score标准化）
#        使用scale()函数进行行标准化：normalized_matrix <- t(scale(t(matrix_df)))
#        同时保留原始矩阵和标准化矩阵
#
#    热图绘制：
#        使用标准化矩阵绘制热图
#        热图标题明确标注为"Z-score"
#        移除了scale = "row"参数（因为我们已经手动标准化）
#
#    结果保存：
#        在Excel文件中同时保存原始矩阵和标准化矩阵
#        工作表命名区分：Matrix_{energy}_raw和Matrix_{energy}_zscore
#
#    日志更新：
#        更新日志信息以反映标准化处理
#        明确标注"生成标准化能量矩阵"和"绘制热图（行标准化）"
#
#标准化处理说明：
#
#行标准化（Z-score标准化）公式：
#
#复制代码
#z = (x - μ) / σ
#
#其中：
#
#    x 是原始值
#    μ 是行的平均值
#    σ 是行的标准差
#
#这种标准化方法将每行的数据转换为均值为0、标准差为1的分布，使得不同抗原之间的比较更加公平。
#输出文件结构：
#
#复制代码
#results/
#├── average_energy.csv                # 平均能量值总表（原始值）
#├── energy_results.xlsx               # Excel汇总文件（包含所有数据）
#│   ├── Average_Energy                # 平均能量值表
#│   ├── Matrix_desolv_raw             # desolv原始矩阵
#│   ├── Matrix_desolv_zscore          # desolv标准化矩阵
#│   ├── Matrix_elec_raw               # elec原始矩阵
#│   ├── Matrix_elec_zscore            # elec标准化矩阵
#│   ├── ...                           # 其他能量项
#│   └── Prism_*                       # Prism数据表
#├── heatmap_desolv.png                # desolv能量热图（标准化）
#├── heatmap_elec.png                  # elec能量热图（标准化）
#├── heatmap_total.png                 # total能量热图（标准化）
#├── heatmap_vdw.png                   # vdw能量热图（标准化）
#├── heatmap_air.png                   # air能量热图（标准化）
#└── prism_data/                       # Prism数据目录
#    ├── prism_3PGW_H.csv              # 3PGW_H抗原的Prism数据
#    ├── prism_3JCR_o.csv              # 3JCR_o抗原的Prism数据
#    └── ...                           # 其他抗原
#
#使用示例：
#
#bash
#
#复制代码
#Rscript analyze_docking_results.R \
#  -i docking_results_with_metadata.csv \
#  -o results
#
#日志输出示例：
#
#复制代码
#对接结果数据分析流程启动
#🔹 [2023-10-02 14:30:00] 开始时间: 2023-10-02 14:30:00
#✅ [2023-10-02 14:30:01] 读取数据: docking_results_with_metadata.csv (1000 行)
#🔹 [2023-10-02 14:30:02] 计算平均能量值...
#✅ [2023-10-02 14:30:03] 创建平均总表: 50 抗原 x 5 抗体
#🔹 [2023-10-02 14:30:04] 生成标准化能量矩阵...
#✅ [2023-10-02 14:30:05] 生成5个标准化能量矩阵
#🔹 [2023-10-02 14:30:06] 绘制热图（行标准化）...
#✅ [2023-10-02 14:30:07] 保存5个热图
#🔹 [2023-10-02 14:30:08] 准备Prism数据...
#✅ [2023-10-02 14:30:09] 保存50个Prism数据文件
#🔹 [2023-10-02 14:30:10] 保存汇总Excel文件
#✅ [2023-10-02 14:30:11] 成功保存所有结果到: results/
#分析完成
#✅ 所有结果已保存至: results/
#
#这个更新版本确保所有值在热图绘制前都进行了行标准化处理，同时保留了原始数据用于其他分析。Excel文件中同时包含原始值和标准化值，方便您进行对比分析。
