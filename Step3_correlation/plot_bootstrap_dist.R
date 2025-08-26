# 1. 加载必要的包
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)

# 2. 设置路径
# 用户指定的包含15个RDS文件的文件夹路径
folder_path <- "D:/datasets/yunfu/interfile_folder/anovaPB_new"
# 设置一个输出文件夹来保存最终的PDF文件
figure_folder <- 'D:/datasets/yunfu/figures/fig2/anovaPB_combined'

# 确保输出文件夹存在
dir.create(figure_folder, showWarnings = FALSE, recursive = TRUE)

# 3. 获取所有目标RDS文件的列表
all_rds_files <- list.files(path = folder_path, pattern = "^anova_simulation_.*\\.rds$", full.names = TRUE)

if (length(all_rds_files) == 0) {
  stop("错误: 在指定文件夹中找不到任何 'anova_simulation' RDS 文件。")
}

# 4. 初始化一个列表来存储所有的图
plot_list <- list()

# 5. 循环处理每个RDS文件并生成一个图
for (file_path in all_rds_files) {
  
  # 从文件名中提取任务变量和心理变量
  file_name <- basename(file_path)
  # 使用正则表达式捕获变量部分
  matches <- str_match(file_name, "anova_simulation_([^_]+)_(dataSDQ_.+)\\.rds$")
  
  if (is.na(matches[1,1])) {
    cat(paste("警告: 文件名", file_name, "不符合预期的格式，跳过。\n"))
    next
  }
  
  task_var <- matches[1, 2] # 提取任务变量, 如: back1, back2, GNGd
  psy_var <- matches[1, 3]  # 提取心理变量, 如: dataSDQ_CP_sum_z
  
  cat(paste("正在处理:", task_var, "-", psy_var, "\n"))
  
  # 读取RDS文件
  sim_data <- readRDS(file_path)
  
  # 提取模拟数据和观测值
  combined_sim_stats <- sim_data$simulation$simulated_stats
  observed_stat <- sim_data$simulation$observed_stat
  
  # 检查数据是否存在
  if (is.null(combined_sim_stats) || is.null(observed_stat)) {
    cat(paste("警告: 在文件", file_name, "中找不到所需数据，跳过。\n"))
    next
  }
  
  # --- 开始绘图 ---
  
  # 1. 准备数据和标签
  sim_df <- data.frame(simulated_stats = combined_sim_stats)
  n_sim <- length(combined_sim_stats)
  final_p_value <- mean(combined_sim_stats >= observed_stat, na.rm = TRUE)
  
  if (final_p_value == 0) {
    p_label <- paste0("p < ", formatC(1/n_sim, format = "e", digits = 0))
  } else {
    p_label <- paste0("p = ", format(final_p_value, digits = 3, nsmall = 3))
  }
  annotation_label <- paste0("Obs. = ", format(round(observed_stat, 2), nsmall = 2), "\n", p_label)
  
  # 2. 统一将注释文本位置设置在右上角
  x_pos <- Inf      
  hjust_val <- 1.5   
  
  # 3. 使用ggplot2绘图
  p <- ggplot(sim_df, aes(x = simulated_stats)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = observed_stat, color = "red", linetype = "dashed", linewidth = 1) +
    geom_point(aes(x = !!observed_stat, y = 0), color = "red", size = 3, shape = 18) +
    labs(
      title = paste(task_var, "-", psy_var),
      x = "Statistic Value",
      y = "Frequency"
    ) +
    annotate(
      "text",
      x = x_pos, y = Inf,
      label = annotation_label,
      hjust = hjust_val, vjust = 1.5,
      color = "black", fontface = "bold", size = 3.5
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
  
  # ================================================================= #
  # --- 【修改点 1】: 创建一个新的排序键 (心理变量_任务变量) --- #
  sort_key <- paste(psy_var, task_var, sep = "_")
  # 将生成的图用新的排序键存入列表中
  plot_list[[sort_key]] <- p
  # ================================================================= #
}

# 6. 检查是否成功生成了任何图
if (length(plot_list) > 0) {
  
  cat(paste("\n成功生成", length(plot_list), "个图，现在按任务列进行排序和拼接...\n"))
  
  # ================================================================= #
  # --- 【修改点 2】: 使用新的排序键对列表进行排序 --- #
  # 1. 获取所有新创建的键并进行字母排序
  sorted_keys <- sort(names(plot_list))
  # 2. 根据排序好的键，重新组织图的列表
  sorted_plot_list <- plot_list[sorted_keys]
  # ================================================================= #
  
  # 7. 使用 patchwork 将所有图拼接成一个 5x3 的网格
  # 由于列表已经按 psy_var -> task_var 排序, wrap_plots 会自动按列填充任务
  combined_plot <- wrap_plots(sorted_plot_list, ncol = 3, nrow = 5)
  
  # 8. 保存拼接后的图为PDF文件
  final_pdf_path <- file.path(figure_folder, "combined_anova_simulation_distributions_by_task.pdf")
  
  # 调整宽度和高度以获得最佳布局 (单位: 英寸)
  ggsave(final_pdf_path, plot = combined_plot, width = 12, height = 15)
  
  cat(paste("所有分布图已成功拼接并保存至:", final_pdf_path, "\n"))
  
} else {
  cat("\n没有生成任何图，请检查文件路径和文件内容。\n")
}

cat("\n脚本执行完毕。\n")