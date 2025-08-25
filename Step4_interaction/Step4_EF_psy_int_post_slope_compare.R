rm(list = ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(cowplot)
# 文件路径设置
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/interfileFolder"
  resultFolder <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/results/Int'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/figure/Int'
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/function"
} else {
  datapath <- ' '
  FigureFolder <- ' '
  interfileFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
}

beta_table <- read.csv("/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/results/corr/corr_results_with_bonf.csv")
# 读取数据
GNGd_data <- readRDS(paste0(datapath, '/Gonogo/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(datapath, '/1-back/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(datapath, '/2-back/back2Acc.deviations.rds'))
GNGd_data$Sex <- as.factor(GNGd_data$Sex)
back1_data$Sex <- as.factor(back1_data$Sex)
back2_data$Sex <- as.factor(back2_data$Sex)

source(paste0(functionFolder, '/gam_varyingcoefficients.R'))
# 定义变量
psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
EFvars.set <- matrix(c("d_prime_deviationZ", "GNGd",
                       "Oneback_acc_deviationZ", "back1",
                       "Twoback_acc_deviationZ", "back2"), byrow=TRUE, ncol=2, dimnames=list(NULL, c("varname", "dataname")))
EFvars.set <- as.data.frame(EFvars.set)
#####
standardize_clean <- function(df, vars) {
  for (var in vars) {
    x <- df[[var]]
    x <- x[!is.na(x)]
    mu <- mean(x)
    sd_val <- sd(x)
    valid_index <- which(abs(df[[var]] - mu) <= 3 * sd_val)
    z_varname <- paste0(var, "_z")
    df[[z_varname]] <- NA
    df[[z_varname]][valid_index] <- scale(df[[var]][valid_index])
  }
  return(df)
}

original_vars <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")

GNGd_data  <- standardize_clean(GNGd_data,  original_vars)
back1_data <- standardize_clean(back1_data, original_vars)
back2_data <- standardize_clean(back2_data, original_vars)

psyc_variables_continous <- paste0(original_vars, "_z")

# 设置参数
knots <- 3
set_fx <- FALSE 
increments <- 1000
draws <- 1000
return_posterior_coefficients <- T

# 初始化结果存储列表
interaction_results <- list()

# 主循环 - 循环计算所有模型并保存结果
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname <- paste0(EFvars.set$dataname[i], "_data")
  
  for (j in 1:length(psyc_variables_continous)) {
    dependentvar <- psyc_variables_continous[j]
    smooth_var <- "Age_year"
    
    cat(paste("Processing:", dataname, "~", int_var, "by", smooth_var, "for dependent var:", dependentvar, "\n"))
    
    # 调用函数，返回后验斜率数据
    result <- gam.varyingcoefficients(
      dependentvar = dependentvar,
      dataname = dataname,
      smooth_var = smooth_var,
      int_var = int_var,
      covariates = "Sex",
      knots = knots,
      set_fx = set_fx,
      increments = increments,
      draws = draws,
      return_posterior_coefficients = return_posterior_coefficients
    )
    
    # 保存结果
    model_name <- paste(EFvars.set$dataname[i], psyc_variables_continous[j], sep = "_")
    interaction_results[[model_name]] <- result
  }
}

# 保存所有结果为RDS文件
saveRDS(interaction_results, file = paste0(resultFolder, "/all_interaction_results.rds"))

#interaction_results <- readRDS(paste0(resultFolder, "/all_interaction_results.rds"))


task_mapping <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
variable_mapping <- c(
  "SDQ_ES_sum_z" = "Emotional Symptoms", "SDQ_PP_sum_z" = "Peer Problems",
  "SDQ_CP_sum_z" = "Conduct Problems", "SDQ_H_sum_z" = "Hyperactivity",
  "SDQ_PB_sum_z" = "Prosocial Behavior"
)

# 定义统一的主题
custom_theme <- theme_minimal() +
  theme(
    axis.line = element_line(size = 0.25, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 9, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.05, "cm")
  )

# 初始化一个空的数据框来存储所有模型的转折点结果
transition_age_summary <- data.frame(
  dataname = character(),
  parcel = character(),
  age_of_transition = numeric(),
  transition_type = character(),
  stringsAsFactors = FALSE
)

# -
y_axis_settings <- list(
  "GNGd_SDQ_ES_sum_z" = list(y_limits = c(-0.08, 0.10), y_breaks = seq(-0.05, 0.10, by = 0.05)),
  "back1_SDQ_ES_sum_z" = list(y_limits = c(-0.08, 0.10), y_breaks = seq(-0.05, 0.10, by = 0.05)),
  "back2_SDQ_ES_sum_z" = list(y_limits = c(-0.08, 0.10), y_breaks = seq(-0.05, 0.10, by = 0.05)),
  
  "GNGd_SDQ_PP_sum_z" = list(y_limits = c(-0.21, -0.0), y_breaks = seq(-0.2, 0, by = 0.05)),
  "back1_SDQ_PP_sum_z" = list(y_limits = c(-0.21, -0.0), y_breaks = seq(-0.2, 0, by = 0.05)),
  "back2_SDQ_PP_sum_z" = list(y_limits = c(-0.21, -0.0), y_breaks = seq(-0.2, 0, by = 0.05)),
  
  "GNGd_SDQ_CP_sum_z" = list(y_limits = c(-0.12, 0.05), y_breaks = seq(-0.10, 0.05, by = 0.05)),
  "back1_SDQ_CP_sum_z" = list(y_limits = c(-0.12, 0.05), y_breaks = seq(-0.10, 0.05, by = 0.05)),
  "back2_SDQ_CP_sum_z" = list(y_limits = c(-0.12, 0.05), y_breaks = seq(-0.10, 0.05, by = 0.05)),
  
  "GNGd_SDQ_H_sum_z" = list(y_limits = c(-0.14, 0.07), y_breaks = seq(-0.1, 0.05, by = 0.05)),
  "back1_SDQ_H_sum_z" = list(y_limits = c(-0.2, 0.05), y_breaks = seq(-0.2, 0, by = 0.05)),
  "back2_SDQ_H_sum_z" = list(y_limits = c(-0.13, 0.07), y_breaks = seq(-0.1, 0.05, by = 0.05)),
  
  "GNGd_SDQ_PB_sum_z" = list(y_limits = c(0, 0.17), y_breaks = seq(0, 0.2, by = 0.05)),
  "back1_SDQ_PB_sum_z" = list(y_limits = c(0, 0.17), y_breaks = seq(0, 0.2, by = 0.05)),
  "back2_SDQ_PB_sum_z" = list(y_limits = c(0, 0.17), y_breaks = seq(0, 0.2, by = 0.05))
)

all_slope_summaries <- list()
for (model_name in names(interaction_results)) {
  result <- interaction_results[[model_name]]
  parts <- strsplit(model_name, "_")[[1]]
  current_task <- parts[1]
  current_variable_raw <- paste(parts[2:length(parts)], collapse = "_")
  current_slope_data <- result[[2]]
  slope_value <- beta_table %>%
    filter(dataname == current_task, parcel == current_variable_raw) %>%
    pull(beta)
  if (length(slope_value) == 0) next
  
  slope_summary_temp <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      p_value = min(mean(slope > slope_value, na.rm = TRUE), mean(slope < slope_value, na.rm = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(
      is_significant = (p_value < 0.025),
      abs_diff = abs(median_slope - slope_value)
    )
  all_slope_summaries[[model_name]] <- slope_summary_temp
}
global_abs_diff_range <- range(do.call(rbind, all_slope_summaries)$abs_diff, na.rm = TRUE)

#######
for (model_name in names(interaction_results)) {
  result <- interaction_results[[model_name]]
  parts <- strsplit(model_name, "_")[[1]]
  current_task <- parts[1]
  current_variable_raw <- paste(parts[2:length(parts)], collapse = "_")

  y_settings <- y_axis_settings[[model_name]]

  current_slope_data <- result[[2]]
  slope_value <- beta_table %>%
    filter(dataname == current_task, parcel == current_variable_raw) %>%
    pull(beta)
  if (length(slope_value) == 0) next
  
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE),
      p_value = min(mean(slope > slope_value, na.rm = TRUE), mean(slope < slope_value, na.rm = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(
      is_significant = (p_value < 0.025),
      abs_diff = abs(median_slope - slope_value)
    )
  
  # 记录转折点
  if (nrow(slope_summary) > 1) {
    for (i in 2:nrow(slope_summary)) {
      prev_sig <- slope_summary$is_significant[i-1]
      curr_sig <- slope_summary$is_significant[i]
      if (prev_sig != curr_sig) {
        transition_age <- slope_summary$Age_year[i]
        type <- if (prev_sig == FALSE && curr_sig == TRUE) "nonsig_to_sig" else "sig_to_nonsig"
        transition_age_summary <- rbind(transition_age_summary,
                                        data.frame(
                                          dataname = current_task, parcel = current_variable_raw,
                                          age_of_transition = transition_age, transition_type = type,
                                          stringsAsFactors = FALSE
                                        )
        )
      }
    }
  }
  
  # 绘制主图
  slope_plot_main <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) +
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") +
    geom_hline(yintercept = slope_value, linetype = "dashed", color = "red", size = 0.4) +
    labs(
      title = paste0(variable_mapping[current_variable_raw], " ~ ", task_mapping[current_task]),
      x = "Age", y = "Slope"
    ) +
    scale_x_continuous(name = "", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
    scale_y_continuous(
      limits = y_settings$y_limits,
      breaks = y_settings$y_breaks,
      labels = label_number(accuracy = 0.01, style_positive = "plus", style_negative = "minus")
    ) +
    custom_theme
  
  # 绘制显著性条形图
  slope_plot_bar <- ggplot(slope_summary, aes(x = Age_year, y = 1)) +
    geom_col(aes(fill = abs_diff), data = filter(slope_summary, is_significant == TRUE),
             width = 1, color = "transparent") +
    geom_col(data = filter(slope_summary, is_significant == FALSE), fill = "white",
             width = 1, color = "transparent") +
    scale_fill_gradient(low = "#B7CDE0", high = "#487BAC", na.value = "white", guide = "none",
                        limits = global_abs_diff_range) +
    geom_rect(aes(xmin = 11, xmax = 18, ymin = -0.0, ymax = 1),
              fill = NA, color = "black", size = 0.25) +
    scale_x_continuous(name = " ", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(t = -15, r = 0, b = 0, l = 0, unit = "pt")
    )
  
  # 将主图和条形图组合
  final_plot <- plot_grid(slope_plot_main, slope_plot_bar, ncol = 1, rel_heights = c(0.95, 0.05), align = "v")
  
  # 保存图片
  ggsave(
    filename = paste0(FigureFolder, "/slope_comparison_", current_task, "_", current_variable_raw, ".pdf"),
    plot = final_plot,
    width = 6, height = 5.5, units = "cm"
  )
}


# 保存转折点总结
write.csv(transition_age_summary, file = paste0(resultFolder, "/transition_age_summary.csv"), row.names = FALSE)
print(transition_age_summary)

