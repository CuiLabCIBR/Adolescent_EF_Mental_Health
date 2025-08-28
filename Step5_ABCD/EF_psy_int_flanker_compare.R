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
library(showtext)
library(patchwork)

# 文件路径设置
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/EF_results'
  demopath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/rawdata_results0616'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/Rcode_EFnorms/functions"
  
} else {
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/interfileFolder/ABCD'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/FigureFolder/figure4'
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/code/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/results/figure4"
}


# 加载自定义函数
source(paste0(functionFolder,"/gamm_varyingcoefficients_new.R"))

# 读取数据
Flanker_data <- read_csv(paste0(datapath, '/Flanker.deviations_addr.csv'))
Flanker_data$Sex <- as.factor(Flanker_data$Sex)
Flanker_data$Sex <- factor(Flanker_data$Sex, levels = c(1, 2), labels = c("M", "F"))

# 定义变量
psyc_variables_continous <- c("cbcl_scr_syn_social_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")
EFvars.set <- matrix(c("nihtbx_flanker_uncorrected_deviationZ", "Flanker"), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("varname", "dataname")))
EFvars.set <- as.data.frame(EFvars.set)

# 数据标准化和异常值处理函数
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

# 对 Flanker 数据进行标准化和异常值处理
Flanker_data <- standardize_clean(Flanker_data, psyc_variables_continous)
psyc_variables_continous <- paste0(psyc_variables_continous, "_z")

# beta table
beta_table <- read_csv("/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/results/figure4/corr_EF_psych_continuous.result_ABCD_withbonf.csv")

# 设置模型参数
knots <- 3
set_fx <- TRUE
increments <- 1000
draws <- 1000
return_posterior_coefficients <- TRUE

# 初始化结果存储列表
interaction_results <- list()
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname <- paste0(EFvars.set$dataname[i], "_data")
  
  for (j in 1:length(psyc_variables_continous)) {
    dependentvar <- psyc_variables_continous[j]
    smooth_var <- "Age_year"
    covariates <- "Sex"
    bestmodel_row <- beta_table %>% 
      filter(parcel == dependentvar & interest.indep.var == int_var)
    bestmodel <- as.character(bestmodel_row$bettermodel[1])
    
    cat(paste("Processing:", dataname, "~", int_var, "by", smooth_var, "for dependent var:", dependentvar, "\n"))
    
    result <- gamm.varyingcoefficients(
      dependentvar = dependentvar,
      dataname = dataname,
      smooth_var = smooth_var,
      int_var = int_var,
      covariates = covariates,
      bestmodel = bestmodel, 
      knots = knots,
      set_fx = set_fx,
      increments = increments,
      draws = draws,
      return_posterior_coefficients = return_posterior_coefficients
    )
    
    # save
    model_name <- paste(EFvars.set$dataname[i], psyc_variables_continous[j], sep = "_")
    interaction_results[[model_name]] <- result
  }
}

# save data as RDS
saveRDS(interaction_results, file = paste0(resultFolder, "/all_interaction_results_ABCD.rds"))

#interaction_results <- readRDS(paste0(resultFolder, "/all_interaction_results_ABCD.rds"))
# plot
ariable_mapping <- c(
  "cbcl_scr_syn_social_r_z" = "Social",
  "cbcl_scr_syn_attention_r_z" = "Attention",
  "cbcl_scr_syn_internal_r_z" = "Internalizing",
  "cbcl_scr_syn_external_r_z" = "Externalizing"
)

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

all_slope_summaries <- list()
for (model_name in names(interaction_results)) {
  result <- interaction_results[[model_name]]
  parts <- strsplit(model_name, "_")[[1]]
  current_task <- parts[1]
  current_variable_raw <- paste(parts[2:length(parts)], collapse = "_")
  current_slope_data <- result[[2]]

  slope_value <- beta_table %>%
    filter(parcel == current_variable_raw) %>%
    pull(slope)
  
  if (length(slope_value) == 0) {
    message(paste("找不到 `slope` 列的值，尝试使用 `correstimate`."))
    slope_value <- beta_table %>%
      filter(parcel == current_variable_raw) %>%
      pull(correstimate)
  }
  
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

# 确定统一的全局颜色范围
global_abs_diff_range <- range(do.call(rbind, all_slope_summaries)$abs_diff, na.rm = TRUE)

for (model_name in names(interaction_results)) {
  result <- interaction_results[[model_name]]
  
  parts <- strsplit(model_name, "_")[[1]]
  current_task <- parts[1]
  current_variable_raw <- paste(parts[2:length(parts)], collapse = "_")
  
  slope_value <- beta_table %>%
    filter(parcel == current_variable_raw) %>%
    pull(slope)
  
  if (length(slope_value) == 0) {
    message(paste("cant find `slope` , try `correstimate`."))
    slope_value <- beta_table %>%
      filter(parcel == current_variable_raw) %>%
      pull(correstimate)
  }
  
  if (length(slope_value) == 0) {
    message(paste("can not find ", current_task, "and ", current_variable_raw, "beta value."))
    next
  }
  
  current_slope_data <- result[[2]]
  
  # 计算每个年龄点的斜率中位数、95% CI，并进行显著性检验
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
  
  # 查找并记录转折点
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
  
  # min_y_local <- -0.09
  # max_y_local <- 0.06
  # y_interval <-  0.05
  # step_expand <- y_interval / 2
  # min_y <- floor(min_y_local / y_interval) * y_interval
  # max_y <- ceiling(max_y_local / y_interval) * y_interval
  # repeat {
  #   y_breaks <- seq(min_y, max_y, by = y_interval)
  #   if (length(y_breaks) >= 3) break
  #   min_y <- min_y - 0.01
  #   max_y <- max_y + 0.01
  # }
  min_y <- - 0.10
  max_y <-  0.06
  # 
  # 绘图
  slope_plot_main <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) +
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") +
    geom_hline(yintercept = slope_value, linetype = "dashed", color = "red", size = 0.4) +
    labs(
      title = paste0(ariable_mapping[current_variable_raw], " ~ ", current_task),
      x = NULL, y = "Slope"
    ) +
    scale_x_continuous(name = "", limits = c(8.8, 16), breaks = seq(9, 16, by = 2)) +
    scale_y_continuous(
      limits = c(min_y, max_y),
      breaks = y_breaks,
      labels = label_number(accuracy = 0.01, style_positive = "plus", style_negative = "minus")
    ) +
    custom_theme
  
  # 绘制显著性条形图（子图）
  slope_plot_bar <- ggplot(slope_summary, aes(x = Age_year, y = 1)) +
    geom_col(aes(fill = abs_diff), data = filter(slope_summary, is_significant == TRUE),
             width = 1, color = "transparent") +
    geom_col(data = filter(slope_summary, is_significant == FALSE), fill = "white",
             width = 1, color = "transparent") +
    scale_fill_gradient(low = "white", high = "#A4C5DF", na.value = "white", guide = "none",
                        limits = global_abs_diff_range) +
    geom_rect(aes(xmin = 8.8, xmax = 16, ymin = -0.0, ymax = 1),
              fill = NA, color = "black", size = 0.25) +
    scale_x_continuous(name = " ", limits = c(8.8, 16), breaks = seq(9, 16, by = 2)) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(t = -15, r = 0, b = 0, l = 0, unit = "pt")
    )
  
  # 合并主图和子图并保存
  final_plot <- plot_grid(slope_plot_main, slope_plot_bar, ncol = 1, rel_heights = c(0.95, 0.05), align = "v")
  final_plot
  ggsave(
    filename = paste0(FigureFolder, "/slope_comparison_", current_task, "_", current_variable_raw, ".pdf"),
    plot = final_plot,
    width = 6, height = 5.5, units = "cm"
  )
}
write.csv(transition_age_summary, file = paste0(resultFolder, "/transition_age_summary_ABCD.csv"), row.names = FALSE)
print(transition_age_summary)


