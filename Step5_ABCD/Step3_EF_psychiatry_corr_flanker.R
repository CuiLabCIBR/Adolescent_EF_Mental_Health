library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(writexl)
library(dplyr)

# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/Rcode/functions_corr"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/interfileFolder/corr"
}else{
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/interfileFolder/ABCD'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/FigureFolder/figure4'
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/code/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy_2508/results/figure4"
}


# source functions
source(paste0(functionFolder, "/gamm_factor_interaction_deviation.R"))

Flanker_data <- read_csv(paste0(datapath, '/Flanker.deviations_addr.csv'))
Flanker_data <- as.data.frame(Flanker_data)

## 1) set up variables
psyc_variables_continous <- c("cbcl_scr_syn_internal_r","cbcl_scr_syn_social_r",
                              "cbcl_scr_syn_external_r","cbcl_scr_syn_attention_r")
# EF vars
EFvar <- "nihtbx_flanker_uncorrected_deviationZ"
## 2) convert variables class & describe variables
Flanker_data[, c(psyc_variables_continous, EFvar)] <- lapply(Flanker_data[, c(psyc_variables_continous, EFvar)], as.numeric)
site_id_ltab <- unique(Flanker_data$site_id_l)
Flanker_data$site_id_l_fac <- factor(Flanker_data$site_id_l, levels=site_id_ltab, labels=paste0("site_id_l", 1:length(site_id_ltab)))
###normalize
# Flanker_data[, paste0(psyc_variables_continous, "_z")] <- scale(Flanker_data[, psyc_variables_continous])

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

original_vars <- c("cbcl_scr_syn_internal_r","cbcl_scr_syn_social_r",
                   "cbcl_scr_syn_external_r","cbcl_scr_syn_attention_r")

Flanker_data  <- standardize_clean(Flanker_data,  original_vars)

psyc_variables_continous <- paste0(original_vars, "_z")

###
describe_tab_Flanker <- CreateTableOne(c(EFvar, psyc_variables_continous, "Age_year", "Sex"), 
                                       data = Flanker_data, testNonNormal = TRUE)
describe_tab_Flanker.continous <- as.data.frame(describe_tab_Flanker$ContTable[["Overall"]])
# save out
write.xlsx(list(Flanker_con=describe_tab_Flanker.continous), paste0(resultFolder, "/description_interest_vars.xlsx"), rowNames=T)

## 3) Correlations in separate age periods
# continuos variables
corr.result <- list()
for (psyvar.tmp in psyc_variables_continous) {
  dependentvar <- psyvar.tmp
  interest.indep.var <- EFvar
  covariates <- "Sex"
  smoothvar <- "Age_year"
  
  # Perform analysis using gamm.smooth.predict.interaction
  result.tmp <- gamm.smooth.predict.interaction(
    dependentvar = dependentvar,
    dataname = "Flanker_data",
    smoothvar = smoothvar,
    interest.indep.var = interest.indep.var,
    covariates = covariates
  )
  
  result.tmp <- as.data.frame(result.tmp)
  corr.result[[psyvar.tmp]] <- result.tmp
}

# Combine results
corr.result.df <- do.call(rbind, corr.result)
write.csv(corr.result.df, paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv"), row.names = FALSE)

#corr.result.df <- read_csv(paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv"))
corr.result.df$correstimate <- as.numeric(corr.result.df$correstimate)
corr.result.df$anovap.bonf <- p.adjust(corr.result.df$boots.pvalues, method = "bonferroni")
corr.result.df$sig <- (corr.result.df$anovap.bonf < 0.05)
corr.result.df$significance <- ifelse(corr.result.df$sig, "*", "")

write.csv(corr.result.df, paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD_withbonf.csv"), row.names = FALSE)



#corr.result.df <- read.csv(paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD_withbonf.csv"))
# Define labels and limits
psy_labels <- c("cbcl_scr_syn_internal_r_z" = "Internalizing", 
                "cbcl_scr_syn_social_r_z" = "Social",
                "cbcl_scr_syn_external_r_z" = "Externalizing",
                "cbcl_scr_syn_attention_r_z" = "Attention")


corr.result.df$parcel <- factor(corr.result.df$parcel,
                                levels = psyc_variables_continous,
                                labels = psy_labels)


y_limits <- c(-0.06, 0.02)
corr.result.df$label_y <- ifelse(
  corr.result.df$slope >= 0, 
  corr.result.df$slope + 0.01, 
  corr.result.df$slope - 0.01
)
task_colors <- c("Flanker" = "#c6d6ea")
vline_positions <- seq(1.5, length(unique(corr.result.df$parcel)) - 0.5, by = 1)

Fig <- ggplot(data = corr.result.df, aes(x = parcel, y = slope, color = "#c6d6ea", fill = "#c6d6ea")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5, size = 1, 
           color = NA, show.legend = T) +  # 去除图例，使用 Task 填充和边框
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.25) +  # 添加 y = 0 的水平线 浅灰色水平辅助线
  geom_text(aes(label = significance, y = label_y), 
            position = position_dodge(width = 0.7), size = 5, color = "black") + 
  scale_fill_manual(values = "#c6d6ea") +  # 设置任务的填充颜色
  scale_color_manual(values = "#c6d6ea") +  # 设置任务的边框颜色
  scale_y_continuous(limits = y_limits, breaks = seq(-0.1, 0.05, by = 0.02), labels = scales::number_format()) +  # 设置 y 轴的上下限和刻度
  labs(x = NULL,
       y = "beta",
       color = "Tasks") +  # 标题和 y 轴标签
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 0.25),  # Add x-axis line
        axis.line.y = element_line(color = "black", size = 0.25),  # Add y-axis line
        axis.title = element_text(size = 9),
        #axis.text.x = element_text(size = 9, hjust = 1,vjust = 1,color = "black", angle = 45),  # x 轴标签
        axis.text.x = element_text(size = 9, hjust = 0.5, color = "black"),  # x 轴标签
        axis.text.y = element_text(size = 9,color = "black"),  
        axis.ticks.x = element_line(color = "black", size = 0.25),
        axis.ticks.y = element_line(color = "black", size = 0.25),
        axis.ticks.length = unit(0.05, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.position = "none",
        legend.key.size = unit(0.05, "cm"),
        panel.grid.major.y = element_blank(),  # Remove x-axis grid lines
        panel.grid.major.x = element_blank(),  # Remove x-axis grid lines
        panel.grid.minor = element_blank()) 

print(Fig)
ggsave(paste0(FigureFolder, "/correlation_barplot_slope.pdf"), plot = Fig, width = 12, height = 5, units = "cm")


# #combined_results <- read.csv("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure4_corr_delet_Z/corr_EF_psych_continuous.result_ABCD_withbonf.csv")
# ## 数据准备：将数据转换为适合绘制柱状图的格式
# y_levels <- c("cbcl_scr_syn_internal_r_z","cbcl_scr_syn_social_r_z","cbcl_scr_syn_external_r_z",
#               "cbcl_scr_syn_attention_r_z")
# combined_results$correstimate <- as.numeric(combined_results$correstimate)
# # 自定义心理变量的标签
# psy_labels <- c("cbcl_scr_syn_internal_r_z"= "Internalizing Symptoms", "cbcl_scr_syn_social_r_z"= "Social","cbcl_scr_syn_external_r_z" = "Externalizing Symptoms",
#                 "cbcl_scr_syn_attention_r_z" = "Attention")
# 
# # 应用新的标签
# 
# 
# 
# combined_results$parcel <- factor(combined_results$parcel,
#                                   levels = y_levels,
#                                   labels = psy_labels)
# combined_results$Task <- factor(combined_results$dataname, levels = c("Flanker"),
#                                 labels = c("Flanker"))
# 
# task_colors <- c("Flanker" = "#c6d6ea")
# # 确定 y 轴的上下限，以包含所有数据
# max_correlation <- max(abs(combined_results$correstimate), na.rm = TRUE)
# y_limits <- c(-0.125, 0.05)  # 对称上下限
# 
# # 创建绘图
# combined_results$significance <- ifelse(combined_results$sig, "*", "")  # 标记星号
# combined_results$label_y <- ifelse(
#   combined_results$correstimate >= 0,
#   combined_results$correstimate + 0.01,  # 正值柱子顶部上方标星
#   combined_results$correstimate - 0.02  # 负值柱子底部下方标星
# )
# 
# vline_positions <- seq(1.5, length(unique(combined_results$parcel)) - 0.5, by = 1)
# 
# Fig <- ggplot(data = combined_results, aes(x = parcel, y = correstimate, color = Task, fill = Task)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5, size = 1,
#            color = NA, show.legend = T) +  # 去除图例，使用 Task 填充和边框
#   geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.4) +  # 添加 y = 0 的水平线 浅灰色水平辅助线
#   geom_text(aes(label = significance, y = label_y),
#             position = position_dodge(width = 0.2), size = 5, color = "black") +
#   scale_fill_manual(values = task_colors) +  
#   scale_color_manual(values = task_colors) + 
#   scale_y_continuous(limits = y_limits, breaks = seq(-0.1, 0.05, by = 0.05), labels = scales::number_format()) +  # 设置 y 轴的上下限和刻度
#   labs(title = "Correlation between Flanker and Mental Health",
#        x = "",
#        y = "correlation coefficient",
#        color = "Tasks") +  # 标题和 y 轴标签
#   theme_minimal() +
#   theme(axis.line.y = element_blank(),  # 去掉 y 轴线
#         axis.title = element_text(size = 9),
#         axis.text.x = element_text(size = 9, hjust = 0.5,color = "black"),  # x 轴标签
#         axis.text.y = element_text(size = 9,color = "black"),
#         plot.title = element_text(size = 9, hjust = 0.5),
#         legend.title = element_text(size = 9),
#         legend.text = element_text(size = 9),
#         legend.position = "bottom",
#         legend.key.size = unit(0.4, "cm"),
#         panel.grid.major.y = element_line(color = "gray90", linetype = "solid", size = 0.2),  # 浅灰色水平辅助线
#         panel.grid.major.x = element_blank(),  # 移除 x 轴的网格线
#         panel.grid.minor = element_blank()) + # 移除次网格线
#   annotate("segment", x = vline_positions, xend = vline_positions, y = -0.005, yend = 0, color = "black", size = 0.4)
# 
# print(Fig)
# ggsave(paste0(FigureFolder, "/correlation_barplot.pdf"), plot = Fig, width = 10, height = 7, units = "cm")








####new plot
psy_labels <- c("cbcl_scr_syn_attention_r_z" = "Attention", 
                "cbcl_scr_syn_external_r_z" = "Externalizing",
                "cbcl_scr_syn_social_r_z" = "Social",
                "cbcl_scr_syn_internal_r_z" = "Internalizing")


psyc_variables_continous <- names(psy_labels)

corr.result.df <- corr.result.df %>%
  arrange(slope) %>%
  mutate(parcel = factor(parcel, levels = parcel))  # 重新设定 factor 顺序

# 设置 bar 的颜色（或使用原有颜色）
corr.result.df$fill_color <- "#c6d6ea"  # 或自定义每个 bar 的颜色

# 显著性星号位置
corr.result.df$label_x <- ifelse(
  corr.result.df$slope >= 0, 
  corr.result.df$slope + 0.004, 
  corr.result.df$slope - 0.004
)

# 图中标注的心理维度名称（贴在 bar 上）
corr.result.df$label_text <- as.character(corr.result.df$parcel)

# 绘图
Fig <- ggplot(data = corr.result.df, aes(y = parcel, x = slope)) +
  geom_col(aes(fill = fill_color), width = 0.5, color = "white", show.legend = FALSE) +
  geom_text(aes(label = significance, x = label_x), size = 5, color = "black", hjust = ifelse(corr.result.df$slope >= 0, 0, 1)) +
  geom_text(aes(label = label_text), x = 0.001, hjust = 0, size = 3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.25) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.05, 0.02), breaks = seq(-0.04, 0.02, 0.02)) +
  labs(x = "beta", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.ticks.x = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.05, "cm"),
    axis.title.x = element_text(size = 9)
  )

print(Fig)
# 保存图片
ggsave(paste0(FigureFolder, "/correlation_barplot_slope.pdf"), plot = Fig, width = 5, height = 6, units = "cm")

