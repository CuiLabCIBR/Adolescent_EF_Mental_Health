rm(list=ls())

# 加载包
library(tidyverse)
library(mgcv)
library(psych)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(ggplot2)
library(parallel)
library(patchwork) 

wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder'
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_2508/function"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_2508/results/figure2_newgng"
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_2508/figure/figure2_newgng"
} else {
  datapath <- '//Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation/data'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z'
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z"
}


source(paste0(functionFolder, "/gamcog_withsmoothvar_deviation.R"))

#head(switch_data)
GNGd_data <- read_rds( '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_2508/results/Gonogo/GNGd_prime.deviations.rds')
back1_data <- read_rds(paste0(datapath, '/1-back/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(datapath, '/2-back/back2Acc.deviations.rds'))


## 1) set up variables
psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
# EF vars
EFvars.set <- matrix(c("d_prime_deviationZ", "GNGd",
                       "Oneback_acc_deviationZ", "back1",
                       "Twoback_acc_deviationZ", "back2"), byrow=TRUE,ncol=2,dimnames=list(NULL,c("varname","dataname")))
EFvars.set <- as.data.frame(EFvars.set)
## 2) convert variables class & describe variables
GNGd_data[,psyc_variables_continous] <- lapply(GNGd_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(GNGd_data$School)
GNGd_data$school_fac <- factor(GNGd_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back1_data[,psyc_variables_continous] <- lapply(back1_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(back1_data$School)
back1_data$school_fac <- factor(back1_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back2_data[,psyc_variables_continous] <- lapply(back2_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(back2_data$School)
back2_data$school_fac <- factor(back2_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

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



# describe
# GNGd
describe_tab_GNGd <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="GNGd"], psyc_variables_continous,  "Age_year", "Sex"), data=GNGd_data,testNonNormal=T)
describe_tab_GNGd.continous <- as.data.frame(describe_tab_GNGd$ContTable[["Overall"]])
describe_tab_GNGd.discrete <- do.call(rbind, lapply(describe_tab_GNGd$CatTable[["Overall"]], as.data.frame))
# back1
describe_tab_back1 <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="back1"], psyc_variables_continous,  "Age_year", "Sex"), data=back1_data, testNonNormal=T)
describe_tab_back1.continous <- as.data.frame(describe_tab_back1$ContTable[["Overall"]])
describe_tab_back1.discrete <- do.call(rbind, lapply(describe_tab_back1$CatTable[["Overall"]], as.data.frame))
# back2
describe_tab_back2 <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="back2"], psyc_variables_continous,  "Age_year", "Sex"), data=back2_data, testNonNormal=T)
describe_tab_back2.continous <- as.data.frame(describe_tab_back2$ContTable[["Overall"]])
describe_tab_back2.discrete <- do.call(rbind, lapply(describe_tab_back2$CatTable[["Overall"]], as.data.frame))

# # save out
write.xlsx(list(GNGd_con=describe_tab_GNGd.continous,GNGd_dis=describe_tab_GNGd.discrete,
                back1_con=describe_tab_back1.continous,back1_dis=describe_tab_back1.discrete,
                back2_con=describe_tab_back2.continous,back2_dis=describe_tab_back2.discrete),
           paste0(resultFolder, "/description_interest_vars.xlsx"), rowNames=T)


## 3) Correlations in separate age periods
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 2)) ###
cl <- makeCluster(n_cores) ###
cat(paste("This run will allocate ", n_cores, "cores.\n")) ###
# continuos variables
corr.result.period <- mclapply(1:nrow(EFvars.set), function(i){
  EFvar.tmp <- EFvars.set$varname[i]
  dataname0 <- paste0(EFvars.set$dataname[i], "_data")
  data.tmp <- get(dataname0)
  
  period <- "EF"
  data_segment <- data.tmp
  
  if (sum(!is.na(data_segment)) < 30) { 
    return(data.frame(
      period = period, 
      dataname = EFvars.set$dataname[i],
      parcel = EFvar.tmp,
      gam.indep.t = NA, gam.indep.pvalue = NA, 
      anova.pvalues = NA, partialRsq = NA, 
      corrmethod = NA, correstimate = NA, 
      corrp = NA, samplesize = NA, beta = NA
    ))
  }
  
  corr.result <- list()
  for (x in seq_along(psyc_variables_continous)) {
    psyvar.tmp <- psyc_variables_continous[x]
    dependentvar <- psyvar.tmp
    interest.indep.var <- EFvar.tmp
    smoothvar <- "Age_year"
    covariates <- "Sex"
    knots = 3
    
    result.tmp <- gam.fit.Independent.var(
      dependentvar = dependentvar,
      dataname = dataname0,
      smoothvar = smoothvar,
      interest.indep.var = interest.indep.var,
      covariates = covariates,
      stats_only = TRUE,
      cl = cl
    )
    result.tmp <- as.data.frame(result.tmp)
    result.tmp$dataname <- EFvars.set$dataname[i]
    result.tmp$period <- period
    corr.result[[x]] <- result.tmp
  }
  corr.result.df <- do.call(rbind, corr.result)
  return(corr.result.df)
})
stopCluster(cl) ###
corr.result.period.df.con <- do.call(rbind, corr.result.period)
write.csv(corr.result.period.df.con, paste0(resultFolder, "/corr_EF_psych_continuous_result_slope.csv"), row.names = FALSE)

corr.result.period.df <- rbind(corr.result.period.df.con)
## plot
corr.result.period.df$correstimate <- as.numeric(corr.result.period.df$correstimate)
lwth <- min(corr.result.period.df$correstimate, na.rm = TRUE)
upth <- max(corr.result.period.df$correstimate, na.rm = TRUE)
y_levels <- c("SDQ_PB_sum_z", "SDQ_H_sum_z", "SDQ_CP_sum_z", "SDQ_PP_sum_z",  "SDQ_ES_sum_z")
# Initialize the result list
updated_results <- list()

# Loop through each EFvar
for (i in 1:nrow(EFvars.set)) {
  # EFvar.tmp <- EFvars.set$varname[i]
  # dataname <- EFvars.set$dataname[i]
  EFvar.tmp <- EFvars.set[i, "varname"]
  dataname <- EFvars.set[i, "dataname"]
  
  # Extract data for the corresponding variable
  corr.result.tmp <- corr.result.period.df[which(corr.result.period.df$interest.indep.var == EFvar.tmp & corr.result.period.df$dataname == dataname), ]
  corr.result.tmp$anova.pvalues <- as.numeric(corr.result.tmp$anova.pvalues)
  
  corr.result.tmp$corrp <- as.numeric(corr.result.tmp$corrp)
  
  # Calculate FDR
  #corr.result.tmp$anovap.fdr <- p.adjust(corr.result.tmp$anova.pvalues, method = "fdr")
  #corr.result.tmp$sig <- (corr.result.tmp$anovap.fdr < 0.05)
  # Calculate bonferroni
  corr.result.tmp$anovap.bonf <- p.adjust(corr.result.tmp$anova.pvalues, method = "bonferroni")
  corr.result.tmp$sig <- (corr.result.tmp$anovap.bonf < 0.05)
  
  # Update the main data frame
  updated_results[[i]] <- corr.result.tmp
  
  # Set y-axis order
  corr.result.tmp$parcel <- factor(corr.result.tmp$parcel, levels = y_levels)
  
  # Plot the figure
  Fig <- ggplot() +
    geom_tile(data = corr.result.tmp, aes(x = period, y = parcel, fill = correstimate), color = "white") +
    geom_text(data = corr.result.tmp[corr.result.tmp$sig == TRUE, ], aes(x = period, y = parcel, label = "*"), vjust = 0.7, hjust = 0.5, size = 6) +
    scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(lwth, upth), direction = -1) +
    scale_y_discrete(limits = y_levels,
                     labels = c("SDQ_PB_sum" = "Prosocial Behavior","SDQ_H_sum" = "Hyperactivity",
                                "SDQ_CP_sum" = "Conduct Problems","SDQ_PP_sum" = "Peer Problems",
                                "SDQ_ES_sum" = "Emotional Symptoms")) +
    labs(title = paste0("Correlation between Executive functions and ", EFvar.tmp, " in ", dataname),
         x = "Periods of adolescence", y = "Psychiatric scores") +
    theme(axis.line = element_blank(),
          aspect.ratio = 1.2,
          axis.text.x = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 1))
  print(Fig)
  #ggsave(paste0(FigureFolder, "/", dataname, "/corr_continuouspsych_", EFvar.tmp, "_3periods.pdf"), plot = Fig, width = 14, height = 18, units = "cm")
}

# Combine all updated results
updated_results_df <- do.call(rbind, updated_results)

# Save as Excel or CSV file
write.csv(updated_results_df, file = paste0(resultFolder,"/corr_results_with_bonf.csv"), row.names = FALSE)


combined_results <- data.frame()
for (i in 1:nrow(EFvars.set)) {
  EFvar.tmp <- EFvars.set[i, "varname"]
  dataname <- EFvars.set[i, "dataname"]
  corr.result.tmp <- corr.result.period.df[which(corr.result.period.df$interest.indep.var == EFvar.tmp & corr.result.period.df$dataname == dataname), ]
  corr.result.tmp$anova.pvalues <- as.numeric(corr.result.tmp$anova.pvalues)
  corr.result.tmp$correstimate <- as.numeric(corr.result.tmp$correstimate)
  corr.result.tmp$anovap.bonf <- p.adjust(corr.result.tmp$anova.pvalues, method = "bonferroni")
  corr.result.tmp$sig <- (corr.result.tmp$anovap.bonf < 0.05)
  corr.result.tmp$corrp.bonf <- p.adjust(corr.result.tmp$corrp, method = "bonferroni")
  corr.result.tmp$period <- factor(corr.result.tmp$period, levels = c("EF"))
  combined_results <- rbind(combined_results, corr.result.tmp)
}

## Data preparation
y_levels <- c("SDQ_ES_sum_z", "SDQ_PP_sum_z", "SDQ_CP_sum_z", "SDQ_H_sum_z", "SDQ_PB_sum_z")
combined_results$correstimate <- as.numeric(combined_results$correstimate)

# Custom labels for psychiatric variables
psy_labels <- c("SDQ_ES_sum_z" = "Emotional
Symptoms",
                "SDQ_PP_sum_z" = "Peer
Problems", 
                "SDQ_CP_sum_z" = "Conduct
Problems",
                "SDQ_H_sum_z" = "Hyperactivity", 
                "SDQ_PB_sum_z" = "Prosocial
Behavior")

# Apply new labels
combined_results$parcel <- factor(combined_results$parcel,
                                  levels = y_levels,
                                  labels = psy_labels)
combined_results$Task <- factor(combined_results$dataname, levels = c("GNGd", "back1", "back2"),
                                labels = c("Go/No-go", "1-back", "2-back"))

# Manually define task colors
task_colors <- c("Go/No-go" = "#E5E9F2",
                 "1-back" = "#A4C5DF",
                 "2-back" = "#4980B5")

# Calculate the position for significance markers
combined_results$significance <- ifelse(combined_results$sig, "*", "")  # Mark with a star
combined_results$beta <- as.numeric(combined_results$beta)
combined_results$label_y <- ifelse(combined_results$beta > 0, 
                                   combined_results$beta + 0.01,  # Place the star above the bar top
                                   combined_results$beta - 0.02)  # Place the star below the bar bottom

## Plotting the figure
y_limits <- c(-0.15, 0.12)  # Symmetric y-limits
vline_positions <- seq(1.5, length(unique(combined_results$parcel)) - 0.5, by = 1)

Fig <- ggplot(data = combined_results, aes(x = parcel, y = beta, fill = Task, color = Task)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75, size = 0) +  # Fill colors and borders
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.25) +  # Add horizontal line at y = 0
  geom_text(aes(label = significance, y = label_y), 
            position = position_dodge(width = 0.8), size = 5, color = "black") +  # Add significance stars
  scale_fill_manual(values = task_colors) +  # Manually set task fill colors
  scale_color_manual(values = task_colors) +  # Manually set task border colors
  scale_y_continuous(limits = y_limits, breaks = seq(-0.15, 0.1, by = 0.05), labels = scales::number_format()) +  # Set y-axis limits and ticks
  labs(title = "",
       x = "",
       y = "β",
       color = "Tasks",
       fill = "Tasks") +  # Legend for tasks
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 0.25),  # Add x-axis line
        axis.line.y = element_line(color = "black", size = 0.25),  # Add y-axis line
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 9, hjust = 0.5, color = "black"),  # Ensure x-axis labels are centered
        axis.text.y = element_text(size = 9, color = "black"),  
        axis.ticks.x = element_line(color = "black", size = 0.25),
        axis.ticks.y = element_line(color = "black", size = 0.25),
        axis.ticks.length = unit(0.05, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        #legend.position = c(0.95, 0.15),
        legend.position =c(0.25, 0.9),
        legend.direction = "horizontal" ,
        legend.margin = margin(t = -10),
        legend.key.size = unit(0.3, "cm"), 
        panel.grid.major.y = element_blank(),  # Remove x-axis grid lines
        panel.grid.major.x = element_blank(),  # Remove x-axis grid lines
        panel.grid.minor = element_blank()) 

print(Fig)
ggsave(paste0(FigureFolder, "/correlation_beta_sdq5.pdf"), plot = Fig, width = 15, height = 8, units = "cm")
