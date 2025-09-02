rm(list=ls())
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
  datapath <- ' '
  FigureFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
} else {
  datapath <- ' '
  FigureFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
}

## ========= source function =========
source(paste0(functionFolder, "/gamcog_withsmoothvar_deviation.R"))

#load data
GNGd_data <- read_rds(paste0(datapath, '/Gonogo/GNGd_prime.deviations.rds'))
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
back1_data[,psyc_variables_continous] <- lapply(back1_data[,psyc_variables_continous], as.numeric)
back2_data[,psyc_variables_continous] <- lapply(back2_data[,psyc_variables_continous], as.numeric)

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
corr.result.all <- mclapply(1:nrow(EFvars.set), function(i){
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
corr.result.all.df <- do.call(rbind, corr.result.all)
write.csv(corr.result.all.df, paste0(resultFolder, "/corr_EF_psych_continuous_result_slope.csv"), row.names = FALSE)

corr.result.all.df$anovap.bonf <- p.adjust(corr.result.all.df$anova.pvalues, method = "bonferroni")
corr.result.all.df$sig <- (corr.result.all.df$anovap.bonf < 0.05)
write.csv(corr.result.all.df, file = paste0(resultFolder,"/corr_results_with_bonf.csv"), row.names = FALSE)

## plot
corr.result.all.df$correstimate <- as.numeric(corr.result.all.df$correstimate)
lwth <- min(corr.result.all.df$correstimate, na.rm = TRUE)
upth <- max(corr.result.all.df$correstimate, na.rm = TRUE)
y_levels <- c("SDQ_ES_sum_z", "SDQ_PP_sum_z", "SDQ_CP_sum_z", "SDQ_H_sum_z", "SDQ_PB_sum_z")
corr.result.all.df$correstimate <- as.numeric(corr.result.all.df$correstimate)

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
corr.result.all.df$parcel <- factor(corr.result.all.df$parcel,
                                  levels = y_levels,
                                  labels = psy_labels)
corr.result.all.df$Task <- factor(corr.result.all.df$dataname, levels = c("GNGd", "back1", "back2"),
                                labels = c("Go/No-go", "1-back", "2-back"))

# Manually define task colors
task_colors <- c("Go/No-go" = "#E5E9F2",
                 "1-back" = "#A4C5DF",
                 "2-back" = "#4980B5")

# Calculate the position for significance markers
corr.result.all.df$significance <- ifelse(corr.result.all.df$sig, "*", "")  # Mark with a star
corr.result.all.df$beta <- as.numeric(corr.result.all.df$beta)
corr.result.all.df$label_y <- ifelse(corr.result.all.df$beta > 0, 
                                   corr.result.all.df$beta + 0.01,  # Place the star above the bar top
                                   corr.result.all.df$beta - 0.02)  # Place the star below the bar bottom

## Plotting the figure
y_limits <- c(-0.15, 0.12)  # Symmetric y-limits
vline_positions <- seq(1.5, length(unique(corr.result.all.df$parcel)) - 0.5, by = 1)

Fig <- ggplot(data = corr.result.all.df, aes(x = parcel, y = beta, fill = Task, color = Task)) +
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
  theme(axis.line.x = element_line(color = "black", size = 0.25),  
        axis.line.y = element_line(color = "black", size = 0.25), 
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 9, hjust = 0.5, color = "black"),  # Ensure x-axis labels are centered
        axis.text.y = element_text(size = 9, color = "black"),  
        axis.ticks.x = element_line(color = "black", size = 0.25),
        axis.ticks.y = element_line(color = "black", size = 0.25),
        axis.ticks.length = unit(0.05, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position =c(0.25, 0.9),
        legend.direction = "horizontal" ,
        legend.margin = margin(t = -10),
        legend.key.size = unit(0.3, "cm"), 
        panel.grid.major.y = element_blank(),  
        panel.grid.major.x = element_blank(),  
        panel.grid.minor = element_blank()) 
print(Fig)
ggsave(paste0(FigureFolder, "/correlation_beta_sdq5.pdf"), plot = Fig, width = 15, height = 8, units = "cm")
