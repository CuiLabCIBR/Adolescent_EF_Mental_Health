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
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- ' '
  functionFolder <- ' '
  FigureFolder <- ' '
  resultFolder <- ' '
}else{
  datapath <- ' '
  FigureFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
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
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))
if (is.na(n_cores) || n_cores < 1) n_cores <- 2

cl <- parallel::makeCluster(n_cores, type = "PSOCK")
on.exit(parallel::stopCluster(cl), add = TRUE)

cat(sprintf("This run will allocate %d cores.\n", n_cores))

parallel::clusterEvalQ(cl, { library(lme4); library(pbkrtest); NULL })
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
stopCluster(cl) ###
# Combine results
corr.result.df <- do.call(rbind, corr.result)
write.csv(corr.result.df, paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv"), row.names = FALSE)

#corr.result.df <- read_csv(paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv"))
corr.result.df$correstimate <- as.numeric(corr.result.df$correstimate)
corr.result.df$anovap.bonf <- p.adjust(corr.result.df$boots.pvalues, method = "bonferroni")
corr.result.df$sig <- (corr.result.df$anovap.bonf < 0.05)
corr.result.df$significance <- ifelse(corr.result.df$sig, "*", "")

write.csv(corr.result.df, paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD_withbonf.csv"), row.names = FALSE)


####new plot
psy_labels <- c("cbcl_scr_syn_attention_r_z" = "Attention", 
                "cbcl_scr_syn_external_r_z" = "Externalizing",
                "cbcl_scr_syn_social_r_z" = "Social",
                "cbcl_scr_syn_internal_r_z" = "Internalizing")


psyc_variables_continous <- names(psy_labels)

corr.result.df <- corr.result.df %>%
  arrange(slope) %>%
  mutate(parcel = factor(parcel, levels = parcel))  

corr.result.df$fill_color <- "#c6d6ea"  
corr.result.df$label_x <- ifelse(
  corr.result.df$slope >= 0, 
  corr.result.df$slope + 0.004, 
  corr.result.df$slope - 0.004
)

corr.result.df$label_text <- as.character(corr.result.df$parcel)

# plot
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
ggsave(paste0(FigureFolder, "/correlation_barplot_slope_ABCD.pdf"), plot = Fig, width = 5, height = 6, units = "cm")

