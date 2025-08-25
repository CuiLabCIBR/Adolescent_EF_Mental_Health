rm(list=ls())
#library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")) {
  datapath        <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/data"
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/interfileFolder"
  functionFolder  <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/function"
  resultFolder    <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/results"
} else {
  datapath        <- " "   # <- 本地路径记得填
  interfileFolder <- " "
  functionFolder  <- " "
  resultFolder    <- " "
}

## ==== Source helpers ====
source(file.path(functionFolder, "Compare_distributions_gamlss_new.R"))

## ==== TASK ====
# GNG
# TASK <- list(
#   datafile         = "Q_GNG.xlsx",
#   dataname         = "GNG_data",
#   smoothvar        = "Age_year",
#   covariates_hint  = "Sex",            
#   modelsum_file    = "modelsum_GNGd_prime_bsdf.rds",
#   performance_file = "performance_GNGd_prime_bsdf.csv",
#   efvar_label      = "GNGd",
#   out_subdir       = "GNGd_prime"     
# )

# 1-back 示例（如需切换，把上面 TASK 替换为这段，但 out_subdir 如需沿用原有目录名也可不变）
TASK <- list(
  datafile         = "Q_1back.xlsx",
  dataname         = "back1_data",
  smoothvar        = "Age_year",
  covariates_hint  = "Sex",
  modelsum_file    = "modelsum_1backAcc_bsdf.rds",
  performance_file = "performance_1backAcc_bsdf.csv",
  efvar_label      = "1back",
  out_subdir       = "1-back"
)

# # 2-back 示例
# TASK <- list(
#   datafile         = "Q_2back.xlsx",
#   dataname         = "back2_data",
#   smoothvar        = "Age_year",
#   covariates_hint  = "Sex",
#   modelsum_file    = "modelsum_2backAcc_bsdf.rds",
#   performance_file = "performance_2backAcc_bsdf.csv",
#   efvar_label      = "2back",
#   out_subdir       = "2-back" 
# )

## =================== Load data ===================
assign(TASK$dataname, read_xlsx(file.path(datapath, TASK$datafile)))
dat <- get(TASK$dataname)
if ("Sex" %in% names(dat)) dat$Sex <- as.factor(dat$Sex)   # ensure factor for Sex
assign(TASK$dataname, dat)

## =================== Load Step2 results & pick best (by min BIC, no checks) ===================
perf <- read.csv(file.path(resultFolder, TASK$performance_file))
modelsum_list <- readRDS(file.path(resultFolder, TASK$modelsum_file))
best_idx <- which.min(perf$BIC)
model_best <- modelsum_list[[best_idx]]
model_obj  <- model_best

## =================== CLI arg for bootstrap replicate ===================
n <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(n) == 0 || is.na(n)) n <- 1
cat("Bootstrap replicate:", n, "\n")

## =================== Run bootstrap ===================
Base.Seed <- 925
dataname  <- TASK$dataname
smoothvar <- TASK$smoothvar
stratify  <- "Sex"   # stratified bootstrap by Sex

bootstrap.out   <- Boot.Function(n, Base.Seed, dataname, smoothvar, model_obj, stratify)
mod.tmp         <- bootstrap.out$mod.tmp
gam.data.subset <- as.data.frame(bootstrap.out$gam.data.subset)

## =================== Estimate centiles (F/M & overall) ===================
quantiles   <- c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99)
n_quantiles <- length(quantiles)
n_points    <- 1000
x_grid      <- seq(min(dat[[smoothvar]]), max(dat[[smoothvar]]), length.out = n_points)

Centiles_male   <- array(NA, dim = c(n_quantiles, n_points)); rownames(Centiles_male)   <- paste0("quantiles", quantiles)
Centiles_female <- array(NA, dim = c(n_quantiles, n_points)); rownames(Centiles_female) <- paste0("quantiles", quantiles)

for (i in seq_len(n_quantiles)) {
  QuaF <- getQuantile(mod.tmp, quantile = quantiles[i], term = "age",
                      fixed.at = list(Sex = "F"), n.points = n_points)
  Centiles_female[i, ] <- QuaF(x_grid)
  
  QuaM <- getQuantile(mod.tmp, quantile = quantiles[i], term = "age",
                      fixed.at = list(Sex = "M"), n.points = n_points)
  Centiles_male[i, ] <- QuaM(x_grid)
}
Centiles <- (Centiles_female + Centiles_male) / 2

## =================== Derivatives (finite difference) ===================
eps  <- 1e-7
x0   <- x_grid
x1   <- x0 + eps
P50M <- getQuantile(mod.tmp, quantile = 0.5, term = "age", fixed.at = list(Sex = "M"), n.points = n_points)
P50F <- getQuantile(mod.tmp, quantile = 0.5, term = "age", fixed.at = list(Sex = "F"), n.points = n_points)

P50_derivative_M <- (P50M(x1) - P50M(x0)) / eps
P50_derivative_F <- (P50F(x1) - P50F(x0)) / eps
P50_derivative   <- (P50_derivative_M + P50_derivative_F) / 2

lvl <- levels(dat$Sex)
new_m0 <- data.frame(age = x0, Sex = factor("M", levels = lvl))
new_m1 <- data.frame(age = x1, Sex = factor("M", levels = lvl))
new_f0 <- data.frame(age = x0, Sex = factor("F", levels = lvl))
new_f1 <- data.frame(age = x1, Sex = factor("F", levels = lvl))

sigma_predictions_male0   <- predict(mod.tmp, newdata = new_m0, what = "sigma", type = "response")
sigma_predictions_male1   <- predict(mod.tmp, newdata = new_m1, what = "sigma", type = "response")
sigma_predictions_female0 <- predict(mod.tmp, newdata = new_f0, what = "sigma", type = "response")
sigma_predictions_female1 <- predict(mod.tmp, newdata = new_f1, what = "sigma", type = "response")

sigma_derivative_male    <- (sigma_predictions_male1 - sigma_predictions_male0) / eps
sigma_derivative_female  <- (sigma_predictions_female1 - sigma_predictions_female0) / eps
sigma_derivative         <- (sigma_derivative_male + sigma_derivative_female) / 2

## =================== Predictions (keep your original object names) ===================
predict_sigma_F <- data.frame(age = x0, Sex = factor("F", levels = lvl))
predict_sigma_M <- data.frame(age = x0, Sex = factor("M", levels = lvl))
sigma_pred_F    <- predict(mod.tmp, what = "sigma", newdata = predict_sigma_F, type = "response")
sigma_pred_M    <- predict(mod.tmp, what = "sigma", newdata = predict_sigma_M, type = "response")
sigma_pred      <- (sigma_pred_F + sigma_pred_M) / 2

predict_F <- data.frame(age = x0, Sex = factor("F", levels = lvl))
predict_M <- data.frame(age = x0, Sex = factor("M", levels = lvl))
mu_pred_F <- predict(mod.tmp, what = "mu", newdata = predict_F, type = "response")
mu_pred_M <- predict(mod.tmp, what = "mu", newdata = predict_M, type = "response")
mu_pred   <- (mu_pred_F + mu_pred_M) / 2

## =================== Collect results (names unchanged) ===================
bootstrap_centiles <- list(
  Centiles_female         = Centiles_female,
  Centiles_male           = Centiles_male,
  Centiles_overall        = Centiles,
  Age_points              = x0,
  EFvar                   = TASK$efvar_label,
  bootstrap_time          = n,
  converged               = mod.tmp$converged,
  mu_pred                 = mu_pred,
  mu_derivative           = bootstrap.out$mu_derivative,
  sigma_derivative        = bootstrap.out$sigma_derivative,
  P50_derivative          = P50_derivative,
  sigma_derivative_all    = sigma_derivative,
  sigma_derivative_male   = sigma_derivative_male,
  sigma_derivative_female = sigma_derivative_female,
  sigma_pred              = sigma_pred,
  sigma_pred_F            = sigma_pred_F,
  sigma_pred_M            = sigma_pred_M
)

## =================== Save (path & filename unchanged) ===================
outdir <- file.path(interfileFolder, "bootstrap", TASK$out_subdir)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
saveRDS(bootstrap_centiles, file.path(outdir, sprintf("centile_bootstrap_%s.rds", n)))

cat("Bootstrap", n, "finished!\n")