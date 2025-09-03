rm(list=ls())
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- ' '
  interfileFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
}else{
  datapath <-  ''
  FigureFolder <- ' '
  interfileFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
  
}

# load data
Flanker_data <- read_xlsx(paste0(datapath, "/ABCD_Flanker.xlsx"))
Flanker_data$Sex <- factor(Flanker_data$Sex, levels = c("1", "2"))
modelsum_Flanker_bsdf <- readRDS(paste0(resultFolder, "/modelsum_Flanker_bsdf.rds"))
performance_Flanker_bsdf <- read.csv(paste0(resultFolder, "/performance_Flanker_bsdf.csv"))
model.Flanker <- modelsum_Flanker_bsdf[[which.min(performance_Flanker_bsdf$BIC)]]
# source functions
source(paste0(functionFolder, "/Compare_distributions_gamlss.R"))
# get input
n <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(n)
# run bootstrap
Flanker_data$Sex <- as.factor(Flanker_data$Sex)
print(n)
Flanker_data$Age_year <- as.numeric(Flanker_data$Age_year)
Base.Seed <- 925
dataname <- "Flanker_data"
smoothvar <- "Age_year"
model_obj <- model.Flanker
stratify <- "Sex"

bootstrap.out <- Boot.Function(n, Base.Seed, dataname,smoothvar, model_obj, stratify)
mod.tmp <- bootstrap.out$mod.tmp
gam.data.subset <- bootstrap.out$gam.data.subset

# estimate quantiles
quantiles<-c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99)
n_quantiles <- length(quantiles)
n_points <- 10000
Centiles_male <- array(NA, dim = c(n_quantiles, n_points))
rownames(Centiles_male) <- paste0("quantiles", quantiles)
Centiles_female <- array(NA, dim = c(n_quantiles, n_points))
rownames(Centiles_female) <- paste0("quantiles", quantiles)
x_female <- seq(min(Flanker_data[[smoothvar]]), max(Flanker_data[[smoothvar]]), length.out = n_points)
x_male <- seq(min(Flanker_data[[smoothvar]]), max(Flanker_data[[smoothvar]]), length.out = n_points)

gam.data.subset <- as.data.frame(gam.data.subset)
for (i in 1:n_quantiles){
  command <- "Qua <- getQuantile(mod.tmp, quantile=quantiles[i], term = 'age', fixed.at = list(Sex='2'), n.points = n_points)"
  eval(parse(text=command))
  Centiles_female[i,] <- Qua(x_female)
  
  command <- "Qua <- getQuantile(mod.tmp, quantile=quantiles[i], term = 'age', fixed.at = list(Sex='1'), n.points = n_points)"
  eval(parse(text=command))
  Centiles_male[i,] <- Qua(x_male)
}

Centiles <- (Centiles_female + Centiles_male) /2


eps <- 1e-07 
x.tmp0 <- seq(min(Flanker_data[[smoothvar]]), max(Flanker_data[[smoothvar]]), length.out = n_points)
x.tmp1 <- x.tmp0 + eps
P50_function_M <- getQuantile(mod.tmp, quantile = 0.5, term = "age", fixed.at = list(Sex = '1'), n.points = n_points)
P50_function_F <- getQuantile(mod.tmp, quantile = 0.5, term = "age", fixed.at = list(Sex = '2'), n.points = n_points)

P50_values_tmp0_M <- P50_function_M(x.tmp0)
P50_values_tmp1_M <- P50_function_M(x.tmp1)
P50_derivative_M <- (P50_values_tmp1_M - P50_values_tmp0_M) / eps

P50_values_tmp0_F <- P50_function_F(x.tmp0)
P50_values_tmp1_F <- P50_function_F(x.tmp1)
P50_derivative_F <- (P50_values_tmp1_F - P50_values_tmp0_F) / eps
P50_derivative <- (P50_derivative_M + P50_derivative_F) / 2

newdata_male0 <- data.frame("age" = x.tmp0, Sex = factor("1", levels = c("1", "2")))
newdata_male1 <- data.frame("age" = x.tmp1, Sex = factor("1", levels = c("1", "2")))
newdata_female0 <- data.frame("age" = x.tmp0, Sex = factor("2", levels = c("1", "2")))
newdata_female1 <- data.frame("age" = x.tmp1, Sex = factor("2", levels = c("1", "2")))

sigma_predictions_male0 <- predict(mod.tmp, newdata = newdata_male0, what = "sigma", type = "response")
sigma_predictions_male1 <- predict(mod.tmp, newdata = newdata_male1, what = "sigma", type = "response")
sigma_predictions_female0 <- predict(mod.tmp, newdata = newdata_female0, what = "sigma", type = "response")
sigma_predictions_female1 <- predict(mod.tmp, newdata = newdata_female1, what = "sigma", type = "response")

sigma_derivative_male <- (sigma_predictions_male1 - sigma_predictions_male0) / eps
sigma_derivative_female <- (sigma_predictions_female1 - sigma_predictions_female0) / eps
sigma_derivative <- (sigma_derivative_male + sigma_derivative_female) / 2

predict_sigma_F <- data.frame("age" = x.tmp0, Sex = '2')
predict_sigma_M <- data.frame("age" = x.tmp0, Sex = '1')

sigma_pred_F <- predict(mod.tmp, what = "sigma", newdata = predict_sigma_F, type = "response",data = gam.data.subset)
sigma_pred_M <- predict(mod.tmp, what = "sigma", newdata = predict_sigma_M, type = "response",data = gam.data.subset)
sigma_pred <- (sigma_pred_F + sigma_pred_M) / 2
# mu predict
predict_F <- data.frame("age" = x.tmp0, Sex = '2')
predict_M <- data.frame("age" = x.tmp0, Sex = '1')

mu_pred_F <- predict(mod.tmp, what = "mu", newdata = predict_F, type = "response",data = gam.data.subset)
mu_pred_M <- predict(mod.tmp, what = "mu", newdata = predict_M, type = "response",data = gam.data.subset)
mu_pred <- (mu_pred_F + mu_pred_M) / 2


has_nu <- "nu.formula" %in% names(mod.tmp)

# Check if the model has a 'tau' parameter
has_tau <- "tau.formula" %in% names(mod.tmp)

# Nu predictions and derivatives
if (has_nu) {
  nu_predictions_male0 <- predict(mod.tmp, newdata = newdata_male0, what = "nu", type = "response")
  nu_predictions_male1 <- predict(mod.tmp, newdata = newdata_male1, what = "nu", type = "response")
  nu_predictions_female0 <- predict(mod.tmp, newdata = newdata_female0, what = "nu", type = "response")
  nu_predictions_female1 <- predict(mod.tmp, newdata = newdata_female1, what = "nu", type = "response")
  
  nu_derivative_male <- (nu_predictions_male1 - nu_predictions_male0) / eps
  nu_derivative_female <- (nu_predictions_female1 - nu_predictions_female0) / eps
  nu_derivative <- (nu_derivative_male + nu_derivative_female) / 2
  
  nu_pred_F <- predict(mod.tmp, what = "nu", newdata = predict_F, type = "response", data = gam.data.subset)
  nu_pred_M <- predict(mod.tmp, what = "nu", newdata = predict_M, type = "response", data = gam.data.subset)
  nu_pred <- (nu_pred_F + nu_pred_M) / 2
} else {
  nu_pred <- NA
  nu_derivative <- NA
  nu_derivative_male <- NA
  nu_derivative_female <- NA
}

# Tau predictions and derivatives (conditional)
if (has_tau) {
  tau_predictions_male0 <- predict(mod.tmp, newdata = newdata_male0, what = "tau", type = "response")
  tau_predictions_male1 <- predict(mod.tmp, newdata = newdata_male1, what = "tau", type = "response")
  tau_predictions_female0 <- predict(mod.tmp, newdata = newdata_female0, what = "tau", type = "response")
  tau_predictions_female1 <- predict(mod.tmp, newdata = newdata_female1, what = "tau", type = "response")
  
  tau_derivative_male <- (tau_predictions_male1 - tau_predictions_male0) / eps
  tau_derivative_female <- (tau_predictions_female1 - tau_predictions_female0) / eps
  tau_derivative <- (tau_derivative_male + tau_derivative_female) / 2
  
  tau_pred_F <- predict(mod.tmp, what = "tau", newdata = predict_F, type = "response", data = gam.data.subset)
  tau_pred_M <- predict(mod.tmp, what = "tau", newdata = predict_M, type = "response", data = gam.data.subset)
  tau_pred <- (tau_pred_F + tau_pred_M) / 2
} else {
  tau_pred <- NA
  tau_derivative <- NA
  tau_derivative_male <- NA
  tau_derivative_female <- NA
}

# Extract results
bootstrap_centiles <- list(
  Centiles_female = Centiles_female,
  Centiles_male = Centiles_male,
  Centiles_overall = Centiles,
  Age_points = x.tmp0,
  EFvar = "Flanker",
  bootstrap_time = n,
  converged = mod.tmp$converged,
  mu_pred = mu_pred,
  mu_derivative = bootstrap.out$mu_derivative,
  sigma_derivative = bootstrap.out$sigma_derivative,
  P50_derivative = P50_derivative,
  sigma_derivative_all = sigma_derivative,
  sigma_derivative_male = sigma_derivative_male,
  sigma_derivative_female = sigma_derivative_female,
  sigma_pred = sigma_pred,
  sigma_pred_F = sigma_pred_F,
  sigma_pred_M = sigma_pred_M,
  nu_pred = nu_pred,
  nu_derivative = nu_derivative,
  nu_derivative_male = nu_derivative_male,
  nu_derivative_female = nu_derivative_female,
  tau_pred = tau_pred,
  tau_derivative = tau_derivative,
  tau_derivative_male = tau_derivative_male,
  tau_derivative_female = tau_derivative_female
)

# save out
if (! dir.exists(paste0(interfileFolder, "/bootstrap/Flanker"))){
  dir.create(paste0(interfileFolder, "/bootstrap/Flanker"))
}

saveRDS(bootstrap_centiles, paste0(interfileFolder, "/bootstrap/Flanker/centile_bootstrap_", n, ".rds"))

print(paste("Bootsrap", n, "finished!"))
