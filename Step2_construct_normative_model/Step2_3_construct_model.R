rm(list=ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath        <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/data"
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/interfileFolder"
  functionFolder  <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/function"
  resultFolder    <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/results"
  FigureFolder    <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_final/figure"
} else {
  datapath        <- " "
  FigureFolder    <- " "
  interfileFolder <- " "
  functionFolder  <- " "
  resultFolder    <- " "
}

## =================== soure function ===================
source(file.path(functionFolder, "Construct_gamlss_set_new.R"))

# Map distribution family string to the correct CDF for centile computation
get_pfun <- function(fam) {
  switch(toupper(fam),
         "SEP3" = pSEP3,
         "SEP2" = pSEP2,
         stop(sprintf("No CDF mapped for family '%s'. Extend get_pfun() if needed.", fam))
  )
}

## =================== TASK: ===================
# ---- GNG ----
TASK <- list(
  datafile        = "Q_GNG.xlsx",       # Input data file
  dataname        = "GNG_data",         # R object name for dataset
  smoothterm      = "Age_year",         # Smooth term in GAMLSS
  dependentvar    = "d_prime",          # Dependent variable
  covariates      = "Sex",              # Covariates in the model
  IDvar           = "ID",               # Subject ID column
  distribution    = "SEP3",             # Distribution family
  mu.df           = 2,                  # df for mu
  sigma.df        = 2,                  # df for sigma
  degree          = 2,                  # spline degree
  stratify        = "Sex",              # stratification var for CV
  quantiles       = c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99),
  out_folder_name = "Gonogo",           # output subfolder
  out_file_base   = "GNGd_prime"        # output filename prefix
)

# # ---- 1-back ----
# TASK <- list(
#   datafile        = "Q_1back.xlsx",
#   dataname        = "back1_data",
#   smoothterm      = "Age_year",
#   dependentvar    = "Oneback_acc",
#   covariates      = "Sex",
#   IDvar           = "ID",
#   distribution    = "SEP2",
#   mu.df           = 2,
#   sigma.df        = 2,
#   degree          = 2,
#   stratify        = "Sex",
#   quantiles       = c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99),
#   out_folder_name = "1-back",
#   out_file_base   = "Oneback_acc"
# )

# # ---- 2-back ----
# TASK <- list(
#   datafile        = "Q_2back.xlsx",
#   dataname        = "back2_data",
#   smoothterm      = "Age_year",
#   dependentvar    = "Twoback_acc",
#   covariates      = "Sex",
#   IDvar           = "ID",
#   distribution    = "SEP3",
#   mu.df           = 2,
#   sigma.df        = 2,
#   degree          = 2,
#   stratify        = "Sex",
#   quantiles       = c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99),
#   out_folder_name = "2-back",
#   out_file_base   = "Twoback_acc"
# )

## =================== Load data ===================
assign(TASK$dataname, read_xlsx(file.path(datapath, TASK$datafile)))
dat <- get(TASK$dataname)

## =================== Two-fold stratified split (by Sex) ===================
set.seed(925)
stratify_var <- dat[[TASK$stratify]]
SPLIT <- split(seq_len(nrow(dat)), stratify_var)
idx_set1 <- unlist(lapply(SPLIT, function(X) sample(x = X, size = floor(length(X) * 1/2), replace = FALSE)))
set1 <- dat[idx_set1, ]
set2 <- dat[-idx_set1, ]

## =================== Select modeling variables + clean ===================
keep_cols <- c(TASK$smoothterm, TASK$covariates, "School", TASK$IDvar, TASK$dependentvar)
set1 <- set1 %>% dplyr::select(any_of(keep_cols)) %>% drop_na()
set2 <- set2 %>% dplyr::select(any_of(keep_cols)) %>% drop_na()
if ("Sex" %in% names(set1)) set1$Sex <- as.factor(set1$Sex)
if ("Sex" %in% names(set2)) set2$Sex <- as.factor(set2$Sex)

# Assign to symbols expected by construct_gamlss()
assign(paste0(TASK$dataname, ".set1"), set1)
assign(paste0(TASK$dataname, ".set2"), set2)
dataname_set1 <- paste0(TASK$dataname, ".set1")
dataname_set2 <- paste0(TASK$dataname, ".set2")

## =================== Fit on set1, predict on set2 ===================
if (!dir.exists(file.path(interfileFolder, TASK$out_folder_name))) {
  dir.create(file.path(interfileFolder, TASK$out_folder_name), recursive = TRUE)
}

if (str_detect(wd, "cuizaixu_lab")){
  mod.set1 <- construct_gamlss(
    dataname        = dataname_set1,
    dependentvar    = TASK$dependentvar,
    smoothterm      = TASK$smoothterm,
    covariates      = TASK$covariates,
    mu.df           = TASK$mu.df,
    sigma.df        = TASK$sigma.df,
    degree          = TASK$degree,
    distribution.fam= TASK$distribution,
    IDvar           = TASK$IDvar,
    quantile.vec    = TASK$quantiles,
    stratify        = TASK$stratify,
    randomvar       = NA
  )
  saveRDS(mod.set1, file.path(interfileFolder, TASK$out_folder_name,
                              sprintf("GAMLSS.%sset1.sum.rds", TASK$out_file_base)))
} else {
  mod.set1 <- readRDS(file.path(interfileFolder, TASK$out_folder_name,
                                sprintf("GAMLSS.%sset1.sum.rds", TASK$out_file_base)))
}

modelperformance.set1 <- mod.set1$performance.tb
cat(sum(modelperformance.set1$converged), "models converged (set1).\n")

# Predict to set2 and compute centiles/Z
mod.tmp   <- mod.set1$mod.tmp
pfun      <- get_pfun(TASK$distribution)
mu_pred   <- predict(mod.tmp, newdata = set2, what = "mu",    type = "response")
sigma_pred<- predict(mod.tmp, newdata = set2, what = "sigma", type = "response")
nu_pred   <- predict(mod.tmp, newdata = set2, what = "nu",    type = "response")
tau_pred  <- predict(mod.tmp, newdata = set2, what = "tau",   type = "response")

obs2 <- set2[[TASK$dependentvar]]
dev_set2 <- data.frame(ID = set2[[TASK$IDvar]])
centile2 <- pfun(obs2, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
dev_set2[[paste0(TASK$dependentvar, "_centile")]]    <- centile2
dev_set2[[paste0(TASK$dependentvar, "_deviationZ")]] <- qnorm(centile2)
dev_set2[[paste0(TASK$dependentvar, "_sigma")]]      <- sigma_pred
dev_set2[[paste0(TASK$dependentvar, "_mu")]]         <- mu_pred
dev_set2[[paste0(TASK$dependentvar, "_nu")]]         <- nu_pred
dev_set2[[paste0(TASK$dependentvar, "_tau")]]        <- tau_pred

saveRDS(dev_set2, file.path(interfileFolder, TASK$out_folder_name,
                            sprintf("EF_%s.set2_deviation.rds", TASK$out_file_base)))

## =================== Fit on set2, predict on set1 ===================
if (str_detect(wd, "cuizaixu_lab")){
  mod.set2 <- construct_gamlss(
    dataname        = dataname_set2,
    dependentvar    = TASK$dependentvar,
    smoothterm      = TASK$smoothterm,
    covariates      = TASK$covariates,
    mu.df           = TASK$mu.df,
    sigma.df        = TASK$sigma.df,
    degree          = TASK$degree,
    distribution.fam= TASK$distribution,
    IDvar           = TASK$IDvar,
    quantile.vec    = TASK$quantiles,
    stratify        = TASK$stratify,
    randomvar       = NA
  )
  saveRDS(mod.set2, file.path(interfileFolder, TASK$out_folder_name,
                              sprintf("GAMLSS.%sset2.sum.rds", TASK$out_file_base)))
} else {
  mod.set2 <- readRDS(file.path(interfileFolder, TASK$out_folder_name,
                                sprintf("GAMLSS.%sset2.sum.rds", TASK$out_file_base)))
}

modelperformance.set2 <- mod.set2$performance.tb
cat(sum(modelperformance.set2$converged), "models converged (set2).\n")

mod.tmp   <- mod.set2$mod.tmp
mu_pred   <- predict(mod.tmp, newdata = set1, what = "mu",    type = "response")
sigma_pred<- predict(mod.tmp, newdata = set1, what = "sigma", type = "response")
nu_pred   <- predict(mod.tmp, newdata = set1, what = "nu",    type = "response")
tau_pred  <- predict(mod.tmp, newdata = set1, what = "tau",   type = "response")

obs1 <- set1[[TASK$dependentvar]]
dev_set1 <- data.frame(ID = set1[[TASK$IDvar]])
centile1 <- pfun(obs1, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
dev_set1[[paste0(TASK$dependentvar, "_centile")]]    <- centile1
dev_set1[[paste0(TASK$dependentvar, "_deviationZ")]] <- qnorm(centile1)
dev_set1[[paste0(TASK$dependentvar, "_sigma")]]      <- sigma_pred
dev_set1[[paste0(TASK$dependentvar, "_mu")]]         <- mu_pred
dev_set1[[paste0(TASK$dependentvar, "_nu")]]         <- nu_pred
dev_set1[[paste0(TASK$dependentvar, "_tau")]]        <- tau_pred

saveRDS(dev_set1, file.path(interfileFolder, TASK$out_folder_name,
                            sprintf("EF_%s.set1_deviation.rds", TASK$out_file_base)))

## =================== Export split datasets (for reproducibility) ===================
write.csv(set1, file.path(interfileFolder, TASK$out_folder_name,
                          sprintf("%s.data1.csv", TASK$out_file_base)), row.names = FALSE)
write.csv(set2, file.path(interfileFolder, TASK$out_folder_name,
                          sprintf("%s.data2.csv", TASK$out_file_base)), row.names = FALSE)

## =================== Merge deviations and save (keeps your naming) ===================
dev_merged <- dplyr::bind_rows(dev_set1, dev_set2) %>%
  # Join back any columns you want to carry through (ID + meta columns)
  dplyr::left_join(
    dat %>% dplyr::select(any_of(unique(c(TASK$IDvar, keep_cols)))),
    by = setNames(TASK$IDvar, "ID")
  )

saveRDS(dev_merged, file.path(interfileFolder, TASK$out_folder_name,
                              sprintf("%s.deviations.rds", TASK$out_file_base)))
write.csv(dev_merged, file.path(interfileFolder, TASK$out_folder_name,
                                sprintf("%s.deviations.csv", TASK$out_file_base)), row.names = FALSE)

## =================== Fit on full data and export ===================
data_all <- dat %>% dplyr::select(any_of(keep_cols)) %>% drop_na()
if ("Sex" %in% names(data_all)) data_all$Sex <- as.factor(data_all$Sex)
assign(paste0(TASK$dataname, "_all"), data_all)

if (str_detect(wd, "cuizaixu_lab")){
  mod.all <- construct_gamlss(
    dataname        = paste0(TASK$dataname, "_all"),
    dependentvar    = TASK$dependentvar,
    smoothterm      = TASK$smoothterm,
    covariates      = TASK$covariates,
    mu.df           = TASK$mu.df,
    sigma.df        = TASK$sigma.df,
    degree          = TASK$degree,
    distribution.fam= TASK$distribution,
    IDvar           = TASK$IDvar,
    quantile.vec    = TASK$quantiles,
    stratify        = TASK$stratify,
    randomvar       = NA
  )
  saveRDS(mod.all, file.path(interfileFolder, TASK$out_folder_name,
                             sprintf("GAMLSS.%s.all.sum.rds", TASK$out_file_base)))
} else {
  mod.all <- readRDS(file.path(interfileFolder, TASK$out_folder_name,
                               sprintf("GAMLSS.%s.all.sum.rds", TASK$out_file_base)))
}

mod.tmp   <- mod.all$mod.tmp
mu_pred   <- predict(mod.tmp, newdata = data_all, what = "mu",    type = "response")
sigma_pred<- predict(mod.tmp, newdata = data_all, what = "sigma", type = "response")
nu_pred   <- predict(mod.tmp, newdata = data_all, what = "nu",    type = "response")
tau_pred  <- predict(mod.tmp, newdata = data_all, what = "tau",   type = "response")

obs_all <- data_all[[TASK$dependentvar]]
dev_all <- data.frame(ID = data_all[[TASK$IDvar]])
cent_all <- get_pfun(TASK$distribution)(obs_all, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
dev_all[[paste0(TASK$dependentvar, "_centile")]]    <- cent_all
dev_all[[paste0(TASK$dependentvar, "_deviationZ")]] <- qnorm(cent_all)
dev_all[[paste0(TASK$dependentvar, "_sigma")]]      <- sigma_pred
dev_all[[paste0(TASK$dependentvar, "_mu")]]         <- mu_pred
dev_all[[paste0(TASK$dependentvar, "_nu")]]         <- nu_pred
dev_all[[paste0(TASK$dependentvar, "_tau")]]        <- tau_pred

dev_all <- merge(dev_all, data_all, by = "ID")
saveRDS(dev_all, file.path(interfileFolder, TASK$out_folder_name,
                           sprintf("%s.deviation_all.rds", TASK$out_file_base)))
write.csv(dev_all, file.path(interfileFolder, TASK$out_folder_name,
                             sprintf("%s.deviation_all.csv", TASK$out_file_base)), row.names = FALSE)

cat("All done for task:", TASK$out_file_base, "\n")