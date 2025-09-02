rm(list=ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(parallel)
library(gamlss)
library(scales)
library(readxl)
library(dplyr)

wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")) {
  datapath        <- " "
  interfileFolder <- " "
  functionFolder  <- " "
  resultFolder    <- " "
} else {
  datapath        <- " "
  interfileFolder <- " "
  functionFolder  <- " "
  resultFolder    <- " "
}

## ---- Source function ----
source(file.path(functionFolder, "Compare_distributions_gamlss.R"))

## ---- Load data ----
GNG_data   <- read_xlsx(file.path(datapath, "Q_GNG.xlsx"))
back1_data <- read_xlsx(file.path(datapath, "Q_1back.xlsx"))
back2_data <- read_xlsx(file.path(datapath, "Q_2back.xlsx"))

## ---- Common settings ----
smoothvar  <- "Age_year"
IDvar      <- "ID"
bs.df      <- 3
covariates <- "Sex"
randomvar  <- NA

GNG_data$Sex      <- as.factor(GNG_data$Sex)
back1_data$Sex    <- as.factor(back1_data$Sex)
back2_data$Sex    <- as.factor(back2_data$Sex)

# =====================================================================================
# Step 1. Select best distribution family (three tasks written separately)
# =====================================================================================

## ---- Task 1: Go/No-Go (d_prime) ----
out_GNG <- gamlss_comparedistribution(
  dataname     = "GNG_data",
  dependentvar = "d_prime",
  smoothvar    = smoothvar,
  IDvar        = IDvar,
  bs.df        = bs.df,
  covariates   = covariates,
  randomvar    = randomvar
)
names(out_GNG) <- c("modelsum","performance")
write.csv(out_GNG$performance, file.path(resultFolder, "performance_GNGd_prime_distribution.csv"), row.names = FALSE)
saveRDS(out_GNG$modelsum,   file.path(resultFolder, "modelsum_GNGd_prime_distribution.rds"))

## ---- Task 2: 1-back (Oneback_acc) ----
out_back1 <- gamlss_comparedistribution(
  dataname     = "back1_data",
  dependentvar = "Oneback_acc",
  smoothvar    = smoothvar,
  IDvar        = IDvar,
  bs.df        = bs.df,
  covariates   = covariates,
  randomvar    = randomvar
)
names(out_back1) <- c("modelsum","performance")
write.csv(out_back1$performance, file.path(resultFolder, "performance_1backAcc_distribution.csv"), row.names = FALSE)
saveRDS(out_back1$modelsum,   file.path(resultFolder, "modelsum_1backAcc_distribution.rds"))

## ---- Task 3: 2-back (Twoback_acc) ----
out_back2 <- gamlss_comparedistribution(
  dataname     = "back2_data",
  dependentvar = "Twoback_acc",
  smoothvar    = smoothvar,
  IDvar        = IDvar,
  bs.df        = bs.df,
  covariates   = covariates,
  randomvar    = randomvar
)
names(out_back2) <- c("modelsum","performance")
write.csv(out_back2$performance, file.path(resultFolder, "performance_2backAcc_distribution.csv"), row.names = FALSE)
saveRDS(out_back2$modelsum,   file.path(resultFolder, "modelsum_2backAcc_distribution.rds"))

# =====================================================================================
# Step 2. Grid search of degrees of freedom for chosen family (three tasks separately)
# =====================================================================================

# ## Candidate df grid (mu, sigma, degree)
# con<-gamlss.control(n.cyc=200)
# bs.df.set <- matrix(c(3,3,3,
#                       3,4,3,
#                       3,5,3,
#                       3,6,3,
#                       4,3,3,
#                       4,4,3,
#                       4,5,3,
#                       4,6,3,
#                       5,3,3,
#                       5,4,3,
#                       5,5,3,
#                       5,6,3,
#                       6,3,3,
#                       6,4,3,
#                       6,5,3,
#                       6,6,3,
#                       2,2,2,
#                       2,3,2,
#                       2,4,2,
#                       2,5,2,
#                       2,6,2,
#                       3,2,2,
#                       3,3,2,
#                       3,4,2,
#                       3,5,2,
#                       3,6,2,
#                       4,2,2,
#                       4,3,2,
#                       4,4,2,
#                       4,5,2,
#                       4,6,2,
#                       5,2,2,
#                       5,3,2,
#                       5,4,2,
#                       5,5,2,
#                       5,6,2,
#                       6,2,2,
#                       6,3,2,
#                       6,4,2,
#                       6,5,2,
#                       6,6,2),
#                     byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma", "degree")))
# 
# ## ---- Task 1: Go/No-Go (d_prime) df grid ----
# distribution.fam <- "SEP3"
# out_GNG_bsdf <- gamlss_compare.bs.df(
#   dataname         = "GNG_data",
#   dependentvar     = "d_prime",
#   smoothvar        = smoothvar,
#   IDvar            = IDvar,
#   bs.df.set        = bs.df.set,          # grid of (mu, sigma, degree)
#   covariates       = covariates,
#   distribution.fam = distribution.fam,
#   randomvar        = randomvar
# )
# names(out_GNG_bsdf) <- c("modelsum","performance")
# write.csv(out_GNG_bsdf$performance, file.path(resultFolder, "performance_GNGd_prime_bsdf.csv"), row.names = FALSE)
# saveRDS(out_GNG_bsdf$modelsum,   file.path(resultFolder, "modelsum_GNGd_prime_bsdf.rds"))
# 
# ## ---- Task 2: 1-back (Oneback_acc) df grid ----
# distribution.fam <- "SEP2"
# out_back1_bsdf <- gamlss_compare.bs.df(
#   dataname         = "back1_data",
#   dependentvar     = "Oneback_acc",
#   smoothvar        = smoothvar,
#   IDvar            = IDvar,
#   bs.df.set        = bs.df.set,
#   covariates       = covariates,
#   distribution.fam = distribution.fam,
#   randomvar        = randomvar
# )
# names(out_back1_bsdf) <- c("modelsum","performance")
# write.csv(out_back1_bsdf$performance, file.path(resultFolder, "performance_1backAcc_bsdf.csv"), row.names = FALSE)
# saveRDS(out_back1_bsdf$modelsum,   file.path(resultFolder, "modelsum_1backAcc_bsdf.rds"))
# 
# ## ---- Task 3: 2-back (Twoback_acc) df grid ----
# distribution.fam <- "SEP3"
# out_back2_bsdf <- gamlss_compare.bs.df(
#   dataname         = "back2_data",
#   dependentvar     = "Twoback_acc",
#   smoothvar        = smoothvar,
#   IDvar            = IDvar,
#   bs.df.set        = bs.df.set,
#   covariates       = covariates,
#   distribution.fam = distribution.fam,
#   randomvar        = randomvar
# )
# names(out_back2_bsdf) <- c("modelsum","performance")
# write.csv(out_back2_bsdf$performance, file.path(resultFolder, "performance_2backAcc_bsdf.csv"), row.names = FALSE)
# saveRDS(out_back2_bsdf$modelsum,   file.path(resultFolder, "modelsum_2backAcc_bsdf.rds"))

