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
  datapath <- ' '
  interfileFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
}else{
  datapath <- ' '
  FigureFolder <- ' '
  interfileFolder <- ' '
  functionFolder <- ' '
  resultFolder <- ' '
}

# load data
source(paste0(functionFolder, "/Construct_gamlss_set.R"))
Flanker_data <- read_xlsx(paste0(datapath, "/ABCD_Flanker.xlsx"))


# full data model
Flanker_data_all <- Flanker_data %>% select(c("Age_year", "Sex", "site_id_l","ID","nihtbx_flanker_uncorrected")) %>% drop_na()
Flanker_data_all$Sex <- as.factor(Flanker_data_all$Sex)

dataname <- "Flanker_data_all"
dependentvar <- "nihtbx_flanker_uncorrected"
smoothterm <- "Age_year"
covariates <- "Sex"
randomvar <- NA
mu.df <- sigma.df <- 4
degree <- 2
distribution.fam <- "GG"
IDvar <- "ID"
stratify <- "Sex"
quantile.vec <- c(0.01,0.025, 0.05, 0.25, 0.5, 0.75, 0.95,0.975, 0.99)

# full model
if (str_detect(wd, "cuizaixu_lab")) {
  mod.all <- construct_gamlss(dataname, dependentvar, smoothterm, covariates,
                              mu.df, sigma.df, degree, distribution.fam,
                              IDvar, quantile.vec, stratify, randomvar)
  saveRDS(mod.all, paste0(interfileFolder, "/Flanker/GAMLSS.Flanker.all.sum.rds"))
} else {
  mod.all <- readRDS(paste0(interfileFolder, "/Flanker/GAMLSS.Flanker.all.sum.rds"))
}

# compute deviation
gam.data2 <- Flanker_data_all
mod.tmp <- mod.all$mod.tmp
mu_pred <- predict(mod.tmp, newdata = Flanker_data_all, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = Flanker_data_all, what = "sigma", type = "response")
nu_pred <- predict(mod.tmp, newdata = Flanker_data_all, what = "nu", type = "response")

# centile and deviationZ ===
observation <- Flanker_data_all[[dependentvar]]
deviation.df <- data.frame(ID = Flanker_data_all$ID)
deviation.df[[paste0(dependentvar, "_centile")]] <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
deviation.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(deviation.df[[paste0(dependentvar, "_centile")]])
deviation.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.df[[paste0(dependentvar, "_mu")]] <- mu_pred
deviation.df[[paste0(dependentvar, "_nu")]] <- nu_pred

#combine data and save
deviation.df <- merge(deviation.df, Flanker_data_all, by = "ID")
saveRDS(deviation.df, paste0(interfileFolder, "/Flanker/Flanker.deviation_all.rds"))
write.csv(deviation.df, paste0(interfileFolder, "/Flanker/Flanker.deviation_all.csv"), row.names = FALSE)



## 1st. 2 fold CV, sampling, set1: 1/2; set2: 1/2, statified by Sex
set.seed(925)
stratify <- "Sex"
stratify.var <- Flanker_data[ ,stratify]
SPLIT <- split(1:NROW(Flanker_data), stratify.var)
LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X)*1/2,replace=F)})
index_set1_b <- unlist(LAPPLY)
Flanker_data.set1 <- Flanker_data[index_set1_b,]
Flanker_data.set2 <- Flanker_data %>% filter(!(row_number() %in% index_set1_b))#Subtract set1 from Flanker_data
## 2nd, construct models
Flanker_data.set1 <- Flanker_data.set1 %>% select(c("Age_year", "Sex", "site_id_l","ID","nihtbx_flanker_uncorrected")) %>% drop_na()
Flanker_data.set2 <- Flanker_data.set2 %>% select(c("Age_year", "Sex", "site_id_l","ID","nihtbx_flanker_uncorrected")) %>% drop_na()
Flanker_data.set1$Sex <- as.factor(Flanker_data.set1$Sex)
Flanker_data.set2$Sex <- as.factor(Flanker_data.set2$Sex)
### set1
dataname <- "Flanker_data.set1"
Flanker_data.set1 <- as.data.frame(Flanker_data.set1)
smoothterm <- "Age_year"
dependentvar <- "nihtbx_flanker_uncorrected"
randomvar <- NA
mu.df <- sigma.df <- 4
degree <- 2
distribution.fam <- "GG"
IDvar <- "ID"
covariates <- "Sex"
stratify <- "Sex"
quantile.vec <- c(0.01,0.025, 0.05, 0.25, 0.5, 0.75, 0.95,0.975, 0.99)

if (str_detect(wd, "cuizaixu_lab")){
  mod.set1 <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify,randomvar=NA)
  saveRDS(mod.set1, paste0(interfileFolder, "/GAMLSS.Flankerset1.sum.rds"))
}else{
  mod.set1 <- readRDS(paste0(interfileFolder, "/GAMLSS.Flankerset1.sum.rds"))
}


# performance
modelperformance.set1 <- mod.set1$performance.tb
print(paste(sum(modelperformance.set1$converged), "models converged.")) 

# compute deviations
gam.data2 <- Flanker_data.set1
mod.tmp <- mod.set1$mod.tmp

mu_pred <- predict(mod.tmp, newdata = Flanker_data.set2, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = Flanker_data.set2, what = "sigma", type = "response")
nu_pred <- predict(mod.tmp, newdata = Flanker_data.set2, what = "nu", type = "response")

dependentvar <- "nihtbx_flanker_uncorrected"
deviation.set2.df <- data.frame(ID=Flanker_data.set2$ID)
observation <- Flanker_data.set2[[dependentvar]]
centile <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
deviation.set2.df[[paste0(dependentvar, "_centile")]] <- centile
deviation.set2.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(centile)
deviation.set2.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.set2.df[[paste0(dependentvar, "_mu")]] <- mu_pred
deviation.set2.df[[paste0(dependentvar, "_nu")]] <- nu_pred

#saveRDS(deviation.set2.df, paste0(interfileFolder, "/EF_Flanker.set2_deviation.rds"))


# set2
dataname <- "Flanker_data.set2"
if (str_detect(wd, "cuizaixu_lab")){
  mod.set2 <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify,randomvar=NA)
  saveRDS(mod.set2, paste0(interfileFolder, "/GAMLSS.Flankerset2.sum.rds"))
}else{
  mod.set2 <- readRDS(paste0(interfileFolder, "/GAMLSS.Flankerset2.sum.rds"))
}

# performance
modelperformance.set2 <- mod.set2$performance.tb
print(paste(sum(modelperformance.set2$converged), "models converged.")) 

# compute deviations
gam.data2 <- Flanker_data.set2
mod.tmp <- mod.set2$mod.tmp
mu_pred <- predict(mod.tmp, newdata = Flanker_data.set1, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = Flanker_data.set1, what = "sigma", type = "response")
nu_pred <- predict(mod.tmp, newdata = Flanker_data.set1, what = "nu", type = "response")

dependentvar <- "nihtbx_flanker_uncorrected"
deviation.set1.df <- data.frame(ID=Flanker_data.set1$ID)
observation <- Flanker_data.set1[[dependentvar]]
centile <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
deviation.set1.df[[paste0(dependentvar, "_centile")]] <- centile
deviation.set1.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(centile)
deviation.set1.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.set1.df[[paste0(dependentvar, "_mu")]] <- mu_pred
deviation.set1.df[[paste0(dependentvar, "_nu")]] <- nu_pred

#saveRDS(deviation.set1.df, paste0(interfileFolder, "/EF_Flanker.set1_deviation.rds"))

## 3rd. merge datasets
# The deviations will be averaged across the two models.
Flanker_deviation.set1 <- deviation.set1.df %>% select(c("ID", ends_with("_centile"), ends_with("_deviationZ"), ends_with("_sigma"), ends_with("_mu"), ends_with("_nu")))
Flanker_deviation.set2 <- deviation.set2.df %>% select(c("ID", ends_with("_centile"), ends_with("_deviationZ"), ends_with("_sigma"), ends_with("_mu"), ends_with("_nu")))

# Create an empty dataframe for averaged deviations
Flanker_deviation <- data.frame(ID = Flanker_deviation.set1$ID)

# Calculate the average centile and deviationZ for the variable
centilename <- "Flanker_centile"
deviationname <- "Flanker_deviationZ"

Flanker_deviation <- rbind(Flanker_deviation.set1, Flanker_deviation.set2) 

# Merge with original data if necessary
Flanker_deviation <- merge(Flanker_deviation, Flanker_data, by="ID")

# Save the result
write.csv(Flanker_data.set1, paste0(interfileFolder,"/Flanker_data.set1.csv"))
write.csv(Flanker_data.set2, paste0(interfileFolder,"/Flanker_data.set2.csv"))
saveRDS(Flanker_deviation, paste0(interfileFolder,"/Flanker.deviations.rds"))
write.csv(Flanker_deviation, paste0(interfileFolder,"/Flanker.deviations.csv"))
