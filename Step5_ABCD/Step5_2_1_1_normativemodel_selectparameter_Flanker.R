rm(list=ls())
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(readxl)
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

# source functions
source(paste0(functionFolder, "/Compare_distributions_gamlss.R"))

Flanker_data <- read_xlsx(paste0(datapath, "/ABCD_Flanker.xlsx"))
##### Step1. select the best distribution
dataname <- "Flanker_data"
Flanker_data$Sex <- as.factor(Flanker_data$Sex)
dependentvar <- "nihtbx_flanker_uncorrected"
smoothvar <- "Age_year"
IDvar <- "ID"
bs.df <- 3
covariates <- "Sex"
randomvar=NA

out_Flanker_distribution <- gamlss_comparedistribution(dataname, dependentvar,
                                                     smoothvar, IDvar, bs.df,
                                                     covariates, randomvar=NA)
names(out_Flanker_distribution) <- c("modelsum", "performance")
performance_Flanker_distribution <- out_Flanker_distribution$performance
modelsum_Flanker_distribution <- out_Flanker_distribution$modelsum
write.csv(performance_Flanker_distribution, paste0(resultFolder, "/performance_Flanker_distribution.csv"), row.names = F)
saveRDS(modelsum_Flanker_distribution, paste0(resultFolder, "/modelsum_Flanker_distribution.rds"))

