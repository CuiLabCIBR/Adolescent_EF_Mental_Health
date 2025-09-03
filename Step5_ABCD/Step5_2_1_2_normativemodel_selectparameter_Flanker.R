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

##### 1.Flanker
##### Step2. select proper model details (degree of freedom -- df)
con<-gamlss.control(n.cyc=200)
df.set <- matrix(c(3,3,3,
                   3,4,3,
                   3,5,3,
                   3,6,3,
                   4,3,3,
                   4,4,3,
                   4,5,3,
                   4,6,3,
                   5,3,3,
                   5,4,3,
                   5,5,3,
                   5,6,3,
                   6,3,3,
                   6,4,3,
                   6,5,3,
                   6,6,3,
                   2,2,2,
                   2,3,2,
                   2,4,2,
                   2,5,2,
                   2,6,2,
                   3,2,2,
                   3,3,2,
                   3,4,2,
                   3,5,2,
                   3,6,2,
                   4,2,2,
                   4,3,2,
                   4,4,2,
                   4,5,2,
                   4,6,2,
                   5,3,2,
                   5,4,2,
                   5,5,2,
                   5,6,2,
                   6,3,2,
                   6,4,2,
                   6,5,2,
                   6,6,2),
                 byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","degree")))
dataname <- "Flanker_data"
Flanker_data$Sex <- as.factor(Flanker_data$Sex)
Flanker_data$site_id_l <- as.factor(Flanker_data$site_id_l)
dependentvar <- "nihtbx_flanker_uncorrected"
smoothvar <- "Age_year"
IDvar <- "ID"
bs.df <- 3
saveout_dir <- interfileFolder
bs.df.set <- df.set
covariates <- "Sex"
distribution.fam <- "GG"
out_Flanker_bsdf <- gamlss_compare.bs.df(dataname, dependentvar, smoothvar, IDvar,
                                       bs.df.set, covariates,
                                       distribution.fam, randomvar=NA)
names(out_Flanker_bsdf) <- c("modelsum", "performance")
performance_Flanker_bsdf <- out_Flanker_bsdf$performance
modelsum_Flanker_bsdf <- out_Flanker_bsdf$modelsum
write.csv(performance_Flanker_bsdf, paste0(resultFolder, "/performance_Flanker_bsdf.csv"), row.names = F)
saveRDS(modelsum_Flanker_bsdf, paste0(resultFolder, "/modelsum_Flanker_bsdf.rds"))


