library(tidyverse)
library(gamlss)
## Function 1. construct GAMLSS & return model objects and performance.
## This function is used to fit a GAMLSS containing a smooth term.
construct_gamlss <- function(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify=NULL,randomvar=NA){
  
  # get data
  gam.data <- get(dataname)
  covariates <- gsub(" ", "", covariates)
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(all_of(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=T)), randomvar, IDvar))) %>% drop_na()
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=T)), IDvar)) %>% drop_na()
  }
  
  con<-gamlss.control(n.cyc=200)
  gam.data2 <- as.data.frame(gam.data2)
  assign("gam.data2", gam.data2, envir = .GlobalEnv)
  # construct model
  if (! is.na(randomvar)){
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree=%s) + %s + random(%s)", dependentvar, smoothterm, mu.df, degree, covariates, randomvar))
    mod.mu.null.formula <- as.formula(sprintf("%s ~ %s + random(%s)", dependentvar, covariates, randomvar))
  }else{
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree=%s) + %s", dependentvar, smoothterm, mu.df, degree, covariates))
    mod.mu.null.formula <- as.formula(sprintf("%s ~ %s", dependentvar, covariates))
  }
  
  command <- paste0("mod.tmp <- gamlss(mod.mu.formula, sigma.formula =~ bs(", smoothterm, ", df = ", sigma.df, ", degree = ", degree, ") + ", 
                    covariates, 
                    ", nu.formula = ~1,family=", distribution.fam,", data=gam.data2, control=con)")
  
  command.null <- paste0("mod.null.tmp <- gamlss(mod.mu.null.formula, sigma.formula =~  ", covariates, 
                         ", nu.formula = ~1,family=", distribution.fam,", data=gam.data2, control=con)")
  
  eval(parse(text = command))
  eval(parse(text = command.null))
  
  # performance
  performance.tb <- data.frame(BIC = mod.tmp$sbc, converged = mod.tmp$converged, partialRsq = NA, Rsq = Rsq(mod.tmp))
  
  sse.model <- sum((mod.tmp$y - fitted(mod.tmp, what = "mu"))^2)
  sse.nullmodel <- sum((mod.null.tmp$y - fitted(mod.null.tmp, what = "mu"))^2)
  partialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
  
  # first derivative
  PEF <- getPEF(mod.tmp, term = smoothterm, n.points = 1000, parameter = "mu", type = "response", plot = FALSE)
  x_values <- seq(min(gam.data2[[smoothterm]]), max(gam.data2[[smoothterm]]), length.out = 1000)
  PEF_test <- PEF(x_values, deriv = 1)
  direction <- sum(PEF_test) / abs(sum(PEF_test))
  performance.tb$partialRsq <- partialRsq * direction
  
  
  # quantiles
  n_quantiles <- length(quantile.vec)
  n_points <- 1000
  x.tmp <- seq(min(gam.data2[[smoothterm]]), max(gam.data2[[smoothterm]]), length.out=n_points)
  centiles_strat <- list()
  for (l in 1:nlevels(gam.data2[[stratify]])){
    centile.tmp <- array(NA, dim=c(n_quantiles, n_points))
    for (q in 1:n_quantiles){
      command <- paste0("Qua <- getQuantile(mod.tmp, quantile=quantile.vec[q], term = smoothterm, fixed.at = list(", stratify, "='", levels(gam.data2[[stratify]])[l], "'), n.points = 1000)")
      eval(parse(text = command))
      centile.tmp[q, ] <- Qua(x.tmp)
    }
    centiles_strat[[l]] <- centile.tmp
  }
  
  
  sumlist <- list(performance.tb=performance.tb, mod.tmp=mod.tmp, centiles_strat=centiles_strat)

  return(sumlist)
}









