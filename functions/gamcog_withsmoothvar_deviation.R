# Referring to Sydnor et al. (2023)
library(tidyr)
library(mgcv)
library(psych)
library(dplyr)
library(ecostats)
library(pbkrtest)

#### Fit the linear effects of a continuous variable with a smooth item as covariates ####
## Function to fit evaluate the linear effects of a continuous variable on dependent variables per edge.
gam.fit.Independent.var <- function(dependentvar, dataname, smoothvar, interest.indep.var, covariates, stats_only = FALSE, cl = cl){  
  #Fit the gam
  gam.data <- get(dataname)
  gam.data <- gam.data %>%
    filter(!is.na(.data[[interest.indep.var]]), !is.na(.data[[dependentvar]]))

  Independent.var <- gam.data[[interest.indep.var]]
  indep_outlier_index <- which(Independent.var < mean(Independent.var) - 3 * sd(Independent.var) |
                                 Independent.var > mean(Independent.var) + 3 * sd(Independent.var))
  if (length(indep_outlier_index) > 0) {
    gam.data <- gam.data[-indep_outlier_index, ]
  }

  parcel <- as.character(dependentvar)
  modelformula <- as.formula(sprintf("%s ~ %s + %s + s(%s, k=3, fx=F)",dependentvar, interest.indep.var, covariates, smoothvar))
  modelformula.null<-as.formula(sprintf("%s ~ %s + s(%s, k=3, fx=F)",dependentvar, covariates, smoothvar))
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  
  #gam statistics
  #t value for the linear term (interest.indep.var)
  gam.indep.t <- gam.results$p.table[2,3]
  gam.indep.pvalue <- gam.results$p.table[2,4]
  
  # Full versus reduced model anova p-value
  if (stats_only){
    # Perform ANOVA simulation
    anova_results <- anovaPB(gam.model.null, gam.model, n.sim = 10000, test = 'Chisq', cl = cl,ncpus = 1)# number of parametric bootstrap simulations
    anova.pvalues <- anova_results$`Pr(>Chi)`[2]  # Get p-values
    simulated_chisq <- anova_results$sim.chisq
    
    # Calculate real (observed) chi-squared value using logLik
    real_chisq <- 2 * (as.numeric(logLik(gam.model)) - as.numeric(logLik(gam.model.null)))
    df_diff <- attr(logLik(gam.model), "df") - attr(logLik(gam.model.null), "df")
    #Save simulated chi-squared distribution
    chisq_save_dir <- file.path(resultFolder, "simulated_chisq")
    if (!dir.exists(chisq_save_dir)) dir.create(chisq_save_dir, recursive = TRUE)
    saveRDS(simulated_chisq, file = file.path(chisq_save_dir, paste0("simulated_chisq_", interest.indep.var, "_", dependentvar, ".rds")))
    
  }else{
    anova.pvalues <- NA
    real_chisq <- NA
    df_diff <- NA
    simulated_chisq <- NA
  }
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  if(gam.indep.t < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~  %s+s(%s, k=3, fx=F)", dependentvar, covariates, smoothvar))
  res1<-residuals(gam(varcorformula1,method="REML", data = gam.data))
  res2<-gam.data[[interest.indep.var]] # Deviation Z scores already adjusted for age/sex; no further residualization
  res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
  res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
  
  if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
    corrmethod <- "pearson"
  }else{
    corrmethod <- "spearman"
  }
  
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-as.numeric(PCorr_Test$r)
  corrp <- as.numeric(PCorr_Test$p)
  
  samplesize <- nrow(gam.data)
  beta = gam.results$p.table[2,1]
  stats.results <- cbind(parcel, interest.indep.var, gam.indep.t, gam.indep.pvalue, anova.pvalues,real_chisq, df_diff, partialRsq, corrmethod, correstimate, corrp, samplesize,beta)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
}
