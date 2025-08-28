library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)
library(psych)
library(cAIC4)

pbootint <- function(modelobj, int_var=NA){
  numsims <- 1000
  set.seed(123)
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  if (!is.na(int_var)){
    indexcomma <- gregexpr(",", theseVars[1])[[1]]
    addsmooth <- theseVars[1]
    addsmooth <- paste0(substr(addsmooth, 1, indexcomma[1]), substr(addsmooth, indexcomma[2]+1, nchar(addsmooth)))
    theseVars <-  c(theseVars, addsmooth)
    }
  if (sum(str_detect(theseVars, "by ="))==1){
    int_var <- str_split_i(theseVars[1], "by = ", 2)
    int_var <- str_split_i(int_var, ", ", 1)
    f2 <- reformulate(c(theseVars[2:(length(theseVars))], int_var),response = thisResp)
  }else{
    f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  }
    
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  src_subject_id<-df$src_subject_id
  y <- df[,thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1|src_subject_id))
  m2 <- lmer(y ~ -1 + mat2 + (1|src_subject_id))
  refdist <- PBrefdist(m1, m2, nsim=numsims)
  pb <- PBmodcomp(m1, m2, ref = refdist)
  int_pval <- pb$test["PBtest","p.value"]
  
  return(int_pval)
}


gamm.smooth.predict.interaction <- function(dependentvar, dataname, smoothvar, interest.indep.var, covariates){
  set.seed(123)
  gam.data <- get(dataname)
  gam.data <- gam.data %>%
    filter(!is.na(.data[[dependentvar]]), !is.na(.data[[interest.indep.var]]))
  
  parcel <- as.character(dependentvar)
  # Process outliers in the independent variable
  Independent.var <- gam.data[[interest.indep.var]]
  outlierindx1 <- which(Independent.var < mean(Independent.var) - 3 * sd(Independent.var) |
                          Independent.var > mean(Independent.var) + 3 * sd(Independent.var))
  if (length(outlierindx1) > 0) {
    gam.data <- gam.data[-outlierindx1, ]
  }
  
  
  #Fit the gam
  if (is.na(covariates) & is.na(smoothvar)){
    modelformula <- as.formula(sprintf("%s ~ %s", dependentvar, interest.indep.var))
    modelformula.null <- as.formula(sprintf("%s ~ 1", dependentvar))
  }else{
    modelformula.nonlinear <- as.formula(sprintf("%s ~ %s + %s + s(%s, k=3, fx=T)", dependentvar, interest.indep.var, covariates, smoothvar))
    modelformula.nonlinear.null <- as.formula(sprintf("%s ~ %s + s(%s, k=3, fx=T)", dependentvar, covariates, smoothvar))
    
    modelformula.linear <- as.formula(sprintf("%s ~ %s + %s + %s + (1|src_subject_id)", dependentvar, interest.indep.var, covariates, smoothvar))
    modelformula.linear.null <- as.formula(sprintf("%s ~ %s + %s + (1|src_subject_id)", dependentvar, covariates, smoothvar))
  }
  
  gamm.model <- gamm4(modelformula.nonlinear, random=~(1|src_subject_id), REML=T, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  lm.model <- lmer(modelformula.linear, REML = T, data=gam.data)
  lm.results <- summary(lm.model)
  
  gamm.model.null <- gamm4(modelformula.nonlinear.null, random=~(1|src_subject_id), REML=T, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  lm.model.null <- lmer(modelformula.linear.null, REML = T, data=gam.data)
  gamm_BIC <- BIC(gamm.model$mer)
  #gamm_AIC <- AIC(gamm.model$mer)
  #lm_AIC   <- AIC(lm.model)
  lm_BIC <- BIC(lm.model)
  if (gamm_BIC < lm_BIC){
    bettermodel <- "GAMM"
  }else{
    bettermodel <- "LM"
  }
  
  # gamm_cAIC <- cAIC(gamm.model$mer)$caic
  # lm_cAIC   <- cAIC(lm.model)$caic
  # if (gamm_cAIC < lm_cAIC){
  #   bettermodel <- "GAMM"
  # }else{
  #   bettermodel <- "LM"
  # }
  
  
  if (bettermodel == "GAMM"){
    # Extract F-value and p-value for the independent variable from the GAMM results
    gamm.independent.t <- gamm.results$p.table[2,3]
    gamm.independent.pvalue <- gamm.results$p.table[2,4]
    slope <- gamm.results$p.table[2,1]  
    boots.pvalues <- pbootint(gamm.model)
    # interaction effect size
    sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
    sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
    partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
    
    # Calculate residuals and check for normality for correlation analysis
    varcorformula1 <- as.formula(sprintf("%s ~  %s+ s(%s, k=3, fx=T)", dependentvar, covariates, smoothvar))
    varcorformula2 <- as.formula(sprintf("%s ~  %s+ s(%s, k=3, fx=T)", interest.indep.var, covariates, smoothvar))
    res1<-residuals(gamm4(varcorformula1,random=~(1|src_subject_id), REML=T, data = gam.data)$gam)
    res2<-gam.data[[interest.indep.var]]
    # Determine correlation method based on normality of residuals
    res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
    res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
  }else{
    gamm.independent.t <- lm.results$coefficients[2,3]
    gamm.independent.pvalue <- NA
    slope <- lm.results$coefficients[2,1]
    set.seed(123)
    boots.pvalues <- PBmodcomp(lm.model, lm.model.null, nsim=1000)$test["PBtest","p.value"]
    # interaction effect size
    sse.model <- sum(residuals(lm.model)^2)
    sse.nullmodel <- sum(residuals(lm.model.null)^2)
    sse.total <- sum((gam.data[[dependentvar]]-mean(gam.data[[dependentvar]]))^2)
    partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
    
    # Calculate residuals and check for normality for correlation analysis
    varcorformula1 <- as.formula(sprintf("%s ~  %s + %s + (1|src_subject_id)", dependentvar, covariates, smoothvar))
    varcorformula2 <- as.formula(sprintf("%s ~  %s + %s + (1|src_subject_id)", interest.indep.var, covariates, smoothvar))
    res1<-residuals(lmer(varcorformula1,REML=T, data = gam.data))
    res2<-gam.data[[interest.indep.var]]
    # Determine correlation method based on normality of residuals
    res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
    res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
    
  }
  
  if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
    corrmethod <- "pearson"
  }else{
    corrmethod <- "spearman"
  }
  # Compute correlation between residuals
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-as.numeric(PCorr_Test$r)
  corrp <- as.numeric(PCorr_Test$p)
  # Return relevant statistics
  stats.results <- data.frame(
    parcel = dependentvar,
    interest.indep.var = interest.indep.var,
    partialRsq = partialRsq,
    bettermodel=bettermodel,
    gamm.independent.t = gamm.independent.t,
    gamm.independent.pvalue = gamm.independent.pvalue,
    boots.pvalues = boots.pvalues,
    correstimate = correstimate,
    corrp = corrp,
    slope = slope
  )
  
  
  return(stats.results)
}
