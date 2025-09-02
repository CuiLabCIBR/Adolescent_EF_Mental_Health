# Referring to Sydnor et al. (2023)
library(tidyr)
library(mgcv)
library(psych)
library(dplyr)
library(ecostats)
gam.varyingcoefficients <- function(dependentvar, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, increments, draws, return_posterior_coefficients = FALSE){
  set.seed(123)
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  gam.data <- get(dataname)
  smoothmin <- min(gam.data[ , smooth_var])
  smoothmax <- max(gam.data[ , smooth_var])
  
  NonNANIndex <- which(
    !is.na(gam.data[[int_var]]) & 
      !is.na(gam.data[[dependentvar]])
  )
  gam.data <- gam.data[NonNANIndex, ]
  cognition <- as.numeric(gam.data[[int_var]])
  outlierindx1 <- which(
    cognition < mean(cognition, na.rm = TRUE) - 3 * sd(cognition, na.rm = TRUE) |
      cognition > mean(cognition, na.rm = TRUE) + 3 * sd(cognition, na.rm = TRUE)
  )
  if (length(outlierindx1) > 0) {
    gam.data <- gam.data[-outlierindx1, ]
  }
  
  parcel <- dependentvar
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ %6$s", dependentvar, smooth_var, knots, set_fx, int_var, covariates))
  gam.model <- gam(modelformula, method ="REML", data = gam.data)
  
  gam.summary <- summary(gam.model)
  modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) +%5$s + %6$s", dependentvar, smooth_var, knots, set_fx, int_var, covariates))
  gam.model.null <- gam(modelformula.null, method = "REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  # Full versus reduced model anova p-value
  set.seed(123)
  anova.pvalues <- anovaPB(gam.model.null, gam.model, n.sim = 10000, test = 'Chisq')$`Pr(>Chi)`[2]
  
  # Calculate partial R² as effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)  # SSE for full model
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)  # SSE for null model
  IntpartialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel  # Partial R²
  
  # Calculate Z-value for ANOVA p-value
  anova.int.zvalue <- qnorm(anova.pvalues / 2, lower.tail = FALSE)
  full.results<-cbind(int_var, parcel, anova.pvalues, IntpartialRsq, anova.int.pvalue, anova.int.zvalue)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate (posterior) model slopes (varying covariate coefficients)
  theseVars <- attr(gam.model$terms,"term.labels") 
  varClasses <- attr(gam.model$terms,"dataClasses") 
  
  #prediction df for int_var P10th at np smooth_var increments
  pred.low <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == int_var) { 
      #pred.low[,int_var] <- (min(df[,int_var],na.rm = T))
      pred.low[,int_var] <- (quantile(df[,int_var], probs=0.1, na.rm = T))
    } else if (thisVar == smooth_var) {
      pred.low[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.low[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {pred.low[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {pred.low[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )}}
  
  #prediction df for int_var P90th at np smooth_var increments
  pred.high <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == int_var) { 
      #pred.high[,int_var] <- (max(df[,int_var],na.rm = T))
      pred.high[,int_var] <- (quantile(df[,int_var], probs=0.9, na.rm = T))
    } else if (thisVar == smooth_var) {
      pred.high[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.high[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {pred.high[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {pred.high[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )}}
  
  pred <- rbind(pred.low, pred.high) #complete pred df 
  pred <- pred %>% dplyr::select(-init)
  
  varyingcoeff.slopes <- NULL
  varyingcoeff.slopes.CI <- NULL
  #Get effects (slopes) along the smooth function for the true model
  if(return_posterior_coefficients == FALSE){
    #varying coefficient slopes
    predicted.values <- fitted_values(object = gam.model, data = pred) #predict y at min and mix int_var along the smooth function
    predicted.values <- predicted.values %>% dplyr::select(.fitted, all_of(smooth_var), all_of(int_var))
    colnames(predicted.values) <- c(".fitted", "smooth.var", "int.var")
    predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
    
    varyingcoeff.slopes <-  predicted.values %>% #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function 
      group_by(smooth.var) %>%
      do(slope = diff(.$.fitted)/diff(.$int.var)) %>%
      unnest(cols = c(slope))
    colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")  
  }
  data.results1 <- list(full.results, varyingcoeff.slopes, pred)
  
  #Estimate posterior distribution of effects (slopes) from simulated GAM posterior distribution
  if(return_posterior_coefficients == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix of fitted gam coefficients
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulated model coefficients for npd draws from the posterior
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") # Compute the linear predictor matrix for the prediction dataset
    predicted.values.posterior <- X0 %*% t(sims) #Generate posterior predicted values (Y_hat) for each draw
    
    predicted.values.posterior <- as.data.frame(predicted.values.posterior)
    colnames(predicted.values.posterior) <- sprintf("draw%s",seq(from = 1, to = npd))
    predicted.values.posterior <- cbind(as.numeric(pred[,smooth_var]), as.numeric(pred[,int_var]), predicted.values.posterior)
    colnames(predicted.values.posterior)[1] <- c("smooth.var")
    colnames(predicted.values.posterior)[2] <- c("int.var")
    predicted.values.posterior$smooth.var <- round(predicted.values.posterior$smooth.var, 3)
    #Pivot wider to get predicted values for EF P90 and P10, then compute slope = ΔY / ΔEF
    varyingcoeff.slopes.CI = predicted.values.posterior %>% pivot_longer(cols = contains("draw"), names_to = "draw",values_to = "posterior.fitted") 
    varyingcoeff.slopes.CI$int.var[varyingcoeff.slopes.CI$int.var == quantile(df[,int_var], 0.1)] <- c("low")
    varyingcoeff.slopes.CI$int.var[varyingcoeff.slopes.CI$int.var == quantile(df[,int_var], 0.9)] <- c("high")
    varyingcoeff.slopes.CI <- varyingcoeff.slopes.CI %>% pivot_wider(names_from = "int.var", values_from = "posterior.fitted") %>% 
      mutate(slope = (high-low)/(quantile(df[,int_var], 0.9) - quantile(df[,int_var], 0.1))) 
    #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function for all draws
    
    varyingcoeff.slopes.CI <- varyingcoeff.slopes.CI %>% dplyr::select(smooth.var, draw, slope)
    varyingcoeff.slopes.CI <- cbind(as.character(parcel), varyingcoeff.slopes.CI)
    colnames(varyingcoeff.slopes.CI) <- c("label", smooth_var, "draw", "slope")
  }
  
  data.results2 <- list(full.results, varyingcoeff.slopes.CI, pred)
  if(return_posterior_coefficients == FALSE)
    return(data.results1)
  if(return_posterior_coefficients == TRUE)
    return(data.results2)
  
}
