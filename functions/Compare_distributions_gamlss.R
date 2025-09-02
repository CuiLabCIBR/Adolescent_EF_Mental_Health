library(tidyverse)
library(gamlss)
library(scales)
library(parallel)

### Function 1
# This function is used to compare GAMLSS with different distribution.
# Referring to Bethlehem et. al., Nature, 2022 & Sun et. al., bioArxiv, 2023,
# we firstly compared the fitting performance of different distributions to select the best one.
# Only continuous distributions with >= 3 parameters were compared.

gamlss_comparedistribution <- function(dataname, dependentvar, smoothvar,IDvar, bs.df, covariates, randomvar=NA){
  
  # get data
  gam.data <- get(dataname)
  gam.data$age <- gam.data[[smoothvar]]
  
  covariates <- gsub(" ", "", covariates)
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)), randomvar, IDvar)) %>% drop_na()
    gam.data1$dependent <- gam.data1[[dependentvar]]
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    gam.data1 <- gam.data1 %>% filter((dependent>mean(dependent)-3*sd(dependent)) & (dependent<mean(dependent)+3*sd(dependent)))
    gam.data1$dependent <- NULL
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s) + %s + random(as.factor(%s))", dependentvar, "age", bs.df, covariates, "randomvar"))
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)), IDvar)) %>% drop_na()
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s) + %s", dependentvar, "age", bs.df, covariates))
  }
  
  # fit models. 26 distributions were adopted.
  # R names	No of parameters	continuous	y_transform	model_ID
  # BEINF()	4	1	rescale to 0.1~0.9	mod1
  # GG()	3	1	rescale to 0.1~0.9	mod10
  # GIG()	3	1	rescale to 0.1~0.9	mod11
  # GT()	4	1	rescale to 0.1~0.9	mod12
  # PE()	3	1	rescale to 0.1~0.9	mod13
  # PE2()	3	1	rescale to 0.1~0.9	mod14
  # SEP1()	4	1	rescale to 0.1~0.9	mod15
  # SEP2()	4	1	rescale to 0.1~0.9	mod16
  # SEP3()	4	1	rescale to 0.1~0.9	mod17
  # SEP4()	4	1	rescale to 0.1~0.9	mod18
  # SHASH()	4	1	rescale to 0.1~0.9	mod19
  # JSU()	4	1	rescale to 0.1~0.9	mod2
  # SHASHo()	4	1	rescale to 0.1~0.9	mod20
  # ST1()	3	1	rescale to 0.1~0.9	mod21
  # ST2()	3	1	rescale to 0.1~0.9	mod22
  # ST3()	3	1	rescale to 0.1~0.9	mod23
  # ST4()	3	1	rescale to 0.1~0.9	mod24
  # ST5()	3	1	rescale to 0.1~0.9	mod25
  # TF()	3	1	rescale to 0.1~0.9	mod26
  # BCCG()	3	1	rescale to 0.1~0.9	mod3
  # BCPE()	4	1	rescale to 0.1~0.9	mod4
  # BCT()	4	1	rescale to 0.1~0.9	mod5
  # exGAUS()	3	1	rescale to 0.1~0.9	mod6
  # EGB2()	4	1	rescale to 0.1~0.9	mod7
  # GB1()	4	1	rescale to 0.1~0.9	mod8
  # GB2()	4	1	rescale to 0.1~0.9	mod9
  
  familylist <- c("JSU", "BCCG", "BCPE", "BCT", "exGAUS", "EGB2", 
                  "GB2", "GG", "GIG", "GT", "PE", "PE2", "SEP1", "SEP2", 
                  "SEP3", "SEP4", "SHASH", "SHASHo", "ST1", "ST2", "ST3", 
                  "ST4", "ST5", "TF")
  con<-gamlss.control(n.cyc=200)
  
  mod.sum <- mclapply(1:length(familylist), function(i){
    command <- paste0("mod.tmp <- gamlss(mod.mu.formula, sigma.formula =~ bs(age) + ", 
                      covariates, 
                      ", nu.formula = ~1,family=", familylist[i],", data=gam.data2, control=con)")
    
    eval(parse(text = command))
    mod.tmp$ID <- gam.data2[[IDvar]]
    return(mod.tmp)
  }, mc.cores = 50)
  
  
  performance <- data.frame(matrix(NA, length(familylist), 4))
  names(performance) <- c("model_ID", "distribution", "converged", "BIC")
    
  for (i in 1:length(mod.sum)){
    if (length(mod.sum[[i]])>30){
      performance$model_ID[i] <- i
      performance$distribution[i] <- mod.sum[[i]]$family[[1]]
      performance$converged[i] <- mod.sum[[i]]$converged
      performance$BIC[i] <- mod.sum[[i]]$sbc
    }
  }
  familylist <- familylist[performance$converged==T]
  print(paste0("The best distribution for ", dependentvar, " is ", familylist[which.min(performance$BIC[performance$converged==T])], "."))
  
  performance$dependentvar <- dependentvar
  
  outlist <- list(mod.sum, performance=performance)
  
  return(outlist)
  
}


### Function 2
# This function is used to compare GAMLSS with different bs degree of freedom.
# Referring to Sun et. al., bioArxiv, 2023, df = 3~6 for mu & sigma were tested.
# bs.df.mat should be a n*2 matrix with the first column is the df of mu. 
# the second column is df for sigma is the df of sigma.
# the third column is the spline degree.

gamlss_compare.bs.df <- function(dataname, dependentvar, smoothvar, IDvar, bs.df.set, covariates, distribution.fam, randomvar=NA){
  
  # get data
  gam.data <- get(dataname)
  gam.data$age <- gam.data[[smoothvar]]
  covariates <- gsub(" ", "", covariates)
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)), IDvar, randomvar)) %>% drop_na()
    gam.data1$dependent <- gam.data1[[dependentvar]]
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    gam.data1 <- gam.data1 %>% filter((dependent>mean(dependent)-3*sd(dependent)) & (dependent<mean(dependent)+3*sd(dependent)))
    gam.data1$dependent <- NULL
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)),IDvar)) %>% drop_na()
  }
  
  con<-gamlss.control(n.cyc=200)
  mod.sum <- mclapply(1:nrow(bs.df.set), function(i){
    # df for mu // df for sigma
    mu.df <- bs.df.set[i, 1]
    sigma.df <- bs.df.set[i, 2]
    degree <- bs.df.set[i, 3]
    
    if (! is.na(randomvar)){
      mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree = %s) + %s + random(as.factor(%s))", dependentvar, "age", mu.df, degree, covariates, "randomvar"))
    }else{
      mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree = %s) + %s", dependentvar, "age", mu.df, degree, covariates))
    }
    # Dynamically construct and evaluate the gamlss call
    # Note: nu.formula = ~1 means the nu parameter is constant
    # The sigma.formula string is constructed directly here
    command <- paste0("mod.tmp <- gamlss(mod.mu.formula, sigma.formula =~ bs(age, df = ", sigma.df, ", degree = ", degree,") + ", 
                      covariates, 
                      ", nu.formula = ~1,family=", distribution.fam,", data=gam.data2, control=con)")
    
    eval(parse(text = command))
    mod.tmp$ID <- gam.data2[[IDvar]]
    return(mod.tmp)
  }, mc.cores = 50)
  
  performance <- data.frame(matrix(NA, nrow(bs.df.set), 7))
  names(performance) <- c("model_ID", "distribution", "converged", "BIC", "mu.df", "sigma.df", "degree")
  
  for (i in 1:nrow(bs.df.set)){
    performance$model_ID[i] <- i
    performance$distribution[i] <- mod.sum[[i]]$family[[1]]
    performance$converged[i] <- mod.sum[[i]]$converged
    performance$BIC[i] <- mod.sum[[i]]$sbc
    performance$mu.df[i] <- bs.df.set[i,1]
    performance$sigma.df[i] <- bs.df.set[i,2]
    performance$degree[i] <- bs.df.set[i,3]
  }
  
  bs.df.set <- bs.df.set[performance$converged==T, ]
  print(paste0("The best distribution for ", dependentvar, " is mu.df = ", bs.df.set[which.min(performance$BIC[performance$converged==T]), 1], 
               "sigma.df = ", bs.df.set[which.min(performance$BIC[performance$converged==T]), 2], ", degree = ", bs.df.set[which.min(performance$BIC[performance$converged==T]), 3], "."))
  
  performance$dependentvar <- dependentvar
  outlist <- list(mod.sum, performance=performance)
  
  return(outlist)

}


### Function 3: Boot.Function
# This function performs bootstrap resampling and fits GAMLSS models
# to calculate the first derivatives (Partial Effect Functions) of mu and sigma parameters.
Boot.Function <- function( n, Base.Seed, dataname,smoothvar, model_obj, stratify=NULL,randomvar=NA) {
  # extract variables
  mod.mu.formula <- model_obj$mu.formula
  formula.vars <- as.character(mod.mu.formula)
  dependentvar <- formula.vars[2]
  
  covariates <- str_split(formula.vars[3], "\\+")
  if (! is.na(randomvar)){
    covariates <- covariates[[1]][2:(length(covariates[[1]])-1)]
  }else{
    covariates <- covariates[[1]][2:length(covariates[[1]])]
  }
  covariates <- gsub(" ", "", covariates)
  # get data
  gam.data <- get(dataname)
  gam.data$age <- gam.data[[smoothvar]]
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(unique(c(dependentvar, smoothvar, unlist(strsplit(covariates, "+", fixed=T)), randomvar, stratify))) %>% drop_na()
    gam.data1$dependent <- gam.data1[[dependentvar]]
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    gam.data1 <- gam.data1 %>% filter((dependent>mean(dependent)-3*sd(dependent)) & (dependent<mean(dependent)+3*sd(dependent)))
    gam.data1$dependent <- NULL
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
    
    gam.data2[[randomvar]] <- droplevels(gam.data2[[randomvar]])
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, "age", covariates)) %>% drop_na()
  }
  
  # gam.data2 is a data frame and assign to global environment
  gam.data2 <- as.data.frame(gam.data2)
  assign("gam.data2", gam.data2, envir = .GlobalEnv)
  
  set.seed( seed=Base.Seed + n )
  if (length(stratify) == 1) {
    stratify.var <- gam.data2[[stratify]]
  } else{
    stratify.var <- interaction(gam.data2[ ,stratify])
  }
  
  
  if (length(stratify)==0){
    INDEX <- sample(1:NROW(gam.data2),NROW(gam.data2),replace=TRUE) ## unstratified bootstrap
  } else {
    SPLIT <- split(1:NROW(gam.data2), stratify.var)
    LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X),replace=TRUE)})
    INDEX <- unsplit(LAPPLY, stratify.var)
  }
  
  gam.data.subset <- gam.data2[INDEX, ] ## generate SUBSET with bootstrapped-SUBSET
  
  # fit gamlss
  sigma.formula <- as.character(model_obj$sigma.formula)[2]
  mod.mu.formula <- paste(formula.vars[2], "~", formula.vars[3])
  distribution.fam <- model_obj$family[1]
  con <- gamlss.control(n.cyc = 200)
  command <- paste0("mod.tmp <- gamlss(", mod.mu.formula, ", sigma.formula =~", sigma.formula,
                    ", nu.formula = ~1,family=", distribution.fam,", data=gam.data.subset, control=con)")
  
  eval(parse(text = command))
  
  assign("gam.data.subset", gam.data.subset, envir = .GlobalEnv)
  # Calculate Partial Effect Functions (PEF) for mu and sigma
  # Note: getPEF is assumed to be a custom function, not from the GAMLSS package
  # And its logic is assumed to run correctly in the current environment
  # Calculate partial effect functions for mu and sigma derivatives over range of smoothvar
  x_values <- seq(min(gam.data.subset$age), max(gam.data.subset$age), length.out = 1000)
  mu.pef <- getPEF(obj = mod.tmp, term = "age", data = gam.data.subset, parameter = "mu", plot = FALSE)
  sigma.pef <- getPEF(obj = mod.tmp, term = "age", data = gam.data.subset, parameter = "sigma", plot = FALSE)
  
  # Calculate mu and sigma first derivatives over x_values
  mu_derivative_values <- sapply(x_values, function(value) mu.pef(value, deriv = 1))
  sigma_derivative_values <- sapply(x_values, function(value) sigma.pef(value, deriv = 1))
  
  # Store results
  bootstrap_list <- list(
    gam.data.subset = gam.data.subset,
    mod.tmp = mod.tmp,
    mu_derivative = mu_derivative_values,
    sigma_derivative = sigma_derivative_values,
    x_values = x_values
  )
  
  return(bootstrap_list)
}







