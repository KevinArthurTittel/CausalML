# IMPORT PACKAGES NECESSARY FOR CAUSAL FOREST
install.packages("grf")
library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
install.packages("standardize")
library(standardize)

# IMPORT STATA DATA FILE
library(foreign)
library(haven)
Grace_Period_Data <- read_dta("Downloads/112672-V1/Grace-Period-Data.dta")
Grace_Period_Data <- as.data.frame(Grace_Period_Data)

# Results matrix
resultsTable5OriginalPaper <- matrix(data = 0, nrow = 5, ncol = 6)
resultsTable5OriginalPaper[1,] <- t(c("dependent variable", "savings", "risk loving", "wage earner", "household member chronically ill",
                   "impatient"))
resultsTable5OriginalPaper[2:5,1] <- c("CF", "Cluster-robust CF", "LLCF", "observations")

for (i in 1:6) {
  # Appoint treatment assignment and outcome variables
  W <- Grace_Period_Data$sec_treat
  W <- as.vector(W)
  
  # Standardize the continuous variables for improved performance of RF
  Grace_Period_Data$Years_Education_C <- scale(Grace_Period_Data$Years_Education_C)
  Grace_Period_Data$Age_C <- scale(Grace_Period_Data$Age_C)
  Grace_Period_Data$SEI <- scale(Grace_Period_Data$SEI)
  Grace_Period_Data$HH_Size_C <- scale(Grace_Period_Data$HH_Size_C)
  
  # Appoint the control variables matrix
  X <- Grace_Period_Data[,8:18]
  
  # Matrix with characteristics for Table 5
  characteristics <- Grace_Period_Data[,(74+i)]
  
  # Create a numeric vector of the character group name vector
  loangroups <- as.numeric(factor(Grace_Period_Data$sec_group_name))
  
  # Create loan group dummies (as in original analysis) to be added to control variables matrix
  loansize1 <- as.integer(c(Grace_Period_Data$sec_loanamount == 4000))
  loansize2 <- as.integer(c(Grace_Period_Data$sec_loanamount == 5000))
  loansize3 <- as.integer(c(Grace_Period_Data$sec_loanamount == 6000))
  loansize4 <- as.integer(c(Grace_Period_Data$sec_loanamount == 8000))
  loansize5 <- as.integer(c(Grace_Period_Data$sec_loanamount == 9000))
  loansize6 <- as.integer(c(Grace_Period_Data$sec_loanamount == 10000))
  loansizematrix <- cbind(loansize1, loansize2, loansize3, loansize4,
                          loansize5, loansize6)
  
  X <- as.matrix(cbind(X, loansizematrix))
  Y <- Grace_Period_Data[,(50)]
  Y <- as.vector(Y)
  
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristics))
  missingvalues[missingvalues == 2] <- 1
  
  resultsTable5OriginalPaper[5,(i+1)] <- (845 - sum(missingvalues))
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  X <- X[!missingvalues,]
  W <- W[!missingvalues]
  loangroups <- loangroups[!missingvalues]
  characteristics <- characteristics[!missingvalues]
  Y <- Y[!missingvalues]
  
  # TRAIN FORESTS FOR W AND Y SEPARATELY
  forest.W <- regression_forest(X, W)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- regression_forest(X, Y)
  Y.hat <- predict(forest.Y)$predictions
  
  # Orthogonalization for the orthogonalized CF
  # orthog.W <- W - W.hat 
  # orthog.Y <- Y - Y.hat
  
  # Implement preliminary standard Causal Forest to select important covariates
  prelim.CF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, 
                             num.trees = 4000)
  prelim.CF.varimp <- variable_importance(prelim.CF)
  selected.vars <- which(prelim.CF.varimp / mean(prelim.CF.varimp) > 0.2)
  
  # Implement final version of CF with selected covariates
  CF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat,
                      num.trees = 8000, tune.parameters = "all")
  tau.hat.CF <- predict(CF)$predictions
  hist(tau.hat.CF)

  # TEST FOR HETEROGENEITY along set of characteristics
  ate.subgroup <- average_treatment_effect(CF, target.sample = "all", subset = (characteristics == 1))
  ate.rest <- average_treatment_effect(CF, target.sample = "all", subset = !(characteristics == 1))
  resultsTable5OriginalPaper[2,(i+1)] <- paste(round(ate.subgroup[1] - ate.rest[1], 3), "(", sqrt(ate.subgroup[2]^2 + ate.rest[2]^2), ")")
  
  # TEST FOR HETEROGENEITY: Approach 2: Best Linear Predictor Analysis
  # results[(i+1),5] <- test_calibration(CF)[1]
  # results[(i+1),6] <- test_calibration(CF)[2]
  
  # Implement preliminary cluster robust Causal Forest to select important covariates
  prelim.CF.CR <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, 
                                num.trees = 4000)
  prelim.CF.CR.varimp <- variable_importance(prelim.CF.CR)
  selected.vars <- which(prelim.CF.CR.varimp / mean(prelim.CF.CR.varimp) > 0.2)
  
  # Implement final version of cluster-robust CF with selected covariates
  CF.CR <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat,
                         clusters = loangroups, 
                         num.trees = 8000, tune.parameters = "all")
  tau.hat.CF.CR <- predict(CF.CR)$predictions
  
  # UNCOVER HETEROGENEITY IN THESE VARIABLES
  
  # TEST FOR HETEROGENEITY: Approach 1: Differential ATE
  ate.subgroup <- average_treatment_effect(CF.CR, target.sample = "all", subset = (characteristics == 1))
  ate.rest <- average_treatment_effect(CF.CR, target.sample = "all", subset = !(characteristics == 1))
  resultsTable5OriginalPaper[3,(i+1)] <- paste(round(ate.subgroup[1] - ate.rest[1], 3), "(", sqrt(ate.subgroup[2]^2 + ate.rest[2]^2), ")")
  
  # TEST FOR HETEROGENEITY: Approach 2: Best Linear Predictor Analysis
  # results[(i+1),8] <- test_calibration(CF.CR)[1]
  # results[(i+1),9] <- test_calibration(CF.CR)[2]
  
  # RUN LLCF
  # Firstly, run the LLF on both (X,W) and (X,Y)
  forest.W <- ll_regression_forest(X, W, honesty = TRUE)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- ll_regression_forest(X, Y, honesty = TRUE)
  Y.hat <- predict(forest.Y)$predictions
  
  # Choose which variables to use in LLCF algorithm, and obtain estimates
  # lasso.mod <- cv.glmnet(X, Y, alpha = 1)
  # selected <- which(coef(lasso.mod) != 0)
  # if(length(selected) < 2) {
  #  selected <- 1:ncol(X)
  # } else {
  #  selected <- selected[-1] - 1 # Remove intercept
  #}
  
  LLCF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, num.trees = 8000, tune.parameters = "all")
  # LLCF.pred <- predict(LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
  # LLCF.CATE <- mean(LLCF.pred$predictions)
  # LLCF.CATE.SE <- mean((LLCF.pred$predictions - mean(LLCF.pred$predictions))^2)
  
  ate.subgroup <- average_treatment_effect(LLCF, target.sample = "all", subset = (characteristics == 1))
  ate.rest <- average_treatment_effect(LLCF, target.sample = "all", subset = !(characteristics == 1))
  resultsTable5OriginalPaper[4,(i+1)] <- paste(round(ate.subgroup[1] - ate.rest[1], 3), "(", sqrt(ate.subgroup[2]^2 + ate.rest[2]^2), ")")
  

 
}


# TO ADD BACK LATER MAYBE
# Predict: tuning done using set of lambdas
llcf.mse.old <- +Inf
for (l in length(lambdas)) {
  llcf.pred.old <- predict(LLCF, newdata = X[(characteristics == 1),], linear.correction.variables = 1:ncol(X), 
                           ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
  predictions <- llcf.pred.old$predictions
  llcf.mse.new <- mean((predictions - mean(predictions))**2)
  if (llcf.mse.new < llcf.mse.old) {
    llcf.mse.old <- llcf.mse.new
    # LLCF.CATE.Variance <- mean(llcf.pred.old$variance.estimates)
    predictions.new <- predictions
    predictions.variances.new <- llcf.pred.old$variance.estimates
  }
}

ate.subgroup <- mean(predictions.new)
ate.subgroup.VAR <- mean(predictions.variances.new)

# Predict: tuning done using set of lambdas
llcf.mse.old <- +Inf
for (l in length(lambdas)) {
  llcf.pred.old <- predict(LLCF, newdata = X[!(characteristics == 1),], linear.correction.variables = 1:ncol(X), 
                           ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
  predictions <- llcf.pred.old$predictions
  llcf.mse.new <- mean((predictions - mean(predictions))**2)
  if (llcf.mse.new < llcf.mse.old) {
    llcf.mse.old <- llcf.mse.new
    # LLCF.CATE.Variance <- mean(llcf.pred.old$variance.estimates)
    predictions.new <- predictions
    predictions.variances.new <- llcf.pred.old$variance.estimates
  }
}

ate.rest <- mean(predictions.new)
ate.rest.VAR <- mean(predictions.variances.new)

# ate.subgroup <- mean(predictions.new[(characteristics == 1)])
# ate.rest <- mean(predictions.new[!(characteristics == 1)])
# ate.subgroup.VAR <- mean(predictions.variances.new[(characteristics == 1)])
# ate.rest.VAR <- mean(predictions.variances.new[!(characteristics == 1)])