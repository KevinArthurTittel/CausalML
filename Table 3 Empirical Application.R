library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(stats)
library(haven)
library(ggplot2)
library(glmnet)
set.seed(123)

# Import and prepare Microfinance data set (Field et al., 2013)
Grace_Period_Data <- read_dta("Downloads/112672-V1/Grace-Period-Data.dta")
Grace_Period_Data <- as.data.frame(Grace_Period_Data)

# Appoint treatment assignment and outcome variables
W <- Grace_Period_Data$sec_treat
W <- as.vector(W)

# Create a numerical vector of the character group name vector
loangroups <- as.numeric(factor(Grace_Period_Data$sec_loanamount))

# Appoint the control variables matrix
X <- Grace_Period_Data[,c(2,6,8:18)]
colnames(X) <- c("Loan.Officer", "Stratification", "Age", "Married", "Literate", "Muslim", "HH.Size", "Years.Education", "Shock", "Has.Business",
                 "Financial.Control", "Home.Owner", "No.Drain")

# Create loan group dummies (as in original analysis) to be added to control variables matrix
loansize1 <- as.integer(c(Grace_Period_Data$sec_loanamount == 4000))
loansize2 <- as.integer(c(Grace_Period_Data$sec_loanamount == 5000))
loansize3 <- as.integer(c(Grace_Period_Data$sec_loanamount == 6000))
loansize4 <- as.integer(c(Grace_Period_Data$sec_loanamount == 8000))
loansize5 <- as.integer(c(Grace_Period_Data$sec_loanamount == 9000))
loansize6 <- as.integer(c(Grace_Period_Data$sec_loanamount == 10000))
loansizematrix <- cbind(loansize1, loansize2, loansize3, loansize4,
                        loansize5, loansize6)

# Initialize parameters
numtrees <- 4000 # Set to 1000 or 5000 to perform sensitivity analysis.
index <- c(1:5) # Concerns 5 characteristics; do not adjust.
dep.var <- 50 # Concerns dependent variable; set to 50 for monhtly profit, 51 for log of monthly HH income, or 52 for capital.
boolean.lambdas <- FALSE # Set to TRUE to use lambdas instead of automatic penalty tuning.
boolean.plot <- FALSE # Set to TRUE to make various plots of interest.
results_DiffATE <- c()

# Estimation procedure
for (i in index) {
  Y <- Grace_Period_Data[,(dep.var)]
  Y <- as.vector(Y)
  
  # Appoint characteristic for Table 5 (Field et al., 2013)
  characteristic <- Grace_Period_Data[,(74+i)]
  
  ###########################################
  ############ CR.LLCF Version 1 ############
  ###########################################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, characteristic))
  
  # Determine missing values
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristic))
  missingvalues[missingvalues == 2] <- 1
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  current.X <- as.matrix(current.X[!missingvalues,])
  current.W <- as.matrix(W[!missingvalues])
  current.loangroups <- loangroups[!missingvalues]
  current.characteristic <- characteristic[!missingvalues]
  current.Y <- Y[!missingvalues]
  
  # Grow preliminary forests for (W, X) and (Y, X) separately
  forest.W <- ll_regression_forest(current.X, current.W, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                   enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- ll_regression_forest(current.X, current.Y, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                   enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  Y.hat <- predict(forest.Y)$predictions
  
  # Select variables to include using Lasso feature selection
  lasso.mod <- cv.glmnet(current.X, current.Y, alpha = 1)
  selected <- which(coef(lasso.mod) != 0)
  if(length(selected) < 2) {
    selected <- 1:ncol(current.X)
  } else {
    selected <- selected[-1] - 1 # Remove intercept
  }
  
  # Implement LLCF
  CR.LLCF1 <- causal_forest(current.X, current.Y, current.W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE, 
                           clusters = current.loangroups, num.trees = numtrees)
  
  # Compute the variable importance
  CR.LLCF1.varimp <- variable_importance(CR.LLCF1) 
  CR.LLCF1.mostimportant <- colnames(current.X)[order(CR.LLCF1.varimp, decreasing = TRUE)]
  
  # Compute HTE with corresponding 95% confidence intervals 
  if (boolean.lambdas == FALSE) {
    # Predict: tuning without grid search over lambdas
    CR.LLCF1.pred <- predict(CR.LLCF1, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
    CR.LLCF1.CATE <- CR.LLCF1.pred$predictions
    CR.LLCF1.CATE.SE <- sqrt(CR.LLCF1.pred$variance.estimates)
  } else {
    # Predict: tuning done using set of lambdas
    CR.LLCF.mse.old <- +Inf
    for (l in length(lambdas)) {
      CR.LLCF.CATE.old <- predict(CR.LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
      predictions <- CR.LLCF.CATE.old$predictions
      CR.LLCF.mse.new <- mean((predictions - mean(predictions))**2)
      if (CR.LLCF.mse.new < CR.LLCF.mse.old) {
        CR.LLCF.mse.old <- CR.LLCF.mse.new
        CR.LLCF.CATE.SE <- sqrt(CR.LLCF.CATE.old$variance.estimates)
        CR.LLCF.CATE <- predictions
      }
    }
  }
  lower.CR.LLCF1 <- CR.LLCF1.CATE - qnorm(0.975)*CR.LLCF1.CATE.SE
  upper.CR.LLCF1 <- CR.LLCF1.CATE + qnorm(0.975)*CR.LLCF1.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  CR.LLCF1.ATE.charact <- average_treatment_effect(CR.LLCF1, target.sample = "all", subset = (current.characteristic == 1))
  CR.LLCF1.ATE.not.charact <- average_treatment_effect(CR.LLCF1, target.sample = "all", subset = !(current.characteristic == 1))
  CR.LLCF1.AIPW.charact <- get_scores(CR.LLCF1, subset = (current.characteristic == 1))
  CR.LLCF1.AIPW.not.charact <- get_scores(CR.LLCF1, subset = !(current.characteristic == 1))
  DiffATE.CR.LLCF1.test <- t.test(CR.LLCF1.AIPW.charact, CR.LLCF1.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ###########################################
  ############ CR.LLCF Version 1 ############
  ###########################################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, characteristic))
  
  # Determine missing values
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristic))
  missingvalues[missingvalues == 2] <- 1
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  current.X <- as.matrix(current.X[!missingvalues,])
  current.W <- as.matrix(W[!missingvalues])
  current.loangroups <- loangroups[!missingvalues]
  current.characteristic <- characteristic[!missingvalues]
  current.Y <- Y[!missingvalues]
  
  # Grow preliminary forests for (W, X) and (Y, X) separately
  forest.W <- ll_regression_forest(current.X, current.W, honesty = TRUE, honesty.fraction = 0.8, honesty.prune.leaves = FALSE,
                                   enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- ll_regression_forest(current.X, current.Y, honesty = TRUE, honesty.fraction = 0.8, honesty.prune.leaves = FALSE,
                                   enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  Y.hat <- predict(forest.Y)$predictions
  
  # Select variables to include using Lasso feature selection
  lasso.mod <- cv.glmnet(current.X, current.Y, alpha = 1)
  selected <- which(coef(lasso.mod) != 0)
  if(length(selected) < 2) {
    selected <- 1:ncol(current.X)
  } else {
    selected <- selected[-1] - 1 # Remove intercept
  }
  
  # Implement LLCF
  CR.LLCF2 <- causal_forest(current.X, current.Y, current.W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, honesty.fraction = 0.8, honesty.prune.leaves = FALSE, 
                            clusters = current.loangroups, num.trees = numtrees)
  
  # Compute the variable importance
  CR.LLCF2.varimp <- variable_importance(CR.LLCF2) 
  CR.LLCF2.mostimportant <- colnames(current.X)[order(CR.LLCF2.varimp, decreasing = TRUE)]
  
  # Compute HTE with corresponding 95% confidence intervals 
  if (boolean.lambdas == FALSE) {
    # Predict: tuning without grid search over lambdas
    CR.LLCF2.pred <- predict(CR.LLCF2, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
    CR.LLCF2.CATE <- CR.LLCF2.pred$predictions
    CR.LLCF2.CATE.SE <- sqrt(CR.LLCF2.pred$variance.estimates)
  } else {
    # Predict: tuning done using set of lambdas
    CR.LLCF.mse.old <- +Inf
    for (l in length(lambdas)) {
      CR.LLCF.CATE.old <- predict(CR.LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
      predictions <- CR.LLCF.CATE.old$predictions
      CR.LLCF.mse.new <- mean((predictions - mean(predictions))**2)
      if (CR.LLCF.mse.new < CR.LLCF.mse.old) {
        CR.LLCF.mse.old <- CR.LLCF.mse.new
        CR.LLCF.CATE.SE <- sqrt(CR.LLCF.CATE.old$variance.estimates)
        CR.LLCF.CATE <- predictions
      }
    }
  }
  lower.CR.LLCF2 <- CR.LLCF2.CATE - qnorm(0.975)*CR.LLCF2.CATE.SE
  upper.CR.LLCF2 <- CR.LLCF2.CATE + qnorm(0.975)*CR.LLCF2.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  CR.LLCF2.ATE.charact <- average_treatment_effect(CR.LLCF2, target.sample = "all", subset = (current.characteristic == 1))
  CR.LLCF2.ATE.not.charact <- average_treatment_effect(CR.LLCF2, target.sample = "all", subset = !(current.characteristic == 1))
  CR.LLCF2.AIPW.charact <- get_scores(CR.LLCF2, subset = (current.characteristic == 1))
  CR.LLCF2.AIPW.not.charact <- get_scores(CR.LLCF2, subset = !(current.characteristic == 1))
  DiffATE.CR.LLCF2.test <- t.test(CR.LLCF2.AIPW.charact, CR.LLCF2.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ###########################################
  ############ CR.LLCF Version 3 ############
  ###########################################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, characteristic))
  current.X.LLF <- as.matrix(cbind(X, loansizematrix, characteristic))
  
  # Determine missing values
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristic))
  missingvalues[missingvalues == 2] <- 1
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  current.X <- as.matrix(current.X[!missingvalues,])
  current.X.LLF <- as.matrix(current.X.LLF[!missingvalues,])
  current.W <- as.matrix(W[!missingvalues])
  current.loangroups <- loangroups[!missingvalues]
  current.characteristic <- characteristic[!missingvalues]
  current.Y <- Y[!missingvalues]
  
  # Grow preliminary forests for (W, X) and (Y, X) separately
  forest.W <- ll_regression_forest(current.X.LLF, current.W, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                   enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- ll_regression_forest(current.X.LLF, current.Y, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                   enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  Y.hat <- predict(forest.Y)$predictions
  
  # Select variables to include using Lasso feature selection
  lasso.mod <- cv.glmnet(current.X, current.Y, alpha = 1)
  selected <- which(coef(lasso.mod) != 0)
  if(length(selected) < 2) {
    selected <- 1:ncol(current.X)
  } else {
    selected <- selected[-1] - 1 # Remove intercept
  }
  
  # Implement LLCF
  CR.LLCF3 <- causal_forest(current.X, current.Y, current.W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE, 
                           clusters = current.loangroups, num.trees = numtrees)
  
  # Compute the variable importance
  CR.LLCF3.varimp <- variable_importance(CR.LLCF3) 
  CR.LLCF3.mostimportant <- colnames(current.X)[order(CR.LLCF3.varimp, decreasing = TRUE)]
  
  # Compute HTE with corresponding 95% confidence intervals 
  if (boolean.lambdas == FALSE) {
    # Predict: tuning without grid search over lambdas
    CR.LLCF3.pred <- predict(CR.LLCF3, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
    CR.LLCF3.CATE <- CR.LLCF.pred$predictions
    CR.LLCF3.CATE.SE <- sqrt(CR.LLCF3.pred$variance.estimates)
  } else {
    # Predict: tuning done using set of lambdas
    CR.LLCF.mse.old <- +Inf
    for (l in length(lambdas)) {
      CR.LLCF.CATE.old <- predict(CR.LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
      predictions <- CR.LLCF.CATE.old$predictions
      CR.LLCF.mse.new <- mean((predictions - mean(predictions))**2)
      if (CR.LLCF.mse.new < CR.LLCF.mse.old) {
        CR.LLCF.mse.old <- CR.LLCF.mse.new
        CR.LLCF.CATE.SE <- sqrt(CR.LLCF.CATE.old$variance.estimates)
        CR.LLCF.CATE <- predictions
      }
    }
  }
  lower.CR.LLCF3 <- CR.LLCF3.CATE - qnorm(0.975)*CR.LLCF3.CATE.SE
  upper.CR.LLCF3 <- CR.LLCF3.CATE + qnorm(0.975)*CR.LLCF3.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  CR.LLCF3.ATE.charact <- average_treatment_effect(CR.LLCF3, target.sample = "all", subset = (current.characteristic == 1))
  CR.LLCF3.ATE.not.charact <- average_treatment_effect(CR.LLCF3, target.sample = "all", subset = !(current.characteristic == 1))
  CR.LLCF3.AIPW.charact <- get_scores(CR.LLCF3, subset = (current.characteristic == 1))
  CR.LLCF3.AIPW.not.charact <- get_scores(CR.LLCF3, subset = !(current.characteristic == 1))
  DiffATE.CR.LLCF3.test <- t.test(CR.LLCF3.AIPW.charact, CR.LLCF3.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ###########################################
  ############ CR.LLCF Version 4 ############
  ###########################################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, characteristic))
  
  # Determine missing values
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristic))
  missingvalues[missingvalues == 2] <- 1
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  current.X <- as.matrix(current.X[!missingvalues,])
  current.W <- as.matrix(W[!missingvalues])
  current.loangroups <- loangroups[!missingvalues]
  current.characteristic <- characteristic[!missingvalues]
  current.Y <- Y[!missingvalues]
  
  # Grow preliminary forests for (W, X) and (Y, X) separately
  forest.W <- ll_regression_forest(current.X, current.W, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                   clusters = current.loangroups, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- ll_regression_forest(current.X, current.Y, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                   clusters = current.loangroups, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
  Y.hat <- predict(forest.Y)$predictions
  
  # Select variables to include using Lasso feature selection
  lasso.mod <- cv.glmnet(current.X, current.Y, alpha = 1)
  selected <- which(coef(lasso.mod) != 0)
  if(length(selected) < 2) {
    selected <- 1:ncol(current.X)
  } else {
    selected <- selected[-1] - 1 # Remove intercept
  }
  
  # Implement LLCF
  CR.LLCF4 <- causal_forest(current.X, current.Y, current.W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE, 
                            clusters = current.loangroups, num.trees = numtrees)
  
  # Compute the variable importance
  CR.LLCF4.varimp <- variable_importance(CR.LLCF4) 
  CR.LLCF4.mostimportant <- colnames(current.X)[order(CR.LLCF4.varimp, decreasing = TRUE)]
  
  # Compute HTE with corresponding 95% confidence intervals 
  if (boolean.lambdas == FALSE) {
    # Predict: tuning without grid search over lambdas
    CR.LLCF4.pred <- predict(CR.LLCF4, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
    CR.LLCF4.CATE <- CR.LLCF4.pred$predictions
    CR.LLCF4.CATE.SE <- sqrt(CR.LLCF4.pred$variance.estimates)
  } else {
    # Predict: tuning done using set of lambdas
    CR.LLCF.mse.old <- +Inf
    for (l in length(lambdas)) {
      CR.LLCF.CATE.old <- predict(CR.LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
      predictions <- CR.LLCF.CATE.old$predictions
      CR.LLCF.mse.new <- mean((predictions - mean(predictions))**2)
      if (CR.LLCF.mse.new < CR.LLCF.mse.old) {
        CR.LLCF.mse.old <- CR.LLCF.mse.new
        CR.LLCF.CATE.SE <- sqrt(CR.LLCF.CATE.old$variance.estimates)
        CR.LLCF.CATE <- predictions
      }
    }
  }
  lower.CR.LLCF4 <- CR.LLCF4.CATE - qnorm(0.975)*CR.LLCF4.CATE.SE
  upper.CR.LLCF4 <- CR.LLCF4.CATE + qnorm(0.975)*CR.LLCF4.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  CR.LLCF4.ATE.charact <- average_treatment_effect(CR.LLCF4, target.sample = "all", subset = (current.characteristic == 1))
  CR.LLCF4.ATE.not.charact <- average_treatment_effect(CR.LLCF4, target.sample = "all", subset = !(current.characteristic == 1))
  CR.LLCF4.AIPW.charact <- get_scores(CR.LLCF4, subset = (current.characteristic == 1))
  CR.LLCF4.AIPW.not.charact <- get_scores(CR.LLCF4, subset = !(current.characteristic == 1))
  DiffATE.CR.LLCF4.test <- t.test(CR.LLCF4.AIPW.charact, CR.LLCF4.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ############################
  ###### Update results ######
  ############################
  
  results_DiffATE <- rbind(results_DiffATE, 
                           data.frame(t(c(paste(round(CR.LLCF1.ATE.charact[1], 3), "(", round(CR.LLCF1.ATE.charact[2], 3), ")"),
                                          paste(round(CR.LLCF1.ATE.not.charact[1], 3), "(", round(CR.LLCF1.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.CR.LLCF1.test$p.value, 3)),
                                          paste(round(CR.LLCF2.ATE.charact[1], 3), "(", round(CR.LLCF2.ATE.charact[2], 3), ")"),
                                          paste(round(CR.LLCF2.ATE.not.charact[1], 3), "(", round(CR.LLCF2.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.CR.LLCF2.test$p.value, 3)),
                                          paste(round(CR.LLCF3.ATE.charact[1], 3), "(", round(CR.LLCF3.ATE.charact[2], 3), ")"),
                                          paste(round(CR.LLCF3.ATE.not.charact[1], 3), "(", round(CR.LLCF3.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.CR.LLCF3.test$p.value, 3)),
                                          paste(round(CR.LLCF4.ATE.charact[1], 3), "(", round(CR.LLCF4.ATE.charact[2], 3), ")"),
                                          paste(round(CR.LLCF4.ATE.not.charact[1], 3), "(", round(CR.LLCF4.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.CR.LLCF4.test$p.value, 3)),
                                          (845 - sum(missingvalues))))))
}

colnames(results_DiffATE) <- c("CR.LLCF 1 CATE subgroup", "CR.LLCF 1 CATE rest", "CR.LLCF 1 p-value difference", 
                               "CR.LLCF 2 CATE subgroup", "CR.LLCF 2 CATE rest", "CR.LLCF 2 p-value difference", 
                               "CR.LLCF 3 CATE subgroup", "CR.LLCF 3 CATE rest", "CR.LLCF 3 p-value difference", 
                               "CR.LLCF 4 CATE subgroup", "CR.LLCF 4 CATE rest", "CR.LLCF 4 p-value difference", "Number of observations")
rownames(results_DiffATE) <- c("savings", "risk loving", "wage earner", "household member chronically ill", "impatient")
results_DiffATE[c(2,1,4,5,3),]

results_table3 <- list('Sheet1' = results_DiffATE)
write.xlsx(results_table5, file = "Table 3.xlsx")
