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
dep.var <- 52 # Concerns dependent variable; set to 50 for monhtly profit, 51 for log of monthly HH income, or 52 for capital.
lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5) # Concerns ridge penalty parameters; do not adjust.
boolean.lambdas <- FALSE # Set to TRUE to use lambdas instead of automatic penalty tuning.
boolean.plot <- FALSE # Set to TRUE to make various plots of interest.
results_DiffATE <- c()

# Estimation procedure
for (i in index) {
  Y <- Grace_Period_Data[,(dep.var)]
  Y <- as.vector(Y)
  
  # Appoint characteristic for Table 5 (Field et al., 2013)
  characteristic <- Grace_Period_Data[,(74+i)]
  
  ###########################
  ########### GRF ###########
  ###########################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, loansizematrix, characteristic))
  
  # Determine missing values
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristic))
  missingvalues[missingvalues == 2] <- 1
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  current.X <- current.X[!missingvalues,]
  current.W <- W[!missingvalues]
  current.loangroups <- loangroups[!missingvalues]
  current.characteristic <- characteristic[!missingvalues]
  current.Y <- Y[!missingvalues]
  
  # Grow preliminary forests for (Y, X) 
  forest.Y <- regression_forest(current.X, current.Y, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
  
  # Compute the variable importance
  GRF.varimp <- variable_importance(forest.Y) 
  GRF.mostimportant <- colnames(current.X)[order(GRF.varimp)[1:4]] # 4 most important variables for splitting
  GRF.mostimportant <- c(GRF.mostimportant, colnames(characteristic))
  
  # Select variables to include using preliminary GRF
  prelim.GRF <- causal_forest(current.X, current.Y, current.W, num.trees = numtrees, honesty = TRUE)
  prelim.GRF.varimp <- variable_importance(prelim.GRF)
  selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)
  
  # Implement GRF
  GRF <- causal_forest(current.X[,selected.vars], current.Y, current.W, num.trees = numtrees, 
                       honesty = TRUE, tune.parameters = "all")
  
  # Compute ATE
  GRF.ATE <- average_treatment_effect(GRF, target.sample = "all")
  
  # Compute HTE with corresponding 95% confidence intervals  
  GRF.pred <- predict(GRF, estimate.variance = TRUE)
  GRF.CATE <- GRF.pred$predictions
  GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
  lower.GRF <- GRF.CATE - qnorm(0.975)*GRF.CATE.SE
  upper.GRF <- GRF.CATE + qnorm(0.975)*GRF.CATE.SE

  # Test of heterogeneity using Differential ATE
  GRF.ATE.charact <- average_treatment_effect(GRF, target.sample = "all", subset = (current.characteristic == 1))
  GRF.ATE.not.charact <- average_treatment_effect(GRF, target.sample = "all", subset = !(current.characteristic == 1))
  GRF.AIPW.charact <- get_scores(GRF, subset = (current.characteristic == 1))
  GRF.AIPW.not.charact <- get_scores(GRF, subset = !(current.characteristic == 1))
  DiffATE.GRF.test <- t.test(GRF.AIPW.charact, GRF.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ############################
  #### Cluster-Robust GRF ####
  ############################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, characteristic))
  
  # Determine missing values
  missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristic))
  missingvalues[missingvalues == 2] <- 1
  
  # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
  current.X <- current.X[!missingvalues,]
  current.W <- W[!missingvalues]
  current.loangroups <- loangroups[!missingvalues]
  current.characteristic <- characteristic[!missingvalues]
  current.Y <- Y[!missingvalues]
  
  # Grow preliminary forests for (Y, X) 
  forest.Y <- regression_forest(current.X, current.Y, num.trees = numtrees, honesty = TRUE, clusters = current.loangroups)
  
  # Select variables to include using preliminary Cluster-Robust GRF
  prelim.CR.GRF <- causal_forest(current.X, current.Y, current.W, honesty = TRUE, clusters = current.loangroups, num.trees = numtrees)
  prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
  selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)
  
  # Compute the variable importance
  CR.GRF.varimp <- variable_importance(forest.Y) 
  CR.GRF.mostimportant <- colnames(current.X)[order(CR.GRF.varimp)[1:4]] # 4 most important variables for splitting
  CR.GRF.mostimportant <- c(CR.GRF.mostimportant, colnames(characteristic))
  
  # Implement Cluster-Robust GRF
  CR.GRF <- causal_forest(current.X[,selected.vars], current.Y, current.W, clusters = current.loangroups, honesty = TRUE, num.trees = numtrees)
  
  # Compute ATE 
  CR.GRF.ATE <- average_treatment_effect(CR.GRF, target.sample = "all")
  
  # Compute HTE with corresponding 95% confidence intervals
  CR.GRF.pred <- predict(CR.GRF, estimate.variance = TRUE)
  CR.GRF.CATE <- CR.GRF.pred$predictions
  CR.GRF.CATE.SE <- sqrt(CR.GRF.pred$variance.estimates)
  lower.CR.GRF <- CR.GRF.CATE - qnorm(0.975)*CR.GRF.CATE.SE
  upper.CR.GRF <- CR.GRF.CATE + qnorm(0.975)*CR.GRF.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  CR.GRF.ATE.charact <- average_treatment_effect(CR.GRF, target.sample = "all", subset = (current.characteristic == 1))
  CR.GRF.ATE.not.charact <- average_treatment_effect(CR.GRF, target.sample = "all", subset = !(current.characteristic == 1))
  CR.GRF.AIPW.charact <- get_scores(CR.GRF, subset = (current.characteristic == 1))
  CR.GRF.AIPW.not.charact <- get_scores(CR.GRF, subset = !(current.characteristic == 1))
  DiffATE.CR.GRF.test <- t.test(CR.GRF.AIPW.charact, CR.GRF.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ############################
  ########### LLCF ###########
  ############################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, loansizematrix, characteristic))
  
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
  forest.W <- ll_regression_forest(current.X, current.W, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                   num.trees = 1000)
  W.hat <- predict(forest.W)$predictions
  forest.Y <- ll_regression_forest(current.X, current.Y, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                   num.trees = 1000)
  Y.hat <- predict(forest.Y)$predictions
  
  # Compute the variable importance
  LLCF.varimp <- variable_importance(forest.Y) 
  LLCF.mostimportant <- colnames(current.X)[order(LLCF.varimp)[1:4]] # 4 most important variables for splitting
  LLCF.mostimportant <- c(LLCF.mostimportant, colnames(characteristic))
  
  # Select variables to include using Lasso feature selection
  lasso.mod <- cv.glmnet(current.X, current.Y, alpha = 1)
  selected <- which(coef(lasso.mod) != 0)
  if(length(selected) < 2) {
    selected <- 1:ncol(current.X)
  } else {
    selected <- selected[-1] - 1 # Remove intercept
  }
  
  # Implement LLCF
  LLCF <- causal_forest(current.X, current.Y, current.W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE,  
                        num.trees = numtrees, tune.parameters = "all")
  
  # Compute ATE 
  LLCF.ATE <- average_treatment_effect(LLCF, target.sample = "all")
  
  # Compute HTE with corresponding 95% confidence intervals 
  if (boolean.lambdas == FALSE) {
    # Predict: tuning without grid search over lambdas
    LLCF.pred <- predict(LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
    LLCF.CATE <- LLCF.pred$predictions
    LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
  } else {
    # Predict: tuning done using set of lambdas
    LLCF.mse.old <- +Inf
    for (l in length(lambdas)) {
      LLCF.CATE.old <- predict(LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
      predictions <- LLCF.CATE.old$predictions
      LLCF.mse.new <- mean((predictions - mean(predictions))**2)
      if (LLCF.mse.new < LLCF.mse.old) {
        LLCF.mse.old <- LLCF.mse.new
        LLCF.CATE.SE <- sqrt(LLCF.CATE.old$variance.estimates)
        LLCF.CATE <- predictions
      }
    }
  }
  lower.LLCF <- LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE
  upper.LLCF <- LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  LLCF.ATE.charact <- average_treatment_effect(LLCF, target.sample = "all", subset = (current.characteristic == 1))
  LLCF.ATE.not.charact <- average_treatment_effect(LLCF, target.sample = "all", subset = !(current.characteristic == 1))
  LLCF.AIPW.charact <- get_scores(LLCF, subset = (current.characteristic == 1))
  LLCF.AIPW.not.charact <- get_scores(LLCF, subset = !(current.characteristic == 1))
  DiffATE.LLCF.test <- t.test(LLCF.AIPW.charact, LLCF.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ###########################################
  ######### LLCF Standard splitting #########
  ###########################################
  
  # Appoint final X-matrix
  current.X <- as.matrix(cbind(X, loansizematrix, characteristic))
  
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
  
  # Compute the variable importance
  CR.LLCF.varimp <- variable_importance(forest.Y) 
  CR.LLCF.mostimportant <- colnames(current.X)[order(CR.LLCF.varimp)[1:4]] # 4 most important variables for splitting
  CR.LLCF.mostimportant <- c(CR.LLCF.mostimportant, colnames(characteristic))
  
  # Select variables to include using Lasso feature selection
  lasso.mod <- cv.glmnet(current.X, current.Y, alpha = 1)
  selected <- which(coef(lasso.mod) != 0)
  if(length(selected) < 2) {
    selected <- 1:ncol(current.X)
  } else {
    selected <- selected[-1] - 1 # Remove intercept
  }
  
  # Implement LLCF
  CR.LLCF <- causal_forest(current.X, current.Y, current.W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, num.trees = numtrees, tune.num.trees = 50,
            honesty.fraction = 0.7, honesty.prune.leaves = FALSE, tune.parameters = c("min.node.size", "sample.fraction", "mtry", "alpha", "imbalance.penalty"))
  
  # Compute ATE 
  CR.LLCF.ATE <- average_treatment_effect(CR.LLCF, target.sample = "all")
  
  # Compute HTE with corresponding 95% confidence intervals 
  if (boolean.lambdas == FALSE) {
    # Predict: tuning without grid search over lambdas
    CR.LLCF.pred <- predict(CR.LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
    CR.LLCF.CATE <- CR.LLCF.pred$predictions
    CR.LLCF.CATE.SE <- sqrt(CR.LLCF.pred$variance.estimates)
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
  lower.CR.LLCF <- CR.LLCF.CATE - qnorm(0.975)*CR.LLCF.CATE.SE
  upper.CR.LLCF <- CR.LLCF.CATE + qnorm(0.975)*CR.LLCF.CATE.SE
  
  # Test of heterogeneity using Differential ATE
  CR.LLCF.ATE.charact <- average_treatment_effect(CR.LLCF, target.sample = "all", subset = (current.characteristic == 1))
  CR.LLCF.ATE.not.charact <- average_treatment_effect(CR.LLCF, target.sample = "all", subset = !(current.characteristic == 1))
  CR.LLCF.AIPW.charact <- get_scores(CR.LLCF, subset = (current.characteristic == 1))
  CR.LLCF.AIPW.not.charact <- get_scores(CR.LLCF, subset = !(current.characteristic == 1))
  DiffATE.CR.LLCF.test <- t.test(CR.LLCF.AIPW.charact, CR.LLCF.AIPW.not.charact, alternative = "two.sided", var.equal = FALSE)
  
  ############################
  ###### Update results ######
  ############################
  
  results_BLP <- rbind(results_BLP, 
                       data.frame(t(c(paste(round(GRF.BLP[1,1], 3), "(", round(GRF.BLP[1,2], 3), ")"),
                                      paste(round(GRF.BLP[2,1], 3), "(", round(GRF.BLP[2,2], 3), ")"),
                                      paste(round(CR.GRF.BLP[1,1], 3), "(", round(CR.GRF.BLP[1,2], 3), ")"),
                                      paste(round(CR.GRF.BLP[2,1], 3), "(", round(CR.GRF.BLP[2,2], 3), ")"),
                                      paste(round(LLCF.BLP[1,1], 3), "(", round(LLCF.BLP[1,2], 3), ")"),
                                      paste(round(LLCF.BLP[2,1], 3), "(", round(LLCF.BLP[2,2], 3), ")"),
                                      paste(round(CR.LLCF.BLP[1,1], 3), "(", round(CR.LLCF.BLP[1,2], 3), ")"),
                                      paste(round(CR.LLCF.BLP[2,1], 3), "(", round(CR.LLCF.BLP[2,2], 3), ")")))))
  
  results_DiffATE <- rbind(results_DiffATE, 
                           data.frame(t(c(paste(round(GRF.ATE.charact[1], 3), "(", round(GRF.ATE.charact[2], 3), ")"),
                                          paste(round(GRF.ATE.not.charact[1], 3), "(", round(GRF.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.GRF.test$p.value, 3)),
                                          paste(round(CR.GRF.ATE.charact[1], 3), "(", round(CR.GRF.ATE.charact[2], 3), ")"),
                                          paste(round(CR.GRF.ATE.not.charact[1], 3), "(", round(CR.GRF.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.CR.GRF.test$p.value, 3)),
                                          paste(round(LLCF.ATE.charact[1], 3), "(", round(LLCF.ATE.charact[2], 3), ")"),
                                          paste(round(LLCF.ATE.not.charact[1], 3), "(", round(LLCF.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.LLCF.test$p.value, 3)),
                                          paste(round(CR.LLCF.ATE.charact[1], 3), "(", round(CR.LLCF.ATE.charact[2], 3), ")"),
                                          paste(round(CR.LLCF.ATE.not.charact[1], 3), "(", round(CR.LLCF.ATE.not.charact[2], 3), ")"),
                                          paste(round(DiffATE.CR.LLCF.test$p.value, 3)),
                                          (845 - sum(missingvalues))))))
}

colnames(results_DiffATE) <- c("GRF CATE subgroup", "GRF CATE rest", "GRF p-value difference", 
                               "CR.GRF CATE subgroup", "CR.GRF CATE rest", "CR.GRF p-value difference",
                               "LLCF CATE subgroup", "LLCF CATE rest", "LLCF p-value difference", 
                               "LLCF Other splitting CATE subgroup", "LLCF Other splitting CATE rest", "LLCF Other splitting p-value difference","Number of observations")
rownames(results_DiffATE) <- c("savings", "risk loving", "wage earner", "household member chronically ill", "impatient")

results_table5 <- list('Sheet1' = results_DiffATE)
write.xlsx(results_table5, file = "Table 5 (Field et al., 2013).xlsx")


