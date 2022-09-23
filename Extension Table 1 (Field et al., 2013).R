library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(haven)

# Import and prepare Microfinance data set (Field et al., 2013)
  Grace_Period_Data <- read_dta("Downloads/112672-V1/Grace-Period-Data.dta")
  Grace_Period_Data <- as.data.frame(Grace_Period_Data)

  # Appoint treatment assignment and outcome variables
    W <- Grace_Period_Data$sec_treat
    W <- as.vector(W)

  # Create a numerical vector of the character group name vector
    loangroups <- as.numeric(factor(Grace_Period_Data$sec_group_name))

  # Standardize the continuous variables
    Grace_Period_Data$Years_Education_C <- scale(Grace_Period_Data$Years_Education_C)
    Grace_Period_Data$Age_C <- scale(Grace_Period_Data$Age_C)
    Grace_Period_Data$SEI <- scale(Grace_Period_Data$SEI)
    Grace_Period_Data$HH_Size_C <- scale(Grace_Period_Data$HH_Size_C)

  # Appoint the control variables matrix
    X <- Grace_Period_Data[,8:18]

  # Create loan group dummies (as in original analysis) to be added to control variables matrix
    loansize1 <- as.integer(c(Grace_Period_Data$sec_loanamount == 4000))
    loansize2 <- as.integer(c(Grace_Period_Data$sec_loanamount == 5000))
    loansize3 <- as.integer(c(Grace_Period_Data$sec_loanamount == 6000))
    loansize4 <- as.integer(c(Grace_Period_Data$sec_loanamount == 8000))
    loansize5 <- as.integer(c(Grace_Period_Data$sec_loanamount == 9000))
    loansize6 <- as.integer(c(Grace_Period_Data$sec_loanamount == 10000))
    loansizematrix <- cbind(loansize1, loansize2, loansize3, loansize4,
                        loansize5, loansize6)

  # Create stratification dummies (fixed effects) to be added to control variables matrix
    stratifgroup1 <- as.integer(c(Grace_Period_Data$Stratification_Dummies == 1))
    stratifgroup2 <- as.integer(c(Grace_Period_Data$Stratification_Dummies == 2))
    stratifgroup3 <- as.integer(c(Grace_Period_Data$Stratification_Dummies == 3))
    stratifgroup4 <- as.integer(c(Grace_Period_Data$Stratification_Dummies == 4))
    stratifgroup5 <- as.integer(c(Grace_Period_Data$Stratification_Dummies == 5))
    stratifgroup6 <- as.integer(c(Grace_Period_Data$Stratification_Dummies == 6))
    stratifmatrix <- cbind(stratifgroup1, stratifgroup2, stratifgroup3, 
                       stratifgroup4, stratifgroup5, stratifgroup6)

  # Create loan officer dummies (fixed effects) to be added to control variables matrix
    loanofficergroup1 <- as.integer(c(Grace_Period_Data$sec_loan_officer == 1))
    loanofficergroup2 <- as.integer(c(Grace_Period_Data$sec_loan_officer == 3))
    loanofficergroup3 <- as.integer(c(Grace_Period_Data$sec_loan_officer == 4))
    loanofficergroup4 <- as.integer(c(Grace_Period_Data$sec_loan_officer == 7))
    loanofficergroup5 <- as.integer(c(Grace_Period_Data$sec_loan_officer == 8))
    loanofficermatrix <- cbind(loanofficergroup1, loanofficergroup2, loanofficergroup3, 
                           loanofficergroup4, loanofficergroup5)

  # Combine all the control covariates in one large matrix
    # X <- as.matrix(cbind(X, stratifmatrix, loansizematrix, loanofficermatrix))
    X <- as.matrix(cbind(X, loansizematrix))
    # X[X == "NA"] <- 0 # Set all missing values to 0 as done in original paper

# Initialize results matrix
  resultsTable1OriginalPaper <- matrix(data = 0, nrow = 13, ncol = 4)
  resultsTable1OriginalPaper[2:13,1] <- c("Total business spending", "Inventory and raw materials", "Business equipment",
                                        "Operating costs", "Total non-business spending", "Home repairs", "Utilities, taxes and rent",
                                        "Human capital", "Money for relending", "Savings", "Food and durable consumption", "New businesses")
  resultsTable1OriginalPaper[1,] <- c("Dependent variable", "CF", "Cluster-robust CF", "LLCF")

# Estimation procedure
  for (i in 1:12) {
    Y <- Grace_Period_Data[,(37+i)]
    Y <- as.vector(Y)
  
    # Causal Forest estimation:
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- regression_forest(X, W)
        W.hat <- predict(forest.W)$predictions
        forest.Y <- regression_forest(X, Y)
        Y.hat <- predict(forest.Y)$predictions

      # Select variables to include using preliminary CF
        prelim.CF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, 
                              num.trees = 4000)
        prelim.CF.varimp <- variable_importance(prelim.CF)
        selected.vars <- which(prelim.CF.varimp / mean(prelim.CF.varimp) > 0.2)

      # Implement CF
        CF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat,
                       num.trees = 8000, tune.parameters = "all")
        tauhat <- predict(CF)$predictions
  
      # Graph the predicted CF heterogeneous treatment effect estimates
        hist(tauhat)

      # Compute CATE with corresponding 95% confidence intervals
        CATE <- average_treatment_effect(CF, target.sample = "all")
        resultsTable1OriginalPaper[(i+1),2] <- paste(round(CATE[1], 3), "(", round(CATE[2], 3), ")")
        # results[i,2] <- paste("[", (round(CATE[1], 3) - (round(qnorm(0.975) * CATE[2], 3))), ",",  
          #                    round(CATE[1], 3) + (round(qnorm(0.975) * CATE[2], 3)), "]")
        # paste("95% CI for the CATE:", (round(CATE.CR[1], 3) - (round(qnorm(0.975) * CATE.CR[2], 3))),
        # round(CATE.CR[1], 3) + (round(qnorm(0.975) * CATE.CR[2], 3)))
  
    # Cluster-Robust Causal Forest estimation:
      # Select variables to include using preliminary Cluster-Robust CF
        prelim.CF.CR <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, 
                           num.trees = 4000)
        prelim.CF.CR.varimp <- variable_importance(prelim.CF.CR)
        selected.vars <- which(prelim.CF.CR.varimp / mean(prelim.CF.CR.varimp) > 0.2)

      # Implement Cluster-Robust CF
        CF.CR <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat,
                    clusters = loangroups, 
                    num.trees = 8000, tune.parameters = "all")
        tauhat <- predict(CF.CR)$predictions
    
      # Graph the predicted Cluster-Robust CF heterogeneous treatment effect estimates
        hist(tauhat)
    
      # Compute CATE with corresponding 95% confidence intervals
        CATE.CR <- average_treatment_effect(CF.CR, target.sample = "all")
        resultsTable1OriginalPaper[(i+1),3] <- paste(round(CATE.CR[1], 3), "(", round(CATE.CR[2], 3), ")")
        # results[i,4] <- paste("[", (round(CATE.CR[1], 3) - (round(qnorm(0.975) * CATE.CR[2], 3))), ",",  
          #                   round(CATE.CR[1], 3) + (round(qnorm(0.975) * CATE.CR[2], 3)), "]")
        #paste("95% CI for the CATE:", (round(CATE.CR[1], 3) - (round(qnorm(0.975) * CATE.CR[2], 3))),
                                                   #round(CATE.CR[1], 3) + (round(qnorm(0.975) * CATE.CR[2], 3)))
  
    # Linear Local Causal Forest estimation:
       # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- ll_regression_forest(X, W, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 2000, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- ll_regression_forest(X, Y, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 2000, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions
  
      # Select variables to include using preliminary LLCF
        lasso.mod <- cv.glmnet(X, Y, alpha = 1)
        selected <- which(coef(lasso.mod) != 0)
        if(length(selected) < 2) {
          selected <- 1:ncol(X)
        } else {
          selected <- selected[-1] - 1 # Remove intercept
        }
  
      # Implement LLCF
        LLCF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE,
                              num.trees = 2000, tune.parameters = "all")
        LLCF.pred <- predict(LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
        LLCF.ATE <- mean(LLCF.pred$predictions)
        LLCF.ATE.SE <- mean((LLCF.pred$predictions - mean(LLCF.pred$predictions))^2)
  
      # Predict: tuning done using set of lambdas
        llcf.mse.old <- +Inf
        for (l in length(lambdas)) {
          llcf.pred.old <- predict(LLCF, linear.correction.variables = 1:ncol(X), ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
          predictions <- llcf.pred.old$predictions
          llcf.mse.new <- mean((predictions - mean(predictions))**2)
          if (llcf.mse.new < llcf.mse.old) {
            llcf.mse.old <- llcf.mse.new
            LLCF.CATE.SE <- sqrt(mean(llcf.pred.old$variance.estimates))
            predictions.new <- predictions
          }
        }
  
        LLCF.ATE <- mean(predictions.new)
        resultsTable1OriginalPaper[(i+1),4] <- paste(LLCF.ATE, "(", LLCF.ATE.SE, ")")
}
