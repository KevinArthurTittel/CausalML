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
 
# Initialize results matrix
  resultsTable5OriginalPaper <- matrix(data = 0, nrow = 5, ncol = 6)
  resultsTable5OriginalPaper[1,] <- t(c("dependent variable", "savings", "risk loving", "wage earner", "household member chronically ill",
                   "impatient"))
  resultsTable5OriginalPaper[2:5,1] <- c("CF", "Cluster-robust CF", "LLCF", "observations")

# Initialize parameters
  numtrees <- 2000
  index <- c(1:6)
  dep.var <- 50 # 50 = monhtly profit, 51 = log of monthly HH income, 52 = capital
  lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5)
  boolean.lambdas <- FALSE
  boolean.plot <- FALSE
  set.seed(123)

# Estimation procedure
run_method = function(numtrees, index, lambdas, boolean.plot, boolean.lambdas) {
  basic.results = sapply(index, function(i) {
    Y <- Grace_Period_Data[,(dep.var)]
    Y <- as.vector(Y)
  
    # Appoint characteristic for Table 5 (Field et al., 2013)
      characteristics <- Grace_Period_Data[,(74+i)]
  
    # Determine missing values
      missingvalues <- as.integer(is.na(Y)) + as.integer(is.na(characteristics))
      missingvalues[missingvalues == 2] <- 1
  
    # Remove the observations for which Y and/or the characteristic for Table 5 has an NA value
      X <- X[!missingvalues,]
      W <- W[!missingvalues]
      loangroups <- loangroups[!missingvalues]
      characteristics <- characteristics[!missingvalues]
      Y <- Y[!missingvalues]
  
    ###########################
    ########### GRF ###########
    ###########################
    
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- regression_forest(X, W, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- regression_forest(X, Y, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions
        
      # Compute the variable importance
        GRF.varimp <- variable_importance(forest.Y) 
    
      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, honesty = TRUE)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)
  
      # Implement GRF
        GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, 
                             honesty = TRUE, tune.parameters = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction"))
        GRF.pred <- predict(GRF, estimate.variance = TRUE)
        GRF.CATE <- GRF.pred$predictions
        GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
  
      # Find lower and upper bounds for 95% confidence intervals
        lower.GRF <- GRF.CATE - qnorm(0.975)*GRF.CATE.SE
        upper.GRF <- GRF.CATE + qnorm(0.975)*GRF.CATE.SE
    
      # Graph the predicted GRF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(GRF.CATE)
        }
    
      # Compute ATE with corresponding 95% confidence intervals
        GRF.ATE <- average_treatment_effect(GRF, target.sample = "all")
    
      # See if the GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        GRF.rate <- rank_average_treatment_effect(GRF, GRF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(GRF.rate)
        }

      # Assessing GRF fit and heterogeneity using the Best Linear Predictor Approach [BLP]
        GRF.BLP <- test_calibration(GRF)
        mean.forest.pred.GRF <- c(GRF.BLP[1,1], GRF.BLP[1,2], (GRF.BLP[1,4] < 0.10))
        diff.forest.pred.GRF <- c(GRF.BLP[2,1], GRF.BLP[2,2], (GRF.BLP[2,4] < 0.10))
    
      # Test of heterogeneity using Differential ATE
        GRF.ATE.charact <- average_treatment_effect(GRF, target.sample = "all", subset = (characteristics == 1))
        GRF.ATE.not.charact <- average_treatment_effect(GRF, target.sample = "all", subset = !(characteristics == 1))
        DiffATE.GRF.mean <- GRF.ATE.charact[1] - GRF.ATE.not.charact[1]
        DiffATE.GRF.SE <- sqrt(GRF.ATE.charact[2]^2 + GRF.ATE.not.charact[2]^2)
        lower.DiffATE.GRF <- (DiffATE.GRF.mean - (qnorm(0.975) * DiffATE.GRF.SE))
        upper.DiffATE.GRF <- (DiffATE.GRF.mean + (qnorm(0.975) * DiffATE.GRF.SE))
  
    ############################
    #### Cluster-Robust GRF ####
    ############################
      
      # Select variables to include using preliminary Cluster-Robust GRF
        prelim.CR.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, clusters = loangroups, num.trees = numtrees)
        prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
        selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)
  
      # Compute the variable importance
        CR.GRF.varimp <- variable_importance(forest.Y) 
              
      # Implement Cluster-Robust GRF
        CR.GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, honesty = TRUE, num.trees = numtrees, 
                                tune.parameters = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction"))
        CR.GRF.pred <- predict(CR.GRF, estimate.variance = TRUE)
        CR.GRF.CATE <- CR.GRF.pred$predictions
        CR.GRF.CATE.SE <- sqrt(CR.GRF.pred$variance.estimates)
    
      # Find lower and upper bounds for 95% confidence intervals
        lower.CR.GRF <- CR.GRF.CATE - qnorm(0.975)*CR.GRF.CATE.SE
        upper.CR.GRF <- CR.GRF.CATE + qnorm(0.975)*CR.GRF.CATE.SE
              
      # Graph the predicted Cluster-Robust GRF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(CR.GRF.CATE)
        }
    
      # Compute ATE with corresponding 95% confidence intervals
        CR.GRF.ATE <- average_treatment_effect(CR.GRF, target.sample = "all")
    
      # See if the Cluster-Robust GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        CR.GRF.rate <- rank_average_treatment_effect(CR.GRF, CR.GRF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(CR.GRF.rate)
        }
  
      # Assessing Cluster-Robust GRF fit using the Best Linear Predictor Approach [BLP]
        CR.GRF.BLP <- test_calibration(CR.GRF)
        mean.forest.pred.CR.GRF <- c(CR.GRF.BLP[1,1], CR.GRF.BLP[1,2], (CR.GRF.BLP[1,4] < 0.10))
        diff.forest.pred.CR.GRF <- c(CR.GRF.BLP[2,1], CR.GRF.BLP[2,2], (CR.GRF.BLP[2,4] < 0.10))
    
      # Test of heterogeneity using Differential ATE
        CR.GRF.ATE.charact <- average_treatment_effect(CR.GRF, target.sample = "all", subset = (characteristics == 1))
        CR.GRF.ATE.not.charact <- average_treatment_effect(CR.GRF, target.sample = "all", subset = !(characteristics == 1))
        DiffATE.CR.GRF.mean <- CR.GRF.ATE.charact[1] - CR.GRF.ATE.not.charact[1]
        DiffATE.CR.GRF.SE <- sqrt(CR.GRF.ATE.charact[2]^2 + CR.GRF.ATE.not.charact[2]^2)
        lower.DiffATE.CR.GRF <- (DiffATE.CR.GRF.mean - (qnorm(0.975) * DiffATE.CR.GRF.SE))
        upper.DiffATE.CR.GRF <- (DiffATE.CR.GRF.mean + (qnorm(0.975) * DiffATE.CR.GRF.SE))
  
    ############################
    ########### LLCF ###########
    ############################
    
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- ll_regression_forest(X, W, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                         num.trees = numtrees, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- ll_regression_forest(X, Y, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                         num.trees = numtrees, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions
    
      # Compute the variable importance
        LLCF.varimp <- variable_importance(forest.Y) 
              
      # Select variables to include using Lasso feature selection
        lasso.mod <- cv.glmnet(X, Y, alpha = 1)
        selected <- which(coef(lasso.mod) != 0)
        if(length(selected) < 2) {
          selected <- 1:ncol(X)
        } else {
          selected <- selected[-1] - 1 # Remove intercept
        }
  
      # Implement LLCF
        LLCF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                              num.trees = numtrees, tune.parameters = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction"))
        LLCF.ATE <- average_treatment_effect(LLCF, target.sample = "all")
    
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
                LLCF.CATE.SE <- sqrt(LLCF.CATE.old$variance.estimates))
                LLCF.CATE <- predictions
              }
            }
        }
      
      # Test of heterogeneity using Differential ATE
        LLCF.ATE.charact <- average_treatment_effect(LLCF, target.sample = "all", subset = (characteristics == 1))
        LLCF.ATE.not.charact <- average_treatment_effect(LLCF, target.sample = "all", subset = !(characteristics == 1))
        DiffATE.LLCF.mean <- LLCF.ATE.charact[1] - LLCF.ATE.not.charact[1]
        DiffATE.LLCF.SE <- sqrt(LLCF.ATE.charact[2]^2 + LLCF.ATE.not.charact[2]^2)
        lower.DiffATE.LLCF <- (DiffATE.LLCF.mean - (qnorm(0.975) * DiffATE.LLCF.SE))
        upper.DiffATE.LLCF <- (DiffATE.LLCF.mean + (qnorm(0.975) * DiffATE.LLCF.SE))

        results_BLP <- data.frame(t(c(mean.forest.pred.GRF, diff.forest.pred.GRF, 
                         mean.forest.pred.CR.GRF, diff.forest.pred.CR.GRF, 
                         mean.forest.pred.LLCF, diff.forest.pred.LLCF)))
    
        results_DiffATE <- data.frame(t(c(paste(round(DiffATE.GRF.mean, 3), "(", round(DiffATE.GRF.SE, 3), ")"),
                             paste(round(DiffATE.CR.GRF.mean, 3), "(", round(DiffATE.CR.GRF.SE, 3), ")"),
                             paste(round(DiffATE.LLCF.mean, 3), "(", round(DiffATE.LLCF.SE, 3), ")")), (845 - sum(missingvalues))))   
    
        data.frame(cbind(results_ATE, results_AUTOC, results_BLP, results_DiffATE)
    })
    results_BLP <- basic.results[,1:6]        
    colnames(results_BLP) <- c("BLP[1] GRF", "BLP[2] GRF", "BLP[1] CR.GRF", "BLP[2] CR.GRF", "BLP[1] LLCF", "BLP[2] LLCF")
    rownames(results_BLP) <- c("savings", "risk loving", "wage earner", "household member chronically ill", "impatient")   
    
    results_DiffATE <- basic.results[,7:10]                                    
    colnames(results_DiffATE) <- c("Differential ATE GRF", "Differential ATE CR.GRF", "Differential ATE LLCF", "Number of observations")
    rownames(results_DiffATE) <- c("savings", "risk loving", "wage earner", "household member chronically ill", "impatient")
    
    results = list("results_BLP" = results_BLP, 
                         "results_DiffATE" = results_DiffATE)                              
    return(results)
}

results = run_method(numtrees, index, lambdas, boolean.plot, boolean.lambdas)
results_table5 <- list('Sheet1' = results[["results_BLP"]], 'Sheet2' = results[["results_DiffATE"]])
write.xlsx(results_table5, file = "Table 5 (Field et al., 2013).xlsx")
