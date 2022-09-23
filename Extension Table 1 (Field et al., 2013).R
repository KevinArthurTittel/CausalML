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
  
    ###########################
    ########### GRF ###########
    ###########################
    
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- regression_forest(X, W, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- regression_forest(X, Y, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions

      # Compute the variable importance
        GRF.varimp <- variable_importance(forest.Y) 
    
      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = 2000)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)

      # Implement GRF
        GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = 2000, tune.parameters = "all")
        GRF.pred <- predict(GRF, estimate.variance = TRUE)
        GRF.CATE <- GRF.pred$predictions
        GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
  
      # Find lower and upper bounds for 95% confidence intervals
        lower.GRF <- GRF.CATE - 1.96*GRF.CATE.SE
        upper.GRF <- GRF.CATE + 1.96*GRF.CATE.SE
    
      # Graph the predicted GRF heterogeneous treatment effect estimates
        hist(GRF.CATE)

      # Compute ATE with corresponding 95% confidence intervals
        GRF.ATE <- average_treatment_effect(GRF, target.sample = "all")
        resultsTable1OriginalPaper[(i+1),2] <- paste(round(GRF.ATE[1], 3), "(", round(GRF.ATE[2], 3), ")")
    
      # See if the GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        GRF.rate <- rank_average_treatment_effect(GRF, GRF.CATE)
        plot(GRF.rate)
        paste("AUTOC:", round(GRF.rate$estimate, 2), "+/", round(1.96 * GRF.rate$std.err, 2)
  
    ############################
    #### Cluster-Robust GRF ####
    ############################
    
      # Select variables to include using preliminary Cluster-Robust GRF
        prelim.CR.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, 
                           num.trees = 4000)
        prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
        selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)

      # Compute the variable importance
        CR.GRF.varimp <- variable_importance(forest.Y) 
              
      # Implement Cluster-Robust GRF
        CR.GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, num.trees = 8000, tune.parameters = "all")
        CR.GRF.pred <- predict(CR.GRF, estimate.variance = TRUE)
        CR.GRF.CATE <- CR.GRF.pred$predictions
        CR.GRF.CATE.SE <- sqrt(CR.GRF.pred$variance.estimates)
    
      # Find lower and upper bounds for 95% confidence intervals
        lower.CR.GRF <- CR.GRF.CATE - 1.96*CR.GRF.CATE.SE
        upper.CR.GRF <- CR.GRF.CATE + 1.96*CR.GRF.CATE.SE
              
      # Graph the predicted Cluster-Robust GRF heterogeneous treatment effect estimates
        hist(CR.GRF.CATE)
    
      # Compute ATE with corresponding 95% confidence intervals
        CR.GRF.ATE <- average_treatment_effect(CR.GRF, target.sample = "all")
        resultsTable1OriginalPaper[(i+1),3] <- paste(round(CR.GRF.ATE[1], 3), "(", round(CR.GRF.ATE[2], 3), ")")
    
      # See if the Cluster-Robust GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        CR.GRF.rate <- rank_average_treatment_effect(CR.GRF, CR.GRF.CATE)
        plot(CR.GRF.rate)
        paste("AUTOC:", round(CR.GRF.rate$estimate, 2), "+/", round(1.96 * CR.GRF.rate$std.err, 2)
  
    ############################
    ########### LLCF ###########
    ############################
    
       # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- ll_regression_forest(X, W, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 2000, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- ll_regression_forest(X, Y, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 2000, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions
  
      # Compute the variable importance
        LLCF.varimp <- variable_importance(forest.Y) 
              
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
        LLCF.ATE <- average_treatment_effect(LLCF, target.sample = "all")
    
      # Predict: tuning without grid search over lambdas
        LLCF.pred <- predict(LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
        LLCF.CATE <- LLCF.pred$predictions
        LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
  
      # Predict: tuning done using set of lambdas
        LLCF.mse.old <- +Inf
        for (l in length(lambdas)) {
          LLCF.CATE.old <- predict(LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
          predictions <- LLCF.CATE.old$predictions
          LLCF.mse.new <- mean((predictions - mean(predictions))**2)
          if (LLCF.mse.new < LLCF.mse.old) {
            LLCF.mse.old <- LLCF.mse.new
            LLCF.CATE.SE <- sqrt(LLCF.CATE.old$variance.estimates))
            predictions.new <- predictions
          }
        }
       
      # Find lower and upper bounds for 95% confidence intervals
        lower.LLCF <- LLCF.CATE - 1.96*LLCF.CATE.SE
        upper.LLCF <- LLCF.CATE + 1.96*LLCF.CATE.SE
              
      # See if the LLCF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        LLCF.rate <- rank_average_treatment_effect(LLCF, LLCF.CATE)
        plot(LLCF.rate)
        paste("AUTOC:", round(LLCF.rate$estimate, 2), "+/", round(1.96 * LLCF.rate$std.err, 2)
    
        resultsTable1OriginalPaper[(i+1),4] <- paste(LLCF.ATE, "(", LLCF.ATE.SE, ")")
}
