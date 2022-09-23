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

# Initialize parameters
  numtrees <- 2000
  index <- c(1:12)
  lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5)
  boolean.lambdas <- FALSE
  boolean.plot <- FALSE

# Estimation procedure
run_method = function(index, lambdas, boolean.plot, boolean.lambdas) {
  basic.results = sapply(index, function(i) {
    Y <- Grace_Period_Data[,(37+i)]
    Y <- as.vector(Y)
  
    ###########################
    ########### GRF ###########
    ###########################
    
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- regression_forest(X, W, num.trees = numtrees, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- regression_forest(X, Y, num.trees = numtrees, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions

      # Compute the variable importance
        GRF.varimp <- variable_importance(forest.Y) 
    
      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)

      # Implement GRF
        GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, tune.parameters = "all")
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
              
      # Assessing GRF heterogeneity using Differential ATE Approach
        GRF.high.effect <- (GRF.CATE > median(GRF.CATE))
        GRF.ATE.high <- average_treatment_effect(GRF, subset = GRF.high.effect)
        GRF.ATE.low <- average_treatment_effect(GRF, subset = !GRF.high.effect)
        DiffATE.GRF.mean <- GRF.ATE.high[1] - GRF.ATE.low[1]
        DiffATE.GRF.SE <- sqrt(GRF.ATE.high[2]^2 + GRF.ATE.low[2]^2)
        lower.DiffATE.GRF <- (DiffATE.GRF.mean - (qnorm(0.975) * DiffATE.GRF.SE))
        upper.DiffATE.GRF <- (DiffATE.GRF.mean + (qnorm(0.975) * DiffATE.GRF.SE))
              
  
    ############################
    #### Cluster-Robust GRF ####
    ############################
    
      # Select variables to include using preliminary Cluster-Robust GRF
        prelim.CR.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, num.trees = numtrees)
        prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
        selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)

      # Compute the variable importance
        CR.GRF.varimp <- variable_importance(forest.Y) 
              
      # Implement Cluster-Robust GRF
        CR.GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = loangroups, num.trees = numtrees, tune.parameters = "all")
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
        resultsTable1OriginalPaper[(i+1),3] <- paste(round(CR.GRF.ATE[1], 3), "(", round(CR.GRF.ATE[2], 3), ")")
    
      # See if the Cluster-Robust GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        CR.GRF.rate <- rank_average_treatment_effect(CR.GRF, CR.GRF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(CR.GRF.rate)
        }
  
      # Assessing Cluster-Robust GRF fit using the Best Linear Predictor Approach [BLP]
        CR.GRF.BLP <- test_calibration(CR.GRF)
        mean.forest.pred.CR.GRF <- c(CR.GRF.BLP[1,1], CR.GRF.BLP[1,2], (CR.GRF.BLP[1,4] < 0.10))
        diff.forest.pred.CR.GRF <- c(CR.GRF.BLP[2,1], CR.GRF.BLP[2,2], (CR.GRF.BLP[2,4] < 0.10))
              
      # Assessing Cluster-Robust GRF heterogeneity using Differential ATE Approach
        CR.GRF.high.effect <- (CR.GRF.CATE > median(CR.GRF.CATE))
        CR.GRF.ATE.high <- average_treatment_effect(CR.GRF, subset = CR.GRF.high.effect)
        CR.GRF.ATE.low <- average_treatment_effect(CR.GRF, subset = !CR.GRF.high.effect)
        DiffATE.CR.GRF.mean <- CR.GRF.ATE.high[1] - CR.GRF.ATE.low[1]
        DiffATE.CR.GRF.SE <- sqrt(CR.GRF.ATE.high[2]^2 + CR.GRF.ATE.low[2]^2)
        lower.DiffATE.CR.GRF <- (DiffATE.CR.GRF.mean - (qnorm(0.975) * DiffATE.CR.GRF.SE))
        upper.DiffATE.CR.GRF <- (DiffATE.CR.GRF.mean + (qnorm(0.975) * DiffATE.CR.GRF.SE))
              
    ############################
    ########### LLCF ###########
    ############################
    
       # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- ll_regression_forest(X, W, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = numtrees, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- ll_regression_forest(X, Y, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = numtrees, tune.parameters = "all")
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
        LLCF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, enable.ll.split = TRUE, 
                              ll.split.weight.penalty = TRUE, num.trees = numtrees, tune.parameters = "all")
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
                predictions.new <- predictions
              }
            }
       
      # Find lower and upper bounds for 95% confidence intervals
        lower.LLCF <- LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE
        upper.LLCF <- LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE
    
      # Graph the predicted LLCF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(LLCF.CATE)
        }
              
      # See if the LLCF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        LLCF.rate <- rank_average_treatment_effect(LLCF, LLCF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(LLCF.rate)
        }
              
      # Assessing LLCF fit using the Best Linear Predictor Approach [BLP]
        LLCF.BLP <- test_calibration(LLCF)
        mean.forest.pred.LLCF <- c(LLCF.BLP[1,1], LLCF.BLP[1,2], (LLCF.BLP[1,4] < 0.10))
        diff.forest.pred.LLCF <- c(LLCF.BLP[2,1], LLCF.BLP[2,2], (LLCF.BLP[2,4] < 0.10))
      
      # Assessing LLCF heterogeneity using Differential ATE Approach
        LLCF.high.effect <- (LLCF.CATE > median(LLCF.CATE))
        LLCF.ATE.high <- average_treatment_effect(LLCF, subset = LLCF.high.effect)
        LLCF.ATE.low <- average_treatment_effect(LLCF, subset = !LLCF.high.effect)
        DiffATE.LLCF.mean <- LLCF.ATE.high[1] - LLCF.ATE.low[1]
        DiffATE.LLCF.SE <- sqrt(LLCF.ATE.high[2]^2 + LLCF.ATE.low[2]^2)
        lower.DiffATE.LLCF <- (DiffATE.LLCF.mean - (qnorm(0.975) * DiffATE.LLCF.SE))
        upper.DiffATE.LLCF <- (DiffATE.LLCF.mean + (qnorm(0.975) * DiffATE.LLCF.SE))
              
        results_ATE <- data.frame(t(c(paste(round(GRF.ATE[1], 3), "(", round(GRF.ATE[2], 3), ")"), 
                         paste(round(CR.GRF.ATE[1], 3), "(", round(CR.GRF.ATE[2], 3), ")"),
                         paste(round(LLCF.ATE[1], 3), "(", round(LLCF.ATE[2], 3), ")"))))
    
        results_AUTOC <- data.frame(t(c(paste(round(GRF.rate$estimate, 2), "+/", round(qnorm(0.975) * GRF.rate$std.err, 2),
                           paste(round(CR.GRF.rate$estimate, 2), "+/", round(qnorm(0.975) * CR.GRF.rate$std.err, 2),
                           paste(round(LLCF.rate$estimate, 2), "+/", round(qnorm(0.975) * LLCF.rate$std.err, 2))))
                                        
        results_BLP <- data.frame(t(c(mean.forest.pred.GRF, diff.forest.pred.GRF, 
                         mean.forest.pred.CR.GRF, diff.forest.pred.CR.GRF, 
                         mean.forest.pred.LLCF, diff.forest.pred.LLCF)))
                                        
        results_DiffATE <- data.frame(t(c(paste(round(DiffATE.GRF.mean, 3), "(", round(DiffATE.GRF.SE, 3), ")"),
                             paste(round(DiffATE.CR.GRF.mean, 3), "(", round(DiffATE.CR.GRF.SE, 3), ")"),
                             paste(round(DiffATE.LLCF.mean, 3), "(", round(DiffATE.LLCF.SE, 3), ")"))))                                
        
        data.frame(cbind(results_ATE, results_AUTOC, results_BLP, results_DiffATE)
    })
    results_ATE <- basic.results[,1:3]
    colnames(results_ATE) <- c("ATE GRF", "ATE CR.GRF", "ATE LLCF")    
    rownames(results_ATE) <- c("Total business spending", "Inventory and raw materials", "Business equipment",
                                        "Operating costs", "Total non-business spending", "Home repairs", "Utilities, taxes and rent",
                                        "Human capital", "Money for relending", "Savings", "Food and durable consumption", "New businesses")
                                        
    results_AUTOC <- basic.results[,4:6]                                    
    colnames(results_AUTOC) <- c("AUTOC GRF", "AUTOC CR.GRF", "AUTOC LLCF")    
    rownames(results_AUTOC) <- c("Total business spending", "Inventory and raw materials", "Business equipment",
                                        "Operating costs", "Total non-business spending", "Home repairs", "Utilities, taxes and rent",
                                        "Human capital", "Money for relending", "Savings", "Food and durable consumption", "New businesses")
    
    results_BLP <- basic.results[,7:12]        
    colnames(results_BLP) <- c("BLP[1] GRF", "BLP[2] GRF", "BLP[1] CR.GRF", "BLP[2] CR.GRF", "BLP[1] LLCF", "BLP[2] LLCF")
    rownames(results_BLP) <- c("Total business spending", "Inventory and raw materials", "Business equipment",
                                        "Operating costs", "Total non-business spending", "Home repairs", "Utilities, taxes and rent",
                                        "Human capital", "Money for relending", "Savings", "Food and durable consumption", "New businesses")
                                        
    results_DiffATE <- basic.results[,13:15]                                    
    colnames(results_DiffATE) <- c("Differential ATE GRF", "Differential ATE CR.GRF", "Differential ATE LLCF")
    rownames(results_DiffATE) <- c("Total business spending", "Inventory and raw materials", "Business equipment",
                                        "Operating costs", "Total non-business spending", "Home repairs", "Utilities, taxes and rent",
                                        "Human capital", "Money for relending", "Savings", "Food and durable consumption", "New businesses")                                    
                                                                            
    results = list("results_ATE" = results_ATE, 
                         "results_AUTOC" = results_AUTOC, 
                         "results_BLP" = results_BLP, 
                         "results_DiffATE" = results_DiffATE)                              
    return(results)
 }

results = run_method(index, lambdas, boolean.plot, boolean.lambdas)
results_table1 <- list('Sheet1' = results[["results_ATE"]], 'Sheet2' = results[["results_AUTOC"]], 
                       'Sheet3' = results[["results_BLP"]], 'Sheet4' = results[["results_DiffATE"]])
write.xlsx(results_table1, file = "Table 1 (Field et al., 2013).xlsx")
