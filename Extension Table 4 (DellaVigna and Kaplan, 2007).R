library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(haven)
library(glmnet)
set.seed(123)

# Initialize parameters
  numtrees <- 2000 # Set to 1000 or 5000 to perform sensitivity analysis.
  index <- c(1:12) # Concerns 12 dependent variables; do not adjust.
  lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5) # Concerns ridge penalty parameters; do not adjust.
  boolean.lambdas <- FALSE # Set to TRUE to use lambdas instead of automatic penalty tuning.
  boolean.plot <- FALSE # Set to TRUE to make various plots of interest.
  filename.plot.GRF.CATE <- "GRF CATE .pdf"
  filename.plot.CR.GRF.CATE <- "CR.GRF CATE .pdf"
  filename.plot.LLCF.CATE <- "LLCF CATE .pdf"

# Estimation procedure

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
        GRF.varimp.ordered <- order(GRF.varimp)
        GRF.mostimportant <- colnames(X)[(GRF.varimp.ordered[1:4])] # 4 most important variables for splitting
        
      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, honesty = TRUE)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)
  
      # Implement GRF
        GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, 
                             honesty = TRUE, tune.parameters = "all")
    
      # Compute ATE 
        GRF.ATE <- average_treatment_effect(GRF, target.sample = "all")
    
      # Compute HTE with corresponding 95% confidence intervals
        GRF.pred <- predict(GRF, estimate.variance = TRUE)
        GRF.CATE <- GRF.pred$predictions
        GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
        lower.GRF <- GRF.CATE - qnorm(0.975)*GRF.CATE.SE
        upper.GRF <- GRF.CATE + qnorm(0.975)*GRF.CATE.SE
    
      # Graph the predicted GRF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(GRF.CATE)
        }
    
      # See if the GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        GRF.rate <- rank_average_treatment_effect(GRF, GRF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(GRF.rate)
        }

      # Assessing GRF fit and heterogeneity using the Best Linear Predictor Approach [BLPredictor]
        GRF.BLPredictor <- test_calibration(GRF)
                                         
      # Assessing general GRF heterogeneity using Differential ATE Approach
        GRF.high.effect <- (GRF.CATE > median(GRF.CATE))
        GRF.ATE.high <- average_treatment_effect(GRF, subset = GRF.high.effect)
        GRF.ATE.low <- average_treatment_effect(GRF, subset = !GRF.high.effect)
        DiffATE.GRF.mean <- GRF.ATE.high[1] - GRF.ATE.low[1]
        DiffATE.GRF.SE <- sqrt(GRF.ATE.high[2]^2 + GRF.ATE.low[2]^2)
        lower.DiffATE.GRF <- (DiffATE.GRF.mean - (qnorm(0.975) * DiffATE.GRF.SE))
        upper.DiffATE.GRF <- (DiffATE.GRF.mean + (qnorm(0.975) * DiffATE.GRF.SE))
    
      # Test of heterogeneity using Differential ATE along each and every variable
        combined.X.median <- apply(combined.X, 2, median)
        results_DiffATE_GRF = sapply(c(1:ncol(combined.X)), function(k) {
          GRF.ATE.abovemedian <- average_treatment_effect(GRF, target.sample = "all", subset = (combined.X[,k] >= combined.X.median[k]))
          GRF.ATE.belowmedian <- average_treatment_effect(GRF, target.sample = "all", subset = (combined.X[,k] < combined.X.median[k]))
          GRF.AIPW.abovemedian <- get_scores(GRF, subset = (combined.X[,k] >= combined.X.median[k]))
          GRF.AIPW.belowmedian <- get_scores(GRF, subset = (combined.X[,k] < combined.X.median[k]))
          DiffATE.GRF.test <- t.test(GRF.AIPW.abovemedian, GRF.AIPW.belowmedian, alternative = "two.sided", var.equal = FALSE)
          data.frame(t(c(paste(round(GRF.ATE.belowmedian[1], 3), "(", round(GRF.ATE.belowmedian[2], 3), ")"),
                         paste(round(GRF.ATE.abovemedian[1], 3), "(", round(GRF.ATE.abovemedian[2], 3), ")"),
                         paste(round(DiffATE.GRF.test$p.value, 3)))))
        }
        colnames(results_DiffATE_GRF) <- c("GRF CATE below median", "GRF CATE above median", "GRF p-value difference")
        rownames(results_DiffATE_GRF) <- colnames(combined.X)                           
        sign.var.DiffATE_GRF <- rownames((results_DiffATE_GRF[,3] < 0.10))
                                     
      # Assessing GRF fit and heterogeneity using the Best Linear Projection Approach [BLProjection]
        full.GRF.BLProjection <- best_linear_projection(GRF, X)
        mostimportant.GRF.BLProjection <- best_linear_projection(GRF, X[,GRF.mostimportant])       
        sign.var.DiffATE_GRF.BLProjection <- best_linear_projection(GRF, X[,sign.var.DiffATE_GRF])                           
        BLProjection_GRF <- texreg(list(full.GRF.BLProjection, mostimportant.GRF.BLProjection, sign.var.DiffATE_GRF.BLProjection),
          custom.model.names = c("All covariates", "Most important covariates", "Significant Differential ATE covariates"),
          table = FALSE,
          use.packages = FALSE,
          dcolumn = TRUE,
          single.row = TRUE,
          custom.coef.names = colnames(combined.X)
         )
                                   
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
         pdf("GRF CATE plots.pdf")
        if (boolean.plot == TRUE) {
          for (k in 1:length(sign.var.DiffATE_GRF)) {
            # Set all variables at their median values
              X.median <- apply(X, 2, median)
              
            # Create ordered vector of important variable
              important.var.test = seq(min(X$sign.var.DiffATE_GRF[k]), max(X$sign.var.DiffATE_GRF[k]))
              
            # Create test set
              X.test <- matrix(rep(X.median, length(important.var.test)), length(important.var.test), byrow = TRUE)
              X.test[,(sign.var.DiffATE_GRF[k])] = important.var.test
              
            # Predict new CATE estimates
              GRF.pred.test <- predict(GRF, X.test, estimate.variance = TRUE)
              GRF.CATE.test <- GRF.pred.test$predictions
              GRF.CATE.SE.test <- sqrt(GRF.pred.test$variance.estimates)
              lower.GRF.test <- GRF.CATE.test - qnorm(0.975)*GRF.CATE.SE.test
              upper.GRF.test <- GRF.CATE.test + qnorm(0.975)*GRF.CATE.SE.test
            
            # Make plot and save
              plot(X.test[,(sign.var.DiffATE_GRF[k])], GRF.CATE.test, type = "l", ylim = range(min(lower.GRF.test), max(upper.GRF.test)), xlab = sign.var.DiffATE_GRF[k], ylab = "CATE")
              lines(X.test[,(sign.var.DiffATE_GRF[k])], upper.GRF.test, col = 1, lty = 2)
              lines(X.test[,(sign.var.DiffATE_GRF[k])], lower.GRF.test, col = 1, lty = 2)
              grid()
          }
         }
         dev.off()
    ############################
    #### Cluster-Robust GRF ####
    ############################
      
      # Select variables to include using preliminary Cluster-Robust GRF
        prelim.CR.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, clusters = district.clusters, num.trees = numtrees)
        prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
        selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)
  
      # Compute the variable importance
        CR.GRF.varimp <- variable_importance(forest.Y) 
        CR.GRF.varimp.ordered <- order(CR.GRF.varimp)                             
        CR.GRF.mostimportant <- colnames(X)[CR.GRF.varimp.ordered[1:4]] # 4 most important variables for splitting
        
      # Implement Cluster-Robust GRF
        CR.GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, clusters = district.clusters, honesty = TRUE, num.trees = numtrees, 
                                tune.parameters = "all")
                                     
      # Compute ATE 
        CR.GRF.ATE <- average_treatment_effect(CR.GRF, target.sample = "all")
                                     
      # Compute HTE with corresponding 95% confidence intervals                               
        CR.GRF.pred <- predict(CR.GRF, estimate.variance = TRUE)
        CR.GRF.CATE <- CR.GRF.pred$predictions
        CR.GRF.CATE.SE <- sqrt(CR.GRF.pred$variance.estimates)
        lower.CR.GRF <- CR.GRF.CATE - qnorm(0.975)*CR.GRF.CATE.SE
        upper.CR.GRF <- CR.GRF.CATE + qnorm(0.975)*CR.GRF.CATE.SE
              
      # Graph the predicted Cluster-Robust GRF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(CR.GRF.CATE)
        }
    
      # See if the Cluster-Robust GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        CR.GRF.rate <- rank_average_treatment_effect(CR.GRF, CR.GRF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(CR.GRF.rate)
        }
  
      # Assessing Cluster-Robust GRF fit using the Best Linear Predictor Approach [BLPredictor]
        CR.GRF.BLPredictor <- test_calibration(CR.GRF)
    
      # Assessing general CR.GRF heterogeneity using Differential ATE Approach
        CR.GRF.high.effect <- (CR.GRF.CATE > median(CR.GRF.CATE))
        CR.GRF.ATE.high <- average_treatment_effect(CR.GRF, subset = CR.GRF.high.effect)
        CR.GRF.ATE.low <- average_treatment_effect(CR.GRF, subset = !CR.GRF.high.effect)
        DiffATE.CR.GRF.mean <- CR.GRF.ATE.high[1] - CR.GRF.ATE.low[1]
        DiffATE.CR.GRF.SE <- sqrt(CR.GRF.ATE.high[2]^2 + CR.GRF.ATE.low[2]^2)
        lower.DiffATE.CR.GRF <- (DiffATE.CR.GRF.mean - (qnorm(0.975) * DiffATE.CR.GRF.SE))
        upper.DiffATE.CR.GRF <- (DiffATE.CR.GRF.mean + (qnorm(0.975) * DiffATE.CR.GRF.SE))
                                   
      # Test of heterogeneity using Differential ATE along each and every variable
        combined.X.median <- apply(combined.X, 2, median)
        results_DiffATE_CR.GRF = sapply(c(1:ncol(combined.X)), function(k) {
          CR.GRF.ATE.abovemedian <- average_treatment_effect(CR.GRF, target.sample = "all", subset = (combined.X[,k] >= combined.X.median[k]))
          CR.GRF.ATE.belowmedian <- average_treatment_effect(CR.GRF, target.sample = "all", subset = (combined.X[,k] < combined.X.median[k]))
          CR.GRF.AIPW.abovemedian <- get_scores(CR.GRF, subset = (combined.X[,k] >= combined.X.median[k]))
          CR.GRF.AIPW.belowmedian <- get_scores(CR.GRF, subset = (combined.X[,k] < combined.X.median[k]))
          DiffATE.CR.GRF.test <- t.test(CR.GRF.AIPW.abovemedian, CR.GRF.AIPW.belowmedian, alternative = "two.sided", var.equal = FALSE)
          data.frame(t(c(paste(round(CR.GRF.ATE.belowmedian[1], 3), "(", round(CR.GRF.ATE.belowmedian[2], 3), ")"),
                         paste(round(CR.GRF.ATE.abovemedian[1], 3), "(", round(CR.GRF.ATE.abovemedian[2], 3), ")"),
                         paste(round(DiffATE.CR.GRF.test$p.value, 3)))))
        }
        colnames(results_DiffATE_CR.GRF) <- c("CR.GRF CATE below median", "CR.GRF CATE above median", "CR.GRF p-value difference")
        rownames(results_DiffATE_CR.GRF) <- colnames(combined.X)                
        sign.var.DiffATE_CR.GRF <- rownames((results_DiffATE_CR.GRF[,3] < 0.10))
      
      # Assessing GRF fit and heterogeneity using the Best Linear Projection Approach [BLProjection]
        full.CR.GRF.BLProjection <- best_linear_projection(CR.GRF, X)
        mostimportant.CR.GRF.BLProjection <- best_linear_projection(CR.GRF, X[,CR.GRF.mostimportant])       
        sign.var.DiffATE_CR.GRF.BLProjection <- best_linear_projection(CR.GRF, X[,sign.var.DiffATE_CR.GRF])                           
        BLProjection_CR.GRF <- texreg(list(full.CR.GRF.BLProjection, mostimportant.CR.GRF.BLProjection, CR.sign.var.DiffATE_GRF.BLProjection),
          custom.model.names = c("All covariates", "Most important covariates", "Significant Differential ATE covariates"),
          table = FALSE,
          use.packages = FALSE,
          dcolumn = TRUE,
          single.row = TRUE,
          custom.coef.names = colnames(combined.X)
         )
                                        
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
        pdf("CR.GRF CATE plots.pdf")
        if (boolean.plot == TRUE) {
          for (k in 1:length(sign.var.DiffATE_CR.GRF)) {
            # Set all variables at their median values
              X.median <- apply(X, 2, median)
              
            # Create ordered vector of important variable
              important.var.test = seq(min(X$sign.var.DiffATE_CR.GRF[k]), max(X$sign.var.DiffATE_CR.GRF[k]))
              
            # Create test set
              X.test <- matrix(rep(X.median, length(important.var.test)), length(important.var.test), byrow = TRUE)
              X.test[,(sign.var.DiffATE_CR.GRF[k])] = important.var.test
              
            # Predict new CATE estimates
              CR.GRF.pred.test <- predict(CR.GRF, X.test, estimate.variance = TRUE)
              CR.GRF.CATE.test <- CR.GRF.pred.test$predictions
              CR.GRF.CATE.SE.test <- sqrt(CR.GRF.pred.test$variance.estimates)
              lower.CR.GRF.test <- CR.GRF.CATE.test - qnorm(0.975)*CR.GRF.CATE.SE.test
              upper.CR.GRF.test <- CR.GRF.CATE.test + qnorm(0.975)*CR.GRF.CATE.SE.test
            
            # Make plot and save
              plot(X.test[,(sign.var.DiffATE_CR.GRF[k])], CR.GRF.CATE.test, type = "l", ylim = range(min(lower.CR.GRF.test), max(upper.CR.GRF.test)), xlab = sign.var.DiffATE_CR.GRF[k], ylab = "CATE")
              lines(X.test[,(sign.var.DiffATE_CR.GRF[k])], upper.CR.GRF.test, col = 1, lty = 2)
              lines(X.test[,(sign.var.DiffATE_CR.GRF[k])], lower.CR.GRF.test, col = 1, lty = 2)
              grid()        
          }
         }
         dev.off()
    ############################
    ########### LLCF ###########
    ############################
    
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- ll_regression_forest(X, W, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                         num.trees = numtrees)
        W.hat <- predict(forest.W)$predictions
        forest.Y <- ll_regression_forest(X, Y, honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                         num.trees = numtrees)
        Y.hat <- predict(forest.Y)$predictions
    
      # Compute the variable importance
        LLCF.varimp <- variable_importance(forest.Y) 
        LLCF.varimp.ordered <- order(LLCF.varimp)                                
        LLCF.mostimportant <- colnames(X)[LLCF.varimp.ordered[1:4]] # 4 most important variables for splitting
        
      # Select variables to include using Lasso feature selection
        lasso.mod <- cv.glmnet(X, Y, alpha = 1)
        selected <- which(coef(lasso.mod) != 0)
        if(length(selected) < 2) {
          selected <- 1:ncol(X)
        } else {
          selected <- selected[-1] - 1 # Remove intercept
        }
  
      # Implement LLCF
        LLCF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE,  
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
      
      # Graph the predicted Cluster-Robust GRF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(LLCF.CATE)
        }
    
      # See if the Cluster-Robust GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        LLCF.rate <- rank_average_treatment_effect(LLCF, LLCF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(LLCF.rate)
        }
  
      # Assessing Cluster-Robust GRF fit using the Best Linear Predictor Approach [BLPredictor]
        LLCF.BLPredictor <- test_calibration(LLCF)
      
      # Assessing general LLCF heterogeneity using Differential ATE Approach
        LLCF.high.effect <- (LLCF.CATE > median(LLCF.CATE))
        LLCF.ATE.high <- average_treatment_effect(LLCF, subset = LLCF.high.effect)
        LLCF.ATE.low <- average_treatment_effect(LLCF, subset = !LLCF.high.effect)
        DiffATE.LLCF.mean <- LLCF.ATE.high[1] - LLCF.ATE.low[1]
        DiffATE.LLCF.SE <- sqrt(LLCF.ATE.high[2]^2 + LLCF.ATE.low[2]^2)
        lower.DiffATE.LLCF <- (DiffATE.LLCF.mean - (qnorm(0.975) * DiffATE.LLCF.SE))
        upper.DiffATE.LLCF <- (DiffATE.LLCF.mean + (qnorm(0.975) * DiffATE.LLCF.SE))
                                                  
      # Test of heterogeneity using Differential ATE
        combined.X.median <- apply(combined.X, 2, median)
        results_DiffATE_LLCF = sapply(c(1:ncol(combined.X)), function(k) {
          LLCF.ATE.abovemedian <- average_treatment_effect(LLCF, target.sample = "all", subset = (combined.X[,k] >= combined.X.median[k]))
          LLCF.ATE.belowmedian <- average_treatment_effect(LLCF, target.sample = "all", subset = (combined.X[,k] < combined.X.median[k]))
          LLCF.AIPW.abovemedian <- get_scores(LLCF, subset = (combined.X[,k] >= combined.X.median[k]))
          LLCF.AIPW.belowmedian <- get_scores(LLCF, subset = (combined.X[,k] < combined.X.median[k]))
          DiffATE.LLCF.test <- t.test(LLCF.AIPW.abovemedian, LLCF.AIPW.belowmedian, alternative = "two.sided", var.equal = FALSE)
          data.frame(t(c(paste(round(LLCF.ATE.belowmedian[1], 3), "(", round(LLCF.ATE.belowmedian[2], 3), ")"),
                         paste(round(LLCF.ATE.abovemedian[1], 3), "(", round(LLCF.ATE.abovemedian[2], 3), ")"),
                         paste(round(DiffATE.LLCF.test$p.value, 3)))))
        }
        colnames(results_DiffATE_LLCF) <- c("LLCF CATE below median", "LLCF CATE above median", "LLCF p-value difference")
        rownames(results_DiffATE_LLCF) <- colnames(combined.X)    
        sign.var.DiffATE_LLCF <- rownames((results_DiffATE_LLCF[,3] < 0.10))
      
      # Assessing LLCF fit and heterogeneity using the Best Linear Projection Approach [BLProjection]
        full.LLCF.BLProjection <- best_linear_projection(LLCF, X)
        mostimportant.LLCF.BLProjection <- best_linear_projection(LLCF, X[,LLCF.mostimportant])       
        sign.var.DiffATE_LLCF.BLProjection <- best_linear_projection(LLCF, X[,sign.var.DiffATE_LLCF])                           
        BLProjection_LLCF <- texreg(list(full.LLCF.BLProjection, mostimportant.LLCF.BLProjection, sign.var.DiffATE_LLCF.BLProjection),
          custom.model.names = c("All covariates", "Most important covariates", "Significant Differential ATE covariates"),
          table = FALSE,
          use.packages = FALSE,
          dcolumn = TRUE,
          single.row = TRUE,
          custom.coef.names = colnames(combined.X)
         )
                                      
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
        pdf("LLCF CATE plots.pdf")
        if (boolean.plot == TRUE) {
          for (k in 1:length(sign.var.DiffATE_LLCF)) {
            # Set all variables at their median values
              X.median <- apply(X, 2, median)
              
            # Create ordered vector of important variable
              important.var.test = seq(min(X$sign.var.DiffATE_LLCF[k]), max(X$sign.var.DiffATE_LLCF[k]))
              
            # Create test set
              X.test <- matrix(rep(X.median, length(important.var.test)), length(important.var.test), byrow = TRUE)
              X.test[,(sign.var.DiffATE_LLCF[k])] = important.var.test
              
            # Predict new CATE estimates
              LLCF.pred.test <- predict(LLCF, X.test, estimate.variance = TRUE)
              LLCF.CATE.test <- LLCF.pred.test$predictions
              LLCF.CATE.SE.test <- sqrt(LLCF.pred.test$variance.estimates)
              lower.LLCF.test <- LLCF.CATE.test - qnorm(0.975)*LLCF.CATE.SE.test
              upper.LLCF.test <- LLCF.CATE.test + qnorm(0.975)*LLCF.CATE.SE.test
            
            # Make plot and save
              plot(X.test[,(sign.var.DiffATE_LLCF[k])], LLCF.CATE.test, type = "l", ylim = range(min(lower.LLCF.test), max(upper.LLCF.test)), xlab = sign.var.DiffATE_LLCF[k], ylab = "CATE")
              lines(X.test[,(sign.var.DiffATE_LLCF[k])], upper.LLCF.test, col = 1, lty = 2)
              lines(X.test[,(sign.var.DiffATE_LLCF[k])], lower.LLCF.test, col = 1, lty = 2)
              grid()
          }
         }
         dev.off()
                                      
    ############################
    ###### Update results ######
    ############################
                                      
       results_BLPredictor <-  data.frame(t(c(paste(round(GRF.BLPredictor[1,1], 3), "(", round(GRF.BLPredictor[1,2], 3), ")"),
                                      paste(round(GRF.BLPredictor[2,1], 3), "(", round(GRF.BLPredictor[2,2], 3), ")"),
                                      paste(round(CR.GRF.BLPredictor[1,1], 3), "(", round(CR.GRF.BLPredictor[1,2], 3), ")"),
                                      paste(round(CR.GRF.BLPredictor[2,1], 3), "(", round(CR.GRF.BLPredictor[2,2], 3), ")"),
                                      paste(round(LLCF.BLPredictor[1,1], 3), "(", round(LLCF.BLPredictor[1,2], 3), ")"),
                                      paste(round(LLCF.BLPredictor[2,1], 3), "(", round(LLCF.BLPredictor[2,2], 3), ")"))))
       colnames(results_BLPredictor) <- c("BLPredictor[1] GRF", "BLPredictor[2] GRF", 
                                  "BLPredictor[1] CR.GRF", "BLPredictor[2] CR.GRF", 
                                  "BLPredictor[1] LLCF", "BLPredictor[2] LLCF")
                                      
       results_DiffATE_General <- data.frame(t(c(paste(round(GRF.ATE.high[1], 3), "(", round(GRF.ATE.high[2], 3), ")"),
                                          paste(round(GRF.ATE.low[1], 3), "(", round(GRF.ATE.low[2], 3), ")"),
                                          paste(round(DiffATE.GRF.mean, 3), "(", round(DiffATE.GRF.SE, 3), ")"),
                                          paste(round(CR.GRF.ATE.high[1], 3), "(", round(CR.GRF.ATE.high[2], 3), ")"),
                                          paste(round(CR.GRF.ATE.low[1], 3), "(", round(CR.GRF.ATE.low[2], 3), ")"),
                                          paste(round(DiffATE.CR.GRF.mean, 3), "(", round(DiffATE.CR.GRF.SE, 3), ")"),
                                          paste(round(LLCF.ATE.high[1], 3), "(", round(LLCF.ATE.high[2], 3), ")"),
                                          paste(round(LLCF.ATE.low[1], 3), "(", round(LLCF.ATE.low[2], 3), ")"),
                                          paste(round(DiffATE.LLCF.mean, 3), "(", round(DiffATE.LLCF.SE, 3), ")"))))                          
       colnames(results_DiffATE) <- c("GRF CATE above median", "GRF CATE below median", "GRF CATE mean difference", 
                                      "CR.GRF CATE above median", "CR.GRF CATE below median", "CR.GRF CATE mean difference",
                                      "LLCF CATE above median", "LLCF CATE below median", "LLCF CATE mean difference")
                                                                                     
fun_insert <- function(x, pos, insert) {
  gsub(paste0("^(.{", pos, "})(.*)$"),
       paste0("\\1", insert, "\\2"),
       x)
}

                                      
results_table4 <- list('Sheet1' = results_DiffATE_GRF, 
                       'Sheet2' = results_DiffATE_CR.GRF,
                       'Sheet3' = results_DiffATE_LLCF,
                       'Sheet4' = results_DiffATE_General,
                       'Sheet5' = GRF.varimp.ordered,
                       'Sheet6' = CR.GRF.varimp.ordered,
                       'Sheet7' = LLCF.varimp.ordered,
                       'Sheet8' = results_BLPredictor,
                       'Sheet9' = BLProjection_GRF,
                       'Sheet10' = BLProjection_CR.GRF,
                       'Sheet11' = BLProjection_LLCF)
                                                              
write.xlsx(results_table4, file = "Table 4 (DellaVigna and Kaplan, 2007).xlsx")

               
