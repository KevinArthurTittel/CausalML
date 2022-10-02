library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(haven)
library(glmnet)

# Import and prepare Fox News data set (DellaVigna and Kaplan, 2007)
  Fox_News_Data <- read_stata("C:/Users/481044kt/AppData/Local/Temp/2/FoxNewsFinalDataQJE_zOeIzx")
  Fox_News_Data <- as.data.frame(Fox_News_Data)

  # Appoint main outcome variable
    Y <- as.vector(Fox_News_Data$reppresfv2p00m96)

  # Appoint treatment assignment variable
    W <- as.vector(Fox_News_Data$foxnews2000)

  # Create a numerical vector of the character group name vector
    district.clusters <- as.vector(Fox_News_Data$diststate)

  # Appoint the control variables matrix
    demographic.controls.2000 <- Fox_News_Data[,c(37:49)]
    colnames(demographic.controls.2000) <- c("Population 2000", "Population (over 18) 2000", "Share with high school 2000", "Share with some college 2000",
                                         "Share with college degree 2000", "Share male 2000", "Share African Americans 2000", "Share Hispanics 2000",
                                         "Employment rate 2000", "Unemployment rate 2000", "Share married 2000", "Median income 2000", "Share urban 2000")

    demographic.controls.Diff19962000 <- Fox_News_Data[,c(139:150)]
    colnames(demographic.controls.Diff19962000) <- c("Population, Diff. btwn. 2000 and 1996", "Share with high school, Diff. btwn. 2000 and 1996", "Share with some college, Diff. btwn. 2000 and 1996",
                                                 "Share with college degree, Diff. btwn. 2000 and 1996", "Share male, Diff. btwn. 2000 and 1996", "Share African Americans, Diff. btwn. 2000 and 1996",
                                                 "Share Hispanics, Diff. btwn. 2000 and 1996", "Employment rate, Diff. btwn. 2000 and 1996", "Unemployment rate, Diff. btwn. 2000 and 1996",
                                                 "Share married, Diff. btwn. 2000 and 1996", "Median income, Diff. btwn. 2000 and 1996", "Share urban, Diff. btwn. 2000 and 1996")

    Decile1NumberChannels <- as.integer((rowSums(Fox_News_Data[,c(62:70)]) == 0))
    Decile1PotentialSubscribers <- as.integer(rowSums(Fox_News_Data[,c(71:79)]) == 0)

    cable.controls <- cbind(Decile1NumberChannels, Fox_News_Data[,c(62:70)], Decile1PotentialSubscribers, Fox_News_Data[,c(71:79)])
    colnames(cable.controls) <- c("Decile 1 no. cable channels available", "Decile 2 no. cable channels available", "Decile 3 no. cable channels available", 
                                  "Decile 4 no. cable channels available", "Decile 5 no. cable channels available", "Decile 6 no. cable channels available", 
                                  "Decile 7 no. cable channels available", "Decile 8 no. cable channels available", "Decile 9 no. cable channels available", 
                                  "Decile 10 no. cable channels available", "Decile 1 no. potential subscribers", "Decile 2 no. potential subscribers", 
                                  "Decile 3 no. potential subscribers", "Decile 4 no. potential subscribers", "Decile 5 no. potential subscribers", 
                                  "Decile 6 no. potential subscribers", "Decile 7 no. potential subscribers", "Decile 8 no. potential subscribers", 
                                  "Decile 9 no. potential subscribers", "Decile 10 no. potential subscribers")


    X <- as.matrix(cbind(demographic.controls.2000, demographic.controls.Diff19962000, cable.controls))

    remaining.variables <- Fox_News_Data[,c(14,59, 204, 205)]
    colnames(remaining.variables) <- c("No. cable channels available", "NumberPotentialSubscribers", "Swing district", "Republican district")

    combined.X <- as.matrix(cbind(X, remaining.variables))
    colnames(combined.X) <- c(colnames(X), colnames(remaining.variables))         

  # Clean the data
    indices <- (Fox_News_Data$sample12000 == 1)
    Y <- Y[indices]
    X <- X[indices,]
    W <- W[indices]
    district.clusters <- district.clusters[indices]
    remaining.variables <- remaining.variables[indices,]

# Initialize parameters
  numtrees <- 2000
  lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5)
  boolean.lambdas <- FALSE
  boolean.plot <- FALSE
  set.seed(123)
  filename.plot.GRF.CATE <- "GRF CATE .pdf"
  filename.plot.CR.GRF.CATE <- "CR.GRF CATE .pdf"
  filename.plot.LLCF.CATE <- "LLCF CATE .pdf"

# Estimation procedure
run_method = function(numtrees, index, lambdas, boolean.plot, boolean.lambdas) {
  basic.results = sapply(index, function(i) {
  
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
                                     
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
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
              pdf(fun_insert(x = filename.plot.GRF.CATE, pos = (nchar(filename.plot.GRF.CATE) - 4), insert = (sign.var.DiffATE_GRF[k])))
              plot(X.test[,(sign.var.DiffATE_GRF[k])], GRF.CATE.test, type = "l", ylim = range(min(lower.GRF.test), max(upper.GRF.test)), xlab = sign.var.DiffATE_GRF[k], ylab = "CATE")
              lines(X.test[,(sign.var.DiffATE_GRF[k])], upper.GRF.test, col = 1, lty = 2)
              lines(X.test[,(sign.var.DiffATE_GRF[k])], lower.GRF.test, col = 1, lty = 2)
              grid()
              dev.off()
          }
         }
        
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
                                             
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
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
              pdf(fun_insert(x = filename.plot.CR.GRF.CATE, pos = (nchar(filename.plot.CR.GRF.CATE) - 4), insert = (sign.var.DiffATE_CR.GRF[k])))
              plot(X.test[,(sign.var.DiffATE_CR.GRF[k])], CR.GRF.CATE.test, type = "l", ylim = range(min(lower.CR.GRF.test), max(upper.CR.GRF.test)), xlab = sign.var.DiffATE_CR.GRF[k], ylab = "CATE")
              lines(X.test[,(sign.var.DiffATE_CR.GRF[k])], upper.CR.GRF.test, col = 1, lty = 2)
              lines(X.test[,(sign.var.DiffATE_CR.GRF[k])], lower.CR.GRF.test, col = 1, lty = 2)
              grid()
              dev.off()
          }
         }
         
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
                                        
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
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
              pdf(fun_insert(x = filename.plot.LLCF.CATE, pos = (nchar(filename.plot.LLCF.CATE) - 4), insert = (sign.var.DiffATE_LLCF[k])))
              plot(X.test[,(sign.var.DiffATE_LLCF[k])], LLCF.CATE.test, type = "l", ylim = range(min(lower.LLCF.test), max(upper.LLCF.test)), xlab = sign.var.DiffATE_LLCF[k], ylab = "CATE")
              lines(X.test[,(sign.var.DiffATE_LLCF[k])], upper.LLCF.test, col = 1, lty = 2)
              lines(X.test[,(sign.var.DiffATE_LLCF[k])], lower.LLCF.test, col = 1, lty = 2)
              grid()
              dev.off()
          }
         }
        
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
                                                                 
       list("results_DiffATE_GRF" = results_DiffATE_GRF,
             "results_DiffATE_CR.GRF" = results_DiffATE_CR.GRF,
             "results_DiffATE_LLCF" = results_DiffATE_LLCF,
             "results_DiffATE_General" = results_DiffATE_General,
             "variable.importance.GRF" = GRF.varimp.ordered,
             "variable.importance.CR.GRF" = CR.GRF.varimp.ordered,
             "variable.importance.LLCF" = LLCF.varimp.ordered,
             "results_BLPredictor" = results_BLPredictor)
    
    })
                    
    return(results)
}


fun_insert <- function(x, pos, insert) {
  gsub(paste0("^(.{", pos, "})(.*)$"),
       paste0("\\1", insert, "\\2"),
       x)
}






                   
                   
  
