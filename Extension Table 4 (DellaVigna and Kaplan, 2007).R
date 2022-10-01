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
        GRF.mostimportant <- colnames(X)[GRF.varimp.ordered[1:4]] # 4 most important variables for splitting
        
      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, honesty = TRUE)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)
  
      # Implement GRF
        GRF <- causal_forest(X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, 
                             honesty = TRUE, tune.parameters = "all")
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
        
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
        if (boolean.plot == TRUE) {
          for (k in 1:length(GRF.mostimportant)) {
            # Set all variables at their median values
              X.median <- apply(X, 2, median)
              
            # Create ordered vector of important variable
              important.var.test = seq(min(X$GRF.mostimportant[k]), max(X$GRF.mostimportant[k]))
              
            # Create test set
              X.test <- matrix(rep(X.median, length(important.var.test)), length(important.var.test), byrow = TRUE)
              X.test[,(GRF.mostimportant[k])] = important.var.test
              
            # Predict new CATE estimates
              GRF.pred.test <- predict(GRF, X.test, estimate.variance = TRUE)
              GRF.CATE.test <- GRF.pred.test$predictions
              GRF.CATE.SE.test <- sqrt(GRF.pred.test$variance.estimates)
              lower.GRF.test <- GRF.CATE.test - qnorm(0.975)*GRF.CATE.SE.test
              upper.GRF.test <- GRF.CATE.test + qnorm(0.975)*GRF.CATE.SE.test
            
            # Make plot and save
              pdf(fun_insert(x = filename.plot.GRF.CATE, pos = (nchar(filename.plot.GRF.CATE) - 4), insert = (GRF.mostimportant[k])))
              plot(X.test[,(GRF.mostimportant[k])], GRF.CATE.test, type = "l", ylim = range(min(lower.GRF.test), max(upper.GRF.test)), xlab = GRF.mostimportant[k], ylab = "CATE")
              lines(X.test[,(GRF.mostimportant[k])], upper.GRF.test, col = 1, lty = 2)
              lines(X.test[,(GRF.mostimportant[k])], lower.GRF.test, col = 1, lty = 2)
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
                                     
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
        if (boolean.plot == TRUE) {
          for (k in 1:length(CR.GRF.mostimportant)) {
            # Set all variables at their median values
              X.median <- apply(X, 2, median)
              
            # Create ordered vector of important variable
              important.var.test = seq(min(X$CR.GRF.mostimportant[k]), max(X$CR.GRF.mostimportant[k]))
              
            # Create test set
              X.test <- matrix(rep(X.median, length(important.var.test)), length(important.var.test), byrow = TRUE)
              X.test[,(CR.GRF.mostimportant[k])] = important.var.test
              
            # Predict new CATE estimates
              CR.GRF.pred.test <- predict(CR.GRF, X.test, estimate.variance = TRUE)
              CR.GRF.CATE.test <- CR.GRF.pred.test$predictions
              CR.GRF.CATE.SE.test <- sqrt(CR.GRF.pred.test$variance.estimates)
              lower.CR.GRF.test <- CR.GRF.CATE.test - qnorm(0.975)*CR.GRF.CATE.SE.test
              upper.CR.GRF.test <- CR.GRF.CATE.test + qnorm(0.975)*CR.GRF.CATE.SE.test
            
            # Make plot and save
              pdf(fun_insert(x = filename.plot.CR.GRF.CATE, pos = (nchar(filename.plot.CR.GRF.CATE) - 4), insert = (CR.GRF.mostimportant[k])))
              plot(X.test[,(CR.GRF.mostimportant[k])], CR.GRF.CATE.test, type = "l", ylim = range(min(lower.CR.GRF.test), max(upper.CR.GRF.test)), xlab = CR.GRF.mostimportant[k], ylab = "CATE")
              lines(X.test[,(CR.GRF.mostimportant[k])], upper.CR.GRF.test, col = 1, lty = 2)
              lines(X.test[,(CR.GRF.mostimportant[k])], lower.CR.GRF.test, col = 1, lty = 2)
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
        
      # Find lower and upper bounds for 95% confidence intervals
        lower.LLCF <- LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE
        upper.LLCF <- LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE
      
      # Graph the predicted Cluster-Robust GRF heterogeneous treatment effect estimates
        if (boolean.plot == TRUE) {
          hist(LLCF.CATE)
        }
    
      # Compute ATE with corresponding 95% confidence intervals
        LLCF.ATE <- average_treatment_effect(LLCF, target.sample = "all")
    
      # See if the Cluster-Robust GRF succeeded in capturing heterogeneity by plotting the TOC and calculating the 95% confidence interval for the AUTOC
        LLCF.rate <- rank_average_treatment_effect(LLCF, LLCF.CATE, target = "AUTOC")
        if (boolean.plot == TRUE) {
          plot(LLCF.rate)
        }
  
      # Assessing Cluster-Robust GRF fit using the Best Linear Predictor Approach [BLP]
        LLCF.BLP <- test_calibration(LLCF)
        
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
                                      
      # Plot the estimated CATE against the covariates with the highest variable importance, and the characteristic vector
        if (boolean.plot == TRUE) {
          for (k in 1:length(LLCF.mostimportant)) {
            # Set all variables at their median values
              X.median <- apply(X, 2, median)
              
            # Create ordered vector of important variable
              important.var.test = seq(min(X$LLCF.mostimportant[k]), max(X$LLCF.mostimportant[k]))
              
            # Create test set
              X.test <- matrix(rep(X.median, length(important.var.test)), length(important.var.test), byrow = TRUE)
              X.test[,(LLCF.mostimportant[k])] = important.var.test
              
            # Predict new CATE estimates
              LLCF.pred.test <- predict(LLCF, X.test, estimate.variance = TRUE)
              LLCF.CATE.test <- LLCF.pred.test$predictions
              LLCF.CATE.SE.test <- sqrt(LLCF.pred.test$variance.estimates)
              lower.LLCF.test <- LLCF.CATE.test - qnorm(0.975)*LLCF.CATE.SE.test
              upper.LLCF.test <- LLCF.CATE.test + qnorm(0.975)*LLCF.CATE.SE.test
            
            # Make plot and save
              pdf(fun_insert(x = filename.plot.LLCF.CATE, pos = (nchar(filename.plot.LLCF.CATE) - 4), insert = (LLCF.mostimportant[k])))
              plot(X.test[,(LLCF.mostimportant[k])], LLCF.CATE.test, type = "l", ylim = range(min(lower.LLCF.test), max(upper.LLCF.test)), xlab = LLCF.mostimportant[k], ylab = "CATE")
              lines(X.test[,(LLCF.mostimportant[k])], upper.LLCF.test, col = 1, lty = 2)
              lines(X.test[,(LLCF.mostimportant[k])], lower.LLCF.test, col = 1, lty = 2)
              grid()
              dev.off()
          }
         }
        
                                      
       
                                      
       list("results_DiffATE_GRF" = results_DiffATE_GRF,
             "results_DiffATE_CR.GRF" = results_DiffATE_CR.GRF,
             "results_DiffATE_LLCF" = results_DiffATE_LLCF,
             "variable.importance.GRF" = GRF.varimp.ordered,
             "variable.importance.CR.GRF" = CR.GRF.varimp.ordered,
             "variable.importance.LLCF" = LLCF.varimp.ordered,)
    
    })
    results_BLP <- t(basic.results[,1:6])   
    colnames(results_BLP) <- c("BLP[1] GRF", "BLP[2] GRF", "BLP[1] CR.GRF", "BLP[2] CR.GRF", "BLP[1] LLCF", "BLP[2] LLCF")
    
                       
    return(results)
}


fun_insert <- function(x, pos, insert) {
  gsub(paste0("^(.{", pos, "})(.*)$"),
       paste0("\\1", insert, "\\2"),
       x)
}






                   
                   
  
