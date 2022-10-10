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
  trainfrac <- 0.7 # Size of the training sample
  B <- 100 # Number of replications
  rows <- c(1:nrow(X))
  samplesize <- nrow(X)
  lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5) # Concerns ridge penalty parameters; do not adjust.
  boolean.lambdas <- FALSE # Set to TRUE to use lambdas instead of automatic penalty tuning.
  boolean.plot <- FALSE # Set to TRUE to make various plots of interest.
  results.GRF <- list()
  results.CR.GRF <- list()
  results.LLCF <- list()
  filename.plot.GRF.CATE <- "GRF CATE .pdf"
  filename.plot.CR.GRF.CATE <- "CR.GRF CATE .pdf"
  filename.plot.LLCF.CATE <- "LLCF CATE .pdf"
  
for (b in 1:B) {
  # Randomly shuffle (sample without replacement) the observations 
    split <- sample(rows, replace = FALSE)
      
  # Split the data into a training and test sample
    Xtrain <- X[split[1:ceiling(trainfrac*samplesize)],]
    Xtest <- X[-split[1:ceiling(trainfrac*samplesize)],]
      
    Wtrain <- W[split[1:ceiling(trainfrac*samplesize)]]
    Wtest <- W[-split[1:ceiling(trainfrac*samplesize)]]
      
    Ytrain <- Y[split[1:ceiling(trainfrac*samplesize)]]
    Ytest <- Y[-split[1:ceiling(trainfrac*samplesize)]]
        
    remaining.variables.train <- remaining.variables[split[1:ceiling(trainfrac*samplesize)],]
    remaining.variables.test <- remaining.variables[-split[1:ceiling(trainfrac*samplesize)],]
      
    district.clusters.train <- district.clusters[split[1:ceiling(trainfrac*samplesize)]]
    district.clusters.test <- district.clusters[-split[1:ceiling(trainfrac*samplesize)]]
        
    dummies.district.clusters.train <- dummies.district.clusters[split[1:ceiling(trainfrac*samplesize)]]
    dummies.district.clusters.test <- dummies.district.clusters[-split[1:ceiling(trainfrac*samplesize)]]
      
    ###########################
    ########### GRF ###########
    ###########################
    
      # For GRF we create the X matrix including district-specific dummies
        current.X.train <- cbind(Xtrain, dummies.district.clusters.train)
        current.X.test <- cbind(Xtest, dummies.district.clusters.test)
        
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- regression_forest(current.X.train, Wtrain, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
        W.hat <- predict(forest.W)$predictions
        forest.Y <- regression_forest(current.X.train, Ytrain, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
        Y.hat <- predict(forest.Y)$predictions
        
      # Compute the variable importance
        GRF.varimp <- variable_importance(forest.Y) 
        GRF.varimp.ordered <- order(GRF.varimp)
        GRF.mostimportant <- colnames(X)[(GRF.varimp.ordered[1:4])] # 4 most important variables for splitting
        
      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(current.X.train, Ytrain, Wtrain, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, honesty = TRUE)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)
  
      # Implement GRF
        GRF <- causal_forest(current.X.train[,selected.vars], Ytrain, Wtrain, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, 
                             honesty = TRUE, tune.parameters = "all")
                             
      # Make GRF predictions based on test set  
        GRF.pred <- predict(GRF, newdata = current.X.test[,selected.vars])
        GRF.CATE <- GRF.pred$predictions
          
    ############################
    #### Cluster-Robust GRF ####
    ############################
      
      # Select variables to include using preliminary Cluster-Robust GRF
        prelim.CR.GRF <- causal_forest(Xtrain, Ytrain, Wtrain, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, 
                                        clusters = district.clusters.train, num.trees = numtrees)
        prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
        selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)
  
      # Compute the variable importance
        CR.GRF.varimp <- variable_importance(forest.Y) 
        CR.GRF.varimp.ordered <- order(CR.GRF.varimp)                             
        CR.GRF.mostimportant <- colnames(X)[CR.GRF.varimp.ordered[1:4]] # 4 most important variables for splitting
        
      # Implement Cluster-Robust GRF
        CR.GRF <- causal_forest(Xtrain[,selected.vars], Ytrain, Wtrain, Y.hat = Y.hat, W.hat = W.hat, 
                                        clusters = district.clusters.train, honesty = TRUE, num.trees = numtrees, tune.parameters = "all")
        
      # Make Cluster-Robust GRF predictions based on test set  
        CR.GRF.pred <- predict(CR.GRF, newdata = Xtest[,selected.vars])
        CR.GRF.CATE <- CR.GRF.pred$predictions
          
    ############################
    ########### LLCF ###########
    ############################
    
      # Grow preliminary forests for (W, X) and (Y, X) separately
        forest.W <- ll_regression_forest(current.X.train, Wtrain, honesty = TRUE, ll.split.variables = selected.vars, num.trees = 50)
        W.hat <- predict(forest.W)$predictions
        forest.Y <- ll_regression_forest(current.X.train, Ytrain, honesty = TRUE, ll.split.variables = selected.vars, num.trees = 50)
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
        LLCF <- causal_forest(current.X.train, Ytrain, Wtrain, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE,  
                              num.trees = numtrees, tune.parameters = "all")
                              
      # Compute HTE with corresponding 95% confidence intervals                                  
        if (boolean.lambdas == FALSE) {
          # Predict: tuning without grid search over lambdas
            LLCF.pred <- predict(LLCF, newdata = current.X.test[,selected.vars], linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE, ll.lambda = 0.1)
            LLCF.CATE <- LLCF.pred$predictions
            LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
        } else {
          # Predict: tuning done using set of lambdas
            LLCF.mse.old <- +Inf
            for (l in length(lambdas)) {
              LLCF.CATE.old <- predict(GRF, newdata = current.X.test[,selected.vars], linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
              predictions <- LLCF.CATE.old$predictions
              LLCF.mse.new <- mean((predictions - mean(predictions))**2)
              if (LLCF.mse.new < LLCF.mse.old) {
              LLCF.mse.old <- LLCF.mse.new
              LLCF.CATE.SE <- sqrt(LLCF.CATE.old$variance.estimates)
              LLCF.CATE <- predictions
            }
          }
        }     
        
    ############################
    ###### Update results ######
    ############################
    
        results.GRF[[b]] <- cbind(current.X.test[], GRF.CATE)
        results.CR.GRF[[b]] <- cbind(X.test[], CR.GRF.CATE)
        results.LLCF[[b]] <- cbind(current.X.test[], LLCF.CATE)
        
 }
 
 results.GRF <- do.call(rbind, results.GRF)
 results.CR.GRF <- do.call(rbind, results.CR.GRF)
 results.LLCF <- do.call(rbind, results.LLCF)
 
 boxplot(results.GRF[,2] ~ results.GRF[,1], main = " ", xlab = " ", ylab = " ")
 boxplot(results.CR.GRF[,2] ~ results.CR.GRF[,1], main = " ", xlab = " ", ylab = " ")
 boxplot(results.LLCF[,2] ~ results.LLCF[,1], main = " ", xlab = " ", ylab = " ")
 
 
 
