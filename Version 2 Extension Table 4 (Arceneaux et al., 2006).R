library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(haven)
library(glmnet)
library(stats4)
set.seed(1234)

# Import and prepare Voting Study data set (Arceneaux et al., 2006)
voting.study.data <- as.data.frame(read_dta("Downloads/ArceneauxGerberGreen_PA_2006_IA_MI_merge040504.dta"))

# Clean data of NA values in "contact" and "treat_real" variables
voting.study.data <- voting.study.data[(is.na(voting.study.data$treat_real) == FALSE),]

# Appoint main outcome variable
Y <- as.vector(voting.study.data$vote02)

# Appoint treatment assignment variable
W <- as.vector(voting.study.data$treat_real)

# Appoint the X-variables
X <- voting.study.data[,c("state", "comp_mi", "comp_ia", "persons", "age", "female2", "newreg", "vote00", "vote98", "fem_miss")]

# Create a numerical vector of the county vector
county.clusters <- as.vector(voting.study.data$county)

dummies.county.clusters <- matrix(0, nrow = length(county.clusters), ncol = length(unique(county.clusters)))
colnames(dummies.county.clusters) <- c(unique(county.clusters))
for (i in 1:length(county.clusters)) {
  dummies.county.clusters[i,toString(county.clusters[i])] <- 1
}

# Synthetic treatment assignment
tau <- -1* (X$vote00 / (2 + (100/X$age)))
# tau <- -1* (X$vote00 / (1 + exp(-1*(X$age)/50)))
# tau <- -1* (X$vote00 / (1.5 + (1000/X$age^2)))
R <- rbinom(n = length(tau), size = 1, prob = abs(tau))
Y0 <- matrix(0, nrow = length(Y), ncol = 1)
Y1 <- matrix(0, nrow = length(Y), ncol = 1)
Ynew <- matrix(0, nrow = length(Y), ncol = 1)
Y0[(R == 0)] <- Y[(R == 0)]
Y1[(R == 0)] <- Y[(R == 0)]
Y0[((R == 1) & (tau > 0))] <- 0
Y1[((R == 1) & (tau > 0))] <- 1
Y0[((R == 1) & (tau < 0))] <- 1
Y1[((R == 1) & (tau < 0))] <- 0
Ynew[(W == 0)] <- Y0[(W == 0)]
Ynew[(W == 1)] <- Y1[(W == 1)]
treatment.group <- sample(which(voting.study.data$treat_real == 1), size = length((which(voting.study.data$treat_real == 1))), replace = FALSE)
control.group <- which(voting.study.data$treat_real == 0)
subsample.control.group <- sample(control.group, size = 89958, replace = FALSE)

# Split the data into a balanced test set (2/5th treated, 3/5th control) of 25000 observations and balanced training set of 100000 observations
# indices.test <- c(treatment.group[1:10000], subsample.control.group[1:15000])
indices.test <-  c(treatment.group[1:800], subsample.control.group[1:1200]) 
indices.training <- c(treatment.group[10001:50000], subsample.control.group[15001:75000])
X.test <- X[indices.test,]
W.test <- W[indices.test]
Ynew.test <- Ynew[indices.test]
county.clusters.test <- county.clusters[indices.test]
dummies.county.clusters.test <- dummies.county.clusters[indices.test,]
tau.test <- tau[indices.test]
current.X.test <- cbind(X.test, dummies.county.clusters.test)

# Initialize parameters
num.reps = 50
numtrees <- 2000 # Set to 1000 or 5000 to perform sensitivity analysis.

# Different training sample sizes
training.sample.size <- c(200, 400, 600, 800, 1000, 1200, 1500, 2000) 
# training.sample.size <- c(2000, 5000, 10000, 20000, 30000)

run_method = function(training.sample.size, num.reps) {
  results = sapply(training.sample.size, function(n.train) {  
    basic.results = replicate(num.reps, {
      
      # Shuffle the training set
      subsample.indices.training <- sample(indices.training, size = n.train, replace = FALSE)
      
      # Appoint training and test sets
      X.train <- X[subsample.indices.training,] 
      W.train <- W[subsample.indices.training] 
      Ynew.train <- Ynew[subsample.indices.training]
      county.clusters.train <- county.clusters[subsample.indices.training]
      dummies.county.clusters.train <- dummies.county.clusters[subsample.indices.training,]
      tau.train <- tau[subsample.indices.training]
      
      # Estimation procedure
      
      ###########################
      ########### GRF ###########
      ###########################
      
      # For GRF we create the X matrix including district-specific dummies
      current.X <- cbind(X.train, dummies.county.clusters.train)
      
      # Implement GRF
      GRF <- causal_forest(current.X, Ynew.train, W.train, num.trees = numtrees, honesty = TRUE)
      
      # Compute HTE with corresponding 95% confidence intervals
      GRF.oob <- predict(GRF, estimate.variance = TRUE)
      GRF.pred <- predict(GRF, current.X.test, estimate.variance = TRUE)
      GRF.CATE.oob <- GRF.oob$predictions
      GRF.CATE <- GRF.pred$predictions
      GRF.CATE.SE.oob <- sqrt(GRF.oob$variance.estimates)
      GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
      lower.GRF <- GRF.CATE - qnorm(0.975)*GRF.CATE.SE
      upper.GRF <- GRF.CATE + qnorm(0.975)*GRF.CATE.SE
      
      # Compute measures
      GRF.MSE.oob <- sqrt(mean((GRF.CATE.oob - tau.train)**2))
      GRF.MSE <- sqrt(mean((GRF.CATE - tau.test)**2))
      
      ############################
      #### Cluster-Robust GRF ####
      ############################
      
      # Implement Cluster-Robust GRF
      CR.GRF <- causal_forest(X.train, Ynew.train, W.train, clusters = county.clusters.train, honesty = TRUE, num.trees = numtrees)
      
      # Compute HTE with corresponding 95% confidence intervals                               
      CR.GRF.oob <- predict(CR.GRF, estimate.variance = TRUE)
      CR.GRF.pred <- predict(CR.GRF, X.test, estimate.variance = TRUE)
      CR.GRF.CATE.oob <- CR.GRF.oob$predictions
      CR.GRF.CATE <- CR.GRF.pred$predictions
      CR.GRF.CATE.SE.oob <- sqrt(CR.GRF.oob$variance.estimates)
      CR.GRF.CATE.SE <- sqrt(CR.GRF.pred$variance.estimates)
      lower.CR.GRF <- CR.GRF.CATE - qnorm(0.975)*CR.GRF.CATE.SE
      upper.CR.GRF <- CR.GRF.CATE + qnorm(0.975)*CR.GRF.CATE.SE
      
      # Compute measures
      CR.GRF.MSE.oob <- sqrt(mean((CR.GRF.CATE.oob - tau.train)**2))
      CR.GRF.MSE <- sqrt(mean((CR.GRF.CATE - tau.test)**2))
      
      ############################
      ########### LLCF ###########
      ############################
      
      # Grow preliminary forests for (W, X) and (Y, X) separately
      # forest.W <- ll_regression_forest(as.matrix(X.train), as.vector(W.train), honesty = TRUE, num.trees = 500, tune.parameters = "all")
      # W.hat <- predict(forest.W)$predictions
      # forest.Y <- ll_regression_forest(as.matrix(X.train), as.vector(Ynew.train), honesty = TRUE, num.trees = 500, tune.parameters = "all")
      # Y.hat <- predict(forest.Y)$predictions
      
      # Select variables to include using Lasso feature selection
      lasso.mod <- cv.glmnet(as.matrix(current.X), Ynew.train, alpha = 1)
      selected <- which(coef(lasso.mod) != 0)
      if(length(selected) < 2) {
        selected <- 1:ncol(current.X)
      } else {
        selected <- selected[-1] - 1 # Remove intercept
      }
      
      # Implement LLCF
      # LLCF <- causal_forest(X.train, Ynew.train, W.train, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, num.trees = numtrees)
      
      # Compute HTE with corresponding 95% confidence intervals                                  
      LLCF.oob <- predict(GRF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
      LLCF.pred <- predict(GRF, current.X.test, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
      LLCF.CATE.oob <- LLCF.oob$predictions
      LLCF.CATE <- LLCF.pred$predictions
      LLCF.CATE.SE.oob <- sqrt(LLCF.oob$variance.estimates)
      LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
      
      # Compute measures
      LLCF.MSE.oob <- sqrt(mean((LLCF.CATE.oob - tau.train)**2))
      LLCF.MSE <- sqrt(mean((LLCF.CATE - tau.test)**2))
      
      return(c(GRF.MSE.oob, CR.GRF.MSE.oob, LLCF.MSE.oob, GRF.MSE, CR.GRF.MSE, LLCF.MSE))
    })
    
    basic.results = data.frame(t(basic.results))
    colMeans(basic.results)
    
  })
  
  results = data.frame(t(results))
  colnames(results) = c("In-sample RMSE GRF", "In-sample RMSE CR.GRF", "In-sample RMSE LLCF", 
                        "Out-of-sample RMSE GRF", "Out-of-sample RMSE CR.GRF", "Out-of-sample RMSE LLCF")
  results$n.train = training.sample.size
  results
}

results = run_method(training.sample.size, num.reps)


