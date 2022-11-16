install.packages("grf")
library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(haven)
install.packages("glmnet")
library(glmnet)
install.packages("resample")
library(resample)
library(stats4)
set.seed(123)

# Import and prepare Voting Study data set (Arceneaux et al., 2006)
voting.study.data <- as.data.frame(read_dta("Downloads/ArceneauxGerberGreen_PA_2006_IA_MI_merge040504.dta"))
# voting.study.data <- read_dta("C:/Users/481044kt/Downloads/ArceneauxGerberGreen_PA_2006_IA_MI_merge040504.dta")

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
tau <- -0.5* (X$vote00 / (1 + (50/X$age)))
# tau <- -1* (X$vote00 / (2 + (100/X$age)))
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

# Determine treatment and control groups s.t. 2/5th is treated (59972) and 3/5th is controlled (89958), total 149930
treatment.group <- sample(which(voting.study.data$treat_real == 1), size = length(which(voting.study.data$treat_real == 1)), replace = FALSE) 
subsample.control.group <- sample(which(voting.study.data$treat_real == 0), size = 89958, replace = FALSE)
training.subsample <- c(treatment.group[1:ceiling((0.8*149930*0.4))], subsample.control.group[1:ceiling((0.8*149930*0.6))])
test.subsample <- c(treatment.group[-(1:ceiling((0.8*149930*0.4)))], subsample.control.group[-(1:ceiling((0.8*149930*0.6)))])

# Initialize parameters
num.reps.training = 100
num.reps.test = 10
numtrees <- 2000 # Set to 2000 or 8000 to perform sensitivity analysis.

# Different training sample sizes
training.sample.size <- c(200, 400, 600, 800, 1000, 1200, 2000, 4000) 
n.test <- 2000

run_method = function(training.sample.size, num.reps) {
  results = sapply(training.sample.size, function(n.train) {  
    
    basic.results = replicate(num.reps.training, {
      
      # Create training sample 
      indices.training = sample(training.subsample, size = n.train, replace = FALSE)
      
      # Appoint training sets
      X.train <- X[indices.training,] 
      W.train <- W[indices.training] 
      Ynew.train <- Ynew[indices.training]
      county.clusters.train <- county.clusters[indices.training]
      dummies.county.clusters.train <- dummies.county.clusters[indices.training,]
      tau.train <- tau[indices.training]
      
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
      GRF.CATE.oob <- GRF.oob$predictions
      GRF.CATE.SE.oob <- sqrt(GRF.oob$variance.estimates)
      
      # Compute measures
      GRF.MSE.oob <- sqrt(mean((GRF.CATE.oob - tau.train)**2))
      GRF.coverage.oob <- mean(as.integer((((GRF.CATE.oob - qnorm(0.975)*GRF.CATE.SE.oob) <= tau.train) & 
                                             (tau.train <= (GRF.CATE.oob + qnorm(0.975)*GRF.CATE.SE.oob)))))
      GRF.length.oob <- mean(abs((GRF.CATE.oob + qnorm(0.975)*GRF.CATE.SE.oob) - (GRF.CATE.oob - qnorm(0.975)*GRF.CATE.SE.oob)))
      GRF.abs.bias.oob <- mean((abs(GRF.CATE.oob - tau.train)))
      
      ############################
      #### Cluster-Robust GRF ####
      ############################
      
      # Implement Cluster-Robust GRF
      CR.GRF <- causal_forest(X.train, Ynew.train, W.train, clusters = county.clusters.train, honesty = TRUE, num.trees = numtrees)
      
      # Compute HTE with corresponding 95% confidence intervals                               
      CR.GRF.oob <- predict(CR.GRF, estimate.variance = TRUE)
      CR.GRF.CATE.oob <- CR.GRF.oob$predictions
      CR.GRF.CATE.SE.oob <- sqrt(CR.GRF.oob$variance.estimates)
      
      # Compute measures
      CR.GRF.MSE.oob <- sqrt(mean((CR.GRF.CATE.oob - tau.train)**2))
      CR.GRF.coverage.oob <- mean(as.integer((((CR.GRF.CATE.oob - qnorm(0.975)*CR.GRF.CATE.SE.oob) <= tau.train) & 
                                                (tau.train <= (CR.GRF.CATE.oob + qnorm(0.975)*CR.GRF.CATE.SE.oob)))))
      CR.GRF.length.oob <- mean(abs((CR.GRF.CATE.oob + qnorm(0.975)*CR.GRF.CATE.SE.oob) - (CR.GRF.CATE.oob - qnorm(0.975)*CR.GRF.CATE.SE.oob)))
      CR.GRF.abs.bias.oob <- mean((abs(CR.GRF.CATE.oob - tau.train)))
      
      ############################
      ########### LLCF ###########
      ############################
      
      # Grow preliminary forests for (W, X) and (Y, X) separately
      forest.W <- ll_regression_forest(as.matrix(X.train), as.vector(W.train), honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                       num.trees = 500)
      W.hat <- predict(forest.W)$predictions
      forest.Y <- ll_regression_forest(as.matrix(X.train), as.vector(Ynew.train), honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE,
                                       num.trees = 500)
      Y.hat <- predict(forest.Y)$predictions
      
      # Select variables to include using Lasso feature selection
      lasso.mod <- cv.glmnet(as.matrix(X.train), Ynew.train, alpha = 1)
      selected <- which(coef(lasso.mod) != 0)
      if(length(selected) < 2) {
        selected <- 1:ncol(X.train)
      } else {
        selected <- selected[-1] - 1 # Remove intercept
      }
      
      # Implement LLCF
      LLCF <- causal_forest(X.train, Ynew.train, W.train, Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, num.trees = numtrees, 
                            honesty.fraction = 0.7, honesty.prune.leaves = FALSE, tune.parameters = c("min.node.size", "sample.fraction", "mtry", "alpha", "imbalance.penalty"))
      
      # Compute HTE with corresponding 95% confidence intervals                                  
      LLCF.oob <- predict(LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, ll.lambda = 0.1, estimate.variance = TRUE)
      LLCF.CATE.oob <- LLCF.oob$predictions
      LLCF.CATE.SE.oob <- sqrt(LLCF.oob$variance.estimates)
      
      # Compute measures
      LLCF.MSE.oob <- sqrt(mean((LLCF.CATE.oob - tau.train)**2))
      LLCF.coverage.oob <- mean(as.integer((((LLCF.CATE.oob - qnorm(0.975)*LLCF.CATE.SE.oob) <= tau.train) & 
                                              (tau.train <= (LLCF.CATE.oob + qnorm(0.975)*LLCF.CATE.SE.oob)))))
      LLCF.length.oob <- mean(abs((LLCF.CATE.oob + qnorm(0.975)*LLCF.CATE.SE.oob) - (LLCF.CATE.oob -  qnorm(0.975)*LLCF.CATE.SE.oob)))
      LLCF.abs.bias.oob <- mean((abs(LLCF.CATE.oob - tau.train)))
      
      sub.basic.results = replicate(num.reps.test, {
        # Create test sample 
        indices.test <- sample(test.subsample, size = n.test, replace = FALSE)
        
        # Appoint test sets
        X.test <- X[indices.test,]
        W.test <- W[indices.test]
        Ynew.test <- Ynew[indices.test]
        county.clusters.test <- county.clusters[indices.test]
        dummies.county.clusters.test <- dummies.county.clusters[indices.test,]
        tau.test <- tau[indices.test]
        current.X.test <- cbind(X.test, dummies.county.clusters.test)
        
        # Predictions GRF
        GRF.pred <- predict(GRF, current.X.test, estimate.variance = TRUE)
        GRF.CATE <- GRF.pred$predictions
        GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
        
        # Predictions CR.GRF
        CR.GRF.pred <- predict(CR.GRF, X.test, estimate.variance = TRUE)
        CR.GRF.CATE <- CR.GRF.pred$predictions
        CR.GRF.CATE.SE <- sqrt(CR.GRF.pred$variance.estimates)
        
        # Predictions LLCF
        LLCF.pred <- predict(LLCF, X.test, linear.correction.variables = selected, ll.weight.penalty = TRUE, ll.lambda = 0.1, estimate.variance = TRUE)
        LLCF.CATE <- LLCF.pred$predictions
        LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
        
        # Compute measures
        GRF.MSE <- sqrt(mean((GRF.CATE - tau.test)**2))
        GRF.coverage <- mean(as.integer((((GRF.CATE - qnorm(0.975)*GRF.CATE.SE) <= tau.test) & 
                                           (tau.test <= (GRF.CATE + qnorm(0.975)*GRF.CATE.SE)))))
        GRF.length <- mean(abs((GRF.CATE + qnorm(0.975)*GRF.CATE.SE) - (GRF.CATE - qnorm(0.975)*GRF.CATE.SE)))
        GRF.abs.bias <- mean((abs(GRF.CATE - tau.test)))
        
        CR.GRF.MSE <- sqrt(mean((CR.GRF.CATE - tau.test)**2))
        CR.GRF.coverage <- mean(as.integer((((CR.GRF.CATE - qnorm(0.975)*CR.GRF.CATE.SE) <= tau.test) & 
                                              (tau.test <= (CR.GRF.CATE + qnorm(0.975)*CR.GRF.CATE.SE)))))
        CR.GRF.length <- mean(abs((CR.GRF.CATE + qnorm(0.975)*CR.GRF.CATE.SE) - (CR.GRF.CATE - qnorm(0.975)*CR.GRF.CATE.SE)))
        CR.GRF.abs.bias <- mean((abs(CR.GRF.CATE - tau.test)))
        
        LLCF.MSE <- sqrt(mean((LLCF.CATE - tau.test)**2))
        LLCF.coverage <- mean(as.integer((((LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE) <= tau.test) & 
                                            (tau.test <= (LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE)))))
        LLCF.length <- mean(abs((LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE) - (LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE)))
        LLCF.abs.bias <- mean((abs(LLCF.CATE - tau.test)))
        
        return(c(GRF.MSE, CR.GRF.MSE, LLCF.MSE,
                 GRF.abs.bias, CR.GRF.abs.bias, LLCF.abs.bias,
                 GRF.coverage, CR.GRF.coverage, LLCF.coverage,
                 GRF.length, CR.GRF.length, LLCF.length))
      })  
      
      
      return(c(GRF.MSE.oob, CR.GRF.MSE.oob, LLCF.MSE.oob, 
               GRF.abs.bias.oob, CR.GRF.abs.bias.oob, LLCF.abs.bias.oob,
               GRF.coverage.oob, CR.GRF.coverage.oob, LLCF.coverage.oob,
               GRF.length.oob, CR.GRF.length.oob, LLCF.length.oob,
               colMeans(t(sub.basic.results))))
    })
    
    basic.results = data.frame(t(basic.results))
    standard.deviation.basic.results <- sqrt(colVars(basic.results))
    mean.basic.results <- colMeans(basic.results)
    make.vector <- c()
    for (i in 1:length(mean.basic.results)) {
      make.vector <- cbind(make.vector, paste(round(mean.basic.results[i], 3), "(", round(standard.deviation.basic.results[i], 3), ")"))
    }
    return(make.vector)
  })
  
  results = data.frame(t(results))
  colnames(results) = c("In-sample RMSE GRF", "In-sample RMSE CR.GRF", "In-sample RMSE LLCF", 
                        "In-sample Absolute Bias GRF", "In-sample Absolute Bias CR.GRF", "In-sample Absolute Bias LLCF",
                        "In-sample Coverage GRF", "In-sample Coverage CR.GRF", "In-sample Coverage LLCF",
                        "In-sample Length GRF", "In-sample Length CR.GRF", "In-sample Length LLCF",
                        "Out-of-sample RMSE GRF", "Out-of-sample RMSE CR.GRF", "Out-of-sample RMSE LLCF",
                        "Out-of-sample Absolute Bias GRF", "Out-of-sample Absolute Bias CR.GRF", "Out-of-sample Absolute Bias LLCF",
                        "Out-of-sample Coverage GRF", "Out-of-sample Coverage CR.GRF", "Out-of-sample Coverage LLCF",
                        "Out-of-sample Length GRF", "Out-of-sample Length CR.GRF", "Out-of-sample Length LLCF")
  results$n.train = training.sample.size
  results
}

results.semi.synthetic = run_method(training.sample.size, num.reps)



