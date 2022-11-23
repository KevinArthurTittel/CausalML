library(grf)
library(glmnet)
library(hte)
set.seed(123)
# Simulation setting Nie & Wager (2021), Set-up A & B combined
  # Uniformly distributed X
  # Bernoulli distributed W
  # Standard normal noise
  # Main effect scaled Friedman function
  # Treatment propensity being trimmed sinus function (part of main effect also): confounding
  # Treatment effect function from set-up B, being smooth logistic

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- max(0.1, min(sin(pi*X[,1]*X[,2]), 0.9))
  main_effect <- sin(pi*X[,1]*X[,2]) + 2(X[,3] - 0.5)^2 + X[,4] + 0.5*X[,5]
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  Tau <- X[,1] + log(1 + exp(X[,2]))
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}

# Simulation settings Nie & Wager (2021) (set-up A) and Athey & Wager (2019) combined
  # Uniformly distributed X
  # Bernoulli distributed W
  # Standard normal noise
  # Main effect scaled Friedman function
  # Treatment propensity being trimmed sinus function (part of main effect also): confounding
  # Treatment effect function from Athey & Wager (2019), being smooth exponentiated

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- max(0.1, min(sin(pi*X[,1]*X[,2]), 0.9))
  main_effect <- sin(pi*X[,1]*X[,2]) + 2*(X[,3] - 0.5)^2 + X[,4] + 0.5*X[,5]
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 1 + 1/(1 + exp(-20*(X[,1] - (1/3))))
  zeta2 <- 1 + 1/(1 + exp(-20*(X[,2] - (1/3))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}

# Simulation settings Nie & Wager (2021) (set-up B & C) and Athey & Wager (2019) combined
  # Normally distributed X
  # Bernoulli distributed W
  # Standard normal noise
  # Main effect 
  # Treatment propensity 
  # Treatment effect function 

simulation <- function(n, d, sigma) {
  X <- matrix(rnorm(n*d, mean = 0, sd = 1), nrow = n, ncol = d)
  treatment_propensity <- 
  main_effect <- 
  noise <- rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 1 + 1/(1 + exp(-20*(X[,1] - (1/3))))
  zeta2 <- 1 + 1/(1 + exp(-20*(X[,2] - (1/3))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}

s <- 1
num.reps <- 50
n.training <- c(200, 400, 600, 800, 1000, 1200, 2000)
d <- 25

run_simulation <- function(n.training, num.reps) {
  results <- sapply(n.training, function(n) {
    basic.results <- replicate(num.reps, {
    
          # Create the training and test data sets
            dat <- simulation(n = n, d = d, sigma = s)
            dat.test <- simulation(n = 2000, d = d, sigma = s)
            tau.training <- dat$Tau
            tau.test <- dat.test$Tau # True tau
        
        ###########################
        ########### GRF ###########
        ###########################
        
          # Run GRF algorithm on test data set, obtain predictions, and compute measures
            GRF <- causal_forest(as.matrix(dat$X), as.vector(dat$Y), as.vector(dat$W), honesty = TRUE, num.trees = 2000)
            GRF.pred.oob <- predict(GRF, estimate.variance = TRUE)
            GRF.pred <- predict(GRF, as.matrix(dat.test$X), estimate.variance = TRUE)
            GRF.CATE.oob <- GRF.pred.oob$predictions
            GRF.CATE.SE.oob <- sqrt(GRF.pred.oob$variance.estimates)
            GRF.CATE <- GRF.pred$predictions
            GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
            
          # Compute measures
            GRF.MSE.oob <- sqrt(mean((GRF.CATE.oob - tau.training)**2))
            GRF.MSE <- sqrt(mean((GRF.CATE - tau.test)**2))
            GRF.abs.bias.oob <- mean((abs(GRF.CATE.oob - tau.training)))
            GRF.abs.bias <- mean((abs(GRF.CATE - tau.test)))
            GRF.coverage.oob <- mean(as.integer((((GRF.CATE.oob - qnorm(0.975)*GRF.CATE.SE.oob) <= tau.training) & 
                                               (tau.training <= (GRF.CATE.oob + qnorm(0.975)*GRF.CATE.SE.oob)))))
            GRF.coverage <- mean(as.integer((((GRF.CATE - qnorm(0.975)*GRF.CATE.SE) <= tau.test) & 
                                           (tau.test <= (GRF.CATE + qnorm(0.975)*GRF.CATE.SE)))))
            GRF.length.oob <- mean(abs((GRF.CATE.oob + qnorm(0.975)*GRF.CATE.SE.oob) - (GRF.CATE.oob - qnorm(0.975)*GRF.CATE.SE.oob)))
            GRF.length <- mean(abs((GRF.CATE + qnorm(0.975)*GRF.CATE.SE) - (GRF.CATE - qnorm(0.975)*GRF.CATE.SE)))
  
        ############################
        ########### LLCF ###########
        ############################
          
          # Grow preliminary forests for (W, X) and (Y, X) separately
            forest.W <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$W), honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                             honesty.fraction = 0.7, honesty.prune.leaves = FALSE, ll.split.lambda = 0.1, num.trees = 500)
            W.hat <- predict(forest.W)$predictions
            forest.Y <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$Y), honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                             honesty.fraction = 0.7, honesty.prune.leaves = FALSE, ll.split.lambda = 0.1, num.trees = 500)
            Y.hat <- predict(forest.Y)$predictions
            
          # Implement LLCF
            LLCF <- causal_forest(as.matrix(dat$X), as.numeric(dat$Y), as.numeric(dat$W), Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, honesty.fraction = 0.7, honesty.prune.leaves = FALSE, 
                                  num.trees = numtrees, tune.parameters = c("min.node.size", "sample.fraction", "mtry", "alpha", "imbalance.penalty"))
          
          # Run LLCF algorithm on test data set, obtain predictions, and compute measures
            lasso.mod <- cv.glmnet(as.matrix(dat$X), as.numeric(dat$Y), alpha = 1)
            selected <- which(coef(lasso.mod) != 0)
            if(length(selected) < 2) {
              selected <- 1:ncol(dat$X)
            } else {
              selected <- selected[-1] - 1 # Remove intercept
            }
            
          # Make predictions
            LLCF.pred.oob <- predict(LLCF, linear.correction.variables = selected, ll.lambda = 0.1, ll.weight.penalty = TRUE, estimate.variance = TRUE)
            LLCF.pred <- predict(LLCF, as.matrix(dat.test$X), linear.correction.variables = selected, ll.lambda = 0.1, ll.weight.penalty = TRUE, estimate.variance = TRUE)
            LLCF.CATE.oob <- LLCF.pred.oob$predictions
            LLCF.CATE.SE.oob <- sqrt(LLCF.pred.oob$variance.estimates)
            LLCF.CATE <- LLCF.pred$predictions
            LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
        
          # Compute measures
            LLCF.MSE.oob <- sqrt(mean((LLCF.CATE.oob - tau.training)**2))
            LLCF.MSE <- sqrt(mean((LLCF.CATE - tau.test)**2))
            LLCF.abs.bias.oob <- mean((abs(LLCF.CATE.oob - tau.training)))
            LLCF.abs.bias <- mean((abs(LLCF.CATE - tau.test)))
            LLCF.coverage.oob <- mean(as.integer((((LLCF.CATE.oob - qnorm(0.975)*LLCF.CATE.SE.oob) <= tau.training) & 
                                                   (tau.training <= (LLCF.CATE.oob + qnorm(0.975)*LLCF.CATE.SE.oob)))))
            LLCF.coverage <- mean(as.integer((((LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE) <= tau.test) & 
                                               (tau.test <= (LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE)))))
            LLCF.length.oob <- mean(abs((LLCF.CATE.oob + qnorm(0.975)*LLCF.CATE.SE.oob) - (LLCF.CATE.oob - qnorm(0.975)*LLCF.CATE.SE.oob)))
            LLCF.length <- mean(abs((LLCF.CATE + qnorm(0.975)*LLCF.CATE.SE) - (LLCF.CATE - qnorm(0.975)*LLCF.CATE.SE)))
        
        ############################
        ###### Update results ######
        ############################
        
            return(c(s, d, n, 
                  GRF.MSE.oob, LLCF.MSE.oob, GRF.abs.bias.oob,  LLCF.abs.bias.oob, 
                  GRF.coverage.oob, LLCF.coverage.oob, GRF.length.oob, LLCF.length.oob,
                  GRF.MSE, LLCF.MSE, GRF.abs.bias, LLCF.abs.bias,
                  GRF.coverage, LLCF.coverage, GRF.length, LLCF.length))
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
  print(results)
  colnames(results) <- c("In-sample RMSE GRF", "In-sample RMSE LLCF", "In-sample absolute bias GRF", "In-sample absolute bias LLCF",
                         "In-sample coverage GRF", "In-sample coverage LLCF", "In-sample length GRF", "In-sample length LLCF",
                         "Out-of-sample RMSE GRF", "Out-of-sample RMSE LLCF", "Out-of-sample absolute bias GRF", "Out-of-sample absolute bias LLCF",
                         "Out-of-sample coverage GRF", "Out-of-sample coverage LLCF", "Out-of-sample length GRF", "Out-of-sample length LLCF")
                         
  results$n.training = training.sample.size
  results
}  

results.simulation <- run_simulation(n.training, num.reps)




