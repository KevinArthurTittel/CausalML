num.reps <- 25
n.training <- c(200, 400, 600, 800, 1000, 1200)
d <- 25
numtrees = 4000

run_simulation <- function(n.training, num.reps) {
  results <- sapply(n.training, function(n) {
    print(n)
    basic.results <- replicate(num.reps, {
      
      # Create the training and test data sets
      dat <- simulation(n = n, d = d, sigma = s)
      tau.training <- dat$Tau
      
      ###########################
      ########### GRF ###########
      ###########################
      
      # Run GRF algorithm on test data set, obtain predictions, and compute measures
      GRF <- causal_forest(as.matrix(dat$X), as.vector(dat$Y), as.vector(dat$W), honesty = TRUE, num.trees = numtrees, tune.parameters = "all")
      GRF.pred.oob <- predict(GRF, estimate.variance = TRUE)
      GRF.CATE.oob <- GRF.pred.oob$predictions
      GRF.CATE.SE.oob <- sqrt(GRF.pred.oob$variance.estimates)
      
      # Compute measures
      GRF.MSE.oob <- sqrt(mean((GRF.CATE.oob - tau.training)**2))
      GRF.abs.bias.oob <- mean((abs(GRF.CATE.oob - tau.training)))
      GRF.coverage.oob <- mean(as.integer((((GRF.CATE.oob - qnorm(0.975)*GRF.CATE.SE.oob) <= tau.training) & 
                                             (tau.training <= (GRF.CATE.oob + qnorm(0.975)*GRF.CATE.SE.oob)))))
      GRF.length.oob <- mean(abs((GRF.CATE.oob + qnorm(0.975)*GRF.CATE.SE.oob) - (GRF.CATE.oob - qnorm(0.975)*GRF.CATE.SE.oob)))
      
      ############################
      ########### LLCF ###########
      ############################
      
      # Grow preliminary forests for (W, X) and (Y, X) separately
      forest.W <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$W), honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                       honesty.fraction = 0.7, honesty.prune.leaves = FALSE, ll.split.lambda = 0.1, num.trees = 2000)
      W.hat <- predict(forest.W)$predictions
      forest.Y <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$Y), honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                       honesty.fraction = 0.7, honesty.prune.leaves = FALSE, ll.split.lambda = 0.1, num.trees = 2000)
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
      LLCF.pred.oob <- predict(LLCF, linear.correction.variables = selected, ll.lambda = 1.5, ll.weight.penalty = TRUE, estimate.variance = TRUE)
      LLCF.CATE.oob <- LLCF.pred.oob$predictions
      LLCF.CATE.SE.oob <- sqrt(LLCF.pred.oob$variance.estimates)
      
      # Compute measures
      LLCF.MSE.oob <- sqrt(mean((LLCF.CATE.oob - tau.training)**2))
      LLCF.abs.bias.oob <- mean((abs(LLCF.CATE.oob - tau.training)))
      LLCF.coverage.oob <- mean(as.integer((((LLCF.CATE.oob - qnorm(0.975)*LLCF.CATE.SE.oob) <= tau.training) & 
                                              (tau.training <= (LLCF.CATE.oob + qnorm(0.975)*LLCF.CATE.SE.oob)))))
      LLCF.length.oob <- mean(abs((LLCF.CATE.oob + qnorm(0.975)*LLCF.CATE.SE.oob) - (LLCF.CATE.oob - qnorm(0.975)*LLCF.CATE.SE.oob)))
      
      ############################
      ###### Update results ######
      ############################
      
      return(c(GRF.MSE.oob, LLCF.MSE.oob, GRF.abs.bias.oob,  LLCF.abs.bias.oob, 
               GRF.coverage.oob, LLCF.coverage.oob, GRF.length.oob, LLCF.length.oob))
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
                         "In-sample coverage GRF", "In-sample coverage LLCF", "In-sample length GRF", "In-sample length LLCF")
  
  results$n.training = n.training
  results
}  

sigma = 1

###############
### Table 7 ###
###############

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- 0.5
  main_effect <- 0
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 2/(1 + exp(-12*(X[,1] - (1/2))))
  zeta2 <- 2/(1 + exp(-12*(X[,2] - (1/2))))
  zeta3 <- 2/(1 + exp(-12*(X[,3] - (1/2))))
  zeta4 <- 2/(1 + exp(-12*(X[,4] - (1/2))))
  Tau <- matrix(((zeta1 * zeta2) + (zeta3 * zeta4)), nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}
results.simulation <- run_simulation(n.training, num.reps)

###############
### Table 8 ###
###############

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- pmax(0.1, pmin(sin(pi*X[,1]*X[,2]), 0.9))
  main_effect <- sin(pi*X[,1]*X[,2]) + 2*(X[,3] - 0.5)^2 + X[,4] + 0.5*X[,5]
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  Tau <- matrix(((X[,1] + X[,2])/2), nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}
results.simulation <- run_simulation(n.training, num.reps)

###############
### Table 9 ###
###############

simulation <- function(n, d, sigma) {
  X <- matrix(rnorm(n*d, mean = 0, sd = 1), nrow = n, ncol = d)
  treatment_propensity <- 0.5
  main_effect <- pmax(X[,1] + X[,2] + X[,3],0) + pmax(X[,4] + X[,5],0)
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  Tau <- matrix((X[,1] + log(1 + exp(X[,2]))), nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}
results.simulation <- run_simulation(n.training, num.reps)

sigma = 2

##################
#### Table 10 ####
##################

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- 0.5
  main_effect <- 0
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 2/(1 + exp(-12*(X[,1] - (1/2))))
  zeta2 <- 2/(1 + exp(-12*(X[,2] - (1/2))))
  zeta3 <- 2/(1 + exp(-12*(X[,3] - (1/2))))
  zeta4 <- 2/(1 + exp(-12*(X[,4] - (1/2))))
  Tau <- matrix(((zeta1 * zeta2) + (zeta3 * zeta4)), nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}
results.simulation <- run_simulation(n.training, num.reps)

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- pmax(0.1, pmin(sin(pi*X[,1]*X[,2]), 0.9))
  main_effect <- sin(pi*X[,1]*X[,2]) + 2*(X[,3] - 0.5)^2 + X[,4] + 0.5*X[,5]
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  Tau <- matrix(((X[,1] + X[,2])/2), nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}
results.simulation <- run_simulation(n.training, num.reps)

simulation <- function(n, d, sigma) {
  X <- matrix(rnorm(n*d, mean = 0, sd = 1), nrow = n, ncol = d)
  treatment_propensity <- 0.5
  main_effect <- pmax(X[,1] + X[,2] + X[,3],0) + pmax(X[,4] + X[,5],0)
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  Tau <- matrix((X[,1] + log(1 + exp(X[,2]))), nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}
results.simulation <- run_simulation(n.training, num.reps)



