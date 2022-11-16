library(grf)
library(glmnet)
library(hte)
set.seed(123)

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

s <- 1
num.reps <- 50
n <- 2000
dimensions <- c(6,10,15,20,25,30)

run_simulation <- function(dimensions, num.reps) {
  results <- sapply(dimensions, function(d) {
    print(d)
    basic.results <- replicate(num.reps, {
      
      # Create the training and test data sets
      dat <- simulation(n = n, d = d, sigma = s)
      tau.training <- dat$Tau
      
      ###########################
      ########### GRF ###########
      ###########################
      
      # Run GRF algorithm on test data set, obtain predictions, and compute measures
      GRF <- causal_forest(as.matrix(dat$X), as.vector(dat$Y), as.vector(dat$W), honesty = TRUE, num.trees = 2000)
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
                                       honesty.fraction = 0.8, honesty.prune.leaves = FALSE, ll.split.lambda = 0.1, num.trees = 2000)
      W.hat <- predict(forest.W)$predictions
      forest.Y <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$Y), honesty = TRUE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, 
                                       honesty.fraction = 0.8, honesty.prune.leaves = FALSE, ll.split.lambda = 0.1, num.trees = 2000)
      Y.hat <- predict(forest.Y)$predictions
      
      # Implement LLCF
      LLCF <- causal_forest(as.matrix(dat$X), as.numeric(dat$Y), as.numeric(dat$W),  Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, 
                            honesty.fraction = 0.8, honesty.prune.leaves = FALSE, num.trees = 2000, tune.parameters = c("min.node.size", "sample.fraction", "mtry", "alpha", "imbalance.penalty"))
      
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
      
      return(c(d, GRF.MSE.oob, LLCF.MSE.oob, GRF.abs.bias.oob,  LLCF.abs.bias.oob, 
               GRF.coverage.oob, LLCF.coverage.oob, GRF.length.oob, LLCF.length.oob))
    })
    
    basic.results = data.frame(t(basic.results))
    mean.basic.results <- colMeans(basic.results)
    return(mean.basic.results)
    
  })
  results = data.frame(t(results))
  print(results)
  colnames(results) <- c("d", "In-sample RMSE GRF", "In-sample RMSE LLCF", "In-sample absolute bias GRF", "In-sample absolute bias LLCF",
                         "In-sample coverage GRF", "In-sample coverage LLCF", "In-sample length GRF", "In-sample length LLCF")
  results
}  

results.simulation <- run_simulation(dimensions, num.reps)

plot(x = dimensions, y = results.simulation[,2], type = "o", pch = 1, col = "red", main = "",
     xlim = c(5,31), ylim = c(0.28, 0.48), xlab = "Dimension", ylab = "RMSE")
lines(x = dimensions, y = results.simulation[,3], type = "o", pch = 22, col = "orange")
legend(x = 24, y = 0.33, box.lty = 0, cex = 0.75, pch = c(1,22), legend = c("GRF", "LLCF"),
       lwd = 1.5, col = c("red", "orange"), xpd = TRUE)

plot(x = dimensions, y = results.simulation[,4], type = "o", pch = 1, col = "red", main = "",
     ylim = c(0.23, 0.42), xlab = "Dimension", ylab = "Absolute Bias")
lines(x = dimensions, cex = 0.75, y = results.simulation[,5], type = "o", pch = 22, col = "orange")
legend(x = 24, y = 0.28, box.lty = 0, cex = 0.75, pch = c(1,22), legend = c("GRF", "LLCF"),
       lwd = 1.5, col = c("red", "orange"))

plot(x = dimensions, y = results.simulation[,6], type = "o", pch = 1, col = "red", main = "",
     ylim = c(0.30, 0.95), xlab = "Dimension", ylab = "Coverage rate")
lines(x = dimensions, y = results.simulation[,7], type = "o", pch = 22, col = "orange")
legend(x = 24, y = 0.45, box.lty = 0, cex = 0.75, pch = c(1,22), legend = c("GRF", "LLCF"),
       lwd = 1.5, col = c("red", "orange"))

plot(x = dimensions, y = results.simulation[,8], type = "o", pch = 1, col = "red", main = "",
     ylim = c(0.55, 1.25), xlab = "Dimension", ylab = "Interval length")
lines(x = dimensions, y = results.simulation[,9], type = "o", pch = 22, col = "orange")
legend(x = 24, y = 0.76, box.lty = 0, cex = 0.75, pch = c(1,22), legend = c("GRF", "LLCF"),
       lwd = 1.5, col = c("red", "orange"))
