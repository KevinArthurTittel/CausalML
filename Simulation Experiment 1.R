library(grf)
library(glmnet)
library(hte)

# Simulation setting Athey and Wager (2019), Section 6.2 Table 1

simulation <- function(n = n, d = d) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  noise <- rnorm(n=n)
  main_effect <- 2*X[,3] - 1
  treatment_propensity <- 0.25*(1 + dbeta(X[,3], 2, 4))
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 1 + 1/(1 + exp(-20*(X[,1] - (1/3))))
  zeta2 <- 1 + 1/(1 + exp(-20*(X[,2] - (1/3))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau)
  return(output)
}

num.reps <- 50
n.training <- c(200, 400, 800, 1200, 2000)
dimensions <- c(2, 4, 6, 10, 15, 20, 25, 30)
lambdas = c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5)

run_simulation <- function(n.training, dimensions, lambdas, num.reps) {
  results.per.dimension <- sapply(dimensions, function(d) {
    results.per.trainingsamplesize <- sapply(ns, function(n) {
      basic.results <- replicate(num.reps, {
      
        # Create the training and test data sets
          dat <- simulation(n = n, d = d)
          dat.test <- simulation(n = 2000, d = d)
          truth <- dat.test$Tau # True tau
      
        # Initialize the measures of interest to compare the models' performances
          GRF.coverage <- 0 # Average coverage of true tau in all 95% confidence intervals amongst 2000 test points
          GRF.length <- 0 # Average 95% confidence interval length for tau amongst 2000 test points
          GRF.mse <- 0 # MSE of estimator for tau amongst 2000 test points
          GRF.abs.bias <- 0 # Absolute bias of estimator for tau amongst 2000 test points
          LLCF.coverage <- 0
          LLCF.length <- 0
          LLCF.mse <- 0
          LLCF.abs.bias <- 0
      
    ###########################
    ########### GRF ###########
    ###########################
     
        # Run GRF algorithm on test data set, obtain predictions, and compute measures
          GRF <- causal_forest(as.matrix(dat$X), as.numeric(dat$Y), as.numeric(dat$W), tune.parameters = "all")
          GRF.pred <- predict(GRF, as.matrix(dat.test$X), estimate.variance = TRUE)
          GRF.CATE <- GRF.pred$predictions
          GRF.CATE.SE <- sqrt(GRF.pred$variance.estimates)
          GRF.coverage <- mean(as.integer((((GRF.CATE - qnorm(0.975)*GRF.CATE.SE) <= truth) & 
                                         (truth <= (GRF.CATE + qnorm(0.975)*GRF.CATE.SE)))))
          GRF.length <- mean(abs((GRF.CATE + 1.96*GRF.CATE.SE) - (GRF.CATE - 1.96*GRF.CATE.SE)))
          GRF.mse <- mean((GRF.CATE - truth)**2)
          GRF.abs.bias <- mean((abs(GRF.CATE - truth)))
   
    ############################
    ########### LLCF ###########
    ############################
    
        # Run LLCF algorithm on test data set, obtain predictions, and compute measures
          lasso.mod <- cv.glmnet(as.matrix(dat$X), as.numeric(dat$Y), alpha = 1)
          selected <- which(coef(lasso.mod) != 0)
          if(length(selected) < 2) {
            selected <- 1:ncol(dat$X)
          }
          else {
            selected <- selected[-1] - 1 # Remove intercept
          }
      
          LLCF.mse.old <- +Inf
          for (l in length(lambdas)) {
            LLCF.CATE.old <- predict(LLCF, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
            predictions <- LLCF.CATE.old$predictions
            LLCF.mse.new <- mean((predictions - mean(predictions))**2)
            if (LLCF.mse.new < LLCF.mse.old) {
              LLCF.mse.old <- LLCF.mse.new
              LLCF.CATE <- predictions
              LLCF.CATE.SE <- sqrt(LLCF.CATE.old$variance.estimates)
            }
          }
      
          LLCF.coverage <- mean(as.integer((((LLCF.CATE - 1.96*LLCF.CATE.SE) <= truth) & 
                                          (truth <= (LLCF.CATE + 1.96*LLCF.CATE.SE)))))
          LLCF.length <- mean(abs((LLCF.CATE + 1.96*LLCF.CATE.SE) - (LLCF.CATE - 1.96*LLCF.CATE.SE)))  
          LLCF.mse <- mean((LLCF.CATE - truth)**2)
          LLCF.abs.bias <- mean(abs(LLCF.CATE - truth))
      
    ############################
    ###### Update results ######
    ############################
    
        return(c(d, n, sqrt(GRF.mse), sqrt(LLCF.mse), GRF.abs.bias, LLCF.abs.bias, 
               GRF.length, LLCF.length, GRF.coverage, LLCF.coverage))
      })
      return(colMeans(data.frame(t(basic.results))))
    })
    return(results.per.trainingsamplesize)
  )} 
  
  results <- data.frame(results.per.dimension)
  colnames(results) <- c("d", "n", "CF MSE", "LLCF MSE", "CF Abs Bias", "LLCF Abs Bias", 
                         "CF Length", "LLCF Length", "CF Coverage", "LLCF Coverage")
  return(results)
}  

results <- run_simulation(n.training, dimensions, lambdas, num.reps)




