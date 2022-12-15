###############
### Table 3 ###
###############

s <- 1
num.reps <- 100
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
      GRF <- causal_forest(as.matrix(dat$X), as.vector(dat$Y), as.vector(dat$W), honesty = TRUE, num.trees = 2000, tune.parameters = "all")
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
    standard.deviation.basic.results <- sqrt(colVars(basic.results))
    mean.basic.results <- colMeans(basic.results)
    make.vector <- c()
    for (i in 1:length(mean.basic.results)) {
      make.vector <- cbind(make.vector, paste(round(mean.basic.results[i], 3), "(", round(standard.deviation.basic.results[i], 3), ")"))
    }
    print(make.vector)
    return(make.vector)
  })
  results = data.frame(t(results))
  print(results)
  colnames(results) <- c("d", "In-sample RMSE GRF", "In-sample RMSE LLCF", "In-sample absolute bias GRF", "In-sample absolute bias LLCF",
                         "In-sample coverage GRF", "In-sample coverage LLCF", "In-sample length GRF", "In-sample length LLCF")
  results
}  

results.simulation <- run_simulation(dimensions, num.reps)
