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
n <- 1000
d = 12

# Create the training and test data sets
dat <- simulation(n = n, d = d, sigma = s)
tau.training <- dat$Tau

###########################
########### GRF ###########
###########################

GRF <- causal_forest(as.matrix(dat$X), as.vector(dat$Y), as.vector(dat$W), honesty = TRUE, num.trees = 4000)
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
forest.W <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$W), honesty = TRUE, 
                                 honesty.fraction = 0.7, honesty.prune.leaves = FALSE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
W.hat <- predict(forest.W)$predictions
selected.W <- order(variable_importance(forest.W), decreasing = TRUE)[1:5]
forest.Y <- ll_regression_forest(as.matrix(dat$X), as.numeric(dat$Y), honesty = TRUE, 
                                 honesty.fraction = 0.7, honesty.prune.leaves = FALSE, enable.ll.split = TRUE, ll.split.weight.penalty = TRUE, num.trees = 4000)
Y.hat <- predict(forest.Y)$predictions
selected.Y <- order(variable_importance(forest.Y), decreasing = TRUE)[1:5]

# Implement LLCF
LLCF <- causal_forest(as.matrix(dat$X), as.numeric(dat$Y), as.numeric(dat$W), Y.hat = Y.hat, W.hat = W.hat, honesty = TRUE, 
                      honesty.fraction = 0.7, honesty.prune.leaves = FALSE, num.trees = 4000, tune.parameters = c("min.node.size", "sample.fraction", "mtry", "alpha", "imbalance.penalty"))

# Run LLCF algorithm on test data set, obtain predictions, and compute measures
lasso.mod <- cv.glmnet(as.matrix(dat$X), as.numeric(dat$Y), alpha = 1)
selected <- which(coef(lasso.mod) != 0)
if(length(selected) < 2) {
  selected <- 1:ncol(dat$X)
} else {
  selected <- selected[-1] - 1 # Remove intercept
}

selected <- c(selected, selected.W, selected.Y)
selected <- sort(selected[!duplicated(selected)])

# Make predictions
LLCF.pred.oob <- predict(LLCF, linear.correction.variables = selected, ll.weight.penalty = TRUE, estimate.variance = TRUE)
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

tau.for.plot <- tau.training
X.for.plot <- dat$X
GRF.for.plot <- GRF.CATE.oob
LLCF.for.plot <- LLCF.CATE.oob
measures.for.plot <- c(GRF.MSE.oob, LLCF.MSE.oob, GRF.abs.bias.oob,  LLCF.abs.bias.oob, 
                       GRF.coverage.oob, LLCF.coverage.oob, GRF.length.oob, LLCF.length.oob)

# Heat map

minp = min(tau.for.plot, GRF.for.plot, LLCF.for.plot)
maxp = max(tau.for.plot, GRF.for.plot, LLCF.for.plot)
rngp = maxp - minp

ncol = 100

true.scl = pmax(ceiling(ncol * (tau.for.plot - minp) / rngp), 1)
grf.scl = pmax(ceiling(ncol * (GRF.for.plot - minp) / rngp), 1)
llcf.scl = pmax(ceiling(ncol * (LLCF.for.plot - minp) / rngp), 1)

hc = heat.colors(ncol)

png("Heatmap Simulation Set-Up 1 n = 1000 d = 12 training.png", 
    width = 350, height = 150)
par(mfrow = c(1, 3))

plot(X.for.plot[,1], X.for.plot[,2], pch = 16, col = hc[true.scl], xlab = "X1", ylab = "X2", main = "")

plot(X.for.plot[,1], X.for.plot[,2], pch = 16, col = hc[grf.scl], xlab = "X1", ylab = "X2", main = "")

plot(X.for.plot[,1], X.for.plot[,2], pch = 16, col = hc[llcf.scl], xlab = "X1", ylab = "X2", main = "")
dev.off()
