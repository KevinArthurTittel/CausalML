library(rlearner)
library(glmnet)
library(dbarts)
library(RColorBrewer)
set.seed(123)

# Import data
  data = read.csv("data_clean.csv")
  
# Indices, and training/test/holdout sample sizes
  idx.all = 1:nrow(data)
  n.train = 100000
  n.test = 25000
  n.holdout = length(idx.all) - n.train - n.test
  
# Create synthetic smooth heterogeneous treatment effect function 
  X.all = as.matrix(data[idx.all,1:11])
  Y.obs.all = data[idx.all,12]
  W.all = data[idx.all,13]
  TAU.all = - X.all[,"vote00"] * 0.5 / (1 + 50 / X.all[,"age"])

  FLIP = rbinom(length(idx.all), 1, abs(TAU.all))
  Y.PO = t(sapply(1:length(idx.all), function(ii) {
    if (FLIP[ii] == 0) {
      return(rep(Y.obs.all[ii], 2))
    } else if (TAU.all[ii] > 0) {
      return(c(0, 1))
    } else {
      return(c(1, 0))
    }
  }))
  
  Y.all = sapply(1:length(idx.all), function(ii) {
    Y.PO[ii, W.all[ii] + 1]
  })
  
# Train data
  X = X.all[1:n.train,]
  W = W.all[1:n.train]
  Y = Y.all[1:n.train]

# Test data
  X.test = X.all[n.train + 1:n.test,]
  W.test = W.all[n.train + 1:n.test]
  Y.test = Y.all[n.train + 1:n.test]
  TAU.test = TAU.all[n.train + 1:n.test]

# Holdout data
  X.holdout = X.all[n.train + n.test + 1:n.holdout,]
  W.holdout = W.all[n.train + n.test + 1:n.holdout]
  Y.holdout = Y.all[n.train + n.test + 1:n.holdout]
  
# Fit propensity model
  W.boost = cvboost(X, W, objective = "binary:logistic", nthread = 4)
  W.hat.boost = predict(W.boost)

  W.lasso = cv.glmnet(X, W, family = "binomial", keep = TRUE)
  W.hat.lasso = W.lasso$fit.preval[,!is.na(colSums(W.lasso$fit.preval))][, W.lasso$lambda == W.lasso$lambda.min]

  W.hat = W.hat.boost

# Fit marginal response model
  Y.boost = cvboost(X, Y, objective = "binary:logistic", nthread = 4)
  Y.hat.boost = predict(Y.boost)

  Y.lasso = cv.glmnet(X, Y, keep = TRUE, family = "binomial")
  Y.hat.lasso = Y.lasso$fit.preval[,!is.na(colSums(Y.lasso$fit.preval))][, Y.lasso$lambda == Y.lasso$lambda.min]

  Y.hat = Y.hat.boost
  
# Fit R-learner given chosen nuisance components
  
  tau.boost = rboost(X, W, Y, p_hat = W.hat, m_hat = Y.hat, nthread = 4)
  tau.hat.boost = predict(tau.boost)

  tau.lasso = rlasso(X, W, Y, p_hat = W.hat, m_hat = Y.hat)
  tau.hat.lasso = predict(tau.lasso)

  # lasso wins on holdout
  tau.hat.boost.holdout = predict(tau.boost, X.holdout)
  tau.hat.lasso.holdout = predict(tau.lasso, X.holdout)

  Y.hat.holdout = predict(Y.boost, X.holdout)
  W.hat.holdout = predict(W.boost, X.holdout)

  tau.hat.boost.test = predict(tau.boost, X.test)
  tau.hat.lasso.test = predict(tau.lasso, X.test)
  
  print(round(c(mean((TAU.test - tau.hat.boost.test)^2),
              mean((TAU.test - tau.hat.lasso.test)^2)), 4))

# Make summaries
  cols = brewer.pal(3, "Set2")

  pdf("tau_histogram.pdf")
  pardef = par(mar = c(4.5, 5, 3, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  hist(TAU.all[1:n.train], main = "", xlab = "CATE", col = cols[1])
  dev.off()

  TAU.bucket = as.numeric(cut(TAU.test, breaks = c(quantile(TAU.test[TAU.test < 0], c(0, 0.25, 0.5, 0.75, 1)) - 0.001, 0.001))) - 1
  TAU.mids = sapply(0:4, function(ii) median(TAU.test[TAU.bucket == ii]))
  ATmat = rbind(TAU.mids, TAU.mids)
  ATmat = ATmat + c(-0.005, 0.005)
  ATvec = c(ATmat)

  BPDF = data.frame(CATE=c(2 * TAU.bucket, 1 + 2 * TAU.bucket),
                  prediction=c(tau.hat.lasso.test, tau.hat.boost.test))
 
  pdf("tau_boxplot.pdf")
  pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 8) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  boxplot(prediction ~ CATE, data = BPDF, col = cols[2:3],
        at = ATvec,  pars = list(boxwex = 0.007), xlim = range(ATvec),
        names = rep("", 10), xlab = "CATE", ylab = "prediction", xaxt = "n")
  axis(1, at = round(TAU.mids, 2))
  abline(0, 1, lwd = 2)

  par(xpd = TRUE)
  legend(0.03, 0.1, c("lasso", "boost"), fill = cols[2:3], cex = 1.5)
  par = pardef
  dev.off()

  save.image("analysis_results.RData")







  
