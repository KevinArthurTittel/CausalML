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
  
  
  
  
  
  
  
  
  
