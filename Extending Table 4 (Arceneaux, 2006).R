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
set.seed(123)

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

# Clean the data
  treatment.group <- sample(which(voting.study.data$treat_real == 1), size = length(which((voting.study.data$treat_real == 1))), replace = FALSE)
  control.group <- which(voting.study.data$treat_real == 0)
  subsample.control.group <- sample(control.group, size = 89958, replace = FALSE)
  X.treatment.group <- X[treatment.group,]
  X.control.group <- X[subsample.control.group,]
  W.treatment.group <- W[treatment.group]
  W.control.group <- W[subsample.control.group]
  Ynew.treatment.group <- Ynew[treatment.group]
  Ynew.control.group <- Ynew[subsample.control.group]
  county.clusters.treatment.group <- county.clusters[treatment.group]
  county.clusters.control.group <- county.clusters[subsample.control.group]
  dummies.county.clusters.treatment.group <- dummies.county.clusters[treatment.group,]
  dummies.county.clusters.control.group <- dummies.county.clusters[subsample.control.group,]

# Appoint training and test sets
  X.train <- rbind(X.treatment.group[1:40000,], X.control.group[1:60000,])
  X.test <- rbind(X.treatment.group[40001:50000,], X.control.group[60001:75000,])
  X.holdout <- rbind(X.treatment.group[50001:nrow(X.treatment.group),], X.control.group[75001:nrow(X.control.group),])
  W.train <- c(W.treatment.group[1:40000], W.control.group[1:60000])
  W.test <- c(W.treatment.group[40001:50000], W.control.group[60001:75000])
  W.holdout <- c(W.treatment.group[50001:length(W.treatment.group)], W.control.group[75001:length(W.control.group)])
  Ynew.train <- c(Ynew.treatment.group[1:40000], Ynew.control.group[1:60000])
  Ynew.test <- c(Ynew.treatment.group[40001:50000], Ynew.control.group[60001:75000])
  Ynew.holdout <- c(Ynew.treatment.group[50001:length(Ynew.treatment.group)], Ynew.control.group[75001:length(Ynew.control.group)])
  county.clusters.train <- c(county.clusters.treatment.group[1:40000], county.clusters.control.group[1:60000])
  county.clusters.test <- c(county.clusters.treatment.group[40001:50000], county.clusters.control.group[60001:75000])
  county.clusters.holdout <- c(county.clusters.treatment.group[50001:length(county.clusters.treatment.group)], county.clusters.control.group[75001:length(county.clusters.control.group)])
  dummies.county.clusters.train <- rbind(dummies.county.clusters.treatment.group[1:40000,], dummies.county.clusters.control.group[1:60000,])
  dummies.county.clusters.test <- rbind(dummies.county.clusters.treatment.group[40001:50000,], dummies.county.clusters.control.group[60001:75000,])
  dummies.county.clusters.holdout <- rbind(dummies.county.clusters.treatment.group[50001:nrow(dummies.county.clusters.treatment.group),], dummies.county.clusters.control.group[75001:nrow(dummies.county.clusters.control.group),])

# Initialize parameters
  numtrees <- 2000 # Set to 1000 or 5000 to perform sensitivity analysis.
  lambdas <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5) # Concerns ridge penalty parameters; do not adjust.
  boolean.lambdas <- FALSE # Set to TRUE to use lambdas instead of automatic penalty tuning.
  boolean.plot <- FALSE # Set to TRUE to make various plots of interest.
  filename.plot.GRF.CATE <- "GRF CATE .pdf"
  filename.plot.CR.GRF.CATE <- "CR.GRF CATE .pdf"
  filename.plot.LLCF.CATE <- "LLCF CATE .pdf"

# Estimation procedure

    ###########################
    ########### GRF ###########
    ###########################

      # For GRF we create the X matrix including district-specific dummies
        current.X <- cbind(X.train, dummies.county.clusters.train)
        current.X.test <- cbind(X.test, dummies.county.clusters.test)

      # Select variables to include using preliminary GRF
        prelim.GRF <- causal_forest(current.X, Ynew.train, W.train, num.trees = numtrees, honesty = TRUE)
        prelim.GRF.varimp <- variable_importance(prelim.GRF)
        selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)

      # Implement GRF
        GRF <- causal_forest(current.X[,selected.vars], Ynew.train, W.train, num.trees = numtrees, honesty = TRUE)

    ############################
    #### Cluster-Robust GRF ####
    ############################
      
      # Select variables to include using preliminary Cluster-Robust GRF
        prelim.CR.GRF <- causal_forest(X.train, Ynew.train, W.train, honesty = TRUE, clusters = county.clusters.train, num.trees = numtrees)
        prelim.CR.GRF.varimp <- variable_importance(prelim.CR.GRF)
        selected.vars <- which(prelim.CR.GRF.varimp / mean(prelim.CR.GRF.varimp) > 0.2)
        
      # Implement Cluster-Robust GRF
        CR.GRF <- causal_forest(X.train[,selected.vars], Ynew.train, W.train, clusters = county.clusters.train, honesty = TRUE, num.trees = numtrees)

    ############################
    ########### LLCF ###########
    ############################

      # Select variables to include using Lasso feature selection
        lasso.mod <- cv.glmnet(as.matrix(current.X[,selected.vars]), Ynew.train, alpha = 1)
        selected <- which(coef(lasso.mod) != 0)
        if(length(selected) < 2) {
          selected <- 1:ncol(current.X)
        } else {
          selected <- selected[-1] - 1 # Remove intercept
        }
        LLCF <- GRF

      # Compute HTE with corresponding 95% confidence intervals                                  
        if (boolean.lambdas == FALSE) {
        # Predict: tuning without grid search over lambdas
          LLCF.pred <- predict(LLCF, current.X.test[,selected.vars], linear.correction.variables = selected, ll.weight.penalty = TRUE, ll.lambda = 0.1, estimate.variance = TRUE)
          LLCF.CATE <- LLCF.pred$predictions
          LLCF.CATE.SE <- sqrt(LLCF.pred$variance.estimates)
        } else {
          # Predict: tuning done using set of lambdas
          LLCF.mse.old <- +Inf
          for (l in length(lambdas)) {
            LLCF.CATE.old <- predict(LLCF, X.test, linear.correction.variables = selected, ll.lambda = lambdas[l], ll.weight.penalty = TRUE, estimate.variance = TRUE)
            predictions <- LLCF.CATE.old$predictions
            LLCF.mse.new <- mean((predictions - mean(predictions))**2)
            if (LLCF.mse.new < LLCF.mse.old) {
              LLCF.mse.old <- LLCF.mse.new
              LLCF.CATE.SE <- sqrt(LLCF.CATE.old$variance.estimates)
              LLCF.CATE <- predictions
            }
          }
        }



