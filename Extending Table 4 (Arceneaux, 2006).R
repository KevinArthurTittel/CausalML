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
treatment.group <- (voting.study.data$treat_real == 1)
control.group <- (voting.study.data$treat_real == 0)
subsample.control.group <- sample(control.group, size = 89958, replace = FALSE)

# Synthetic treatment assignment


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
current.X <- cbind(X, dummies.county.clusters)

# Grow preliminary forests for (W, X) and (Y, X) separately
forest.W <- regression_forest(current.X, W, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
W.hat <- predict(forest.W)$predictions
forest.Y <- regression_forest(current.X, Y, num.trees = numtrees, honesty = TRUE, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

# Compute the variable importance
GRF.varimp <- variable_importance(forest.Y) 
GRF.varimp.ordered <- order(GRF.varimp)
GRF.mostimportant <- colnames(X)[(GRF.varimp.ordered[1:4])] # 4 most important variables for splitting

# Select variables to include using preliminary GRF
prelim.GRF <- causal_forest(current.X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, honesty = TRUE)
prelim.GRF.varimp <- variable_importance(prelim.GRF)
selected.vars <- which(prelim.GRF.varimp / mean(prelim.GRF.varimp) > 0.2)

# Implement GRF
GRF <- causal_forest(current.X[,selected.vars], Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = numtrees, 
                     honesty = TRUE, tune.parameters = "all")
