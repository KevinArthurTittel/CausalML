set.seed(123)

# Import and prepare Voting Study data set (Arceneaux et al., 2006)
voting.study.data <- as.data.frame(read_dta("Downloads/ArceneauxGerberGreen_PA_2006_IA_MI_merge040504.dta"))
# voting.study.data <- read_dta("C:/Users/481044kt/Downloads/ArceneauxGerberGreen_PA_2006_IA_MI_merge040504.dta")

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
# tau <- -0.5* (X$vote00 / (1 + (50/X$age)))
# tau <- -1* (X$vote00 / (2 + (100/X$age)))
# tau <- -1* (X$vote00 / (1 + exp(-1*(X$age)/50)))
tau <- -1* (X$vote00 / (1.5 + (1000/X$age^2)))
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

# Determine treatment and control groups s.t. 2/5th is treated (59972) and 3/5th is controlled (89958), total 149930
treatment.group <- sample(which(voting.study.data$treat_real == 1), size = length(which(voting.study.data$treat_real == 1)), replace = FALSE) 
subsample.control.group <- sample(which(voting.study.data$treat_real == 0), size = 89958, replace = FALSE)
combined <- c(treatment.group, subsample.control.group)
training.subsample <- c(treatment.group[1:ceiling((0.8*149930*0.4))], subsample.control.group[1:ceiling((0.8*149930*0.6))])
test.subsample <- c(treatment.group[-(1:ceiling((0.8*149930*0.4)))], subsample.control.group[-(1:ceiling((0.8*149930*0.6)))])

################
### Figure 3 ###
################

hist(X$age, xlab = "Age", ylab = "Number of observations", xlim = c(15, 115), ylim = c(0, 230000), main = "", col = "orange")
age.ordered <- order(X$age)
plot(x = X[age.ordered,"age"], y = tau[age.ordered], type = "p", pch = 16, col = "orange", main = "",
     xlim = c(15,120), ylim = c(-0.70, 0.10), xlab = "Age", ylab = "True treatment effect")

################
### Figure 4 ###
################

Y0before <- Y[subsample.control.group]
Y1before <- Y[treatment.group]
Y0after <- Ynew[subsample.control.group]
Y1after <- Ynew[treatment.group]

binary.data <- matrix(data = 0, nrow = 2, ncol = 4)
rownames(binary.data) <- c("Original voting", "Synthetic voting")
colnames(binary.data) <- c("W = 0, Vote = 0", "W = 0, Vote = 1",
                           "W = 1, Vote = 0", "W = 1, Vote = 1")
binary.data[1,1] <- sum((Y0before == 0))
binary.data[2,1] <- sum((Y0after == 0))
binary.data[1,2] <- sum((Y0before == 1))
binary.data[2,2] <- sum((Y0after == 1))
binary.data[1,3] <- sum((Y1before == 0))
binary.data[2,3] <- sum((Y1after == 0))
binary.data[1,4] <- sum((Y1before == 1))
binary.data[2,4] <- sum((Y1after == 1))
barplot(as.matrix(binary.data), beside = TRUE, main = "", col = c("red", "orange"), xlab = "Votes", ylab = "Number of observations",
        ylim = c(0,70000))
legend(x = 9, y = 60000, c("Original voting", "Synthetic voting"), fill = c("red", "orange"), box.lty = 0, cex = 0.90)

################
### Figure 5 ###
################

hist(tau[combined], xlab = "True treatment effect", ylab = "Number of observations", xlim = c(-0.65,0.10), ylim = c(0,60000), main = "", col = "orange")

