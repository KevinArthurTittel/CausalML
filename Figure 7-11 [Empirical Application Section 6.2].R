library(tidyverse)
install.packages("hexbin")
library(hexbin)
install.packages("car")
library(car)
library(haven)
library(corrplot)
library("Hmisc")

########################
### Data preparation ###
########################

# Import and prepare Microfinance data set (Field et al., 2013)
Grace_Period_Data <- read_dta("Downloads/112672-V1/Grace-Period-Data.dta")
Grace_Period_Data <- as.data.frame(Grace_Period_Data)

# Appoint treatment assignment and outcome variables
W <- Grace_Period_Data$sec_treat
W <- as.vector(W)

# Create a numerical vector of the character group name vector
loangroups <- as.numeric(factor(Grace_Period_Data$sec_loanamount))

# Appoint the control variables matrix
X <- Grace_Period_Data[,c(2,6,8:18)]

# Create loan group dummies (as in original analysis) to be added to control variables matrix
loansize1 <- as.integer(c(Grace_Period_Data$sec_loanamount == 4000))
loansize2 <- as.integer(c(Grace_Period_Data$sec_loanamount == 5000))
loansize3 <- as.integer(c(Grace_Period_Data$sec_loanamount == 6000))
loansize4 <- as.integer(c(Grace_Period_Data$sec_loanamount == 8000))
loansize5 <- as.integer(c(Grace_Period_Data$sec_loanamount == 9000))
loansize6 <- as.integer(c(Grace_Period_Data$sec_loanamount == 10000))
loansizematrix <- cbind(loansize1, loansize2, loansize3, loansize4,
                        loansize5, loansize6)

X <- cbind(X, loansizematrix)
colnames(X) <- c("Loan.Officer", "Stratification", "Age", "Married", "Literate", "Muslim", "HH Size", "Years of Education", "Shock", "Has Business",
                 "Financial Control", "Home Owner", "No Drain", "Rs 4000 loan", "Rs 5000 loan",
                 "Rs 6000 loan", "Rs 8000 loan", "Rs 9000 loan", "Rs 10000 loan")

####################
##### Figure 7 #####
####################

hist(X$Age, xlab = "Age", ylab = "Number of observations", xlim = c(15, 60), ylim = c(0, 200), main = "", col = "orange")
hist(X$`HH Size`, xlab = "Household Size", ylab = "Number of observations", xlim = c(0, 12), ylim = c(0, 350), main = "", col = "orange")
hist(X$`Years of Education`, xlab = "Years of Education", ylab = "Number of observations", xlim = c(0, 15), ylim = c(0, 200), main = "", col = "orange")

####################
##### Figure 8 #####
####################

binary.data <- matrix(data = 0, nrow = 2, ncol = 6)
rownames(binary.data) <- c("No", "Yes")
colnames(binary.data) <- c("Muslim", "Literate", "No drain", "Married", "Financial control", "Home owner")
binary.data[1,1] <- sum(X$Muslim == 0)
binary.data[2,1] <- sum(X$Muslim == 1)
binary.data[1,2] <- sum(X$Literate == 0)
binary.data[2,2] <- sum(X$Literate == 1)
binary.data[1,3] <- sum(X$`No Drain` == 0)
binary.data[2,3] <- sum(X$`No Drain` == 1)
binary.data[1,4] <- sum(X$Married == 0)
binary.data[2,4] <- sum(X$Married == 1)
binary.data[1,5] <- sum(X$`Financial Control` == 0)
binary.data[2,5] <- sum(X$`Financial Control` == 1)
binary.data[1,6] <- sum(X$`Home Owner` == 0)
binary.data[2,6] <- sum(X$`Home Owner` == 1)
barplot(as.matrix(binary.data), beside = TRUE, main = "", col = c("red", "orange"), xlab = "Binary variables", ylab = "Number of observations")
legend(x = 15.5, y = 865, c("Yes", "No"), fill = c("orange", "red"), box.lty = 0, cex = 0.75)

####################
##### Figure 9 #####
####################

ggplot(data = Grace_Period_Data, mapping = aes(x = Years_Education_C, y = Age_C)) +
  geom_hex() +
  guides(fill = guide_colourbar(title = "Obs")) +
  scale_fill_gradient(low="orangered4",high="orangered",trans="log10") +
  xlab("Years of education") +
  ylab("Age") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour="grey"),
        panel.grid.major = element_line(colour = "grey", linetype = "solid", size = 0.10),
        panel.grid.minor = element_line(colour = "grey", linetype = "solid", size = 0.05),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10))

ggplot(data = Grace_Period_Data, mapping = aes(x = Years_Education_C, y = HH_Size_C)) +
  geom_hex() +
  guides(fill = guide_colourbar(title = "Obs")) +
  scale_fill_gradient(low="orangered4",high="orangered",trans="log10") +
  xlab("Years of education") +
  ylab("Household size") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour="grey"),
        panel.grid.major = element_line(colour = "grey", linetype = "solid", size = 0.10),
        panel.grid.minor = element_line(colour = "grey", linetype = "solid", size = 0.05),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10))

ggplot(data = Grace_Period_Data, mapping = aes(x = Years_Education_C, y = Financial_Control_C)) +
  geom_hex() +
  guides(fill = guide_colourbar(title = "Obs")) +
  scale_fill_gradient(low="orangered4",high="orangered",trans="log10") +
  xlab("Years of education") +
  ylab("Financial control") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour="grey"),
        panel.grid.major = element_line(colour = "grey", linetype = "solid", size = 0.10),
        panel.grid.minor = element_line(colour = "grey", linetype = "solid", size = 0.05),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10))

#####################
##### Figure 10 #####
#####################

mydata.cor <- cor(X)
mydata.p <- round(mydata.p, 3)
mydata.p[is.na(mydata.p) == TRUE] <- 0
heatmap(as.matrix(mydata.p), Colv = NA, Rowv = NA, cexRow=0.90, cexCol = 0.60)

#####################
##### Figure 11 #####
#####################

# Create graphs to illustrate smoothness
Grace_Period_Data <- Grace_Period_Data[!(Grace_Period_Data$Age_C == 0),]
Grace_Period_Data <- Grace_Period_Data[!(Grace_Period_Data$Profit > 50000),]
ggplot(data = Grace_Period_Data, mapping = aes(x = Age_C, y = Profit)) +
  geom_point() +
  xlab("Age") +
  ylab("Monthly profit (in Rs)") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour="grey"),
        panel.grid.major = element_line(colour = "grey", linetype = "solid", size = 0.10),
        panel.grid.minor = element_line(colour = "grey", linetype = "solid", size = 0.05),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = "none")

ggplot(data = Grace_Period_Data, mapping = aes(x = HH_Size_C, y = Profit)) +
  geom_point() +
  xlab("Household size") +
  ylab("Monthly profit (in Rs)") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour="grey"),
        panel.grid.major = element_line(colour = "grey", linetype = "solid", size = 0.10),
        panel.grid.minor = element_line(colour = "grey", linetype = "solid", size = 0.05),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = "none")

ggplot(data = Grace_Period_Data, mapping = aes(x = Years_Education_C, y = Profit)) +
  geom_point() +
  xlab("Years of education") +
  ylab("Monthly profit (in Rs)") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour="grey"),
        panel.grid.major = element_line(colour = "grey", linetype = "solid", size = 0.10),
        panel.grid.minor = element_line(colour = "grey", linetype = "solid", size = 0.05),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = "none")

