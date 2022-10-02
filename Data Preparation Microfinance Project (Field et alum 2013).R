# Import and prepare Microfinance data set (Field et al., 2013)
  Grace_Period_Data <- read_dta("Downloads/112672-V1/Grace-Period-Data.dta")
  Grace_Period_Data <- as.data.frame(Grace_Period_Data)

  # Appoint treatment assignment and outcome variables
    W <- Grace_Period_Data$sec_treat
    W <- as.vector(W)

  # Create a numerical vector of the character group name vector
    loangroups <- as.numeric(factor(Grace_Period_Data$sec_group_name))

  # Standardize the continuous variables
    Grace_Period_Data$Years_Education_C <- scale(Grace_Period_Data$Years_Education_C)
    Grace_Period_Data$Age_C <- scale(Grace_Period_Data$Age_C)
    Grace_Period_Data$SEI <- scale(Grace_Period_Data$SEI)
    Grace_Period_Data$HH_Size_C <- scale(Grace_Period_Data$HH_Size_C)

  # Appoint the control variables matrix
    X <- Grace_Period_Data[,c(2,6,8:18)]
    colnames(X) <- c("Loan.Officer", "Stratification", "Age", "Married", "Literate", "Muslim", "HH.Size", "Years.Education", "Shock", "Has.Business",
                 "Financial.Control", "Home.Owner", "No.Drain")
                 
  # Create loan group dummies (as in original analysis) to be added to control variables matrix
    loansize1 <- as.integer(c(Grace_Period_Data$sec_loanamount == 4000))
    loansize2 <- as.integer(c(Grace_Period_Data$sec_loanamount == 5000))
    loansize3 <- as.integer(c(Grace_Period_Data$sec_loanamount == 6000))
    loansize4 <- as.integer(c(Grace_Period_Data$sec_loanamount == 8000))
    loansize5 <- as.integer(c(Grace_Period_Data$sec_loanamount == 9000))
    loansize6 <- as.integer(c(Grace_Period_Data$sec_loanamount == 10000))
    loansizematrix <- cbind(loansize1, loansize2, loansize3, loansize4,
                        loansize5, loansize6)

  # Combine all the control covariates in one large matrix
    X <- as.matrix(cbind(X, loansizematrix))
