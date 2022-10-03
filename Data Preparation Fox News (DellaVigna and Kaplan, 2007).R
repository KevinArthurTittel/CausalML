# Import and prepare Fox News data set (DellaVigna and Kaplan, 2007)
  Fox_News_Data <- read_dta("Downloads/FoxNewsDataQJEMay07/FoxNewsFinalDataQJE.dta")
  Fox_News_Data <- as.data.frame(Fox_News_Data)

  # Appoint main outcome variable
    Y <- as.vector(Fox_News_Data$reppresfv2p00m96)

  # Appoint treatment assignment variable
    W <- as.vector(Fox_News_Data$foxnews2000)

  # Create a numerical vector of the character group name vector
    district.clusters <- as.vector(Fox_News_Data$diststate)

  # Appoint the control variables matrix
    demographic.controls.2000 <- Fox_News_Data[,c(37:49)]
    colnames(demographic.controls.2000) <- c("Population 2000", "Population (over 18) 2000", "Share with high school 2000", "Share with some college 2000",
                                         "Share with college degree 2000", "Share male 2000", "Share African Americans 2000", "Share Hispanics 2000",
                                         "Employment rate 2000", "Unemployment rate 2000", "Share married 2000", "Median income 2000", "Share urban 2000")

    demographic.controls.Diff19962000 <- Fox_News_Data[,c(139:150)]
    colnames(demographic.controls.Diff19962000) <- c("Population, Diff. btwn. 2000 and 1996", "Share with high school, Diff. btwn. 2000 and 1996", "Share with some college, Diff. btwn. 2000 and 1996",
                                                 "Share with college degree, Diff. btwn. 2000 and 1996", "Share male, Diff. btwn. 2000 and 1996", "Share African Americans, Diff. btwn. 2000 and 1996",
                                                 "Share Hispanics, Diff. btwn. 2000 and 1996", "Employment rate, Diff. btwn. 2000 and 1996", "Unemployment rate, Diff. btwn. 2000 and 1996",
                                                 "Share married, Diff. btwn. 2000 and 1996", "Median income, Diff. btwn. 2000 and 1996", "Share urban, Diff. btwn. 2000 and 1996")

    Decile1NumberChannels <- as.integer((rowSums(Fox_News_Data[,c(62:70)]) == 0))
    Decile1PotentialSubscribers <- as.integer(rowSums(Fox_News_Data[,c(71:79)]) == 0)

    cable.controls <- cbind(Decile1NumberChannels, Fox_News_Data[,c(62:70)], Decile1PotentialSubscribers, Fox_News_Data[,c(71:79)])
    colnames(cable.controls) <- c("Decile 1 no. cable channels available", "Decile 2 no. cable channels available", "Decile 3 no. cable channels available", 
                                  "Decile 4 no. cable channels available", "Decile 5 no. cable channels available", "Decile 6 no. cable channels available", 
                                  "Decile 7 no. cable channels available", "Decile 8 no. cable channels available", "Decile 9 no. cable channels available", 
                                  "Decile 10 no. cable channels available", "Decile 1 no. potential subscribers", "Decile 2 no. potential subscribers", 
                                  "Decile 3 no. potential subscribers", "Decile 4 no. potential subscribers", "Decile 5 no. potential subscribers", 
                                  "Decile 6 no. potential subscribers", "Decile 7 no. potential subscribers", "Decile 8 no. potential subscribers", 
                                  "Decile 9 no. potential subscribers", "Decile 10 no. potential subscribers")

    X <- as.matrix(cbind(demographic.controls.2000, demographic.controls.Diff19962000, cable.controls))

  # Create dummy variables for the different districts
    dummies.district.clusters <- matrix(0, nrow = length(district.clusters), ncol = length(unique(district.clusters)))
    colnames(dummies.district.clusters) <- c(unique(district.clusters))
    for (i in 1:length(district.clusters)) {
      dummies.district.clusters[i,toString(district.clusters[i])] <- 1
    }
  
  # Create a matrix for the remaining variables of interest
    remaining.variables <- Fox_News_Data[,c(14,59, 204, 205)]
    colnames(remaining.variables) <- c("No. cable channels available", "NumberPotentialSubscribers", "Swing district", "Republican district")

    combined.X <- as.matrix(cbind(X, remaining.variables))
    colnames(combined.X) <- c(colnames(X), colnames(remaining.variables))         

  # Clean the data
    indices <- (Fox_News_Data$sample12000 == 1)
    Y <- Y[indices]
    X <- X[indices,]
    W <- W[indices]
    district.clusters <- district.clusters[indices]
    dummies.district.clusters <- dummies.district.dlusters[indices,]
    remaining.variables <- remaining.variables[indices,]
