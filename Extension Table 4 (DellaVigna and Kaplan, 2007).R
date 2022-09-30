library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(standardize)
library(foreign)
library(haven)
library(glmnet)

# Import and prepare Fox News data set (DellaVigna and Kaplan, 2007)
  Fox_News_Data <- read_stata("C:/Users/481044kt/AppData/Local/Temp/2/FoxNewsFinalDataQJE_zOeIzx")
  Fox_News_Data <- as.data.frame(Fox_News_Data)

  # Appoint main outcome variable
    Y <- as.vector(Fox_News_Data$reppresfv2p00m96)

  # Appoint treatment assignment variable
    W <- as.vector(foxnews2000)

  # Create a numerical vector of the character group name vector
    district.clusters <- as.vector(Fox_News_Data$diststate)

  # Appoint the control variables matrix
    demographic.controls.2000 <- Fox_News_Data[,c(37:49)]
    colnames(demographic.controls.2000) <- c("Population2000", "Population2000Over18", "HighSchoolFraction2000", "SomeCollegeFraction2000",
                                         "CollegeGraduatesFraction2000", "MalesFraction2000", "BlackFraction2000", "HispanicsFraction2000",
                                         "EmploymentFraction2000", "UnemploymentRate2000", "MarriedFraction2000", "Income2000", "ShareUrban2000")

    demographic.controls.Diff19962000 <- Fox_News_Data[,c(139:150)]
    colnames(demographic.controls.Diff19962000) <- c("PopulationDiff19962000", "SomeCollegeFractionDiff19962000", "SomeCollegeFractionDiff19962000",
                                                 "CollegeGraduatesDiff19962000", "MalesFractionDiff19962000", "BlackFractionDiff19962000",
                                                 "HispanicsFractionDiff19962000", "EmploymentFractionDiff19962000", "UnemploymentRateDiff19962000",
                                                 "MarriedFractionDiff19962000", "IncomeDiff19962000", "ShareUrbanDiff19962000")

    Decile1NumberChannels <- as.integer((rowSums(Fox_News_Data[,c(62:70)]) == 0))
    Decile1PotentialSubscribers <- as.integer(rowSums(Fox_News_Data[,c(71:79)]) == 0)

    cable.controls <- cbind(Decile1NumberChannels, Fox_News_Data[,c(62:70)], Decile1PotentialSubscribers, Fox_News_Data[,c(71:79)])
    colnames(cable.controls) <- c("Decile1NumberChannels, Decile2NumberChannels", "Decile3NumberChannels", "Decile4NumberChannels", 
                              "Decile5NumberChannels", "Decile6NumberChannels", "Decile7NumberChannels", "Decile8NumberChannels", 
                              "Decile9NumberChannels", "Decile10NumberChannels", "Decile1PotentialSubscribers", "Decile2PotentialSubscribers", 
                              "Decile3PotentialSubscribers", "Decile4PotentialSubscribers", "Decile5PotentialSubscribers", 
                              "Decile6PotentialSubscribers", "Decile7PotentialSubscribers", "Decile8PotentialSubscribers", "Decile9PotentialSubscribers", "Decile10PotentialSubscribers")

    X <- as.matrix(cbind(demographic.controls.2000, demographic.controls.Diff19962000, cable.controls))

    remaining.variables <- Fox_News_Data[,c(14,59)]
    colnames(remaining.variables) <- c("NumberCableChannels", "NumberPotentialSubscribers")








                   
                   
  
