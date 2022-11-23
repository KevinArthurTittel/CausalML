# Robust Smooth Heterogeneous Treatment Effect Estimation using Local Linear Causal Forests

For the MSc Econometrics & Management Science, specializing in Business Analytics & Quantitative Marketing [BAQM], at the Erasmus School of Economics [ESE] in Rotterdam, I have written a theoretical research project titled "Robust Smooth Heterogeneous Treatment Effect Estimation using Local Linear Causal Forests". In particular, I have implemented several novel causal machine learning methods in order to assess the added value in empirical economics research. Concretely, I have re-examined two relatively recent observational studies by implementing Local Linear Causal Forests [LLCF] of Friedberg et al. (2021) and the competing Generalized Random Forests with Local Centering [GRF] method by Athey et al. (2019) on these studies. LLCF is expected to outperform competing methods in case the data has strong smooth underlying signals; since LLCF currently has limited applications in empirical economics, I aim to show that this method can be considered important and effective subject to meaningful and careful choice justifications. To complement the empirical results, I add several simulation results to further demonstrate its added value. 

# Replication files

This folder contains replication files for my theoretical MSc thesis titled "Robust Smooth Heterogeneous Treatment Effect Estimation using Local Linear Causal Forests", as well as a general introduction for the reader. All selected settings for the algorithms here are described in my thesis. Contained in this folder are:

The following files replicate the Tables of the paper:
- Data Preparation Microfinance Project (Field et al., 2013).R: prepares and cleans the publicly available Microfinance Project dataset, with which we revisit an empirical application study of Field et al. (2013).
- Table 1 Empirical Application.R: implements GRF, Cluster-Robust GRF, and LLCF to estimate CATEs of a grace period on three long-run profit, income and capital variables along five characteristics ((1) owning a savings account, (2) being risk loving, (3) a chronically ill family member, (4) wage earner, (5) being impatient) (extending Table 5 of the original paper of Field et al. (2013)).
- Table 3 Empirical application.R: implements LLCF with four different ways of achieving cluster-robustness to estimate CATEs of a grace period on long-run monthly profit (extending Table 5 of the original paper of Field et al. (2013)). 
- Table 4 & 5 Semi-synthetic Application.R: implements GRF, Cluster-Robust GRF, and LLCF to estimate the CATEs of brief paid get-out-the-vote phone calls on voter turnout in the United States in 2002 (elaborating on Section 4.1 of the original paper of Nie and Wager (2021)). In particular, in-sample fit and out-of-sample performance measures are provided to examine the methods' performance. 
- Table 6 Simulation.R: implements GRF and LLCF to estimate CATEs in simulation set-ups 1 and 2 and provide in-sample fit measures, where we created a simulated experiment with difficult smooth noise components under confounding, and strong smooth heterogeneous treatment effects given many noisy covariates in small samples.

The following files replicate the Figures of the paper:
- Figure 1-4 Data Analysis Microfinance Project.R: produces (Fig. 1) histograms for the three numerical covariates from the Microfinance Project data, (Fig. 2) barplots of six unevenly distributed binary variables from the Microfinance Project data, (Fig. 3) hexagonal heatmaps for pairs of numerical covariates to illustrate the spread of the data and provide insights into sparse regions in the covariate space, and (Fig. 4) plots of long-run monthly profits (in Rs) against the numerical covariates to indicate smoothness.
- Figure 5 Heatmaps Simulation.R: heatmaps of treatment effects in simulation set-up 1 with true treatment effect, GRF CATE and LLCF CATE in-sample estimates.
- Figure 6 Plots in-sample measures Simulation.R: plots of four in-sample fit measures for GRF and LLCF in simulation set-up 2.

# Required packages

