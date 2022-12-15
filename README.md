# Robust Smooth Heterogeneous Treatment Effect Estimation using Local Linear Causal Forests

For the MSc Econometrics & Management Science, specializing in Business Analytics & Quantitative Marketing [BAQM], at the Erasmus School of Economics [ESE] in Rotterdam, I have written a theoretical research project titled "Robust Smooth Heterogeneous Treatment Effect Estimation using Local Linear Causal Forests". In particular, I have implemented several novel causal machine learning methods in order to assess the added value in empirical economics research. Concretely, I have re-examined two relatively recent observational studies by implementing Local Linear Causal Forests [LLCF] of Friedberg et al. (2021) and the competing Generalized Random Forests with Local Centering [GRF] method by Athey et al. (2019) on these studies. LLCF is expected to outperform competing methods in case the data has strong smooth underlying signals; since LLCF currently has limited applications in empirical economics, I aim to show that this method can be considered important and effective subject to meaningful and careful choice justifications. To complement the empirical results, I add several simulation results to further demonstrate its added value. 

# Replication files

This folder contains replication files for my theoretical MSc thesis titled "Robust Smooth Heterogeneous Treatment Effect Estimation using Local Linear Causal Forests", as well as a general introduction for the reader. All selected settings for the algorithms here are described in my thesis. Contained in this folder are:

The following files replicate the Tables of the paper:
- Table 1-2 [Simulation Experiments Section 4.2].R: applies GRF and LLCF to estimate CATEs in simulation set-ups 1 and 2 and provides in-sample fit measures for varying small training sample sizes given a fixed high dimension, where we created a simulated experiment with difficult smooth noise components under confounding, and strong smooth heterogeneous treatment effects given many noisy covariates in small samples.
- Table 3 [Simulation Experiments Section 4.2].R: applies GRF and LLCF to estimate CATEs in simulation set-up 2 and provides in-sample fit measures for varying dimensions giving a fixed training sample size.
- Table 4-5 [Semi-Synthetic Application Section 5.2].R: applies GRF, Cluster-Robust GRF, LLCF and Cluster-Robust LLCF to estimate the CATEs of brief paid get-out-the-vote phone calls on voter turnout in the United States in 2002 (elaborating on Section 4.1 of the original paper of Nie and Wager (2021)). In particular, in-sample fit and out-of-sample performance measures are provided to examine the methods' performance.
- Table 6 [Empirical Application Section 6.3].R: applies GRF, Cluster-Robust GRF, LLCF and Cluster-Robust LLCF to estimate CATEs of a grace period on three long-run monthly profits, income and capital variables along five characteristics ((1) owning a savings account, (2) being risk loving, (3) a chronically ill family member, (4) wage earner, (5) being impatient) (extending Table 5 of the original paper of Field et al. (2013)).
- Table 7-10 [Simulation Supplement Appendix].R: applies GRF and LLCF to estimate CATEs in supplement simulations 1-4 and provides in-sample fit measures for varying small training sample sizes given a fixed high dimension, in which we adjust several components of the simulation settings investigated in Section 4 by varying the relative degree of smoothness and difficulty of several functions.

The following files replicate the Figures of the paper:
- Figure 1 [Simulation Experiments Section 4.1].R: produces heatmaps of treatment propensities of simulation set-ups 1 and 2, and the true treatment effect spread for a set of 10,000 randomly drawn data points according to simulation set-up 1.
- Figure 2 [Simulation Experiments Section 4.2].R: produces a plot of relative differences in four in-sample fit measures between LLCF and GRF in CATE estimation in small training samples on simulation set-ups 1 and 2.
- Figure 3-5 [Semi-Synthetic Application Section 5.1].R: produces (Fig. 3) a histogram of the variable age spread over the data and the true synthetic CATE function spiked on the empirical voting study data set, (Fig. 4) barplots of the synthetic manipulation of voting outcomes given one's actual treatment assignment spiked on the random subsample, and (Fig. 5) a histogram of the synthetic CATEs spiked on the random subsample.
- Figure 6 [Semi-Synthetic Application Section 5.2].R: produces plots of in-sample and out-of-sample performance of GRF, CR.GRF, LLCF and CR.LLCF in terms of RMSE as a function of the proportion of treated observations on the empirical voting study data set.
- Figure 7-11 [Empirical Application Section 6.1].R: produces (Fig. 7) histograms for the three numerical covariates from the Microfinance Project data, (Fig. 8) barplots of six unevenly distributed binary variables from the Microfinance Project data, (Fig. 9) hexagonal heatmaps for pairs of numerical covariates to illustrate the spread of the data and provide insights into sparse regions in the covariate space, (Fig. 10) a heatmap of the p-values corresponding to the correlation matrix of the covariates, and (Fig. 11) plots of long-run monthly profits (in Rs) against the numerical covariates to indicate smoothness.

