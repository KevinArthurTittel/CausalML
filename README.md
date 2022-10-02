# Robust Smooth Heterogeneous Treatment Effect Estimation using Causal Machine Learning

For the MSc Econometrics & Management Science, specializing in Business Analytics & Quantitative Marketing [BAQM], at the Erasmus School of Economics [ESE] in Rotterdam, I have written a theoretical research project titled "Robust Smooth Heterogeneous Treatment Effect Estimation using Causal Machine Learning". In particular, I have implemented several novel causal machine learning methods in order to assess the added value in empirical economics research. Concretely, I have re-examined two relatively recent observational studies by implementing Local Linear Causal Forests [LLCF] of Friedberg et al. (2021) and Generalized Random Forests with Local Centering [GRF] by Athey et al. (2019) on these studies. LLCF is expected to outperform competing methods in case the data has strong smooth underlying signals; since LLCF currently has limited applications in empirical economics, I aim to show that this method can be considered important and effective subject to meaningful and careful choice justifications. To complement the empirical results, I add several simulation results to further demonstrate its added value. 

# Replication files

This folder contains replication files for my theoretical MSc thesis titled "Robust Smooth Heterogeneous Treatment Effect Estimation using Causal Machine Learning", as well as a general introduction for the reader. All selected settings for the algorithms here are described in my thesis. Contained in this folder are:

The following replication files concern the first empirical application:
- Data Preparation Microfinance Project (Field et al., 2013).R: prepares and cleans the publicly available Microfinance Project dataset, with which we revisit an empirical application study of Field et al. (2013).
- Extension Table 1 (Field et al., 2013).R: implements GRF, Cluster-Robust GRF, and LLCF to estimate the ATE of a grace period on twelve loan use and business formation outcome variables (extending Tables 1 of the original paper of Field et al. (2013)). 
- Extension Table 2 (Field et al., 2013).R: implements GRF, Cluster-Robust GRF, and LLCF to estimate the ATE of a grace period on three long-run profit, income and capital variables (extending Tables 2 of the original paper of Field et al. (2013)). 
- Extension Table 3 (Field et al., 2013).R: implements GRF, Cluster-Robust GRF, and LLCF to estimate the ATE of a grace period on seven default variables (extending Tables 3 of the original paper of Field et al. (2013)). 
- Extension Table 4 (Field et al., 2013).R: implements GRF, Cluster-Robust GRF, and LLCF to estimate the ATE of a grace period on eight business size and business behavior variables (extending Tables 4 of the original paper of Field et al. (2013)).  
- Extension Table 5 (Field et al., 2013).R: implements GRF, Cluster-Robust GRF, and LLCF to estimate ATE and HTE of a grace period on three long-run profit, income and capital variables along five characteristics ((1) owning a savings account, (2) being risk loving, (3) a chronically ill family member, (4) wage earner, (5) being impatient) (extending Table 5 of the original paper of Field et al. (2013)).

The following replication files concern the second empirical application:
- Data Preparation Fox News (DellaVigna and Kaplan, 2007).R: prepares and cleans the publicly available Fox News dataset, with which we revisit an empirical application study of DellaVigna and Kaplan (2007).

# Required packages

