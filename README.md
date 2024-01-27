# RESSM-EEG

R code for the analysis from the paper "A Hierarchical Random Effects State-space Model for Modeling Brain Activities from Electroencephalogram Data" by Xingche Guo, Bin Yang, Ji Meng Loh, Qinxia Wang, and Yuanjia Wang


### algo_init
The folder *"algo_init"* contains functions dedicated to the initialization step of the MCMC algorithm. These functions aim to generate effective initial values for the temporal dynamical and spatial mapping matrices, ensuring a robust starting point for the MCMC iterations.


### algo_mcmc_sign
The folder *"algo_mcmc_sign"* contains functions used for our proposed MCMC algorithm. 
* **algo_mcmc_sign/mcmc_A.R** The .R function for the full conditional posteriors of the temporal dynamical matrices.
