# RESSM-EEG

R code for the analysis from the paper "A Hierarchical Random Effects State-space Model for Modeling Brain Activities from Electroencephalogram Data" by Xingche Guo, Bin Yang, Ji Meng Loh, Qinxia Wang, and Yuanjia Wang

## Overview

### algo_init
The folder *"algo_init"* contains functions dedicated to the initialization step of the MCMC algorithm. These functions aim to generate effective initial values for the temporal dynamical and spatial mapping matrices, ensuring a robust starting point for the MCMC iterations.


### algo_mcmc_sign
The folder *"algo_mcmc_sign"* contains functions used for our proposed MCMC algorithm. 
* **algo_mcmc_sign/mcmc_one_iter.R** The main .R function for the proposed MCMC algorithm.
* **algo_mcmc_sign/mcmc_one_iter.R** The .R function for one iteration of the proposed MCMC algorithm.
* **algo_mcmc_sign/mcmc_A.R** The .R function for the full conditional posteriors of the temporal dynamical matrices in RESSM.
* **algo_mcmc_sign/mcmc_T.R** The .R function for the full conditional posteriors of the spatial mapping matrices in RESSM.
* **algo_mcmc_sign/mcmc_M.R** The .R function for the full conditional posteriors of the latent EEG signals in RESSM.
* **algo_mcmc_sign/mcmc_S.R** The .R function for the full conditional posteriors of the variance components in RESSM.
* **algo_mcmc_sign/rcnorm.R** The .R function for simulating canonical multivariate Gaussian distributions.


