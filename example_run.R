source("algo_init/mcmc_main.R")
source("gener_EEG_funcs.R")
get_pars <- readRDS("get_pars_Q2_g1.rds")

###### Note that in this example, there is only one group, analysis on multiple groups is similar.


########## problem size ##########
## subject number
n <- 10
## trial number
J <- 10
## time grid
K <- 100
tt <- 1:K
## time grid (before starting)
K0 <- 30
## number of EEG channels
P <- nrow(get_pars$theta)
## number of latent channels
Q <- get_pars$Q
## lag for AR(m) model
m <- get_pars$m
## error (first assume equal variance)
sigma2e <- 0.16
## Theta 
Theta <- get_pars$theta
cov_Theta <- get_pars$cov_Theta
cov_Theta_i <- get_pars$cov_Theta_i
## A 
vA <- get_pars$vA
cov_A <- get_pars$cov_A
cov_Ai <- get_pars$cov_Ai

## simulate EEG data (one group)
simu_out <- gener_EEG(n, J, K, K0, P, Q, m, sigma2e, vA, 
                      cov_A, cov_Ai, Theta, cov_Theta, cov_Theta_i, 
                      rseed=121121)




############# initial values ##########
## group
R <- 1
Y <- vector("list", R)
Y[[1]] <- simu_out$Y_list
## Y[[2]] <- .... when there is a second group.

trial_id <- vector("list", R)
trial_id[[1]] <- simu_out$trial_id_list
## trial_id[[2]] <- .... when there is a second group.

NJKP_effective <- simu_out$NJKP_effective
## NJKP_effective <- c(simu_out1$NJKP_effective, simu_out2$NJKP_effective) when there are 2 groups.

nn <- n
## nn <- c(n1, n2) when there are 2 groups.



## initialize A (with bias)
vA_init <- c(1.1, 0, 0, 1.1, -0.5, 0, 0, -0.5)

## initialize Theta (all elements start at 0)
vTheta_init <- rep(0, P*Q)

## initialize sigma2e (with bias)
sigma2e_init <- 0.2




############# run MCMC initialization-step to get better initial values ##########
#### change cores (i.e. CPU cores) for parallel computing
#### change M_iters (i.e. maximum iterations) and burn_in to be larger (e.g. 2500, 1250).

set.seed(123123)
out <- MCMC_multi_group(M_iters = 500, burn_in = 250, thin = 5, Y = Y, R = R, nn = nn, 
                        J = J, K = K, P = P, Q = Q, m = m, trial_id = trial_id, 
                        NJKP_effective = NJKP_effective, 
                        vA_init = vA_init, vTheta_init = vTheta_init, sigma2e_init = sigma2e_init, 
                        cores = 10)


## get initial values
vTheta_init <- apply(out$vTheta_g_samps, 2, mean)
vA_init <- apply(out$vA_samps[[1]], 2, mean)
sigma2e_init <- mean(out$sigma2e_samps)



############# MCMC main-step ##########
#### change cores (i.e. CPU cores) for parallel computing
#### change M_iters (i.e. maximum iterations) and burn_in to be larger (e.g. 7500, 2500).
source("algo_mcmc_sign/mcmc_main.R")

out1 <- MCMC_multi_group(M_iters = 1000, burn_in = 500, thin = 5, Y = Y, R = R, nn = nn, 
                         J = J, K = K, P = P, Q = Q, m = m, trial_id = trial_id, 
                         NJKP_effective = NJKP_effective, 
                         vA_init = vA_init, vTheta_init = vTheta_init, sigma2e_init = sigma2e_init, 
                         cores = 10, sign_thres = -0.2, 
                         kappa_v = 1e-3, kappa_u = 1e-3, 
                         kappa_r = 2e-2, kappa_p = 2e-2,
                         kappa_a = 100, kappa_t = 100)


