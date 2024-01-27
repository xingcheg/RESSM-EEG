source("algo_mcmc_sign/mcmc_M.R")
source("algo_mcmc_sign/mcmc_A.R")
source("algo_mcmc_sign/mcmc_T.R")
source("algo_mcmc_sign/mcmc_S.R")
source("algo_mcmc_sign/rcnorm.R")
source("algo_mcmc_sign/mcmc_one_iter.R")

##### function to provide valid initial values
get_init <- function(R, nn, J, K, P, Q, m, vA_init, vTheta_init, sigma2e_init){
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  vA_g0 <- vA_init
  vTheta_g0 <- vTheta_init
  Sa0 <- diag(df_v) * 1e-3
  St0 <- diag(df_u) * 1e-3
  vA0 <- vector("list", R)
  vTheta0 <- vector("list", R)
  S_r0 <- vector("list", R)
  S_p0 <- vector("list", R)
  vAn0 <- vector("list", R)
  vTheta_n0 <- vector("list", R)
  vA_nJ_0 <- vector("list", R)
  vTheta_nJ_0 <- vector("list", R)
  S_V0 <- vector("list", R)
  S_U0 <- vector("list", R)
  M_nJ_0 <- vector("list", R)
  sigma2e0 <- rep(sigma2e_init, R)
  for (r in 1:R){
    n <- nn[r]
    M_nJ_0[[r]] <- array(0, dim = c(n, J, K, Q))
    vA0[[r]] <- vA_g0
    vTheta0[[r]] <- vTheta_g0
    S_r0[[r]] <- diag(df_v)
    S_p0[[r]] <- diag(df_u)
    S_V0[[r]] <- diag(df_v)
    S_U0[[r]] <- diag(df_u)
    vAn0[[r]] <-  matrix(rep(vA_g0, n), nrow = n, byrow = TRUE)
    vTheta_n0[[r]] <- matrix(rep(vTheta_g0, n), nrow = n, byrow = TRUE)
    vA_nJ_0_r <- array(0, dim = c(n, J, m*(Q^2)))
    vTheta_nJ_0_r <- array(0, dim = c(n, J, P*Q))
    for (i in 1:n){
      vA_nJ_0_r[i,,] <- matrix(rep(vA_g0, J), nrow = J, byrow = TRUE)
      vTheta_nJ_0_r[i,,] <- matrix(rep(vTheta_g0, J), nrow = J, byrow = TRUE)
    }
    vA_nJ_0[[r]] <- vA_nJ_0_r
    vTheta_nJ_0[[r]] <- vTheta_nJ_0_r
  }
  return(list(vA_g0 = vA_g0, vTheta_g0 = vTheta_g0,
              Sa0 = Sa0, St0 = St0, vA0 = vA0, vTheta0 = vTheta0,
              S_r0 = S_r0, S_p0 = S_p0, vAn0 = vAn0, vTheta_n0 = vTheta_n0,
              vA_nJ_0 = vA_nJ_0, vTheta_nJ_0 = vTheta_nJ_0, S_V0 = S_V0,
              S_U0 = S_U0, M_nJ_0 = M_nJ_0, sigma2e0 = sigma2e0))
}



##### function to create empty MCMC samples
get_empty_samps <- function(M_use, R, nn, J, K, P, Q, m){
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  vA_g_samps <- matrix(0, M_use, m*(Q^2))
  vTheta_g_samps <- matrix(0, M_use, P*Q)
  Sa_samps <- array(0, dim = c(M_use, df_v, df_v))
  St_samps <- array(0, dim = c(M_use, df_u, df_u))
  Res_nJ_postmean <- vector("list", R)
  M_nJ_postmean <- vector("list", R)
  vA_nJ_postmean <- vector("list", R)
  vTheta_nJ_postmean <- vector("list", R)
  vAn_samps <- vector("list", R)
  vTheta_n_samps <- vector("list", R)
  vA_samps <- vector("list", R)
  vTheta_samps <- vector("list", R)
  S_v_samps <- vector("list", R)
  S_u_samps <- vector("list", R)
  S_r_samps <- vector("list", R)
  S_p_samps <- vector("list", R)
  sigma2e_samps <- matrix(0, M_use, R)
  deviance_samps0 <- rep(0, M_use)
  deviance_samps <- rep(0, M_use)
  pd_samps <- rep(0, M_use)
  
  for (r in 1:R){
    n <- nn[r] 
    Res_nJ_postmean[[r]] <- array(0, dim = c(n, J, K, P))
    M_nJ_postmean[[r]] <- array(0, dim = c(n, J, K, Q))
    vA_nJ_postmean[[r]] <- array(0, dim = c(n, J, m*(Q^2)))
    vTheta_nJ_postmean[[r]] <- array(0, dim = c(n, J, P*Q))
    vAn_samps[[r]] <- array(0, dim = c(M_use, n, m*(Q^2))) 
    vTheta_n_samps[[r]] <- array(0, dim = c(M_use, n, P*Q)) 
    vA_samps[[r]] <- matrix(0, M_use, m*(Q^2))
    vTheta_samps[[r]] <- matrix(0, M_use, P*Q)
    S_v_samps[[r]] <- array(0, dim = c(M_use, df_v, df_v))
    S_u_samps[[r]] <- array(0, dim = c(M_use, df_u, df_u))
    S_r_samps[[r]] <- array(0, dim = c(M_use, df_v, df_v))
    S_p_samps[[r]] <- array(0, dim = c(M_use, df_u, df_u))
  }
  return(list(vA_g_samps = vA_g_samps, vTheta_g_samps = vTheta_g_samps,
              Sa_samps = Sa_samps, St_samps = St_samps, 
              M_nJ_postmean = M_nJ_postmean, Res_nJ_postmean = Res_nJ_postmean,
              vA_nJ_postmean = vA_nJ_postmean, vTheta_nJ_postmean = vTheta_nJ_postmean,
              vAn_samps = vAn_samps, vTheta_n_samps = vTheta_n_samps,
              vA_samps = vA_samps, vTheta_samps = vTheta_samps,
              S_v_samps = S_v_samps, S_u_samps = S_u_samps,
              S_r_samps = S_r_samps, S_p_samps = S_p_samps,
              sigma2e_samps = sigma2e_samps, 
              deviance_samps0 = deviance_samps0, deviance_samps = deviance_samps,
              pd_samps = pd_samps))
}
  
  
##### MCMC algorithm
MCMC_multi_group <- function(M_iters, burn_in, thin, Y, R, nn, J, K, P, Q, m, 
                             trial_id, NJKP_effective, vA_init, vTheta_init, 
                             sigma2e_init, cores = 1, sign_thres,
                             kappa_v, kappa_u, kappa_r, kappa_p,
                             kappa_a, kappa_t){
  burn_in_1 <- floor(burn_in / 2)
  M_use <- ceiling( (M_iters - burn_in) / thin )
  init_all <- get_init(R, nn, J, K, P, Q, m, vA_init, vTheta_init, sigma2e_init)
  vA_g0 <- init_all$vA_g0
  vTheta_g0 <- init_all$vTheta_g0
  Sa0 <- init_all$Sa0
  St0 <- init_all$St0
  vA0 <- init_all$vA0
  vTheta0 <- init_all$vTheta0
  S_r0 <- init_all$S_r0
  S_p0 <- init_all$S_p0
  vAn0 <- init_all$vAn0
  vTheta_n0 <- init_all$vTheta_n0
  vA_nJ_0 <- init_all$vA_nJ_0
  vTheta_nJ_0 <- init_all$vTheta_nJ_0
  S_V0 <- init_all$S_V0
  S_U0 <- init_all$S_U0
  M_nJ_0 <- init_all$M_nJ_0 
  sigma2e0 <- init_all$sigma2e0
  
  empty_samps <- get_empty_samps(M_use, R, nn, J, K, P, Q, m)
  vA_g_samps <- empty_samps$vA_g_samps
  vTheta_g_samps <- empty_samps$vTheta_g_samps
  Sa_samps <- empty_samps$Sa_samps
  St_samps <- empty_samps$St_samps
  Res_nJ_postmean <- empty_samps$Res_nJ_postmean
  M_nJ_postmean <- empty_samps$M_nJ_postmean
  vA_nJ_postmean <- empty_samps$vA_nJ_postmean
  vTheta_nJ_postmean <- empty_samps$vTheta_nJ_postmean
  vAn_samps <- empty_samps$vAn_samps 
  vTheta_n_samps <- empty_samps$vTheta_n_samps
  vA_samps <- empty_samps$vA_samps
  vTheta_samps <- empty_samps$vTheta_samps
  S_v_samps <- empty_samps$S_v_samps 
  S_u_samps <- empty_samps$S_u_samps
  S_r_samps <- empty_samps$S_r_samps
  S_p_samps <- empty_samps$S_p_samps
  sigma2e_samps <- empty_samps$sigma2e_samps
  deviance_samps0 <- empty_samps$deviance_samps0
  deviance_samps <- empty_samps$deviance_samps
  pd_samps <- empty_samps$pd_samps
  
  u_idx <- upper_idx(P, Q)
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  
  for (kk in 1:M_iters){
    cat("-----------------iter = ", kk, "--------------\n")
    st1 <- Sys.time()
    dev0 <- 0
    dev <- 0
    pd <- 0
    Res_kk <- vector("list", R)
    sign_flag <- ifelse((kk > burn_in_1) & (kk <= burn_in) & ((kk-burn_in_1)%%thin == 1), TRUE, FALSE)
    for (r in 1:R){
      cat("group = ", r, "\n")
      output_kk <- update_group_pars(Y[[r]], M_nJ_0[[r]], vA_nJ_0[[r]], vTheta_nJ_0[[r]], 
                                     S_V0[[r]], S_U0[[r]], vAn0[[r]], vTheta_n0[[r]], 
                                     sigma2e0[r], vA0[[r]], vTheta0[[r]], S_r0[[r]], S_p0[[r]], 
                                     u_idx, trial_id[[r]], 
                                     NJKP_effective[r], K, P, Q=Q, m, vA_g0, vTheta_g0, Sa0, St0, 
                                     cores=cores, sign_thres = sign_thres, sign_flag = sign_flag,
                                     kappa_v = kappa_v, kappa_u = kappa_u, 
                                     kappa_r = kappa_r, kappa_p = kappa_p)

      M_nJ_0[[r]] <- output_kk$M_nJ_1
      vA_nJ_0[[r]] <- output_kk$vA_nJ_1
      vTheta_nJ_0[[r]] <- output_kk$vTheta_nJ_1
      S_V0[[r]] <- output_kk$S_v1
      S_U0[[r]] <- output_kk$S_u1
      vAn0[[r]] <- output_kk$vAn1
      vTheta_n0[[r]] <- output_kk$vTheta_n1
      vA0[[r]] <- output_kk$vA1
      vTheta0[[r]] <- output_kk$vTheta1
      S_r0[[r]] <- output_kk$S_r1
      S_p0[[r]] <- output_kk$S_p1
      sigma2e0[r] <- output_kk$sigma2e1
      dev0 <- dev0 + output_kk$dev0
      dev <- dev + output_kk$dev
      pd <- pd + output_kk$pd
      Res_kk[[r]] <- output_kk$Res
    }
    
    if (R > 1){
      vA0_arr <- do.call("rbind", vA0)
      vTheta0_arr <- do.call("rbind", vTheta0)
      vA_g1 <- q_vA(vA_init, vA0_arr, diag(df_v) * 1e-4, Sa0)
      vTheta_g1 <- q_vTheta(vTheta_init, vTheta0_arr, diag(df_u) * 1e-4, St0, u_idx)
      Sa1 <- q_S1(vA0_arr, vA_g1, rho=kappa_a, df0=df_v+3)
      St1 <- q_S1(vTheta0_arr, vTheta_g1, u_ind=u_idx, rho=kappa_t, df0=df_u+3)
      vA_g0 <- vA_g1
      vTheta_g0 <- vTheta_g1
      Sa0 <- Sa1
      St0 <- St1
    }
    
    st2 <- Sys.time()
    st12 <- difftime(st2, st1, units = "mins")
    cat("one iteration using:", st12, "minutes, done! \n")

    if ( (kk > burn_in) & ((kk-burn_in)%%thin == 1) ){
      kk1 <- (kk - burn_in - 1)/thin + 1
      for (r in 1:R){
        Res_nJ_postmean[[r]] <- Res_nJ_postmean[[r]] + Res_kk[[r]] * (1/M_use)
        M_nJ_postmean[[r]] <- M_nJ_postmean[[r]] + M_nJ_0[[r]] * (1/M_use)
        vA_nJ_postmean[[r]] <- vA_nJ_postmean[[r]] + vA_nJ_0[[r]] * (1/M_use)
        vTheta_nJ_postmean[[r]] <- vTheta_nJ_postmean[[r]] + vTheta_nJ_0[[r]] * (1/M_use)
        vAn_samps[[r]][kk1,,] <- vAn0[[r]]
        vTheta_n_samps[[r]][kk1,,] <- vTheta_n0[[r]]
        vA_samps[[r]][kk1,] <- vA0[[r]]
        vTheta_samps[[r]][kk1,] <- vTheta0[[r]]
        S_v_samps[[r]][kk1,,] <- S_V0[[r]]
        S_u_samps[[r]][kk1,,] <- S_U0[[r]]
        S_r_samps[[r]][kk1,,] <- S_r0[[r]]
        S_p_samps[[r]][kk1,,] <- S_p0[[r]]
      }
      sigma2e_samps[kk1, ] <- sigma2e0
      vA_g_samps[kk1,] <- vA_g0
      vTheta_g_samps[kk1,] <- vTheta_g0
      Sa_samps[kk1,,] <- Sa0
      St_samps[kk1,,] <- St0
      deviance_samps0[kk1] <- dev0
      deviance_samps[kk1] <- dev
      pd_samps[kk1] <- pd
    }
    
    
  }
  
  
  dev_postmean0 <- mean(deviance_samps0)
  dev_insert0 <- 0
  for (r in 1:R){
    dev_insert0 <- dev_insert0 + deviance_Y(sum(Res_nJ_postmean[[r]]^2), NJKP_effective[r], mean(sigma2e_samps[,r]))
  }
  PD0_a <- dev_postmean0 - dev_insert0
  DIC_orig_a <- dev_insert0 + 2 * PD0_a  
  DIC_orig_a2 <- dev_insert0 + 3 * PD0_a 
  PD0_b <- var(deviance_samps0) * 0.5
  DIC_orig_b <- dev_insert0 + 2 * PD0_b  
  
  dev_postmean <- mean(deviance_samps)
  PD <- mean(pd_samps)
  DIC_modify <- dev_postmean + PD


  MCMC_full_samps <- list(Res_nJ_postmean = Res_nJ_postmean,
                          M_nJ_postmean = M_nJ_postmean,
                          vA_nJ_postmean = vA_nJ_postmean, 
                          vTheta_nJ_postmean = vTheta_nJ_postmean, 
                          vAn_samps = vAn_samps, 
                          vTheta_n_samps = vTheta_n_samps, 
                          vA_samps = vA_samps, 
                          vTheta_samps = vTheta_samps,
                          S_v_samps = S_v_samps, 
                          S_u_samps = S_u_samps, 
                          S_r_samps = S_r_samps, 
                          S_p_samps = S_p_samps, 
                          sigma2e_samps = sigma2e_samps,
                          vA_g_samps = vA_g_samps,
                          vTheta_g_samps = vTheta_g_samps,
                          Sa_samps = Sa_samps,
                          St_samps = St_samps,
                          dev_postmean = dev_postmean,
                          PD = PD,
                          DIC_modify = DIC_modify,
                          dev_postmean0 = dev_postmean0,
                          PD0_a = PD0_a,
                          PD0_b = PD0_b,
                          DIC_orig_a = DIC_orig_a,
                          DIC_orig_a2 = DIC_orig_a2,
                          DIC_orig_b = DIC_orig_b)
  
  return(MCMC_full_samps)
  
  
}






