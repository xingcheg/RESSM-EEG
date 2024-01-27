source("algo_init/mcmc_M.R")
source("algo_init/mcmc_A.R")
source("algo_init/mcmc_T_ini.R")
source("algo_init/mcmc_S_init.R")
source("algo_init/rcnorm.R")
source("algo_init/mcmc_one_iter.R")

##### function to provide valid initial values
get_init <- function(R, nn, J, K, P, Q, m, vA_init, vTheta_init, sigma2e_init){
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  vA_g0 <- vA_init
  vTheta_g0 <- vTheta_init
  Sa0 <- diag(df_v) * 1e-3
  vA0 <- vector("list", R)
  S_r0 <- vector("list", R)
  vAn0 <- vector("list", R)
  vA_nJ_0 <- vector("list", R)
  S_V0 <- vector("list", R)
  M_nJ_0 <- vector("list", R)
  sigma2e0 <- rep(sigma2e_init, R)
  for (r in 1:R){
    n <- nn[r]
    M_nJ_0[[r]] <- array(0, dim = c(n, J, K, Q))
    vA0[[r]] <- vA_g0
    S_r0[[r]] <- diag(df_v)
    vAn0[[r]] <-  matrix(rep(vA_g0, n), nrow = n, byrow = TRUE)
    vA_nJ_0_r <- array(0, dim = c(n, J, m*(Q^2)))
    S_V0[[r]] <- diag(df_v)
    for (i in 1:n){
      vA_nJ_0_r[i,,] <- matrix(rep(vA_g0, J), nrow = J, byrow = TRUE)
    }
    vA_nJ_0[[r]] <- vA_nJ_0_r
  }
  return(list(vA_g0 = vA_g0, vTheta_g0 = vTheta_g0,
              Sa0 = Sa0, vA0 = vA0,
              S_r0 = S_r0, vAn0 = vAn0,
              vA_nJ_0 = vA_nJ_0, S_V0 = S_V0,
              M_nJ_0 = M_nJ_0, sigma2e0 = sigma2e0))
}



##### function to create empty MCMC samples
get_empty_samps <- function(M_use, R, nn, J, K, P, Q, m){
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  vA_g_samps <- matrix(0, M_use, m*(Q^2))
  vTheta_g_samps <- matrix(0, M_use, P*Q)
  Sa_samps <- array(0, dim = c(M_use, df_v, df_v))
  Res_nJ_postmean <- vector("list", R)
  M_nJ_postmean <- vector("list", R)
  vA_nJ_postmean <- vector("list", R)
  vAn_samps <- vector("list", R)
  vA_samps <- vector("list", R)
  S_v_samps <- vector("list", R)
  S_r_samps <- vector("list", R)
  sigma2e_samps <- matrix(0, M_use, R)
  
  for (r in 1:R){
    n <- nn[r] 
    Res_nJ_postmean[[r]] <- array(0, dim = c(n, J, K, P))
    M_nJ_postmean[[r]] <- array(0, dim = c(n, J, K, Q))
    vA_nJ_postmean[[r]] <- array(0, dim = c(n, J, m*(Q^2)))
    vAn_samps[[r]] <- array(0, dim = c(M_use, n, m*(Q^2))) 
    vA_samps[[r]] <- matrix(0, M_use, m*(Q^2))
    S_v_samps[[r]] <- array(0, dim = c(M_use, df_v, df_v))
    S_r_samps[[r]] <- array(0, dim = c(M_use, df_v, df_v))
  }
  return(list(vA_g_samps = vA_g_samps, vTheta_g_samps = vTheta_g_samps,
              Sa_samps = Sa_samps, 
              M_nJ_postmean = M_nJ_postmean, Res_nJ_postmean = Res_nJ_postmean,
              vA_nJ_postmean = vA_nJ_postmean, 
              vAn_samps = vAn_samps, 
              vA_samps = vA_samps, 
              S_v_samps = S_v_samps, 
              S_r_samps = S_r_samps, 
              sigma2e_samps = sigma2e_samps))
}
  
  
##### MCMC algorithm
MCMC_multi_group <- function(M_iters, burn_in, thin, Y, R, nn, J, K, P, Q, m, 
                             trial_id, NJKP_effective, vA_init, vTheta_init, 
                             sigma2e_init, cores = 1){
  M_use <- ceiling( (M_iters - burn_in) / thin )
  init_all <- get_init(R, nn, J, K, P, Q, m, vA_init, vTheta_init, sigma2e_init)
  vA_g0 <- init_all$vA_g0
  vTheta_g0 <- init_all$vTheta_g0
  Sa0 <- init_all$Sa0
  vA0 <- init_all$vA0
  S_r0 <- init_all$S_r0
  vAn0 <- init_all$vAn0
  vA_nJ_0 <- init_all$vA_nJ_0
  S_V0 <- init_all$S_V0
  M_nJ_0 <- init_all$M_nJ_0 
  sigma2e0 <- init_all$sigma2e0
  
  empty_samps <- get_empty_samps(M_use, R, nn, J, K, P, Q, m)
  vA_g_samps <- empty_samps$vA_g_samps
  vTheta_g_samps <- empty_samps$vTheta_g_samps
  Sa_samps <- empty_samps$Sa_samps
  Res_nJ_postmean <- empty_samps$Res_nJ_postmean
  M_nJ_postmean <- empty_samps$M_nJ_postmean
  vA_nJ_postmean <- empty_samps$vA_nJ_postmean
  vAn_samps <- empty_samps$vAn_samps 
  vA_samps <- empty_samps$vA_samps
  S_v_samps <- empty_samps$S_v_samps 
  S_r_samps <- empty_samps$S_r_samps
  sigma2e_samps <- empty_samps$sigma2e_samps
  deviance_samps <- empty_samps$deviance_samps
  pd_samps <- empty_samps$pd_samps
  
  u_idx <- upper_idx(P, Q)
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  
  for (kk in 1:M_iters){
    cat("-----------------iter = ", kk, "--------------\n")
    st1 <- Sys.time()
    Res_kk <- vector("list", R)
    Qt <- 0
    bt <- 0
    for (r in 1:R){
      cat("group = ", r, "\n")
      output_kk <- update_group_pars(Y[[r]], M_nJ_0[[r]], vA_nJ_0[[r]], 
                                     S_V0[[r]], vAn0[[r]], 
                                     sigma2e0[r], vA0[[r]], S_r0[[r]],
                                     u_idx, trial_id[[r]], 
                                     NJKP_effective[r], K, P, Q=Q, m, vA_g0, vTheta_g0, Sa0, 
                                     cores=cores)

      M_nJ_0[[r]] <- output_kk$M_nJ_1
      vA_nJ_0[[r]] <- output_kk$vA_nJ_1
      S_V0[[r]] <- output_kk$S_v1
      vAn0[[r]] <- output_kk$vAn1
      vA0[[r]] <- output_kk$vA1
      S_r0[[r]] <- output_kk$S_r1
      sigma2e0[r] <- output_kk$sigma2e1
      Res_kk[[r]] <- output_kk$Res
      Qt <- Qt + output_kk$Qt
      bt <- bt + output_kk$bt
      
      cat("A = ", "\n")
      print(round(matrix(vA0[[r]], Q, m*Q), 3))
    }
    
    vTheta_g1_x <- rcnorm0(bt, Qt)$x
    vTheta_g1 <- rep(0, P*Q)
    vTheta_g1[-u_idx] <- vTheta_g1_x
    vTheta_g0 <- vTheta_g1

    
    if (R > 1){
      vA0_arr <- do.call("rbind", vA0)
      vA_g1 <- q_vA(vA_init, vA0_arr, 1e-4 * diag(df_v), Sa0)
      Sa1 <- q_S1(vA0_arr, vA_g1, rho = 1, df0 = df_v + 3)
      vA_g0 <- vA_g1
      Sa0 <- Sa1
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
        vAn_samps[[r]][kk1,,] <- vAn0[[r]]
        vA_samps[[r]][kk1,] <- vA0[[r]]
        S_v_samps[[r]][kk1,,] <- S_V0[[r]]
        S_r_samps[[r]][kk1,,] <- S_r0[[r]]
      }
      sigma2e_samps[kk1, ] <- sigma2e0
      vA_g_samps[kk1,] <- vA_g0
      vTheta_g_samps[kk1,] <- vTheta_g0
      Sa_samps[kk1,,] <- Sa0
    }
    
    
  }
  


  MCMC_full_samps <- list(Res_nJ_postmean = Res_nJ_postmean,
                          M_nJ_postmean = M_nJ_postmean,
                          vA_nJ_postmean = vA_nJ_postmean, 
                          vAn_samps = vAn_samps, 
                          vA_samps = vA_samps, 
                          S_v_samps = S_v_samps, 
                          S_r_samps = S_r_samps, 
                          sigma2e_samps = sigma2e_samps,
                          vA_g_samps = vA_g_samps,
                          vTheta_g_samps = vTheta_g_samps,
                          Sa_samps = Sa_samps)
  
  return(MCMC_full_samps)
  
  
}






