library(parallel)


vA2A <- function(vA, m, q){
  qsq <- q^2
  A <- vector("list", m)
  for (i in 1:m){
    vAi <- vA[(1:qsq)+(i-1)*qsq]
    A[[i]] <- matrix(vAi, q, q)
  }
  return(A)
}



q_residuals <- function(Y, M, Theta){
  Yhat <- M %*% t(Theta)
  Res <- Y - Yhat
  return(Res)
}


deviance_YM <- function(norm_y, norm_m, size_y, size_m, sigmasq){
  dev_y <- size_y * log(2*pi*sigmasq) + norm_y / sigmasq
  dev_m <- size_m * log(2*pi) + norm_m
  return(dev_y + dev_m)
}


update_subject_pars <- function(Y_array, M_iJ_0, vA_iJ_0, S_v0, vA_i0, 
                                s, vA, vTheta, S_r, u_idx, trial_idx, 
                                P, Q, m){
  Res <- array(0, dim = dim(Y_array))
  M_iJ_1 <- M_iJ_0
  vA_iJ_1 <- vA_iJ_0
  Qt <- 0
  bt <- 0
  for (j in trial_idx){
    Theta <- matrix(vTheta, P, Q) 
    A <- vA2A(vA_iJ_0[j,], m=m, q=Q)
    M_iJ_1[j,,] <- q_M(Y_array[j,,], M_iJ_0[j,,], Theta, A, s)
    vA_list <- q_vA_ij(M_iJ_1[j,,], vA_i0, S_v0, m)
    vTheta_list <- q_vTheta_ij_ini(M_iJ_1[j,,], Y_array[j,,], s, u_ind=u_idx)
    Qt <- Qt + vTheta_list$Qt
    bt <- bt + vTheta_list$bt
    vA_iJ_1[j,] <- vA_list$x
    Res[j,,] <- q_residuals(Y_array[j,,], M_iJ_1[j,,], Theta)
  }
  vA_i1 <- q_vA(vA, vA_iJ_1[trial_idx,], S_r, S_v0)
  ##
  tX_vi <- t(vA_iJ_1[trial_idx,]) - vA_i1
  S_vi <- tcrossprod(tX_vi) 
  ##
  out <- list(M_iJ_1 = M_iJ_1,
              vA_iJ_1 = vA_iJ_1,
              vA_i1 = vA_i1,
              Res = Res,
              Qt = Qt,
              bt = bt,
              S_vi = S_vi)
  return(out)
}






update_group_pars <- function(Y_list, M_nJ_0, vA_nJ_0, S_v0, vAn0, 
                              sigma2e0, vA0, S_r0, u_idx, trial_idx_list, 
                              NJKP_effective, K, P, Q, m, vA_g, vTheta_g, Sa,
                              cores=1){
  n <- length(Y_list)
  s <- 1/sigma2e0
  out <- mclapply(1:n, FUN = function(i){
    update_subject_pars(Y_array = Y_list[[i]], 
                        M_iJ_0 = M_nJ_0[i,,,], 
                        vA_iJ_0 = vA_nJ_0[i,,], 
                        S_v0 = S_v0, 
                        vA_i0 = vAn0[i,], 
                        s = s, 
                        vA = vA0, 
                        vTheta = vTheta_g, 
                        S_r = S_r0, 
                        u_idx = u_idx, 
                        trial_idx = trial_idx_list[[i]], 
                        P = P, 
                        Q = Q, 
                        m = m)}, mc.cores = cores)
  M_nJ_1 <- M_nJ_0
  vA_nJ_1 <- vA_nJ_0
  vAn1 <- vAn0 

  norm_res_y2 <- 0
  Res <- array(0, dim = c(n, dim(Y_list[[1]])))
  Qt <- 0
  bt <- 0
  S_vi_total <- 0
  for (i in 1:n){
    Qt <- Qt + out[[i]]$Qt
    bt <- bt + out[[i]]$bt
    M_nJ_1[i,,,] <- out[[i]]$M_iJ_1
    vA_nJ_1[i,,] <- out[[i]]$vA_iJ_1
    vAn1[i,] <- out[[i]]$vA_i1
    S_vi_total <- S_vi_total + out[[i]]$S_vi
    Res_i <- out[[i]]$Res
    norm_res_y2 <- norm_res_y2 + sum(Res_i^2)
    Res[i,,,] <- Res_i
  }
  vA1 <- q_vA(vA_g, vAn1, Sa, S_r0)
  df_v <- m*(Q^2)
  S_v1 <- rpWishart(df_v, NJKP_effective/(K*P), diag(df_v) * (1e-3), S_vi_total)
  S_r1 <- q_S(vAn1, vA1, rho = 1e-3)
  sigma2e1 <- q_inv_gamma(a0=1, b0=1e-2, size=NJKP_effective, norm_sq=norm_res_y2)

  result <- list(M_nJ_1 = M_nJ_1, 
                 vA_nJ_1 = vA_nJ_1, 
                 vAn1 = vAn1, 
                 vA1 = vA1,
                 S_v1 = S_v1,
                 S_r1 = S_r1,
                 sigma2e1 = sigma2e1,
                 Res = Res,
                 Qt = Qt,
                 bt = bt)
  return(result)
}


