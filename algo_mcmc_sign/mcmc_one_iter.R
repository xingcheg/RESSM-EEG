library(parallel)

norm_dist <- function(x){
  return(sqrt(sum(x^2)))
}

cos_dist <- function(x, y){
  return(sum(x * y) / (norm_dist(x) * norm_dist(y)))
}


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

deviance_Y <- function(norm_y, size_y, sigmasq){
  dev_y <- size_y * log(2*pi*sigmasq) + norm_y / sigmasq
  return(dev_y)
}

update_subject_pars <- function(Y_array, M_iJ_0, vA_iJ_0, vTheta_iJ_0, S_v0, S_u0, vA_i0, vTheta_i0,
                                s, vA, vTheta, S_r, S_p, u_idx, trial_idx, 
                                P, Q, m, sign_thres, sign_flag){
  Theta_i0 <- matrix(vTheta_i0, P, Q)
  Res <- array(0, dim = dim(Y_array))
  Res1 <- array(0, dim = dim(Y_array))
  Res_m <- array(0, dim = dim(M_iJ_0))
  Res_m1 <- array(0, dim = dim(M_iJ_0))
  M_iJ_1 <- M_iJ_0
  vA_iJ_1 <- vA_iJ_0
  vTheta_iJ_1 <- vTheta_iJ_0
  for (j in trial_idx){
    v_Theta_ij <- vTheta_iJ_0[j,]
    Theta <- matrix(v_Theta_ij, P, Q) 
    A <- vA2A(vA_iJ_0[j,], m=m, q=Q)
    M_iJ_1[j,,] <- q_M(Y_array[j,,], M_iJ_0[j,,], Theta, A, s)
    vA_list <- q_vA_ij(M_iJ_1[j,,], vA_i0, S_v0, m)
    vTheta_list <- q_vTheta_ij(M_iJ_1[j,,], Y_array[j,,], vTheta_i0, 
                                   S_u0, s, u_ind=u_idx)
    Theta <- matrix(vTheta_list$x, P, Q) 
    Theta_mu <- matrix(vTheta_list$mean, P, Q) 
    ##### check sign #####
    if (sign_flag){
      for (q in 1:Q){
        check_sign <- cos_dist(Theta_mu[,q], Theta_i0[,q])
        if (check_sign < sign_thres){
          Theta_mu[,q] <- -Theta_mu[,q]
          Theta[,q] <- -Theta[,q]
          M_iJ_1[j,,q] <- -M_iJ_1[j,,q]
        }
      }
    }
    ######################
    A <- vA2A(vA_list$x, m=m, q=Q)
    A_mu <- vA2A(vA_list$mean, m=m, q=Q)
    vA_iJ_1[j,] <- vA_list$x
    vTheta_iJ_1[j,] <- c(Theta)
    Res[j,,] <- q_residuals(Y_array[j,,], M_iJ_1[j,,], Theta)
    Res1[j,,] <- q_residuals(Y_array[j,,], M_iJ_1[j,,], Theta_mu)
    Res_m[j,,] <- q_M_res(M_iJ_1[j,,], A)
    Res_m1[j,,] <- q_M_res(M_iJ_1[j,,], A_mu)
  }
  vA_i1 <- q_vA(vA, vA_iJ_1[trial_idx,], S_r, S_v0)
  vTheta_i1 <- q_vTheta(vTheta, vTheta_iJ_1[trial_idx,], S_p, S_u0, u_idx)
  ##
  tX_vi <- t(vA_iJ_1[trial_idx,]) - vA_i1
  S_vi <- tcrossprod(tX_vi) 
  tX_ui <- t(vTheta_iJ_1[trial_idx, -u_idx]) - vTheta_i1[-u_idx]
  S_ui <- tcrossprod(tX_ui) 
  ##
  out <- list(M_iJ_1 = M_iJ_1,
              vA_iJ_1 = vA_iJ_1,
              vTheta_iJ_1 = vTheta_iJ_1,
              vA_i1 = vA_i1,
              vTheta_i1 = vTheta_i1,
              S_vi = S_vi,
              S_ui = S_ui,
              Res = Res,
              Res1 = Res1,
              Res_m = Res_m,
              Res_m1 = Res_m1)
  return(out)
}






update_group_pars <- function(Y_list, M_nJ_0, vA_nJ_0, vTheta_nJ_0, S_v0, S_u0, vAn0, vTheta_n0, 
                              sigma2e0, vA0, vTheta0, S_r0, S_p0, u_idx, trial_idx_list, 
                              NJKP_effective, K, P, Q, m, vA_g, vTheta_g, Sa, St, 
                              cores=1, sign_thres, sign_flag, kappa_v, kappa_u, 
                              kappa_r, kappa_p){
  Theta0 <- matrix(vTheta0, P, Q)
  n <- length(Y_list)
  s <- 1/sigma2e0
  out <- mclapply(1:n, FUN = function(i){
    update_subject_pars(Y_array = Y_list[[i]], 
                        M_iJ_0 = M_nJ_0[i,,,], 
                        vA_iJ_0 = vA_nJ_0[i,,], 
                        vTheta_iJ_0 = vTheta_nJ_0[i,,], 
                        S_v0 = S_v0, 
                        S_u0 = S_u0, 
                        vA_i0 = vAn0[i,], 
                        vTheta_i0 = vTheta_n0[i,],
                        s = s, 
                        vA = vA0, 
                        vTheta = vTheta0, 
                        S_r = S_r0, 
                        S_p = S_p0, 
                        u_idx = u_idx, 
                        trial_idx = trial_idx_list[[i]], 
                        P = P, 
                        Q = Q, 
                        m = m,
                        sign_thres = sign_thres,
                        sign_flag = sign_flag)}, mc.cores = cores)
  M_nJ_1 <- M_nJ_0
  vA_nJ_1 <- vA_nJ_0
  vTheta_nJ_1 <- vTheta_nJ_0
  vAn1 <- vAn0 
  vTheta_n1 <- vTheta_n0
  
  norm_res_y2 <- 0
  norm_res_y2_1 <- 0
  norm_res_m2 <- 0
  norm_res_m2_1 <- 0
  Res <- array(0, dim = c(n, dim(Y_list[[1]])))
  S_vi_total <- 0
  S_ui_total <- 0
  
  for (i in 1:n){
    ##### check sign #####
    Theta_i <- matrix(out[[i]]$vTheta_i1, P, Q)
    if (sign_flag){
      for (q in 1:Q){
        check_sign <- cos_dist(Theta_i[,q], Theta0[,q])
        if (check_sign < sign_thres){
          Theta_i[,q] <- -Theta_i[,q]
          out[[i]]$M_iJ_1[,,q] <- -out[[i]]$M_iJ_1[,,q]
          out[[i]]$vTheta_iJ_1[,(1:P)+(q-1)*P] <- -out[[i]]$vTheta_iJ_1[,(1:P)+(q-1)*P]
        }
      }
    }
    ######################
    M_nJ_1[i,,,] <- out[[i]]$M_iJ_1
    vA_nJ_1[i,,] <- out[[i]]$vA_iJ_1
    vTheta_nJ_1[i,,] <- out[[i]]$vTheta_iJ_1
    S_vi_total <- S_vi_total + out[[i]]$S_vi
    S_ui_total <- S_ui_total + out[[i]]$S_ui
    vAn1[i,] <- out[[i]]$vA_i1
    vTheta_n1[i,] <- c(Theta_i)
    Res_i <- out[[i]]$Res
    norm_res_y2 <- norm_res_y2 + sum(Res_i^2)
    norm_res_y2_1 <- norm_res_y2_1 + sum((out[[i]]$Res1)^2)
    norm_res_m2 <- norm_res_m2 + sum((out[[i]]$Res_m)^2)
    norm_res_m2_1 <- norm_res_m2_1 + sum((out[[i]]$Res_m1)^2)
    Res[i,,,] <- Res_i
  }
  vA1 <- q_vA(vA_g, vAn1, Sa, S_r0)
  vTheta1 <- q_vTheta(vTheta_g, vTheta_n1, St, S_p0, u_idx)
  ##
  df_v <- m*(Q^2)
  df_u <- P*Q - Q*(Q-1)/2
  S_v1 <- rpWishart(df_v, NJKP_effective/(K*P), diag(df_v) * (kappa_v), S_vi_total)
  S_u1 <- rpWishart(df_u, NJKP_effective/(K*P), diag(df_u) * (kappa_u), S_ui_total)
  ##
  S_r1 <- q_S(vAn1, vA1, rho=kappa_r)
  S_p1 <- q_S(vTheta_n1, vTheta1, u_ind=u_idx, rho=kappa_p)
  sigma2e1 <- q_inv_gamma(a0=1, b0=1e-2, size=NJKP_effective, norm_sq=norm_res_y2)
  dev0 <- deviance_Y(norm_res_y2, NJKP_effective, sigma2e1)
  dev <- deviance_YM(norm_res_y2, norm_res_m2, NJKP_effective, NJKP_effective*Q/P, sigma2e1)
  dev1 <- deviance_YM(norm_res_y2_1, norm_res_m2_1, NJKP_effective, NJKP_effective*Q/P, sigma2e1)
  pd <- dev - dev1
  result <- list(M_nJ_1 = M_nJ_1, 
                 vA_nJ_1 = vA_nJ_1, 
                 vTheta_nJ_1 = vTheta_nJ_1, 
                 S_v1 = S_v1, 
                 S_u1 = S_u1, 
                 vAn1 = vAn1, 
                 vTheta_n1 = vTheta_n1,
                 vA1 = vA1,
                 vTheta1 = vTheta1,
                 S_r1 = S_r1,
                 S_p1 = S_p1,
                 sigma2e1 = sigma2e1,
                 dev0 = dev0,
                 dev = dev,
                 pd = pd,
                 Res = Res)
  return(result)
}


