##### Assist function
vA2A <- function(vA, m, q){
  qsq <- q^2
  A <- vector("list", m)
  for (i in 1:m){
    vAi <- vA[(1:qsq)+(i-1)*qsq]
    A[[i]] <- matrix(vAi, q, q)
  }
  return(A)
}


##### Check stationary 
check_nonstationary <- function(vAnJ, n, J, m, Q){
  nonstationary_id <- vector("list", n)
  for (i in 1:n){
    id <- NULL
    for (j in 1:J){
      v_aij <- vAnJ[i,j,]
      Aijs <- vA2A(v_aij, m, Q)
      if (m==1){
        big_Aij <- Aijs[[1]]
      }
      if (m > 1){
        B <- do.call(cbind, Aijs)
        C <- diag((m-1)*Q)
        D <- matrix(0, (m-1)*Q, Q)
        E <- cbind(C, D)
        big_Aij <- rbind(B, E)
      }
      eig_val <- eigen(big_Aij)$values
      mod_eig_val <- Mod(eig_val)
      if (any(mod_eig_val > 1-0.01)){
        id <- c(id, j)
      }
    }
    if (length(id)>=1){
      nonstationary_id[[i]] <- id
    }
    else {
      nonstationary_id[[i]] <- NA
    }
  }
  return(nonstationary_id)
}





##### generate EEG data  
gener_EEG <- function(n, J, K, K0, P, Q, m, sigma2e, 
                      vA, cov_A, cov_Ai, 
                      Theta, cov_vTheta, cov_vTheta_i, 
                      rseed){
  
  set.seed(rseed)
  
  ########## A matrices ##########
  df_v <- m*(Q^2)
  ## A_i
  randm <- matrix(rnorm(n*df_v), df_v, n)
  vAn <- t(vA + t(chol(cov_A)) %*% randm)
  ## A_ij
  vAnJ <- array(0, dim = c(n, J, df_v))
  for (i in 1:n){
    randm <- matrix(rnorm(J*df_v), df_v, J)
    vAnJ[i,,] <- t(vAn[i,] + t(chol(cov_Ai)) %*% randm)
  }

  ## non-stationary
  nonstationary_id <- check_nonstationary(vAnJ, n, J, m, Q)
  
  trial_id_list <- lapply(1:n, function(i){
    setdiff(1:J, nonstationary_id[[i]])
  })
  
  NJ_effective <- do.call(sum, lapply(trial_id_list, length))
  NJKP_effective <- NJ_effective * K * P
  
  ########## Theta matrices ##########
  vTheta <- c(Theta)
  ## degenerate distribution
  u_idx <- which(c(upper.tri(Theta))==TRUE)
  R <- length(u_idx)
  df_v <- P*Q-R
  ## Theta_i
  randm <- matrix(0, P*Q, n)
  randm0 <- matrix(rnorm(n*df_v), df_v, n)
  randm[-u_idx,] <- t(chol(cov_Theta)) %*% randm0
  vTheta_n <- t(vTheta + randm)
  ## Theta_ij
  vTheta_nJ <- array(0, dim = c(n, J, P*Q))
  for (i in 1:n){
    randm <- matrix(0, P*Q, J)
    randm0 <- matrix(rnorm(J*df_v), df_v, J)
    randm[-u_idx,] <- t(chol(cov_Theta_i)) %*% randm0
    vTheta_nJ[i,,] <- t(vTheta_n[i,] + randm)
  }
  
  ########## Simulate M ##########
  M_list <- vector("list", n)
  for (i in 1:n){
    M_sub <- array(0, dim = c(J, K, Q))
    for (j in 1:J){
      v_aij <- vAnJ[i,j,]
      Aijs <- vA2A(v_aij, m, Q)
      M_chain <- matrix(0, K+K0+m, Q)
      for (s in 1:m){
        M_chain[s,] <- rnorm(Q)
      }
      for (k in 1:(K+K0)){
        M_chain_new <- rnorm(Q)
        for (s in 1:m){
          M_chain_new <- M_chain_new + Aijs[[s]] %*% M_chain[k+m-s,]
        }
        M_chain[k+m,] <- as.numeric(M_chain_new)
      }
      M_sub[j,,] <- M_chain[m+K0+(1:K),]
    }
    M_list[[i]] <- M_sub
  }
  
  ########## Simulate Y ##########
  Y_list <- vector("list", n)
  for (i in 1:n) {
    Y_sub <- array(0, dim = c(J, K, P))
    for (j in 1:J){
      v_theta_ij <- vTheta_nJ[i,j,]
      theta_ij <- matrix(v_theta_ij, P, Q)
      Y_chain <- matrix(0, K, P)
      for (k in 1:K){
        M_ijk <- M_list[[i]][j,k,]
        Y_chain[k,] <- as.numeric(theta_ij %*% M_ijk) + sqrt(sigma2e) * rnorm(P)
      }
      Y_sub[j,,] <- Y_chain
    }
    Y_list[[i]] <- Y_sub
  }
  
  

  
  return(list(n=n, J=J, K=K, K0=K0, P=P, Q=Q, m=m, sigma2e=sigma2e, 
              vA=vA, cov_A=cov_A, cov_Ai=cov_Ai, 
              Theta=Theta, cov_vTheta=cov_vTheta, cov_vTheta_i=cov_vTheta_i, 
              rseed=rseed, vAn = vAn, vAnJ = vAnJ, 
              nonstationary_id = nonstationary_id, trial_id_list = trial_id_list, 
              NJKP_effective = NJKP_effective, vTheta_n = vTheta_n, 
              vTheta_nJ = vTheta_nJ, M_list = M_list, Y_list = Y_list))
  
}
  
  
