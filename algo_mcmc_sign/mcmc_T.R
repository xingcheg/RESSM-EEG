upper_idx <- function(p, q){
  X <- matrix(0, p, q)
  idx <- which(c(upper.tri(X))==TRUE)
  return(idx)
}

q_vTheta_ij <- function(M, Y, vTheta_i, S_u, s, u_ind){
  K <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(M)
  Q0 <- 0
  b0 <- 0
  for (k in 1:K){
    Yk <- Y[k,]
    Mk <- M[k,]
    Q0 <- Q0 + tcrossprod(Mk)
    b0 <- b0 + kronecker(Mk, Yk)
  }
  b0 <- b0[-u_ind] * s + as.numeric(S_u %*% vTheta_i[-u_ind])
  Q0 <-  kronecker(Q0, diag(p))
  Q0 <- Q0[-u_ind,-u_ind] * s + S_u
  out0 <- rcnorm0(b0, Q0)
  out_mean <- rep(0, p*q)
  out_mean[-u_ind] <- out0$mean
  out_x <- rep(0, p*q)
  out_x[-u_ind] <- out0$x
  return(list(mean = out_mean, x = out_x))
}




q_vTheta <- function(vTheta0, vTheta2, S0, S2, u_ind){
  df <- nrow(vTheta2)
  s_vTheta2 <- apply(vTheta2, 2, sum)
  Q0 <- S0 + df * S2
  b0 <- S0 %*% vTheta0[-u_ind] + S2 %*% s_vTheta2[-u_ind]
  out0 <- rcnorm(b0, Q0)
  out <- rep(0, length(vTheta0))
  out[-u_ind] <- out0
  return(out)
}
