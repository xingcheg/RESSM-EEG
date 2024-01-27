upper_idx <- function(p, q){
  X <- matrix(0, p, q)
  idx <- which(c(upper.tri(X))==TRUE)
  return(idx)
}



q_vTheta_ij_ini <- function(M, Y, s, u_ind){
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
  b0 <- b0[-u_ind] * s
  Q0 <-  kronecker(Q0, diag(p))
  Q0 <- Q0[-u_ind,-u_ind] * s
  return(list(Qt = Q0, bt = b0))
}
