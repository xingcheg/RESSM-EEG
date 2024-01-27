q_vA_ij <- function(M, vA_i, S_v, m){
  n <- nrow(M)
  q <- ncol(M)
  Q0 <- 0
  b0 <- as.numeric(S_v %*% vA_i)
  for (k in (m+1):n){
    Mk <- M[k,]
    Mk0 <- M[(k-1):(k-m),]
    Mk_star <- c(t(Mk0))
    Q0 <- Q0 + tcrossprod(Mk_star) 
    b0 <- b0 + kronecker(Mk_star, Mk)
  }
  Q0 <- kronecker(Q0, diag(q)) + S_v
  out <- rcnorm0(b0, Q0)
  return(out)
}



q_vA <- function(vA0, vA2, S0, S2){
  df <- nrow(vA2)
  s_vA2 <- apply(vA2, 2, sum)
  Q0 <- S0 + df * S2
  b0 <- S0 %*% vA0 + S2 %*% s_vA2
  out <- rcnorm(b0, Q0)
  return(out)
}

