rpWishart <- function(df0, df, inv_U0, S){
  df1 <- df0 + df
  inv_U1 <- inv_U0 + S
  U1 <- solve(inv_U1)
  out <- rWishart(1, df1, U1)[,,1]
  return(out)
}



q_S <- function(X, mu, u_ind = NULL, rho = 1e-3){
  if (!is.null(u_ind)){
    X <- X[,-u_ind]
    mu <- mu[-u_ind]
  }
  df <- nrow(X)
  df0 <- ncol(X)
  tXc <- t(X) - mu
  S <- tcrossprod(tXc)
  out <- rpWishart(df0=df0, df=df, inv_U0=diag(df0)*rho, S=S)
  return(out)
}



q_S1 <- function(X, mu, u_ind = NULL, rho, df0){
  if (!is.null(u_ind)){
    X <- X[,-u_ind]
    mu <- mu[-u_ind]
  }
  df_x <- ncol(X)
  df <- nrow(X)
  tXc <- t(X) - mu
  S <- tcrossprod(tXc)
  out <- rpWishart(df0=df0, df=df, inv_U0=diag(df_x)*rho, S=S)
  return(out)
}



q_inv_gamma <- function(a0, b0, size, norm_sq){
  a1 <- a0 + size/2
  b1 <- b0 + norm_sq/2
  return(1/rgamma(1, shape = a1, rate = b1))
}


