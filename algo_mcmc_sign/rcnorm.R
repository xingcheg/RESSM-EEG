rcnorm <- function(b, Q){
  n <- length(b)
  L <- t( chol(Q) )
  z <- rnorm(n)
  x0 <- backsolve(t(L), z)
  u <- forwardsolve(L, b)
  mu <- backsolve(t(L), u)
  x <- x0 + mu
  return(as.numeric(x))
}


rcnorm0 <- function(b, Q){
  n <- length(b)
  L <- t( chol(Q) )
  z <- rnorm(n)
  x0 <- backsolve(t(L), z)
  u <- forwardsolve(L, b)
  mu <- backsolve(t(L), u)
  x <- x0 + mu
  return(list(mean = as.numeric(mu), x = as.numeric(x)))
}




rcnorm1 <- function(b, R1, R2){
  n <- length(b)
  z <- rnorm(n)
  x <- R1%*%z + R2%*%b
  return(as.numeric(x))
}
