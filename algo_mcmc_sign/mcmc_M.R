QM <- function(Theta, A, s){
  m <- length(A)
  q <- nrow(A[[1]])
  Q1 <- t(Theta) %*% Theta * s
  W <- lapply(1:m, function(i){crossprod(A[[i]])})
  Q0 <- Reduce('+', W) + diag(q)
  Q <- Q0 + Q1
  return(Q)
}

invertQ <- function(Q){
  eig_Q <- eigen(Q)
  P <- eig_Q$vectors
  d <- eig_Q$values
  return(list(R1 = P%*%diag(1/sqrt(d))%*%t(P),
              R2 = P%*%diag(1/d)%*%t(P)))
}



q_Mk <- function(b1, Mk, A, R1, R2){
  m <- length(A)
  q <- nrow(A[[1]])
  A1 <- vector("list", m+1)
  A1[[1]] <- - diag(q)
  A1[-1] <- A
  b0 <- 0
  for (i in 0:m){
    d0 <- 0
    for (j in 0:m){
      if (j!=i){
        d0 <- d0 + A1[[j+1]] %*% Mk[m+1+i-j,]
      }
    }
    d0 <- t(A1[[i+1]]) %*% d0
    b0 <- b0 + d0
  }
  b <- b1 - b0
  out <- rcnorm1(b, R1, R2)
  return(out)
}



q_M <- function(Y, M, Theta, A, s){
  B1 <- t(Theta) %*% t(Y) * s
  Q <- QM(Theta, A, s)
  Q_inv <- invertQ(Q)
  R1 <- Q_inv$R1
  R2 <- Q_inv$R2
  K <- nrow(Y)
  q <- ncol(M)
  m <- length(A)
  M1 <- matrix(0, K+2*m, q)
  M1[m+(1:K),] <- M
  for (k in 1:K){
    Mk <- M1[k:(k+2*m),]
    b1 <- B1[,k]
    M1[k+m,] <- q_Mk(b1, Mk, A, R1, R2)
  }
  return(M1[m+(1:K),])
}



q_M_res <- function(M, A){
  K <- nrow(M)
  q <- ncol(M)
  m <- length(A)
  M_res <- matrix(0, K, q)
  for (k in (m+1):K){
    r <- M[k,]
    for (j in 1:m){
      r <- r - as.numeric(A[[j]] %*% M[k-j,])
    }
    M_res[k,] <- r
  }
  return(M_res)
}

