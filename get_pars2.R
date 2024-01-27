theta <- readRDS("theta.rds")
theta <- theta / theta[1,1]

get_pars <- NULL

sd_theta_scale <- 1
scale <- 0.5
P <- 20
Q <- 2
m <- 2
get_pars$P <- P
get_pars$Q <- Q
get_pars$m <- m

df_v <- m * Q^2
df_u <- (2*P-Q+1)*Q/2
get_pars$theta <- theta[1:P,1:Q] * scale

## A
A1 <- diag(c(0.95, 0.9))
A2 <- diag(c(-0.55, -0.5)) 
vA <- c(c(A1), c(A2))

diag_pos <- 1:Q + Q*(0:(Q-1))

sd_A <- rep(0.03, df_v)
sd_A[diag_pos] <- rep(0.10, Q)
sd_A[Q^2+diag_pos] <- rep(0.05, Q)
cov_A <- diag(sd_A^2)

sd_Ai <- rep(0.03, df_v)
sd_Ai[diag_pos] <- rep(0.08, Q)
sd_Ai[Q^2+diag_pos] <- rep(0.04, Q)
cov_Ai <- diag(sd_Ai^2)


## Theta
sd_Theta <- rep(0.15, df_u) * scale * sd_theta_scale
cov_Theta <- diag(sd_Theta^2)
sd_Theta_i <- rep(0.08, df_u) * scale * sd_theta_scale
cov_Theta_i <- diag(sd_Theta_i^2)

get_pars$vA <- vA
get_pars$cov_A <- cov_A 
get_pars$cov_Ai <- cov_Ai
get_pars$cov_Theta <- cov_Theta 
get_pars$cov_Theta_i <- cov_Theta_i

saveRDS(get_pars,  "get_pars_Q2_g2.rds")
