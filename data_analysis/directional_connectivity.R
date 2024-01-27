library(reshape2)
library(ggplot2)
library(igraph)

out <- readRDS("MCMC_samples.rds")
P <- 54
Q <- 5
eeg_coord <- readRDS("eegcoord.rds")
ch_info <- readRDS("ch_names_mcmc.rds")
eeg_coord1 <- eeg_coord[,c("x", "y")] * 10
eeg_coord1$x <- eeg_coord1$x * 1.2
eeg_coord1$electrode <- row.names(eeg_coord1)

idx_use <- NULL
for (p in 1:54){
  idx <- which(eeg_coord1$electrode==ch_info$ch_names_use[p])
  idx_use <- c(idx_use, idx)
}

eeg_coord2 <- eeg_coord1[idx_use,]


############################################################################
H_array <- array(0, dim = c(P, P, 700))
for (i in 1:700){
  theta <- matrix(out$vTheta_samps[[4]][i,], ncol = Q)
  A <- matrix(out$vA_samps[[4]][i,], ncol = Q)
  H <- theta %*% A %*% solve(crossprod(theta)) %*% t(theta)
  H_array[,,i] <- H
}


H_mean1 <- apply(H_array, c(1,2), mean)


##############################################################################
dd_label <- data.frame(x = 1:54, label = eeg_coord2$electrode)

rH1 <- melt(H_mean1)
ggplot(data = rH1) + 
  geom_raster(aes(x = Var2, y = -Var1, fill = value)) + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limits=c(-0.06, 0.17)) + 
  geom_text(data = dd_label, aes(x = x, y = 0.5, label=label), size = 4, angle = 90) + 
  geom_text(data = dd_label, aes(x = -0.5, y = -x, label=label), size = 4) + 
  theme_void(base_size = 25) + 
  theme(legend.position = "none")

