library(Rtsne)
library(ggplot2)
library(dplyr)
library(latex2exp)

MCMC_results <- readRDS("MCMC_samples.rds")
nn <- c(75, 10, 53, 9)
Q <- 5
m <- 1
cbPalette <- c("#00BFFF", "#FF1493", "#FF8C00", "#228B22", 
               "#F0E442", "#0072B2", "#00FF00", "#79028c",
               "#DC143C", "#00FFFF")


my_resize <- function(X){
  n <- dim(X)[1]
  M <- NULL
  for (i in 1:n){
    M <- rbind(M, X[i,,])
  }
  return(M)
}


##
vA_mean_all <- rbind(my_resize(MCMC_results$vA_nJ_postmean[[1]]),
                     my_resize(MCMC_results$vA_nJ_postmean[[2]]),
                     my_resize(MCMC_results$vA_nJ_postmean[[3]]),
                     my_resize(MCMC_results$vA_nJ_postmean[[4]]))


##
q1 <- 2
q2 <- 3
diag_id1 <- (q1-1)*Q + q1
diag_id2 <- (q2-1)*Q + q2

A23_df <- data.frame(x = vA_mean_all[,diag_id1], y = vA_mean_all[,diag_id2], 
                     group = rep(c("MDD", "Control", "MDD", "Control"), nn*40),
                     Subject = as.factor(rep(1:sum(nn), each=40)))

ggplot(data = A23_df) + 
  geom_point(aes(x = x, y = y, colour = group, shape = group), size=3) + 
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(16,2)) +
  xlab(TeX("$A_{rij}(2,2)$")) +
  ylab(TeX("$A_{rij}(3,3)$")) +
  theme_bw(base_size = 24) +
  theme(
    legend.position =  c(0.15, 0.85),
    legend.title= element_blank(),
    legend.box.background = element_rect(colour = "black"),
    axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())









TX_case_idx <- (20*40+1):(30*40)
TX_ctr_idx <- (75*40+1):(85*40)
CU_case_idx <- (90*40+1):(100*40)
CU_ctr_idx <- (138*40+1):(147*40)

ggplot(data = A23_df[TX_case_idx,]) + 
  stat_ellipse(aes(x = x, y = y, fill = Subject), 
               type = "norm", geom = "polygon", alpha = 0.5) + 
  geom_point(aes(x = x, y = y, colour = Subject), size=3.5) + 
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 1)) +
  xlab(TeX("$A_{rij}(2,2)$")) +
  ylab(TeX("$A_{rij}(3,3)$")) +
  theme_bw(base_size = 24) +
  theme(
    legend.position =  c(0.125, 0.77),
    axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())



ggplot(data = A23_df[TX_ctr_idx,]) + 
  stat_ellipse(aes(x = x, y = y, fill = Subject), 
               type = "norm", geom = "polygon", alpha = 0.5) + 
  geom_point(aes(x = x, y = y, colour = Subject), size=3.5) + 
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 1)) +
  xlab(TeX("$A_{rij}(2,2)$")) +
  ylab(TeX("$A_{rij}(3,3)$")) +
  theme_bw(base_size = 24) +
  theme(
    legend.position =  c(0.125, 0.77),
    axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

