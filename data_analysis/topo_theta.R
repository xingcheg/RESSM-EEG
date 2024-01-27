library(eegUtils)
library(cowplot)
library(ggplot2)
library(reshape2)

MCMC_results <- readRDS("MCMC_samples.rds")
eeg_coord <- readRDS("eegcoord.rds")
ch_info <- readRDS("ch_names_mcmc.rds")
eeg_coord1 <- eeg_coord[,c("x", "y")] * 10
eeg_coord1$x <- eeg_coord1$x * 1.2
eeg_coord1$electrode <- row.names(eeg_coord1)

vtheta_mean <- apply(MCMC_results$vTheta_samps[[1]], 2, mean)
theta_mean <- matrix(vtheta_mean, nrow=54)
r_theta_mean <- melt(theta_mean)
r_theta_mean$Var1 <- -r_theta_mean$Var1


#### visualize posterior mean for theta
ggplot(data = r_theta_mean) + 
  geom_rect(aes(xmin = Var2+1, ymin = Var1+1,
                xmax = Var2, ymax = Var1, fill = value)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "value") + 
  theme_void(base_size = 16)


#### mapping posterior mean of theta using Topographies (for each column of theta)
theta_col_dat1 <- data.frame(electrode=ch_info$ch_names_use, theta=theta_mean[,1])
theta_topo_dat1 <- merge(x=eeg_coord1, y=theta_col_dat1, by="electrode")
theta_col_dat2 <- data.frame(electrode=ch_info$ch_names_use, theta=theta_mean[,2])
theta_topo_dat2 <- merge(x=eeg_coord1, y=theta_col_dat2, by="electrode")
theta_col_dat3 <- data.frame(electrode=ch_info$ch_names_use, theta=theta_mean[,3])
theta_topo_dat3 <- merge(x=eeg_coord1, y=theta_col_dat3, by="electrode")
theta_col_dat4 <- data.frame(electrode=ch_info$ch_names_use, theta=theta_mean[,4])
theta_topo_dat4 <- merge(x=eeg_coord1, y=theta_col_dat4, by="electrode")
theta_col_dat5 <- data.frame(electrode=ch_info$ch_names_use, theta=theta_mean[,5])
theta_topo_dat5 <- merge(x=eeg_coord1, y=theta_col_dat5, by="electrode")

## M1
def_plot1 <- topoplot(theta_topo_dat1, chan_marker = "name", r=110, 
                      quantity="theta", interp_limit="head", limits = c(-0.22, 0.28)) 
def_plot1$guides$fill$title <- ""
def_plot1$theme$legend.position <- "none"
def_plot1

## M2
def_plot2 <- topoplot(theta_topo_dat2, chan_marker = "name", r=110, 
                      quantity="theta", interp_limit="head", limits = c(-0.22, 0.28)) 
def_plot1$guides$fill$title <- ""
def_plot2$theme$legend.position <- "none"
def_plot2

## M3
def_plot3 <- topoplot(theta_topo_dat3, chan_marker = "name", r=110, 
                      quantity="theta", interp_limit="head", limits = c(-0.22, 0.28)) 
def_plot3$guides$fill$title <- ""
def_plot3$theme$legend.position <- "none"
def_plot3

## M4
def_plot4 <- topoplot(theta_topo_dat4, chan_marker = "name", r=110, 
                      quantity="theta", interp_limit="head", limits = c(-0.22, 0.28)) 
def_plot4$guides$fill$title <- ""
def_plot4$theme$legend.position <- "none"
def_plot4

## M5
def_plot5 <- topoplot(theta_topo_dat5, chan_marker = "name", r=110, 
                      quantity="theta", interp_limit="head", limits = c(-0.22, 0.28)) 
def_plot5$guides$fill$title <- ""
def_plot5$guides$fill$barwidth=unit(1.2, 'cm')
def_plot5$guides$fill$barheight=unit(5, 'cm')
def_plot5$theme$legend.text = element_text(size=20)
def_plot5



