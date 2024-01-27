library(reshape2)
library(ggplot2)
library(latex2exp)

##### this R file shows the posterior means and 95% credible intervals 
##### for Ar (i.e. group-level temporal dynamical matrices)


out <- readRDS("MCMC_samples.rds")
Q <- 5


#########################################################################
r_vA_samps <- NULL
for (i in 1:4){
  A_mean <- apply(out$vA_samps[[i]], 2, mean)
  A_sd <- apply(out$vA_samps[[i]], 2, sd)
  A_lb <- A_mean - 1.96 * A_sd
  A_ub <- A_mean + 1.96 * A_sd
  df <- data.frame(Row = rep(1:Q, Q),
                   Column = rep(1:Q, each = Q),
                   mean = A_mean,
                   lower = A_lb,
                   upper = A_ub)
  r_vA_samps <- rbind(r_vA_samps, df)
}

r_vA_samps$Groups <- rep(c("MDD (TX)", "Control (TX)", 
                           "MDD (CU)", "Control (CU)"),
                         each = Q^2)

r_vA_samps$Entries <- paste0("entry=(", r_vA_samps$Row, ",", r_vA_samps$Column, ")")

cbPalette_color <- c("#00BFFF", "#FF1493", "#FF8C00", "#228B22")

ggplot(data = r_vA_samps) + 
  geom_errorbar(aes(x = Groups, ymin=lower, ymax=upper, colour = Groups), 
                width = 0.4, linewidth = 1.2) +
  geom_point(aes(x = Groups, y = mean, colour = Groups), size = 3.5) + 
  facet_wrap(~Entries, ncol = 5, scales = "free") + 
  scale_color_manual(values=cbPalette_color)   +
  xlab("") + 
  ylab(TeX("Values of $A_r$")) + 
  theme_bw(base_size = 14) + 
  theme(legend.position =  "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key.width=unit(2,"cm")) 





#########################################################################
r_vA_diff <- NULL
for (i in 1:2){
  vA_diff <- out$vA_samps[[2+2*(i-1)]] - out$vA_samps[[1+2*(i-1)]]
  A_mean <- apply(vA_diff, 2, mean)
  A_sd <- apply(vA_diff, 2, sd)
  A_lb <- A_mean - 1.96 * A_sd
  A_ub <- A_mean + 1.96 * A_sd
  df <- data.frame(Row = rep(1:Q, Q),
                   Column = rep(1:Q, each = Q),
                   mean = A_mean,
                   lower = A_lb,
                   upper = A_ub)
  r_vA_diff <- rbind(r_vA_diff, df)
}

r_vA_diff$Groups <- rep(c("TX", "CU"), each = Q^2)
r_vA_diff$Entries <- paste0("entry=(", r_vA_diff$Row, ",", r_vA_diff$Column, ")")
r_vA_diff$Significant <- (r_vA_diff$lower > 0) | (r_vA_diff$upper < 0)

r_vA_diff_overall <- NULL
vA_diff <- (out$vA_samps[[2]] + out$vA_samps[[4]] - out$vA_samps[[1]] - out$vA_samps[[3]]) / 2
A_mean <- apply(vA_diff, 2, mean)
A_sd <- apply(vA_diff, 2, sd)
A_lb <- A_mean - 1.96 * A_sd
A_ub <- A_mean + 1.96 * A_sd
r_vA_diff_overall <- data.frame(Row = rep(1:Q, Q), 
                                Column = rep(1:Q, each = Q), 
                                mean = A_mean, 
                                lower = A_lb, 
                                upper = A_ub)

r_vA_diff_overall$Groups <- rep("Overall", each = Q^2)
r_vA_diff_overall$Entries <- paste0("entry=(", r_vA_diff_overall$Row, ",", r_vA_diff_overall$Column, ")")
r_vA_diff_overall$Significant <- (r_vA_diff_overall$lower > 0) | (r_vA_diff_overall$upper < 0)
r_vA_diff1 <- rbind(r_vA_diff, r_vA_diff_overall)
r_vA_diff1$Groups <- factor(r_vA_diff1$Groups,
                            levels = c("Overall", "CU", "TX"),
                            ordered = TRUE)



cbPalette_color1 <- c("black", "red")
ggplot(data = r_vA_diff1) + 
  geom_errorbar(aes(x = Groups, ymin=lower, ymax=upper, colour = Significant), 
                width = 0.4, linewidth = 1.2) +
  geom_point(aes(x = Groups, y = mean, shape = Groups,  colour = Significant), size = 4.6) + 
  geom_hline(yintercept = 0, linewidth = 1.2, colour = "blue", linetype = "longdash") + 
  facet_wrap(~Entries, ncol = 5, scales = "free") + 
  scale_color_manual(values=cbPalette_color1)   +
  scale_shape_manual(values=c(7, 1, 2)) + 
  xlab("") + 
  ylab(TeX("Group Difference of $A_r$")) + 
  theme_bw(base_size = 14) + 
  theme(legend.position =  "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key.width=unit(1,"cm")) 

