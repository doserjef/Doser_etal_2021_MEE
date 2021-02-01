# Author: Jeffrey W. Doser
# Code: script to analyze simulations for Figure 4 in Doser et al. (2020).
#       Goal of these simulations is to determine how adding high
#       quality point count data to acoustic recordings influences
#       precision and accuracy of abundance estimates. 

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2020).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for efficient estimation of bird abundance. 

# For variable defintions, see "sim-cov-acoustic.R" script

rm(list = ls())
library(coda)
library(dplyr)
library(jagsUI)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Read in data ------------------------------------------------------------
I.vals <- c(5, 10, 20, 30, 50)
beta.true <- c(2, 0.3)
n.sims <- 100
# Read in data and extract lambda samples 
load("../../results/spatial-sim-results-100-simulations-2021-01-03.R")
lowest.beta.0 <- unlist(lapply(out.model, FUN = function(a) {a$q2.5$beta.0}))
low.beta.0 <- unlist(lapply(out.model, FUN = function(a) {quantile(a$sims.list$beta.0, prob = 0.25)}))
med.beta.0 <- unlist(lapply(out.model, FUN = function(a) {a$q50$beta.0}))
highest.beta.0 <- unlist(lapply(out.model, FUN = function(a) {a$q97.5$beta.0}))
high.beta.0 <- unlist(lapply(out.model, FUN = function(a) {quantile(a$sims.list$beta.0, prob = 0.75)}))
low.beta.1 <- unlist(lapply(out.model, FUN = function(a) {quantile(a$sims.list$beta.1, prob = 0.25)}))
lowest.beta.1 <- unlist(lapply(out.model, FUN = function(a) {a$q2.5$beta.1}))
med.beta.1 <- unlist(lapply(out.model, FUN = function(a) {a$q50$beta.1}))
high.beta.1 <- unlist(lapply(out.model, FUN = function(a) {quantile(a$sims.list$beta.1, prob = 0.75)}))
highest.beta.1 <- unlist(lapply(out.model, FUN = function(a) {a$q97.5$beta.1}))
out.rhat.beta.0 <- unlist(lapply(out.model, FUN = function(a) {a$Rhat$beta.0}))
out.rhat.beta.1 <- unlist(lapply(out.model, FUN = function(a) {a$Rhat$beta.1}))

# Format data for plot and output -----------------------------------------
dat <- data.frame(lowest.beta.0, low.beta.0, med.beta.0, high.beta.0, highest.beta.0, 
		  lowest.beta.1, low.beta.1, med.beta.1, high.beta.1, highest.beta.1)
dat$I <- factor(rep(I.vals, n.sims) / max(I.vals))
dat <- dat %>%
  group_by(I) %>%
  summarize_at(vars(lowest.beta.0:highest.beta.1), median)

# Make plots --------------------------------------------------------------
beta.0.plot <- ggplot(data = dat, aes(x = I)) + 
  geom_boxplot(aes(ymin = lowest.beta.0, lower = low.beta.0, middle = med.beta.0, 
		   upper = high.beta.0, ymax = highest.beta.0), stat = 'identity', 
	       fill = 'lightskyblue1', size = 0.75, alpha = 0.75, width = 0.5) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Proportion of Sites with Count Data', 
       y = 'Intercept') + 
  geom_hline(yintercept = beta.true[1], col = 'gray69', size = 0.75)

beta.1.plot <- ggplot(data = dat, aes(x = I)) + 
  geom_boxplot(aes(ymin = lowest.beta.1, lower = low.beta.1, middle = med.beta.1, 
		   upper = high.beta.1, ymax = highest.beta.1), stat = 'identity', 
	       fill = 'lightskyblue1', size = 0.75, alpha = 0.75, width = 0.5) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Proportion of Sites with Count Data', 
       y = 'Covariate Effect') + 
  geom_hline(yintercept = 0.3, col = 'gray69', size = 0.75)

pdf('../../figures/simulation2Results.pdf', height = 8)
ggarrange(beta.0.plot, beta.1.plot, ncol = 1, nrow = 2, labels = c('A', 'B'))
dev.off()

dat$rel.bias.beta.0 <- (dat$med.beta.0 - beta.true[1]) / beta.true[1]
dat$ci.width.beta.0 <- dat$highest.beta.0 - dat$lowest.beta.0
dat$rel.bias.beta.1 <- (dat$med.beta.1 - beta.true[2]) / beta.true[2]
dat$ci.width.beta.1 <- dat$highest.beta.1 - dat$lowest.beta.1

# Complete simulation study summary
write.csv(dat, "../../results/spatial-simulation-results.csv", row.names = FALSE)
