# Author: Jeffrey W. Doser
# Code: script to obtain simulation results from the simulations
#       discussed in the appendix of Doser et al. (2020). 

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2020).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for efficient estimation of bird abundance. 
#     arXiv preprint arXiv:2011.11047.

# For variable defintions, see "sim-cov-acoustic.R" script

rm(list = ls())
library(coda)
library(dplyr)
library(jagsUI)
library(tidyr)
library(ggplot2)

# Read in data ------------------------------------------------------------
# Define parameters used in simulation
n.site.vals <- list(small = 50, large = 100)
site.type.vals <- c('equal', 'acoustic', 'nmix')
n.count.vals <- list(small = 3, large = 5)
beta.vals <- list(small = c(0.2, 0.3), large = c(1, 0.3))
alpha.vals <- list(small.small = c(-2.19, 1.2, 0.2), small.large = c(-2.19, 3, 0.2))
phi.vals <- list(small.small = c(-0.99, 0.2), small.large = c(0.81, 0.2))
n.scenarios <- length(n.site.vals) * length(site.type.vals) * length(n.count.vals) * 
	length(beta.vals) * length(alpha.vals)
param.vals <- expand.grid(n.site.vals, site.type.vals, n.count.vals, beta.vals, 
			  alpha.vals, KEEP.OUT.ATTRS = FALSE)
names(param.vals) <- c('n.site', 'site.type', 'n.count', 'beta', 'alpha')

# Read in data and extract lambda samples 
# File is not on GitHub because of file size. Contact Jeff Doser (doserjef@msu.edu)
# if file is desired. 
#load("../../results/covariate-model-sim-results-model-1-100-simulations-2021-01-06.R")
out.model.1 <- out.model
m.AV.r.hat <- unlist(lapply(out.model.1, FUN = function(a) {a$Rhat$beta.0}))
m.AV.low.beta.0 <- unlist(lapply(out.model.1, FUN = function(a) {a$q2.5$beta.0}))
m.AV.med.beta.0 <- unlist(lapply(out.model.1, FUN = function(a) {a$q50$beta.0}))
m.AV.high.beta.0 <- unlist(lapply(out.model.1, FUN = function(a) {a$q97.5$beta.0}))
m.AV.low.beta.1 <- unlist(lapply(out.model.1, FUN = function(a) {a$q2.5$beta.1}))
m.AV.med.beta.1 <- unlist(lapply(out.model.1, FUN = function(a) {a$q50$beta.1}))
m.AV.high.beta.1 <- unlist(lapply(out.model.1, FUN = function(a) {a$q97.5$beta.1}))
rm(out.model, out.model.1)
# File is not on GitHub because of file size. Contact Jeff Doser (doserjef@msu.edu)
# if file is desired. 
load("../../results/covariate-model-sim-results-model-2-100-simulations-2020-12-29.R")
out.model.2 <- out.model
m.C.r.hat <- unlist(lapply(out.model.2, FUN = function(a) {a$Rhat$beta.0}))
m.C.low.beta.0 <- unlist(lapply(out.model.2, FUN = function(a) {a$q2.5$beta.0}))
m.C.med.beta.0 <- unlist(lapply(out.model.2, FUN = function(a) {a$q50$beta.0}))
m.C.high.beta.0 <- unlist(lapply(out.model.2, FUN = function(a) {a$q97.5$beta.0}))
m.C.low.beta.1 <- unlist(lapply(out.model.2, FUN = function(a) {a$q2.5$beta.1}))
m.C.med.beta.1 <- unlist(lapply(out.model.2, FUN = function(a) {a$q50$beta.1}))
m.C.high.beta.1 <- unlist(lapply(out.model.2, FUN = function(a) {a$q97.5$beta.1}))
rm(out.model, out.model.2)
# File is not on GitHub because of file size. Contact Jeff Doser (doserjef@msu.edu)
# if file is desired. 
load("../../results/covariate-model-sim-results-model-3-100-simulations-2021-01-03.R")
out.model.3 <- out.model
m.AC.r.hat <- unlist(lapply(out.model.3, FUN = function(a) {a$Rhat$beta.0}))
m.AC.low.beta.0 <- unlist(lapply(out.model.3, FUN = function(a) {a$q2.5$beta.0}))
m.AC.med.beta.0 <- unlist(lapply(out.model.3, FUN = function(a) {a$q50$beta.0}))
m.AC.high.beta.0 <- unlist(lapply(out.model.3, FUN = function(a) {a$q97.5$beta.0}))
m.AC.low.beta.1 <- unlist(lapply(out.model.3, FUN = function(a) {a$q2.5$beta.1}))
m.AC.med.beta.1 <- unlist(lapply(out.model.3, FUN = function(a) {a$q50$beta.1}))
m.AC.high.beta.1 <- unlist(lapply(out.model.3, FUN = function(a) {a$q97.5$beta.1}))
rm(out.model, out.model.3)
# File is not on GitHub because of file size. Contact Jeff Doser (doserjef@msu.edu)
# if file is desired. 
load("../../results/covariate-model-sim-results-model-4-100-simulations-2021-01-06.R")
out.model.4 <- out.model
m.AVC.r.hat <- unlist(lapply(out.model.4, FUN = function(a) {a$Rhat$beta.0}))
m.AVC.low.beta.0 <- unlist(lapply(out.model.4, FUN = function(a) {a$q2.5$beta.0}))
m.AVC.med.beta.0 <- unlist(lapply(out.model.4, FUN = function(a) {a$q50$beta.0}))
m.AVC.high.beta.0 <- unlist(lapply(out.model.4, FUN = function(a) {a$q97.5$beta.0}))
m.AVC.low.beta.1 <- unlist(lapply(out.model.4, FUN = function(a) {a$q2.5$beta.1}))
m.AVC.med.beta.1 <- unlist(lapply(out.model.4, FUN = function(a) {a$q50$beta.1}))
m.AVC.high.beta.1 <- unlist(lapply(out.model.4, FUN = function(a) {a$q97.5$beta.1}))
rm(out.model, out.model.4)

n.sims <- length(m.AVC.med.beta.1) / nrow(param.vals)

# Assess convergence
1 - sum(m.AV.r.hat > 1.1) / length(m.AVC.med.beta.1)
1 - sum(m.C.r.hat > 1.1) / length(m.AVC.med.beta.1)
1 - sum(m.AC.r.hat > 1.1) / length(m.AVC.med.beta.1)
1 - sum(m.AVC.r.hat > 1.1) / length(m.AVC.med.beta.1)

# Format summary for plots and output -------------------------------------
dat <- data.frame(cbind(m.AV.low.beta.0, m.AV.med.beta.0, m.AV.high.beta.0, 
			m.C.low.beta.0, m.C.med.beta.0, m.C.high.beta.0,
			m.AC.low.beta.0, m.AC.med.beta.0, m.AC.high.beta.0, 
			m.AVC.low.beta.0, m.AVC.med.beta.0, m.AVC.high.beta.0))
dat$n.sites <- rep(unlist(param.vals$n.site), n.sims)
dat$site.type <- rep(unlist(param.vals$site.type), n.sims)
dat$n.count <- rep(unlist(param.vals$n.count), n.sims)
dat$beta.0 <- rep(rep(rep(c(0.2, 1), each = n.scenarios / 4), 2), n.sims)
dat$beta.1 <- rep(.3, n.sims * n.scenarios)
dat$alpha.0 <- rep(-2.19, n.sims * n.scenarios)
dat$alpha.1 <- rep(rep(c(1.2, 3), each = n.scenarios / 2), n.sims)

dat.grouped <- dat %>%
  group_by(n.sites, site.type, n.count, beta.0, beta.1, alpha.0, alpha.1) %>%
  summarize_at(vars(m.AV.low.beta.0:m.AVC.high.beta.0), median)
write.csv(dat.grouped, "../../results/covariate-simulation-results-beta-0.csv", 
	  row.names = FALSE)

# Create plot to summarize subset of simulation results -------------------
n.count.curr <- 5
n.sites.curr <- 50
plot.dat <- dat.grouped %>% 
  filter(n.count == n.count.curr, n.sites == n.sites.curr)
plot.dat$lambdaPlot <- factor(ifelse(plot.dat$beta.0 == 0.2, 'Low Abundance', 
			      'High Abundance'), levels = c('Low Abundance', 
			                                    'High Abundance'))
plot.dat$alphaPlot <- factor(ifelse(plot.dat$alpha.1 == 1.2, 'Low Detection', 
			     'High Detection'), levels = c('Low Detection', 
			                                   'High Detection'))
plot.dat$index <- rep(1:3, times = 4)
plot.dat$sites <- rep(c('A = C', 'A > C', 'A < C'), times = 4)
my.colors <- c('Model AV' = 'lightskyblue1', 
	       'Model C' = 'aquamarine4', 'Model AC' = 'firebrick4', 
	      'Model AVC' = 'darkorchid4')
my.shapes <- c('Model AV' = 16,
	       'Model C' = 17, 'Model AC' = 18,
	      'Model AVC' = 15)
plot.dat$beta.0 <- ifelse(plot.dat$lambdaPlot == 'Low Abundance', 0.2, 
			  1)
pdf('../../figures/covariateSimulationsResults-beta0.pdf', width = 10)
ggplot(data = plot.dat, aes(x = index, y = m.AVC.med.beta.0)) +
  geom_point(aes(x = index + 0.125, col = 'Model AVC', shape = 'Model AVC'), 
	     size = 3.5) +
  geom_point(aes(x = index - 0.125, y = m.AC.med.beta.0, col = 'Model AC',
		 shape = 'Model AC'), size = 3.5) +
  geom_point(aes(x = index + 0.375, y =  m.C.med.beta.0, col = 'Model C',
		 shape = 'Model C'), size = 3.5) +
  geom_point(aes(x = index - 0.375, y =  m.AV.med.beta.0, col = 'Model AV',
		 shape = 'Model AV'), size = 3.5) +
  geom_segment(aes(x = index + 0.125, y = m.AVC.low.beta.0, xend = index + 0.125,
		   yend = m.AVC.high.beta.0, col = 'Model AVC'), size = 1) +
  geom_segment(aes(x = index - 0.125, y = m.AC.low.beta.0, xend = index - 0.125,
		   yend = m.AC.high.beta.0, col = 'Model AC'), size = 1) +
  geom_segment(aes(x = index + 0.375, y = m.C.low.beta.0, xend = index + 0.375,
		   yend = m.C.high.beta.0, col = 'Model C'), size = 1) +
  geom_segment(aes(x = index - 0.375, y = m.AV.low.beta.0, xend = index - 0.375,
		   yend = m.AV.high.beta.0, col = 'Model AV'), size = 1) +
  scale_color_manual(values = my.colors) +
  scale_shape_manual(values = my.shapes) + 
  facet_grid(vars(lambdaPlot), vars(alphaPlot), scales = 'free') +
  geom_hline(aes(yintercept = beta.0), linetype = 2) +
  geom_vline(xintercept = 1.5, linetype = 2, col = 'grey') +
  geom_vline(xintercept = 2.5, linetype = 2, col = 'grey')  +
  geom_text(aes(x = index, y = ifelse(beta.0 == 0.2, -0.6, 0.25), label = sites)) +
  theme_bw(base_size = 18) +
  labs(y = 'Intercept', color = 'Model', shape = 'Model') +
  theme(axis.text.x = element_blank(),
	axis.ticks.x = element_blank(),
	axis.title.x = element_blank(),
	legend.position = 'bottom')
dev.off()

# Data Summary ------------------------------------------------------------
dat <- dat %>%
  mutate(rel.bias.m.AV = (m.AV.med.beta.0 - beta.0) / beta.0,
	 rel.bias.m.C = (m.C.med.beta.0 - beta.0) / beta.0,
	 rel.bias.m.AC = (m.AC.med.beta.0 - beta.0) / beta.0,
	 rel.bias.m.AVC = (m.AVC.med.beta.0 - beta.0) / beta.0,
         ci.width.m.AV = m.AV.high.beta.0 - m.AV.low.beta.0,
         ci.width.m.C = m.C.high.beta.0 - m.C.low.beta.0,
         ci.width.m.AC = m.AC.high.beta.0 - m.AC.low.beta.0,
         ci.width.m.AVC = m.AVC.high.beta.0 - m.AVC.low.beta.0)

# Get summary values across all models
apply(select(dat, rel.bias.m.AV:ci.width.m.AVC), 2, median)

# Extract beta.1 samples --------------------------------------------------
dat <- data.frame(cbind(m.AV.low.beta.1, m.AV.med.beta.1, m.AV.high.beta.1, 
			m.C.low.beta.1, m.C.med.beta.1, m.C.high.beta.1,
			m.AC.low.beta.1, m.AC.med.beta.1, m.AC.high.beta.1, 
			m.AVC.low.beta.1, m.AVC.med.beta.1, m.AVC.high.beta.1))
dat$n.sites <- rep(unlist(param.vals$n.site), n.sims)
dat$site.type <- rep(unlist(param.vals$site.type), n.sims)
dat$n.count <- rep(unlist(param.vals$n.count), n.sims)
dat$beta.0 <- rep(rep(rep(c(0.2, 1), each = n.scenarios / 4), 2), n.sims)
dat$beta.1 <- rep(.3, n.sims * n.scenarios)
dat$alpha.0 <- rep(-2.19, n.sims * n.scenarios)
dat$alpha.1 <- rep(rep(c(1.2, 3), each = n.scenarios / 2), n.sims)

dat.grouped <- dat %>%
  group_by(n.sites, site.type, n.count, beta.1, beta.0, alpha.0, alpha.1) %>%
  summarize_at(vars(m.AV.low.beta.1:m.AVC.high.beta.1), median)
write.csv(dat.grouped, "../../results/covariate-simulation-results-beta-1.csv", 
	  row.names = FALSE)

n.count.curr <- 5
n.sites.curr <- 50
plot.dat <- dat.grouped %>% 
  filter(n.count == n.count.curr, n.sites == n.sites.curr)
plot.dat$lambdaPlot <- factor(ifelse(plot.dat$beta.0 == 0.2, 'Low Abundance', 
			      'High Abundance'), levels = c('Low Abundance', 
			                                    'High Abundance'))
plot.dat$alphaPlot <- factor(ifelse(plot.dat$alpha.1 == 1.2, 'Low Detection', 
			     'High Detection'), levels = c('Low Detection', 
			                                   'High Detection'))
plot.dat$index <- rep(1:3, times = 4)
plot.dat$sites <- rep(c('A = C', 'A > C', 'A < C'), times = 4)
my.colors <- c('Model AV' = 'lightskyblue1', 
	       'Model C' = 'aquamarine4', 'Model AC' = 'firebrick4', 
	      'Model AVC' = 'darkorchid4')
my.shapes <- c('Model AV' = 16,
	       'Model C' = 17, 'Model AC' = 18,
	      'Model AVC' = 15)
plot.dat$beta.1 <- ifelse(plot.dat$lambdaPlot == 'Low Abundance', 0.3, 
			  0.3)
pdf('../../figures/covariateSimulationsResults-beta1.pdf', width = 10)
ggplot(data = plot.dat, aes(x = index, y = m.AVC.med.beta.1)) +
  geom_point(aes(x = index + 0.125, col = 'Model AVC', shape = 'Model AVC'), 
	     size = 3.5) +
  geom_point(aes(x = index - 0.125, y = m.AC.med.beta.1, col = 'Model AC',
		 shape = 'Model AC'), size = 3.5) +
  geom_point(aes(x = index + 0.375, y =  m.C.med.beta.1, col = 'Model C',
		 shape = 'Model C'), size = 3.5) +
  geom_point(aes(x = index - 0.375, y =  m.AV.med.beta.1, col = 'Model AV',
		 shape = 'Model AV'), size = 3.5) +
  geom_segment(aes(x = index + 0.125, y = m.AVC.low.beta.1, xend = index + 0.125,
		   yend = m.AVC.high.beta.1, col = 'Model AVC'), size = 1) +
  geom_segment(aes(x = index - 0.125, y = m.AC.low.beta.1, xend = index - 0.125,
		   yend = m.AC.high.beta.1, col = 'Model AC'), size = 1) +
  geom_segment(aes(x = index + 0.375, y = m.C.low.beta.1, xend = index + 0.375,
		   yend = m.C.high.beta.1, col = 'Model C'), size = 1) +
  geom_segment(aes(x = index - 0.375, y = m.AV.low.beta.1, xend = index - 0.375,
		   yend = m.AV.high.beta.1, col = 'Model AV'), size = 1) +
  scale_color_manual(values = my.colors) +
  scale_shape_manual(values = my.shapes) + 
  facet_grid(vars(lambdaPlot), vars(alphaPlot), scales = 'free') +
  geom_hline(aes(yintercept = beta.1), linetype = 2) +
  geom_vline(xintercept = 1.5, linetype = 2, col = 'grey') +
  geom_vline(xintercept = 2.5, linetype = 2, col = 'grey')  +
  geom_text(aes(x = index, y = -0.05, label = sites)) +
  theme_bw(base_size = 18) +
  labs(y = 'Covariate Effect', color = 'Model', shape = 'Model') +
  theme(axis.text.x = element_blank(),
	axis.ticks.x = element_blank(),
	axis.title.x = element_blank(),
	legend.position = 'bottom')
dev.off()

# Data Summary ------------------------------------------------------------
dat <- dat %>%
  mutate(rel.bias.m.AV = (m.AV.med.beta.1 - beta.1) / beta.1,
	 rel.bias.m.C = (m.C.med.beta.1 - beta.1) / beta.1,
	 rel.bias.m.AC = (m.AC.med.beta.1 - beta.1) / beta.1,
	 rel.bias.m.AVC = (m.AVC.med.beta.1 - beta.1) / beta.1,
         ci.width.m.AV = m.AV.high.beta.1 - m.AV.low.beta.1,
         ci.width.m.C = m.C.high.beta.1 - m.C.low.beta.1,
         ci.width.m.AC = m.AC.high.beta.1 - m.AC.low.beta.1,
         ci.width.m.AVC = m.AVC.high.beta.1 - m.AVC.low.beta.1)

# Get summary values across all models
apply(select(dat, rel.bias.m.AV:ci.width.m.AVC), 2, median)


