# Author: Jeffrey W. Doser
# Code: script to obtain simulation results from the basic simulations
#       discussed in Doser et al. (2020). 

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2020).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for efficient estimation of bird abundance. 

# For variable defintions, see "eawp-jags.txt" script

rm(list = ls())
library(jagsUI)
library(coda)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(sp)
library(maps)
library(rgdal)
library(GISTools)
library(dplyr)

load("eawp-data.R")

# J.r contains the number of surveys at each acoustic site that contains
# at least 1 detected vocalization. 
J.r <- apply(v, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)
# A.times is a n.sites x J matrix that contains the indices of the
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[i, j] that
# are used in the zero-truncated Poisson.
A.times <- matrix(NA, R, max(J))
tmp <- apply(v, 1, function(a) {which(a != 0)})
for (i in 1:R) {
  if (length(tmp[[i]]) > 0) {
    A.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}

J.val <- apply(v[sites.a, ], 1, function(a){sum(!is.na(a))})
val.times  <- matrix(NA, R.val, max(J))
tmp <- apply(v[sites.a, ], 1, function(a) {which(!is.na(a))})
for (i in 1:R.val) {
  val.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
}

# For day random effect on cue production rate
n.days <- 14
days <- matrix(NA, R, max(J))
days[76, ] <- rep(1:n.days, each = 3)
days[77, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[78, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)
days[79, ] <- rep(1:n.days, each = 3)
days[80, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[81, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)
days[82, ] <- rep(1:n.days, each = 3)
days[83, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[84, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)
days[85, ] <- rep(1:n.days, each = 3)
days[86, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[87, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)

# For Bayesian P-value
sites.c <- which(n.count != 0)
n.count.c <- n.count[which(n.count != 0)]
R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)
# Bundle data ------------------------------------------------------------
bugs.data.A.C <- list(c = c, R = R, n.count = n.count, J = J, X.lambda = X.lambda, 
		      v = v, y = y, sites.a = sites.a, J.r = J.r, A.times = A.times, 
                      sites.c = sites.c, R.C = length(sites.c), n.count.c = n.count.c, 
                      R.A = R.A, J.A = J.A, sites.a.v = sites.a.v, 
                      n.sites = R / 3, days = days, n.days = max(days, na.rm = TRUE) 
		      )
bugs.data.C <- list(c = c, R = R, n.count = n.count, X.lambda = X.lambda,
                    sites.c = sites.c, R.C = length(sites.c), 
                    n.count.c = n.count.c, n.sites = R / 3)
bugs.data.A.V.C <- list(c = c, R = R, n.count = n.count, J = J, 
                        X.lambda = X.lambda, v = v, y = y, k = k, n = n, 
                        sites.a = sites.a, R.val = R.val, J.val = J.val, 
                        J.r = J.r, A.times = A.times, val.times = val.times, 
			sites.c = sites.c, R.C = length(sites.c), 
			n.count.c = n.count.c, R.A = R.A, J.A = J.A, 
			sites.a.v = sites.a.v, n.sites = R / 3, days = days,
			n.days = max(days, na.rm = TRUE)) 
bugs.data.A <- list(R = R, n.count = n.count, J = J, 
                    X.lambda = X.lambda, v = v, y = y, k = k, n = n, 
                    sites.a = sites.a, R.val = R.val, J.val = J.val, 
                    J.r = J.r, A.times = A.times, val.times = val.times, 
		    R.A = R.A, J.A = J.A, 
		    sites.a.v = sites.a.v, n.sites = R / 3, days = days,
		    n.days = max(days, na.rm = TRUE))
# Initial Values ----------------------------------------------------------
N.init <- rep(1, R)
c.max <- apply(c, 1, max, na.rm = TRUE)
N.init <- ifelse(c.max > 1, c.max + 1, 1)
inits <- function() {
  list(
    N = N.init, 
    beta.0 = rnorm(1),
    beta.1 = rnorm(1), 
    omega = runif(1, 0, 10), 
    tau.day = runif(1, 0.1, 1),
    tau = runif(1, 0.1, 1),
    tau.p = runif(1, 0.1, 1),
    tau.day.p = runif(1, 0.1, 1),
    alpha.1 = runif(1, 0, 1), 
    a.phi = runif(1, 0, 5)
  )
}

# Parameters monitored ----------------------------------------------------
params.A.C <- c('beta.1', 'tau', 'beta.0',
		'alpha.0', 'alpha.1', 
		'tau.day',
		'a.phi', 'omega',
	        'p', 'bp.y', 'bp.c', 'bp.v')
params.C <- c('beta.1', 'tau', 'beta.0',
	      'p', 'bp.c')
params.A.V.C <- c('beta.1', 'tau', 'beta.0', 
		  'alpha.0', 'alpha.1', 
		  'tau.day',
		  'a.phi', 'omega', 
	          'p', 'bp.y', 'bp.c', 'bp.v')
params.A <- c('beta.1', 'tau', 'beta.0', 
              'alpha.0', 'alpha.1', 
	      'tau.day', 
	      'a.phi', 'omega', 'bp.y', 'bp.v')

# MCMC settings -----------------------------------------------------------
n.iter <- 200000
n.thin <- 50
n.burn <- 60000
n.chain <- 3
n.adapt <- 5000
# ## Takes a few minutes to complete
# out.model.A.C <- jags(bugs.data.A.C, inits, params.A.C, 'model-AC-eawp-jags.txt', 
#  	              n.iter = n.iter, n.thin = n.thin, n.burn = n.burn,
#  	              n.chain = n.chain, n.adapt = n.adapt, parallel = TRUE)
# out.model.C <- jags(bugs.data.C, inits, params.C, 'model-C-eawp-jags.txt', 
#  		    n.iter = n.iter, n.thin = n.thin, n.burn = n.burn, 
#  		    n.chain = n.chain, n.adapt = n.adapt, parallel = TRUE)
# out.model.A.V.C <- jags(bugs.data.A.V.C, inits, params.A.V.C, 'model-AVC-eawp-jags.txt', 
#  			n.iter = n.iter, n.thin = n.thin, n.burn = n.burn, 
#  			n.chain = n.chain, n.adapt = n.adapt, parallel = TRUE)
# out.model.A <- jags(bugs.data.A, inits, params.A, 'model-AV-eawp-jags.txt', 
#  		    n.iter = n.iter, n.thin = n.thin, n.burn = n.burn, 
#   		    n.chain = n.chain, n.adapt = n.adapt, parallel = TRUE)

#save(out.model.A, out.model.C, out.model.A.C, out.model.A.V.C, file = '../../results/eawp-results.R')
# Read in results so don't have to run every time. 
load('../../results/eawp-results.R')
# Abundance Intercept -----------------
summary(mcmc(out.model.A$sims.list$beta.0))$quantiles
summary(mcmc(out.model.C$sims.list$beta.0))$quantiles
summary(mcmc(out.model.A.C$sims.list$beta.0))$quantiles
summary(mcmc(out.model.A.V.C$sims.list$beta.0))$quantiles

# Year Effect -------------------------
summary(mcmc(out.model.A$sims.list$beta.1))$quantiles
summary(mcmc(out.model.C$sims.list$beta.1))$quantiles
summary(mcmc(out.model.A.C$sims.list$beta.1))$quantiles
summary(mcmc(out.model.A.V.C$sims.list$beta.1))$quantiles

sum(out.model.A$sims.list$beta.1 < 0) / length(out.model.A.C$sims.list$beta.1)
sum(out.model.C$sims.list$beta.1 < 0) / length(out.model.A.C$sims.list$beta.1)
sum(out.model.A.C$sims.list$beta.1 < 0) / length(out.model.A.C$sims.list$beta.1)
sum(out.model.A.V.C$sims.list$beta.1 < 0) / length(out.model.A.C$sims.list$beta.1)

# Bayesian p-values -------------------------------------------------------
mean(out.model.A$sims.list$bp.v)
mean(out.model.C$sims.list$bp.c)
mean(out.model.A.C$sims.list$bp.v)
mean(out.model.A.C$sims.list$bp.c)
mean(out.model.A.V.C$sims.list$bp.v)
mean(out.model.A.V.C$sims.list$bp.c)

# Figure of trend estimates -----------------------------------------------
a.v.trends <- summary(mcmc(out.model.A$sims.list$beta.1))$quantiles[c(1, 3, 5)]
c.trends <- summary(mcmc(out.model.C$sims.list$beta.1))$quantiles[c(1, 3, 5)]
a.c.trends <- summary(mcmc(out.model.A.C$sims.list$beta.1))$quantiles[c(1, 3, 5)]
a.c.v.trends <- summary(mcmc(out.model.A.V.C$sims.list$beta.1))$quantiles[c(1, 3, 5)]
trend.dat <- as.data.frame(rbind(a.c.trends, a.c.v.trends, 
				 a.v.trends, c.trends))
names(trend.dat) <- c('lowest', 'med', 'highest')
trend.dat$model <- c('AC', 'AVC', 'AV', 'C')


my.colors <- c('AV' = 'lightskyblue1',
	       'C' = 'aquamarine4', 'AC' = 'firebrick4',
	      'AVC' = 'darkorchid4')
trend.dat$x.val <- c(3, 4, 1, 2)

trend.plot <- ggplot(data = trend.dat, aes(x = x.val, y = med, col = model)) + 
  geom_point(size = 3.5) + 
  theme_bw(base_size = 24) + 
  scale_color_manual(values = my.colors) + 
  geom_segment(aes(x = x.val, y = lowest, xend = x.val, yend = highest), size = 1, 
	       lineend = 'butt') + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), 
		     label = c('AV', 'C', 'AC', 'AVC', 'RS')) + 
  labs(x = 'Model', y = expression('Trend Estimate ('~beta[1]~')')) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(col = FALSE)

# Intercept Plot ----------------------------------------------------------
a.v.int <- summary(mcmc(out.model.A$sims.list$beta.0))$quantiles[c(1, 3, 5)]
c.int <- summary(mcmc(out.model.C$sims.list$beta.0))$quantiles[c(1, 3, 5)]
a.c.int <- summary(mcmc(out.model.A.C$sims.list$beta.0))$quantiles[c(1, 3, 5)]
a.c.v.int <- summary(mcmc(out.model.A.V.C$sims.list$beta.0))$quantiles[c(1, 3, 5)]
int.dat <- as.data.frame(rbind(a.c.int, a.c.v.int,
				 a.v.int, c.int))
names(int.dat) <- c('lowest', 'med', 'highest')
int.dat$model <- c('AC', 'AVC', 'AV', 'C')


my.colors <- c('AV' = 'lightskyblue1',
	       'C' = 'aquamarine4', 'AC' = 'firebrick4',
	      'AVC' = 'darkorchid4', 'Removal Sampling' = 'black')
int.dat$x.val <- c(3, 4, 1, 2)

int.plot <- ggplot(data = int.dat, aes(x = x.val, y = med, col = model)) +
  geom_point(size = 3.5) +
  theme_bw(base_size = 24) +
  scale_color_manual(values = my.colors) +
  geom_segment(aes(x = x.val, y = lowest, xend = x.val, yend = highest), size = 1,
	       lineend = 'butt') +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
		     label = c('AV', 'C', 'AC', 'AVC')) +
  labs(x = 'Model', y = expression('Intercept ('~beta[0]~')')) +
  #geom_vline(xintercept = 1.5, linetype = 2, col = 'grey') +
  geom_hline(yintercept = 0, linetype = 2) +
  guides(col = FALSE)

# Figure for results ------------------------------------------------------

pdf('../../figures/eawp-summary-plot.pdf', width = 14)
ggarrange(int.plot, trend.plot, ncol = 2, nrow = 1,
	  labels = c('A', 'B'), font.label = list(size = 24))
dev.off()

