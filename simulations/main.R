# Author: Jeffrey W. Doser
# Code: script to run integrated acoustic and point count data models from 
#       Doser et al. (2020)

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2021).
#     Integrating acoustic recordings and point count surveys for 
#     efficient estimation of bird abundance. 
#     arXiv preprint arXiv:2011.11047.

rm(list = ls())
library(AHMbook)
library(jagsUI)
source("sim-acoustic-cov.R")
# set.seed(109)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulate data -----------------------------------------------------------

n.sites <- 50
R <- 50
I <- 25
n.count <- 4
J <- 10
beta <- c(1, -0.2)
gamma <- c(1.0, 0.5)
omega <- 3
prop <- 0.2
alpha <- c(-2, 3, -0.2)
phi <- c(0, 0.2)

# Use the following parameters for sim.acoustic.count.data() for the 
# following models:

# Model AV: type = 'acoustic', val = TRUE
# Model C: type = 'count', val = FALSE
# Model AC: type = 'integrated', val = FALSE
# Model AVC: type = 'integrated', val = TRUE

# variable definitions
# n.sites = total number of sites
# j = number of acoustic repeat visits
# n.count = number of count repeat visits (l in manuscript)
# beta = covariates on abundance
# omega = mean # of false positive acoustic detections
# prop = proportion (or number) of automated acoustic vocalizations 
#        that are validated
# gamma = covariate effects on true positive vocalization detection rate
#         in acoustic data
# r = number of sites with acoustic data
# i = number of sites with point count data
# alpha = covariates on probability of detecting at least one vocalization
#         in acoustic data
# val = logical value indicating whether or not any validation of acoustic
#       data will occur
# type = type of model to fit. set to 'integrated', 'acoustic', or 'count'
# prop.valid = logical value indicating whether "prop" is a proportion
#              (true) or a number of vocalizations (false)
# print.info = logical value indicating whether or not summar info on the 
#              simulated data set is printed. 

dat <- sim.acoustic.count.data(n.sites = n.sites, 
			       J = J, 
			       n.count = n.count, 
			       beta = beta,
			       omega = omega,
			       prop = prop, 
			       gamma = gamma, 
			       R = R, 
			       I = I, 
			       alpha = alpha,
			       phi = phi,
			       val = TRUE, 
			       type = 'integrated', 
			       prop.valid = TRUE,
			       print.info = TRUE)
# J.r contains the number of surveys at each acoustic site that contains
# at least 1 detected vocalization. 
J.r <- apply(dat$v, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)
# A.times is a n.sites x J matrix that contains the indices of the 
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[i, j] that 
# are used in the zero-truncated Poisson. 
A.times <- matrix(NA, n.sites, J)
tmp <- apply(dat$v, 1, function(a) {which(a != 0)})
for (q in 1:n.sites) {
  if (length(tmp[[q]]) > 0) {
    A.times[q, 1:length(tmp[[q]])] <- tmp[[q]]
  }
}
J.r <- ifelse(dat$J == 0, 0, J.r)
# Bundle data -------------------------------------------------------------
bugs.data <- list(v = dat$v, 
		  y = dat$y, 
		  c = dat$c,
		  n.sites = n.sites, 
		  J = dat$J, 
		  n = dat$n, 
		  k = dat$k, 
		  n.count = dat$n.count, 
		  sites.a = dat$sites.1, 
		  sites.val = dat$sites.val, 
		  X.1 = dat$X.lambda, 
		  A.times = A.times, 
		  J.r = J.r,
		  X.4 = dat$X.phi[2, , ], 
		  X.3 = dat$X.delta[2, , ], 
		  X.2 = dat$X.alpha[3, , ]
)

# Initial values ----------------------------------------------------------
N.init <- rep(1, n.sites)
c.max <- apply(dat$c, 1, max)
N.init <- ifelse(c.max > 1, c.max + 1, 1)
N.init[which(is.na(N.init))] <- 1
inits <- function() {
  list(
    N = N.init, 
    K = dat$k,
    beta.0 = beta[1] + runif(1, -0.3, 0.3),
    beta.1 = beta[2] + runif(1, -0.2, 0.2),
    omega = omega + runif(1, -1, 1), 
    gamma.0 = gamma[1] + runif(1, -0.3, 0.3), 
    gamma.1 = gamma[2] + runif(1, -0.3, 0.3),
    mu.alpha = logit.inv(alpha[1]) + runif(1, 0.1, 0.3), 
    alpha.1 = alpha[2] + runif(1, -0.2, 0.2),
    alpha.2 = alpha[3] + runif(1, -0.2, 0.2), 
    mu.phi = logit.inv(phi[1]) + runif(1, -0.2, 0.2),
    phi.1 = phi[2] + runif(1, -0.2, 0.2)
  )
}

# Parameters monitored ----------------------------------------------------
params <- c('beta.0', 'beta.1', 'alpha.0', 'alpha.1', 'alpha.2', 'phi.0', 'phi.1', 
	   'gamma.0', 'gamma.1', 'omega')

# MCMC settings -----------------------------------------------------------
n.iter <- 5000
n.thin <- 2
n.adapt <- 3000
n.burn <- 2000
n.chain <- 3

out <- jags(bugs.data, inits, params, 'acoustic-cov-abundance-jags.txt', 
	    n.iter = n.iter, n.thin = n.thin, n.adapt = n.adapt,
	    n.burn = n.burn, n.chain = n.chain, parallel = TRUE)

out
#plot(out)
