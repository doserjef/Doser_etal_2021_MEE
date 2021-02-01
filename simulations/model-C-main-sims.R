# Author: Jeffrey W. Doser
# Code: script to run simulations for Model C in Doser et al. (2020). 
#       Model C is the basic N-mixture model. Abundance is varying
#       over space

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2020).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for efficient estimation of bird abundance. 
#     arXiv preprint arXiv:2011.11047.

# For variable defintions, see "sim-cov-acoustic.R" script

rm(list = ls())
library(AHMbook)
library(jagsUI)
source("sim-acoustic-cov.R")
set.seed(109)

# Subroutines -------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}


# MCMC settings -----------------------------------------------------------
n.iter <- 10000
n.thin <- 5
n.burn <- 3000
n.adapt <- 5000
n.chain <- 3

# Simulation setup --------------------------------------------------------
n.sims <- 100
seeds <- sample(1:73629, n.sims)
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

out.model <- list()


for (i in 1:n.sims) {
  print(paste("Current iteration: ", i, " out of ", n.sims, sep = ''))
  set.seed(seeds[i])
  for (j in 1:n.scenarios) {
    print(paste('Current Scenario: ', j, ' out of ', n.scenarios, sep = ''))
    curr.params <- param.vals[j, ] 
    n.sites <- curr.params$n.site[[1]]
    J <- 10
    gamma.true <- c(1.0, 0.5)
    omega.true <- 3
    prop <- 0.2
    site.type <- curr.params$site.type
    R <- ifelse(site.type %in% c('equal', 'acoustic'), n.sites, round(n.sites * 0.5))
    I <- ifelse(site.type %in% c('equal', 'nmix'), n.sites, round(n.sites * 0.5))
    beta.true <- curr.params$beta[[1]]
    alpha.true <- curr.params$alpha[[1]]
    phi.true <- c(alpha.true[1] + alpha.true[2], alpha.true[3])
    n.count <- curr.params$n.count[[1]]
    dat <- sim.acoustic.count.data(n.sites = n.sites,
			       J = J,
			       n.count = n.count,
			       beta = beta.true,
			       omega = omega.true,
			       prop = prop,
			       gamma = gamma.true,
			       R = R,
			       I = I,
			       alpha = alpha.true,
			       phi = phi.true,
			       val = FALSE,
			       type = 'count',
			       prop.valid = TRUE,
			       print.info = FALSE)

    J.r <- apply(dat$v, 1, function(a) {sum(a != 0)})
    J.r <- ifelse(is.na(J.r), 0, J.r)
    A.times <- matrix(NA, R, J)
    tmp <- apply(dat$v, 1, function(a) {which(a != 0)})
    for (q in 1:R) {
      if (length(tmp[[q]]) > 0) {
        A.times[q, 1:length(tmp[[q]])] <- tmp[[q]]
      }
    }
    J.r <- ifelse(dat$J == 0, 0, J.r)
    n.count <- dat$n.count
    sites.val <- dat$sites.val
    J <- dat$J
    # Bundle Data -------------------------------------------------------      
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
		  X.2 = dat$X.alpha[3, , ])

    # Initial Values ----------------------------------------------------
    N.init <- rep(1, n.sites)
    c.max <- apply(dat$c, 1, max)
    N.init <- ifelse(c.max > 1, c.max + 1, 1)
    N.init[which(is.na(N.init))] <- 1
    inits <- function() {
      list(
        N = N.init,
        K = dat$k,
        beta.0 = beta.true[1] + runif(1, -0.3, 0.3),
        beta.1 = beta.true[2] + runif(1, -0.2, 0.2),
        omega = omega.true + runif(1, -1, 1),
        gamma.0 = gamma.true[1] + runif(1, -0.3, 0.3),
        gamma.1 = gamma.true[2] + runif(1, -0.3, 0.3),
        mu.alpha = logit.inv(alpha.true[1]) + runif(1, 0.1, 0.3),
        alpha.1 = alpha.true[2] + runif(1, -0.2, 0.2),
        alpha.2 = alpha.true[3] + runif(1, -0.2, 0.2),
        mu.phi = logit.inv(phi.true[1]) + runif(1, -0.2, 0.2),
        phi.1 = phi.true[2] + runif(1, -0.2, 0.2)
      )
    }
    # Parameters monitored ----------------------------------------------
    params <- c('beta.0', 'beta.1')
    # Run model ---------------------------------------------------------
    out <- try(jags(bugs.data, inits, params, 'acoustic-cov-abundance-jags.txt', 
                n.iter = n.iter, n.thin = n.thin, n.adapt = n.adapt,
                n.burn = n.burn, n.chain = n.chain, parallel = TRUE, 
      	  verbose = FALSE))
    out.model[[(i-1) * n.scenarios + j]] <- out
  } # j
  if (i %in% c(5, 10, 30)) {
    date <- Sys.Date()
    file.name <- paste('../../results/covariate-model-sim-results-model-2-', i, '-simulations-', date, '.R', sep = '')
    save(out.model, file = file.name)
  }
} # i


date <- Sys.Date()
file.name <- paste('../../results/covariate-model-sim-results-model-2-', i, '-simulations-', date, '.R', sep = '')
save(out.model, file = file.name)


