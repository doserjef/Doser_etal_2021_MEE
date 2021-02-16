# Author: Jeffrey W. Doser
# Code: script to run simulations for Figure 4 in Doser et al. (2021).
#       Goal of these simulations is to determine how adding high
#       quality point count data to acoustic recordings influences
#       precision and accuracy of abundance estimates. 

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2021).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for estimation of bird abundance. In press.
#     Methods in Ecology and Evolution. 

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
n.thin <- 2
n.burn <- 3000
n.adapt <- 5000
n.chain <- 3

# Simulation setup --------------------------------------------------------
n.sims <- 100
seeds <- sample(1:73629, n.sims)
I.vals <- c(5, 10, 20, 30, 50)
out.model <- list()
n.scenarios <- length(I.vals)

for (i in 1:n.sims) {
  print(paste("Current iteration: ", i, " out of ", n.sims, sep = ''))
  set.seed(seeds[i])
  # Set parameters
  J <- 10
  gamma.true <- c(1, 0.5)
  omega.true <- 3
  prop <- 0.1
  n.sites <- 50
  R <- 50
  I <- 50
  beta.true <- c(2, 0.3)
  alpha.true <- c(-2.19, 3, 0.2)
  phi.true <- c(alpha.true[1] + alpha.true[2], alpha.true[3])
  n.count <- 4
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
			         type = 'integrated',
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
  for (j in 1:n.scenarios) {
    print(paste('Current Scenario: ', j, ' out of ', n.scenarios, sep = ''))
    I.curr <- I.vals[j]
    sites.2.curr <- sort(sample(dat$sites.2, I.curr, replace = FALSE))
    c.curr <- dat$c
    c.curr[which(!(dat$sites.2 %in% sites.2.curr)), ] <- NA
    n.count <- apply(c.curr, 1, function(a) {sum(!is.na(a))})
    # Bundle Data -------------------------------------------------------      
    bugs.data <- list(v = dat$v, 
    		      y = dat$y, 
      		      c = c.curr,
      		      n.sites = n.sites, 
      		      J = J, 
      		      n = dat$n, 
      		      k = dat$k, 
      		      n.count = n.count, 
      		      sites.a = dat$sites.1, 
      		      sites.val = sites.val, 
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
                n.iter = n.iter, n.thin = n.thin, 
                n.burn = n.burn, n.chain = n.chain, parallel = TRUE, 
      	  verbose = FALSE))
    out.model[[(i-1) * n.scenarios + j]] <- out
  } # j
  if (i %in% c(20)) {
    date <- Sys.Date()
    file.name <- paste('../../results/spatial-sim-results-', i, '-simulations-', date, '.R', sep = '')
    save(out.model, file = file.name)
  }

} # i

# Save and output results -------------------------------------------------
date <- Sys.Date()
file.name <- paste('../../results/spatial-sim-results-', i, '-simulations-', date, '.R', sep = '')
save(out.model, file = file.name)


