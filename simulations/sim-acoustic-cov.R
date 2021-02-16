# Author: Jeffrey W. Doser
# Code: function to simulate acoustic and point count data for use in the 
#       simulation study of Doser et al. (2021)

# This code extends code from the following sources: 
#     1. Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., 
#        and Nichols, J. D. (2018b). A new framework for analysing 
#        automated acoustic species detection data: Occupancy estimation 
#        and optimization of recordings post-processing. 
#        Methods in Ecology and Evolution, 9(3):560–570.
#     2. Kery, M. and Royle, J. A. (2020). Applied hierarchical modeling 
#        in ecology: Analysis of distribution, abundance, and species 
#        richness in R and BUGS: Volume 2: Dynamic and advanced models. 
#        Academic Press.

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2021).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for estimation of bird abundance. In press.
#     Methods in Ecology and Evolution. 

# Variable definitions

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

# sim-raw.R ---------------------------------------------------------------

sim.acoustic.count.data <- function(n.sites, J, n.count, beta, omega, prop, gamma,
				    R, I, alpha, phi, val = TRUE, 
				    type = 'integrated', prop.valid = TRUE, 
				    print.info = TRUE) {

  require(AHMbook)
  
  # Subroutines -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Get sites for each method ---------------------------------------------
  if (R >= I) {
    sites <- 1:n.sites
    sites.1 <- sort(sample(n.sites, R, replace = FALSE))
    sites.tmp <- sites[!(sites %in% sites.1)]
    tmp.R <- I - length(sites.tmp)
    sites.2 <- sort(c(sites.tmp, sample(sites.1, tmp.R, replace = FALSE)))
  } else {
    sites <- 1:n.sites
    sites.2 <- sort(sample(n.sites, I, replace = FALSE))
    sites.tmp <- sites[!(sites %in% sites.2)]
    tmp.R <- R - length(sites.tmp)
    sites.1 <- sort(c(sites.tmp, sample(sites.2, tmp.R, replace = FALSE)))
  }


  # Form covariates -------------------------------------------------------
  
  # Space-varying abundance -----------
  n.beta <- length(beta)
  X.lambda <- matrix(1, nrow = n.sites, ncol = n.beta)
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X.lambda[, i] <- rnorm(n.sites)
    }
  }
  lambda <- exp(X.lambda %*% as.matrix(beta))

  # Space-time-varying cue rate -------
  n.gamma <- length(gamma)
  X.delta <- array(NA, dim = c(n.gamma, n.sites, J))
  X.delta[1,sites.1, ] <- 1
  for (k in 2:n.gamma) {
    for (j in 1:J) {
      X.delta[k, sites.1, j] <- rnorm(R)  
    } # j
  } # k
  delta <- matrix(NA, nrow = n.sites, ncol = J)
  for (j in 1:J) {
    delta[sites.1, j] <- exp(gamma %*% X.delta[, sites.1, j])
  } # j
  
  # Define abundance of all sites -----------------------------------------
  N <- rpois(n.sites, lambda)

  # Define detection process ----------------------------------------------
  # Acoustic data ---------------------
  n.alpha <- length(alpha)
  X.alpha <- array(NA, dim = c(n.alpha, n.sites, J))
  X.alpha[1, sites.1, ] <- 1
  X.alpha[2, sites.1, ] <- N[sites.1]
  for (k in 3:n.alpha) {
    for (j in 1:J) {
      X.alpha[k, sites.1, j] <- rnorm(R)
    } # j
  } # k
  p.a <- matrix(NA, nrow = n.sites, ncol = J)
  for (j in 1:J) {
    p.a[sites.1, j] <- logit.inv(alpha %*% X.alpha[, sites.1, j])
  }
  # Count data ------------------------
  n.phi <- length(phi)
  X.phi <- array(NA, dim = c(n.phi, n.sites, n.count))
  X.phi[1, sites.2, ] <- 1
  p <- matrix(NA, nrow = n.sites, ncol = n.count) 
  for (k in 2:n.phi) {
    for (j in 1:n.count) {
      X.phi[k, sites.2, j] <- rnorm(I)
    } # j
  } # k
  for (j in 1:n.count) {
    p[sites.2, j] <- logit.inv(phi %*% X.phi[, sites.2, j])
  }

  # Define Data -----------------------------------------------------------
  # Binary acoustic data
  y <- matrix(NA, n.sites, J)
  # Acoustic automated vocalization data
  v <- matrix(NA, n.sites, J)
  # True number of vocalizations in acoustic data
  K <- matrix(NA, n.sites, J)
  # False number vocalizations in acoustic data
  Q <- matrix(NA, n.sites, J)
  # True number of verified vocalizations in acoustic data
  k <- matrix(NA, nrow = R, ncol = J)
  # Total number of verified vocalizations in acoustic data
  n <- matrix(NA, nrow = R, ncol = J)
  
  for (i in 1:R) {
    for (j in 1:J) {
      y[sites.1[i], j] <- rbinom(1, 1, p.a[sites.1[i], j])
      K[sites.1[i], j] <- rpois(1, delta[sites.1[i], j] * N[sites.1[i]] * y[sites.1[i], j])
      Q[sites.1[i], j] <- rpois(1, omega * y[sites.1[i], j])
      v[sites.1[i], j] <- K[sites.1[i], j] + Q[sites.1[i], j]
    }
  }

  # Acoustic validation data ----------------------------------------------
  for (j in 1:J) {
    tmp <- valid_data(v[sites.1, j], tp = K[sites.1, j], 
        	      n.valid = prop, prop.valid = prop.valid)
    n[, j] <- tmp$n
    k[, j] <- tmp$k
  }
  # Define point count data -----------------------------------------------
  c <- matrix(NA, n.sites, n.count)
  for (j in 1:n.count) {
    c[sites.2, j] <- rbinom(I, N[sites.2], p[sites.2, j])
  }

  # Re-assign values based on what model will be fit
  sites.val <- ifelse(val, R, 0)
  if (type == 'acoustic') {
    n.count <- rep(0, n.sites)  
  } else {
    n.count <-  apply(c, 1, function(a) {sum(!is.na(a))})
  }
  if (type == 'count') {
    J <- rep(0, n.sites)
  } else {
    J <- apply(v, 1, function(a) {sum(!is.na(a))})
  }
	      
  if (print.info) {
    print(paste("Mean abundance: ", round(mean(lambda), 3), sep = ''))
    print(paste('Mean cue rate: ', round(mean(delta, na.rm = TRUE), 3), sep = ''))
    print(paste('Mean acoustic detection: ', round(mean(p.a, na.rm = TRUE), 3), sep = ''))
    print(paste('Mean count detection: ', round(mean(p, na.rm = TRUE), 3), sep = ''))
  }

  return(
    list(
      v = v, K = K, N = N, Q = Q, k = k, n = n, sites.1 = sites.1, 
      sites.2 = sites.2, c = c, sites = sites, y = y, 
      p.a = p.a, sites.val = sites.val, n.count = n.count, J = J, lambda = lambda, 
      X.lambda = X.lambda, X.alpha = X.alpha, X.phi = X.phi, X.delta = X.delta, 
      delta = delta, p = p, p.a = p.a
    )
  )
}

