# Libraries.
library(mvtnorm)
library(coda)
library(microbenchmark)
library(mixtools)
library(MASS)
library(RColorBrewer)

#
# Welling-Teh functions.
# 

# Bimodal posterior draws.
r.bimodal = function(n, params) {
  rnormmix(n, c(0.5, 0.5), 
           c(params[["theta.1"]], params[["theta.1"]] + params[["theta.2"]]),
           rep(sqrt(params[["sigma.x.2"]]), 2)
           )
}

# Bimodal posterior density.
d.bimodal = function(proposal, data, params) {
  
  theta.1 = proposal[1]
  theta.2 = proposal[2]
  sigma.x.2 = params[["sigma.x.2"]]
  sigma.1.2 = params[["sigma.1.2"]]
  sigma.2.2 = params[["sigma.2.2"]]
  
  theta.cov = diag(c(sigma.1.2, sigma.2.2))
  (sum(log((1/2) * (dnorm(data, theta.1, sqrt(sigma.x.2)) + dnorm(data, theta.1 + theta.2, sqrt(sigma.x.2))))) +
    log(dmvnorm(proposal, c(0, 0), theta.cov)))
}

# Gradient of bimodal posterior density.
grad.d.bimodal = function(proposal, data, params) {
  
  theta.1 = proposal[1]
  theta.2 = proposal[2]
  sigma.x.2 = params[["sigma.x.2"]]
  sigma.1.2 = params[["sigma.1.2"]]
  sigma.2.2 = params[["sigma.2.2"]]
  
  d.1 = dnorm(data, theta.1, sqrt(sigma.x.2))
  d.2 = dnorm(data, theta.1 + theta.2, sqrt(sigma.x.2))
  
  common = ((1 / 2) * (d.1 + d.2))^(-1)
  
  grad.1 = common * (1 / 2) * (d.1 * sigma.x.2^(-1) * (data - theta.1) + d.2 * sigma.x.2^(-1) * (data - theta.1 - theta.2))
  
  grad.2 = common * (1 / 2) * (d.2 * sigma.x.2^(-1) * (data - theta.1 - theta.2))
  
  theta.cov = diag(c(sigma.1.2, sigma.2.2))
  prior = - solve(theta.cov) %*% t(t(c(theta.1, theta.2)))
  
  t(c(sum(grad.1), sum(grad.2)) + prior)
  
}

mc.fisher = function(proposal, params, sample.n, pseudo.n) {
  
  params[["theta.1"]] = proposal[1]
  params[["theta.2"]] = proposal[2]
  
  theta.1 = params[["theta.1"]]
  theta.2 = params[["theta.2"]]
  sigma.x.2 = params[["sigma.x.2"]]
  sigma.1.2 = params[["sigma.1.2"]]
  sigma.2.2 = params[["sigma.2.2"]]
  
  hessian = matrix(rep(0, length(proposal) * 2), length(proposal), length(proposal))
  
  for (i in 1:sample.n) {
    pseudodata = r.bimodal(pseudo.n, params)
    
    perturb = t((rbinom(length(proposal), 1, 0.5) * 2 - 1) * 0.0001)
    
    diff.g = t(grad.d.bimodal(proposal + perturb, pseudodata, params) - grad.d.bimodal(proposal - perturb, pseudodata, params))
    hessian = hessian + (diff.g %*% perturb^(-1) / 2 + t(diff.g %*% perturb^(-1)) / 2) / 2
  }
  
  - hessian / (sample.n * pseudo.n)
}

microbenchmark(mc.fisher(c(1,2), welling.teh, 1000, 100),
               mc.fisher(c(1,2), welling.teh, 100, 100),
               times = 5)


mc.fisher(c(1,2), welling.teh, 50, 100)

#
# MCMC algorithms.
#

# MCMC wrapper function.
mcmc = function(start, iterations, burn, type, ...) {
  
  # Initialize storage.
  draws = array(dim = c(iterations + 1, length(start)))
  draws[1, ] = start
  
  # Iterate.
  for (i in 1:iterations) {
    draws[i + 1, ] = 
      if (type == "mh") {
        mcmc.mh(draws[i, ], ...)
      } else if (type == "hmc") {
        mcmc.hmc(draws[i,], ...)
      }
  }
  
  tail(draws, -burn)
  
}

# Metropolis-Hastings with normal proposals.
mcmc.mh = function(current, data, d.posterior, r.proposal) {
  
    # Propose and find density according to random walk.
    proposal = r.proposal(length(current), current)
    alpha = min(1, exp(d.posterior(proposal) - d.posterior(current)))
    
    # Accept or reject.
    if (runif(1) < alpha) {
      proposal
    } else {
      current
    }
  
}

# Standard HMC.
mcmc.hmc = function(current.q, U, grad.U, epsilon) {
  
  # Initialize state and momentum.
  q = current.q
  p = rnorm(length(q), 0, 1)
  current.p = p
  L = sample(13:17, 1)
  
  # Make a half step for momentum at the beginning.
  p = p - epsilon * grad.U(q) / 2
  
  for (i in 1:L) {
    
    # Make a full step for the position.
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory.
    if (i != L) {
      p = p - epsilon * grad.U(q)
    } 
    
  }
  
  # Make a half step for momentum at the end.
  p = p - epsilon * grad.U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric.
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory.
  current.U = U(current.q) 
  current.K = sum(current.p^2) / 2 
  proposed.U = U(q)
  proposed.K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either 
  # the position at the end of the trajectory or the initial position.
  if (runif(1) < min(1, exp(current.U - proposed.U + current.K - proposed.K))) {
    q
  } else {
    current.q
  }
}

mcmc.rmhmc = function(current.theta, H, rm.U, rm.grad.U, fixed.point.steps) {
  theta = current.theta
  p = rnorm(length(theta), 0, 1)
  current.H = H(theta, p)
  N = sample(10:25, 1)
  
  for (n in 1:N) {
    
    # Update momentum with fixed point iterations.
    hat.p = p
    for (i in 1:fixed.point.steps) {
      hat.p = p - epsilon * rm.grad.U(theta, hat.p) / 2
    }
    p = hat.p
    
    # Update parameters with fixed point iterations.
    hat.theta = theta
    for (i in 1:fixed.point.steps) {
      hat.theta = theta + (epsilon * rm.U(theta, p) + epsilon * rm.U(hat.theta, p)) / 2
    }
    theta = hat.theta
    
    # Update momentum exactly.
    p = p - epsilon * rm.grad.U(theta, p)
  }
  
  proposed.H = H(theta, p)
  
  ratio = -log(proposed.H) + log(current.H)
  
  if (ratio > 0 | ratio > log(runif(1))) {
    theta
  } else{
    current.theta
  }
  
}