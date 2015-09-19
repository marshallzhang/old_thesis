# Libraries.
library(mvtnorm)
library(coda)
library(microbenchmark)
library(mixtools)
library(MASS)
library(RColorBrewer)

#
# Normal functions.
# 

# Normal draws.
r.normal = function(n, params) {
  rnorm(n, mean = params[1], sd = params[2])
}

# Normal density.
d.normal = function(params, data) {
  sum(dnorm(data, mean = params[1], sd = params[2], log = T))
}

# Gradient of bivariate normal density.
grad.d.normal= function(params, data) {
  -c((length(data) * params[1] - sum(data)) / params[2]^2,
     (length(data) / params[2]) - (sum((data - params[1])^2) / params[2]^3))
}

# Compute numerical Fisher information and its derivatives.
mc.fisher = function(proposal, sample.n, pseudo.n, snd = 0) {
  
  c = 0.0001
  
  # Get perturbation for one dimension if taking derivative of Hessian.
  if (snd > 0) {
    perturb.2 = rep(0, length(proposal))
    perturb.2[snd] = c
  }
  
  # Start Hessian.
  hessian = matrix(rep(0, length(proposal) * 2), length(proposal), length(proposal))
  
  for (i in 1:sample.n) {
    
    # Generate vector from proposed parameters.
    pseudodata = r.normal(pseudo.n, proposal)
    
    # Make perturbation vector.
    perturb = ((rbinom(length(proposal), 1, 0.5) * 2 - 1) * c)
    
    # For derivative of Fisher information, or for regular Fisher information, find numerical derivatives.
    if (snd > 0) {
      diff.g.plus = (grad.d.normal(proposal + perturb + perturb.2, pseudodata) - grad.d.normal(proposal + perturb - perturb.2, pseudodata)) * c^(-1)
      diff.g.minus = (grad.d.normal(proposal - perturb + perturb.2, pseudodata) - grad.d.normal(proposal - perturb - perturb.2, pseudodata)) * c^(-1)
      diff.g = (diff.g.plus - diff.g.minus)
    } else {
      diff.g = (grad.d.normal(proposal + perturb, pseudodata) - grad.d.normal(proposal - perturb, pseudodata))
    }
    
    # Get Hessian.
    hessian = hessian + (diff.g %*% t(perturb^(-1)) / 2 + t(diff.g %*% t(perturb^(-1))) / 2) / 2
  }
  
  - hessian / sample.n
}

# Compute RM Hamiltonian.
H = function(theta, p, data, stochastic = F, ...) {
  if (stochastic) {
    G = mc.fisher(theta, ...)
  } else {
    N = length(data)
    G = diag(c(N / theta[2]^2, (2 * N) / theta[2]^2))
  }
  -d.normal(theta, data) + log((2 * pi)^(2) * det(G)) / 2 + p %*% solve(G) %*% t(p)
}

# Compute gradient of RM Hamiltonian wrt theta.
grad.theta.H = function(theta, p, data, stochastic = F, ...) {
  if (stochastic) {
    G = mc.fisher(theta, ...)
    d.G.1 = mc.fisher(theta, ..., snd = 1)
    d.G.2 = mc.fisher(theta, ..., snd = 2)
  } else {
    N = length(data)
    G = diag(c(N / theta[2]^2, (2 * N) / theta[2]^2))
    d.G.1 = diag(c(0, 0))
    d.G.2 = diag(c(- (2 * N) / theta[2]^3, - (4 * N) / theta[2]^3))
  }

  inv.G = solve(G)
  
  -grad.d.normal(theta, data) + 
  c(sum(diag(inv.G %*% d.G.1)) / 2 - p %*% inv.G %*% d.G.1 %*% inv.G %*% t(p) / 2,
    sum(diag(inv.G %*% d.G.2)) / 2 - p %*% inv.G %*% d.G.2 %*% inv.G %*% t(p) / 2)
}

# Compute gradient of RM Hamiltonian wrt p.
grad.p.H = function(theta, p, data, stochastic = F, ...) {
  if (stochastic) {
    G = mc.fisher(theta, ...)
  } else {
    N = length(data)
    G = diag(c(N / theta[2]^2, (2 * N) / theta[2]^2))
  }
  solve(G) %*% t(p)
}

#
# MCMC algorithms.
#

# MCMC wrapper function.
mcmc = function(start, iterations, burn, type, trace, ...) {
  
  # Initialize storage.
  draws = array(dim = c(iterations + 1, length(start)))
  draws[1, ] = start
  
  # Iterate.
  for (i in 1:iterations) {
    if (trace != 0) if (i %% trace == 0) print(paste("Iteration", i, "of", iterations))
    draws[i + 1, ] = 
      if (type == "mh") {
        mcmc.mh(draws[i, ], ...)
      } else if (type == "hmc") {
        mcmc.hmc(draws[i,], ...)
      } else if (type == "srm-hmc") {
        mcmc.srm.hmc(draws[i,], ...)
      }
  }
  
  if (burn == 0) {
    draws
  } else {
    tail(draws, -burn)
  }
  
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
  L = sample(4:6, 1)
  
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

# RM HMC with stochastic metric.
mcmc.srm.hmc = function(current.theta, params, sample.n, pseudo.n, H, grad.theta.H, grad.p.H, fixed.point.steps, epsilon) {
  theta = current.theta
  G = diag(c(30 / theta[2]^2, (2 * 30) / theta[2]^2))
  p = (rmvnorm(1, c(0, 0), diag(c(1/current.theta[2], 1 / (2 * current.theta[2]^2)))))
  current.H = H(theta, p)
  N = sample(4:6, 1)
  
  for (n in 1:N) {
    
    # Update momentum with fixed point iterations.
    hat.p = p
    for (i in 1:fixed.point.steps) {
      hat.p = p - epsilon * grad.theta.H(theta, hat.p) / 2
    }
    p = hat.p
    
    # Update parameters with fixed point iterations.
    hat.theta = theta
    for (i in 1:fixed.point.steps) {
      hat.theta = theta + (epsilon * grad.p.H(theta, p) + epsilon * grad.p.H(hat.theta, p)) / 2
    }
    theta = hat.theta
    
    # Update momentum exactly.
    p = p - epsilon * grad.theta.H(theta, p) / 2
  }
  
  proposed.H = H(t(theta), p)
  
  ratio = -log(proposed.H) + log(current.H)
  
  if (!is.na(ratio) & (ratio > 0 | ratio > log(runif(1)))) {
    print("accept")
    print(theta)
    theta
  } else{
    current.theta
  }
  
}