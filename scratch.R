#
# Reproduce Welling-Teh bimodal posterior.
#

# Parameters from Welling-Teh.
welling.teh = list(theta.1 = 0, 
                   theta.2 = 1, 
                   sigma.x.2 = 2, 
                   sigma.1.2 = 10, 
                   sigma.2.2 = 1)

# Generate data.
data = r.bimodal(100, welling.teh)

# Metropolis-Hastings.
mh.draws = NA
mh.time = microbenchmark(mh.draws <- mcmc(start = c(-5, 2), 
             iterations = 1150, 
             burn = 150,
             type = "mh",
             d.posterior = function(proposal) { d.bimodal(proposal, data, welling.teh)},
             r.proposal = function(n, mean) { rnorm(n, mean = mean, sd = 0.5) }
             ), times = 3)

mh.esps = min(effectiveSize(as.mcmc(mh.draws))) / (mh.time$time[1] / 1e9)
plot(mh.draws, pch = ".", col = rgb(0, 0, 0, 0.7), xlim = c(-1.3, 2.3), ylim = c(-3, 3))

# HMC.
hmc.time = microbenchmark(hmc.draws <- mcmc(start = c(-5, 2),
                               iterations = 1150,
                               burn = 150,
                               type = "hmc",
                               U = function(proposal) { -d.bimodal(proposal, data, welling.teh)},
                               grad.U = function(proposal) { -grad.d.bimodal(proposal, data, welling.teh)},
                               epsilon = 0.15), times = 1)

hmc.esps = min(effectiveSize(as.mcmc(hmc.draws))) / (hmc.time$time[1] / 1e9)
plot(hmc.draws, pch = ".", col = rgb(0, 0, 0, 0.5), xlim = c(-1.3, 2.3), ylim = c(-3, 3))
