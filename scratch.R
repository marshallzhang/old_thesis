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
             iterations = 5550, 
             burn = 550,
             type = "mh",
             trace = 500,
             d.posterior = function(proposal) { d.bimodal(proposal, data, welling.teh)},
             r.proposal = function(n, mean) { rnorm(n, mean = mean, sd = 0.5) }
             ), times = 1)

mh.esps = min(effectiveSize(as.mcmc(mh.draws))) / (mh.time$time[1] / 1e9)
z = kde2d(mh.draws[, 1], mh.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-3, 3), ylim = c(-3, 3))

# HMC.
hmc.time = microbenchmark(hmc.draws <- mcmc(start = c(-5, 2),
                               iterations = 5550,
                               burn = 550,
                               type = "hmc",
                               trace = 500,
                               U = function(proposal) { -d.bimodal(proposal, data, welling.teh)},
                               grad.U = function(proposal) { -grad.d.bimodal(proposal, data, welling.teh)},
                               epsilon = 0.5), times = 1)

hmc.esps = min(effectiveSize(as.mcmc(hmc.draws))) / (hmc.time$time[1] / 1e9)
z = kde2d(hmc.draws[, 1], hmc.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-3, 3), ylim = c(-3, 3))

# SRM-HMC.
sample.n = 1000
pseudo.n = 1000
srm.hmc.time = microbenchmark(srm.hmc.draws <- mcmc(start = c(-5, 2),
                               iterations = 500,
                               burn = 50,
                               type = "srm-hmc",
                               trace = 1,
                               params = welling.teh,
                               sample.n = sample.n,
                               pseudo.n = pseudo.n,
                               H = function(theta, p) { H(theta, p, data, welling.teh, sample.n, pseudo.n)},
                               grad.theta.H = function(theta, p) { grad.theta.H(theta, p, data, welling.teh, sample.n, pseudo.n)},
                               grad.p.H = function(theta, p) { grad.p.H(theta, p, data, welling.teh, sample.n, pseudo.n)},
                               fixed.point.steps = 3,
                               epsilon = 0.5), times = 1)

srm.hmc.esps = min(effectiveSize(as.mcmc(srm.hmc.draws))) / (srm.hmc.time$time[1] / 1e9)
z = kde2d(srm.hmc.draws[, 1], srm.hmc.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-3, 3), ylim = c(-3, 3))
