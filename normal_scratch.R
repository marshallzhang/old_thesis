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
data = r.normal(30, c(0, 10))

# Metropolis-Hastings.
mh.draws = NA
mh.time = microbenchmark(mh.draws <- mcmc(start = c(5, 40), 
             iterations = 200,
             burn = 0,
             type = "mh",
             trace = 500,
             d.posterior = function(proposal) { d.normal(proposal, data)},
             r.proposal = function(n, mean) { rnorm(n, mean = mean, sd = 0.75) }
             ), times = 1)

mh.esps = min(effectiveSize(as.mcmc(mh.draws))) / (mh.time$time[1] / 1e9)
z = kde2d(mh.draws[, 1], mh.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-10, 10), ylim = c(5, 40))
lines(mh.draws[,1], mh.draws[,2])

# HMC.
hmc.time = microbenchmark(hmc.draws <- mcmc(start = c(5, 40),
                               iterations = 200,
                               burn = 0,
                               type = "hmc",
                               trace = 500,
                               U = function(proposal) { -d.normal(proposal, data)},
                               grad.U = function(proposal) { -grad.d.normal(proposal, data)},
                               epsilon = 0.75), times = 1)

hmc.esps = min(effectiveSize(as.mcmc(hmc.draws))) / (hmc.time$time[1] / 1e9)
z = kde2d(hmc.draws[, 1], hmc.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-10, 10), ylim = c(5, 40))
lines(hmc.draws[,1], hmc.draws[,2])

# RM-HMC.
rm.hmc.time = microbenchmark(rm.hmc.draws <- mcmc(start = c(5, 40),
                               iterations = 200,
                               burn = 0,
                               type = "srm-hmc",
                               trace = 1,
                               params = welling.teh,
                               sample.n = sample.n,
                               pseudo.n = pseudo.n,
                               H = function(theta, p) { H(theta, p, data)},
                               grad.theta.H = function(theta, p) { grad.theta.H(theta, p, data)},
                               grad.p.H = function(theta, p) { grad.p.H(theta, p, data)},
                               fixed.point.steps = 6,
                               epsilon = 0.75), times = 1)

rm.hmc.esps = min(effectiveSize(as.mcmc(rm.hmc.draws))) / (rm.hmc.time$time[1] / 1e9)
z = kde2d(rm.hmc.draws[, 1], rm.hmc.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-10, 10), ylim = c(5, 40))
lines(rm.hmc.draws[,1], rm.hmc.draws[,2])

# SRM-HMC.
sample.n = 100
pseudo.n = 100
srm.hmc.time = microbenchmark(srm.hmc.draws <- mcmc(start = c(5, 40),
                               iterations = 200,
                               burn = 0,
                               type = "srm-hmc",
                               trace = 1,
                               params = welling.teh,
                               sample.n = sample.n,
                               pseudo.n = pseudo.n,
                               H = function(theta, p) { H(theta, p, data, stochastic = T, sample.n = sample.n, pseudo.n = pseudo.n)},
                               grad.theta.H = function(theta, p) { grad.theta.H(theta, p, data, stochastic = T, sample.n = sample.n, pseudo.n = pseudo.n)},
                               grad.p.H = function(theta, p) { grad.p.H(theta, p, data, stochastic = T, sample.n = sample.n, pseudo.n = pseudo.n)},
                               fixed.point.steps = 6,
                               epsilon = 0.75), times = 1)

srm.hmc.esps = min(effectiveSize(as.mcmc(srm.hmc.draws))) / (srm.hmc.time$time[1] / 1e9)
z = kde2d(srm.hmc.draws[, 1], srm.hmc.draws[, 2])
contour(z, col = brewer.pal(9, "Blues"), xlim = c(-10, 10), ylim = c(5, 40))
lines(srm.hmc.draws[,1], srm.hmc.draws[,2])