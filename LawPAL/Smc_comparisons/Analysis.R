library(Rcpp)
library(truncnorm)
sourceCpp('SIRover.cpp')
sourceCpp('mcmc/cpp_LawPAL.cpp')
sourceCpp('mcmc/SMC_likelihood.cpp')


T_length <- 100
# set population size
# n <-  
init_dist <- c(0.995,0.005,0)
xi <- rep(1,T_length)
params = c(0.3,0.2)
qpars <- c(0.5,0.1)
par = c(params,qpars)
q = rtruncnorm(T_length,0,1,qpars[1],qpars[2])
sim1 <- SIRsim_over(T_length, n, init_dist, xi, params, (q),1,1)
par(mfrow = c(1,1))
ts.plot(sim1$obs)
sim1$obs

pastel_colors <- hcl(h = seq(15, 375, length = 5)[1:4], c = 35, l = 85)
parnames = c(expression(beta),expression(gamma),expression(mu[q]),expression(sigma[q]^2))

# burnin
mcmc_chain <- LawPAL_mcmc(sim1$obs, n*init_dist, c(0.3,0.2,0.5,0.1), 10000, rw_params = c(0.01,0.01,0.05,0.01))
mcmc_chain$acceptance_ratio


# set covariance proposal
covar <- cov(t(mcmc_chain$param_samples[1:4,]))

# burn in pmmh
pmmh_chain <- pmmh_prop(sim1$obs, n*init_dist, c(0.3,0.2,0.5,0.1), 10000)

#sample
mcmc_chain <- LawPAL_mcmc_prop(sim1$obs, n*init_dist, c(0.3,0.2,0.5,0.1), 50000)

covar <- 1*cov(t(pmmh_chain$param_samples[1:4,1:1000]))
pmmh_chain <- pmmh_prop(sim1$obs, n*init_dist, c(0.3,0.2,0.5,0.1), 50000)




### plots

indices <- sample(5000:50000, 30000)

for (i in 1:4) {
  hist(mcmc_chain$param_samples[i,indices],  main = '', xlab = parnames[i], freq = F, breaks = 20, col = pastel_colors[i])
  #abline(v = par[i], col = 'red', lwd = 3)
}
mtext("LawPALmh", side = 3, outer = TRUE, line = -1, cex = 1.25)

for (i in 1:4) {
  hist(pmmh_chain$param_samples[i,indices],  main = '', xlab = parnames[i], freq = F, breaks = 20, col = pastel_colors[i])
  #abline(v = par[i], col = 'red', lwd = 3)
}
mtext("pmmh", side = 3, outer = TRUE, line = -1, cex = 1.25)

# Diagnostics
par(mfrow= c(4,2))
# Traceplots
for (i in 1:4) {
  ts.plot(mcmc_chain$param_samples[i,],  main = '', ylab = parnames[i], xlab = 'iteration',  col = pastel_colors[i])

}

for (i in 1:4) {
  acf(mcmc_chain$param_samples[i,],  main = '', ylab = parnames[i], xlab = 'iteration',  col = pastel_colors[i])

}

par(mfrow= c(4,2))
# Traceplots
for (i in 1:4) {
  ts.plot(pmmh_chain$param_samples[i,],  main = '', ylab = parnames[i], xlab = 'iteration',  col = pastel_colors[i])

}

for (i in 1:4) {
  acf(pmmh_chain$param_samples[i,],  main = '', ylab = parnames[i], xlab = 'iteration',  col = pastel_colors[i])
}


