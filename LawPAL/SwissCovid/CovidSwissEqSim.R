library(truncnorm)
library(matrixStats)

simulate_SEIR_forced <- function(
    T = 100,
    N = 8570000,
    beta = 1.5,
    rho = 0.38,
    gamma = 0.38,
    I0 = 20,
    E0 = 10,
    R0 = 0,
    q = 0.75,
    sigma_q = 0.2,
    eta = 0.1,
    nu = 3.5,
    xi = 0.6,
    t_q = 24,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  switch_eta <- function(t, eta, nu, xi, t_q) {
    eta + (1 - eta) / (1 + exp(xi * (t - (t_q) - nu)))
  }
  
  S <- numeric(T); E <- numeric(T); I <- numeric(T); R <- numeric(T)
  reported_cases <- numeric(T)
  xi <- xi + 0.5
  S[1] <- rpois(1,N - I0 - E0 - R0)
  E[1] <- rpois(1,E0)
  I[1] <- rpois(1,I0)
  R[1] <- R0
  
  for (t in 2:T) {
    forcing <- switch_eta(t-1, eta, nu, xi, t_q)
    beta_eff <- beta * forcing
    
    lambda <- beta_eff * I[t-1] / N
    new_E <- rbinom(1, S[t-1], 1 - exp(-lambda))
    new_I <- rbinom(1, E[t-1], 1 - exp(-rho))
    new_R <- rbinom(1, I[t-1], 1 - exp(-gamma))
    reported_cases[t-1] <- rbinom(1, new_I, q)
    
    S[t] <- S[t-1] - new_E
    E[t] <- E[t-1] + new_E - new_I
    I[t] <- I[t-1] + new_I - new_R
    R[t] <- R[t-1] + new_R
  }
  
  return(list(
    time = 1:T,
    S = S,
    E = E,
    I = I,
    R = R,
    reported_cases = reported_cases
  ))
}
par(mfrow = c(1,1))

indices <- sample(10000:50000,1000)
i=1
sim <- simulate_SEIR_forced(T=109, N = 8570000 , beta = samples_eq$beta[i],rho = samples_eq$rho[i], gamma = samples_eq$gamma[i],I0 = samples_eq$i0[i],
                            E0 = samples_eq$e0[i],R0 = 0, q = samples_eq$q[i],sigma_q = 2, eta = samples_eq$eta[i],
                            nu = samples_eq$nu[i], xi = samples_eq$xi_raw[i])

nsamples = 1000
group1 <- matrix(nrow = nsamples, ncol = 109)
par(cex.main = 1.5,cex.lab = 1.5)
ts.plot(sim$reported_cases, ylim = c(0,2300), ylab = "Reported Cases" , main = "Posterior predictive plot: equi-dispersed observations", col = alpha('black', 0.2) )
j=0
for (i in indices) {
  j=j+1
  sim <- simulate_SEIR_forced(T=109, N = 8570000 , beta = samples$beta[i],rho = samples$rho[i], gamma = samples$gamma[i],I0 = samples$i0[i],
                              E0 = samples$e0[i],R0 = 0, q = samples$q[i],sigma_q = 2, eta = samples$eta[i],
                              nu = samples$nu[i], xi = samples$xi_raw[i])
  if(i%%100==0){
    lines(sim$reported_cases, col = alpha(i/100, 0.5) )}
  group1[j,] <- sim$reported_cases
  
}
lines(cases, col = 'red')

par(mfrow = c(1,1))
par(cex.main = 1.5,cex.lab = 1.5)
ts.plot(colMeans(group1),ylim = c(0, 2300), ylab = 'New reported cases', main = "Posterior predictive check: equi-dispersed observations")
polygon(c(1:109,rev(1:109)),c(colQuantiles(group1, probs = 0.75),rev(colQuantiles(group1, probs = 0.25))),lty=0,col=rgb(0,0.3,1,0.4))
polygon(c(1:109,rev(1:109)),c(colQuantiles(group1, probs = 0.975),rev(colQuantiles(group1, probs = 0.025))),lty=0,col=rgb(0,0.3,1,0.2))
points(cases, pch = 19)
legend("topright", legend=c(" Data", " Posterior predictive mean", " 50% credible interval", " 90% credible interval"), 
       box.lwd = 1, bty = 'n',bg = "white",
       col=c("black", "blue", rgb(0,0.3,1,0.65), rgb(0,0.3,1,0.2)), lty= c(NA, 1, 1, 1), lwd = c(NA, 2, 8, 10), pch=c(19, NA, NA, NA),
       seg.len=0.25, y.intersp=0.65, x.intersp=0.25, cex=1.4)


