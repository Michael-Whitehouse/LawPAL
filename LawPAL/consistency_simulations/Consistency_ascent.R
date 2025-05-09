library(Rcpp)
library(truncnorm)
sourceCpp('SIRover.cpp')


q_laplace_filter <- function(y,params, qparams, initdist,n){
  
  K <- matrix(data = 0,3,3)  
  K[1,2] = 1 - exp(-params[1]*init_dist[2])
  K[2,3] = 1 - exp(- params[2]) 
  K[1,1] = exp(-params[1]*init_dist[2])
  K[2,2] = exp(- params[2])
  K[3,3] = 1
  K1 <- K
  q_est <- c()
  qvar <- c()
  q1 = 0.5
  state = matrix(0,3, T_length+1)
  state[,1] = n*initdist
  state1 <- state
  PALyq <- 0
  for (t in 1:T_length){
    Lambda_ = (outer(state[,t],c(1,1,1)))*K
    
    q_est[t] = 0.5*(qparams[1] - Lambda_[1,2]*qparams[2]^2 + sqrt((Lambda_[1,2]*qparams[2]^2 - qparams[1])^2 + 4*y[t]*qparams[2]^2))
    
    Lambda = Lambda_
    
    Lambda[1,2] = (1-q_est[t])*Lambda_[1,2] + y[t]
    
    state[,t+1] = colSums(Lambda)
    
    K[1,2] = 1 - exp(-params[1]*state[2,t+1]/sum(state[,t+1]))
    K[1,1] = exp(-params[1]*state[2,t+1]/sum(state[,t+1]))
    qvar[t] = 1/(y[t]/q_est[t]^2 + 1/qparams[2]^2)
    PALyq <- PALyq + dpois(y[t], q_est[t]*Lambda_[1,2], log = T) + log(dtruncnorm(q_est[t],0,1,qparams[1],qparams[2]))  + log(sqrt(qvar[t]*2*pi))
  }
  
  log_lik <- PALyq
  
  return(list(q_est = q_est, qvar = qvar,log_lik=log_lik))
}


Ns <- c(5000,10e3,10e4,10e5)
T_lengths = c(50,100,150,200)
init_dist = c(0.995,0.005,0)
muMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
sigmaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
betaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
gammaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
params = c(0.15,0.1)
qpars = c(0.5,0.1)
n = Ns[3]
T_length = T_lengths[4]
q = rtruncnorm(T_length,0,1,qpars[1],qpars[2])   
xi = rep(1,T_length)
sim1 <- SIRsim_over(T_lengths[4], n, init_dist, xi, params, (q),1,1)
ts.plot(sim1$obs)
sim1$obs[1]
ts.plot(sim1$state[2,])
par <- c(params[1], params[2], qpars[1], qpars[2])
coordinate_ascent_alg <- function(n_steps,n){
  print(par)
  tstart <- Sys.time()
  traj <- matrix(nrow = 4, ncol= n_steps+1)
  lik <- c()
  lik[1] = q_laplace_filter(sim1$obs, c(par[1],par[2]),c(par[3],par[4]), init_dist, n)$log_lik
  for (i in 1:n_steps) {
    if(i%%100==0){print(paste('iteration =',i))}
    t1 <- Sys.time()
    
    for (j in 1:4) {
      # print(j)
      #Delta <- sample(c(-1,1),1)
      Delta <- 1
      par_pos <- par
      par_neg <- par
      par_pos[j] <- par[j]+0.001
      par_neg[j] <- par[j]-0.001
      
      grad_pos <- q_laplace_filter(sim1$obs, c(par_pos[1],par_pos[2]),c(par_pos[3],par_pos[4]), init_dist, n)$log_lik
      grad_neg <- q_laplace_filter(sim1$obs, c(par_neg[1],par_neg[2]),c(par_neg[3],par_neg[4]), init_dist, n)$log_lik
    
      grad_est <- sign(grad_pos - grad_neg)
 
      par[j] <- par[j] + 0.0001*grad_est
    }
    traj[,i+1] = par
    lik[i+1] <- q_laplace_filter(sim1$obs, c(par[1],par[2]),c(par[3],par[4]), init_dist, n)$log_lik
    t2 <- Sys.time()
  }
  tend <- Sys.time()
  time <- tstart - tend
  output <- list(traj = traj, lik =lik, time=time)
  return(output)
}


opt <- coordinate_ascent_alg(1000,n)
ts.plot(opt$lik)
par(mfrow = c(2,2))
ts.plot(opt$traj[1,])
ts.plot(opt$traj[2,])
ts.plot(opt$traj[3,])
ts.plot(opt$traj[4,])
opt$traj[,1001]






# Ns <- c(5000,10e3,10e4,10e5)
# T_lengths = c(50,100,150,200)

Ns <- c(5000,10e3,10e4,10e5)
T_lengths = c(50,100,150,200)
init_dist = c(0.995,0.005,0)
muMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
sigmaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
betaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
gammaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
params = c(0.15,0.1)
qpars = c(0.5,0.1)
init_dist = c(0.995,0.005,0)
samples=100
muMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
sigmaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
betaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
gammaMLEsamples <- array(0, dim = c(length(Ns),length(T_lengths),100))
params = c(0.15,0.1)
qpars = c(0.5,0.1)

count = 0
countT = 0
for (N in Ns) {
  countT = 0
  count = count + 1
  print("count = ")
  print(count)
  for (T_length in T_lengths) {
    countT = countT + 1
    print("countT = ")
    print(countT)
    for (i in 1:100) {
      print(i)
      par <- c(params[1], params[2], qpars[1], qpars[2])
      q = rtruncnorm(T_length,0,1,qpars[1],qpars[2])   
      xi = rep(1,T_length)
      sim1 <- SIRsim_over(T_length, N, init_dist, xi, params, (q),1,1)
      opt <- coordinate_ascent_alg(1000,N)
      
      print(opt$traj[,1000])
      betaMLEsamples[count,countT, i] = opt$traj[1,1000]
      gammaMLEsamples[count,countT, i] = opt$traj[2,1000]
      muMLEsamples[count, countT, i] = opt$traj[3,1000]
      sigmaMLEsamples[count,countT, i] = opt$traj[4,1000]
    }
  }
}



par(mfrow = c(4,4))
for (i in 1:4) {
  for(j in 1:4){
    boxplot(MLEsamples$betaMLEsamples[i,j,], ylim = c(0.1,0.25), main = paste('time ', j, 'n', i,'mean',mean(MLEsamples$betaMLEsamples[i,j,], na.rm = T), 'sd',sd(MLEsamples$betaMLEsamples[i,j,], na.rm = T) ))
    print(c(i,j))
    print("mean = ")
    print(mean(MLEsamples$betaMLEsamples[i,j,], na.rm = T))
    print("sd = ")
    print(sd(MLEsamples$betaMLEsamples[i,j,], na.rm = T))
    abline( h = 0.15, col = 'red')
  }
}

par(mfrow = c(4,4))
for (i in 1:4) {
  for(j in 1:4){
    boxplot(MLEsamples$gammaMLEsamples[i,j,], ylim = c(0.01,0.2), main = paste('time ', j, 'n', i,'mean',mean(MLEsamples$gammaMLEsamples[i,j,], na.rm = T), 'sd',sd(MLEsamples$gammaMLEsamples[i,j,], na.rm = T) ))
    print(c(i,j))
    print("mean = ")
    print(mean(MLEsamples$gammaMLEsamples[i,j,], na.rm = T))
    print("sd = ")
    print(sd(MLEsamples$gammaMLEsamples[i,j,], na.rm = T))
    abline( h = 0.1, col = 'red')
  }
}

par(mfrow = c(4,4))
for (i in 1:4) {
  for(j in 1:4){
    boxplot(MLEsamples$muMLEsamples[i,j,], ylim = c(0.35,0.65), main = paste('time ', j, 'n', i,'mean',mean(MLEsamples$muMLEsamples[i,j,], na.rm = T), 'sd',sd(MLEsamples$muMLEsamples[i,j,], na.rm = T) ))
    print(c(i,j))
    print("mean = ")
    print(mean(MLEsamples$muMLEsamples[i,j,], na.rm = T))
    print("sd = ")
    print(sd(MLEsamples$muMLEsamples[i,j,], na.rm = T))
    abline( h = 0.5, col = 'red')
  }
}

par(mfrow = c(4,4))
for (i in 1:4) {
  for(j in 1:4){
    boxplot(MLEsamples$sigmaMLEsamples[i,j,], ylim = c(0.01,0.2), main = paste('time ', j, 'n', i,'mean',mean(MLEsamples$sigmaMLEsamples[i,j,], na.rm = T), 'sd',1.96*sd(MLEsamples$sigmaMLEsamples[i,j,], na.rm = T) ))
    print(c(i,j))
    print("mean = ")
    print(mean(MLEsamples$sigmaMLEsamples[i,j,], na.rm = T))
    print("sd = ")
    print(sd(MLEsamples$sigmaMLEsamples[i,j,], na.rm = T))
    abline( h = 0.1, col = 'red')
  }
}

MLEsamples = list(betaMLEsamples=betaMLEsamples,gammaMLEsamples=gammaMLEsamples,muMLEsamples=muMLEsamples,sigmaMLEsamples=sigmaMLEsamples)









