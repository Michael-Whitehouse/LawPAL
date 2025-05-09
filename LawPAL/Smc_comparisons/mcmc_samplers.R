library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)

# covar <- cov(t(mcmc_chain$param_samples[1:4,]))
# rmvnorm(1,mean = c(0.3,0.2,0.5,0.1), 2.38^2/10*covar)

logprior_poisson <- function(params){
  return(dmvnorm(params, mean = c(0,0,0.5,0), sigma = diag(c(10,10,10,10), nrow = 4), log = TRUE))
}

proposal <- function(par, var, indicator = 0){
  prop <- rnorm(1, par, var)
  if(prop <= 0){prop = par}
  if(prop>=1 && indicator == 1){prop = par}
  return(prop)
}

joint_proposal <- function(par, covar){
  rmvnorm(1,mean = par, 2.38^2/4*covar)
}



LawPAL_mcmc <- function(y, init_dist, init_params, n_iter, rw_params){
  param_samples <- matrix(nrow = 4, ncol = n_iter+1)
  param_samples[,1] <- init_params
  ind <- c(0,0,1,0)
  old_Loglik <- SIR_approx_lik_over(y, init_dist, params = init_params)
  accepted_params <- c(0,0,0,0)
  for (i in 1:n_iter) {
    print(i)
    sample <- param_samples[,i]
    for (j in 1:4){
      prop <- sample
      prop[j] <- proposal(param_samples[j,i], rw_params[j], ind[j])
      
      new_Loglik <- SIR_approx_lik_over(y, init_dist, params = prop)
      prior_diff <- logprior_poisson(prop) - logprior_poisson(sample)
      Log_lik_diff <- new_Loglik - old_Loglik
      
      u <- runif(1)
      
      if(u < exp(prior_diff + Log_lik_diff)){
        sample <- prop
        old_Loglik <- new_Loglik
        accepted_params[j] = accepted_params[j] + 1 
        #- as.numeric(prop[j]==sample[j]) 
        print('accepted')
        print(sample)
      }
      else{print('rejected')}
    }
    param_samples[,i+1] = sample
  }
  
  acceptance_ratio <- accepted_params/n_iter
  
  out <- list(param_samples = param_samples, acceptance_ratio = acceptance_ratio)
}


LawPAL_mcmc_prop <- function(y, init_dist, init_params, n_iter){
  param_samples <- matrix(nrow = 4, ncol = n_iter+1)
  param_samples[,1] <- init_params
  old_Loglik <- SIR_approx_lik_over(y, init_dist, params = init_params)
  accepted_params <- 0
  for (i in 1:n_iter) {
    print(i)
    sample <- param_samples[,i]
    prop <- joint_proposal(sample,covar)
    
    new_Loglik <- SIR_approx_lik_over(y, init_dist, params = prop)
    prior_diff <- logprior_poisson(prop) - logprior_poisson(sample)
    Log_lik_diff <- new_Loglik - old_Loglik
    
    u <- runif(1)
    condition <- u < exp(prior_diff + Log_lik_diff)
    
    if(is.na(condition)){condition <-FALSE}   
    print(condition)
    if(condition){
      sample <- prop
      old_Loglik <- new_Loglik
      accepted_params = accepted_params +1
      print('accepted')
      print(sample)
    }
    else{print('rejected')}
    param_samples[,i+1] = sample
  }
  
  acceptance_ratio <- accepted_params/n_iter
  
  out <- list(param_samples = param_samples, acceptance_ratio = acceptance_ratio)
}


pmmh_prop <- function(y, init_dist, init_params, n_iter){
  param_samples <- matrix(nrow = 4, ncol = n_iter+1)
  param_samples[,1] <- init_params
  old_Loglik <- Particle_likelihood_SIR(y, init_dist, params = init_params,1000)
  accepted_params <- 0
  for (i in 1:n_iter) {
    print(i)
    sample <- param_samples[,i]
    prop <- joint_proposal(sample,covar)
    
    try(new_Loglik <- Particle_likelihood_SIR(y, init_dist, params = prop,1000))
    prior_diff <- logprior_poisson(prop) - logprior_poisson(sample)
    Log_lik_diff <- new_Loglik - old_Loglik
    
    u <- runif(1)
    condition <- u < exp(prior_diff + Log_lik_diff)
    
    if(is.na(condition)){condition <-FALSE}   
    print(condition)
    if(condition){
      sample <- prop
      old_Loglik <- new_Loglik
      accepted_params = accepted_params +1
      print('accepted')
      print(sample)
    }
    else{print('rejected')}
    param_samples[,i+1] = sample
  }
  
  acceptance_ratio <- accepted_params/n_iter
  
  out <- list(param_samples = param_samples, acceptance_ratio = acceptance_ratio)
}



