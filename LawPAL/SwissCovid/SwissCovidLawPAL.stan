functions{
  real switch_alpha(int t, real alpha, real d, real b) {
    return(alpha + (1 - alpha) / (1 + exp(b * (t - (23) - d))));
  }
  
  real LawPAL(int[] y, real[] theta){
    
    real beta = theta[1];
    real gamma = theta[2];
    real rho = theta[3];
    real alpha = theta[4]; // reduction in transmission rate after quarantine
    real d = theta[5]; // shift of quarantine implementation
    real b = theta[6]; // slope of quarantine implementation
    real i0 = theta[7];
    real e0 = theta[8];
    real muq = theta[9];
    real sigma_q = theta[10];
    real lik = 0;
    
    matrix[4,4] K = [[0,0,0,0],[0,exp(-rho),1-exp(-rho),0],[0,0,exp(-gamma),1-exp(-gamma)],[0,0,0,1]];
    matrix[4,4] Lambda = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    matrix[4,4] barLambda = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    vector[4] lam = [8570000 - i0 - e0 , e0, i0, 0]'; // initialise
    
    for(t in 1:109){
      
      real forcing_function = switch_alpha(t,alpha,d,b); // switch function
      real beta_eff = beta * forcing_function;
      K[1,2] = 1-exp(-lam[3]*beta_eff/8570000);
      K[1,1] = 1 - K[1,2];
      
      Lambda = (lam*[1,1,1,1]).*K;
      
      real qbar = 0.5*(muq-Lambda[2,3]*sigma_q*sigma_q + sqrt((Lambda[2,3]*sigma_q*sigma_q-muq)*(Lambda[2,3]*sigma_q*sigma_q-muq) + 4*y[t]*sigma_q*sigma_q));
      real s = sqrt(1/(y[t]/(qbar*qbar) + 1/(sigma_q*sigma_q)));

      
      lik += poisson_lpmf(y[t] |qbar * Lambda[2,3]) + normal_lpdf(qbar|muq,sigma_q) + log(s*sqrt(2*3.14159));
      
      barLambda = Lambda;
      
      barLambda[2,3] = y[t] +(1-qbar)*Lambda[2,3];
      
      lam = (([1,1,1,1])*barLambda)';
      
  }
  return(lik);
}
}

data{
 int y[109];
}

parameters {
  real<lower=0> beta; // force of infection
  real<lower=0> rho; // exposed period
  real<lower=0> gamma; // recovery rate
  real<lower=0> sigma_q; // 
  real<lower=0,upper=1> alpha; // reduction in transmission due to control measures (in proportion of beta)
  real<lower=0> d; // shift of quarantine implementation (strictly positive as it can only occur after tswitch)
  real<lower=0,upper=1> b_raw; // slope of quarantine implementation (strictly positive as the logistic must be downward)
  real<lower=0, upper=1> q; // proportion of infected (symptomatic) people reported
  real<lower=0> i0; // number of infected people inititally
  real<lower=0> e0;
}


transformed parameters{
  real b = b_raw + 0.5;
  real theta[10];
  theta = {beta, gamma, rho, alpha, d, b, i0, e0,q,sigma_q};
}

model{
  beta ~ normal(2, 0.5);
  gamma ~ normal(0.2, 0.1);
  rho ~ normal(0.2, 0.1);
  sigma_q ~ exponential(0.1);
  i0 ~ normal(0, 20);
  e0 ~ normal(0, 20);
  q ~ beta(1, 2);
  alpha ~ beta(2.5, 4);
  d ~ exponential(1./5);
  b_raw ~ beta(1, 1);
  
  target += LawPAL(y, theta);
}

