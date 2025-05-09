// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]
using namespace Rcpp;
using namespace arma;

void propagate_particles(mat resamp_parts,
                         mat& parts,
                         vec& new_inf,
                         vec& q,
                         vec& mu,
                         vec& sig,
                         vec pars,
                         double y,
                         int n_parts,
                         int pop){
 
  int new_rem;
  
  for(int j = 0; j< n_parts; j++){
    new_inf(j) = R::rbinom(resamp_parts(0,j), 1 - exp(-pars(0)*resamp_parts(1,j)/pop));
    new_rem = R::rbinom(resamp_parts(1,j), 1 - exp(-pars(1)));
    
    parts(0,j) = resamp_parts(0,j) - new_inf(j);
    parts(1,j) = resamp_parts(1,j) + new_inf(j) - new_rem;
    parts(2,j) = resamp_parts(2,j) + new_rem;
    
    mu(j) = (0.5)*((pars(2) -   new_inf(j)*pars(3)*pars(3)) + sqrt((new_inf(j)*pars(3)*pars(3) - pars(2))*(new_inf(j)*pars(3)*pars(3) - pars(2)) + 4*y*pars(3)*pars(3)));
;
    sig(j) = sqrt(1/(y/(mu(j)*mu(j)) + 1/(pars(3)*pars(3))));
    
    q(j) = r_truncnorm( mu(j), sig(j),  0, 1);
    
  }
}


void weight_particles( vec& weights,
                       double obs,
                       vec& mu,
                       vec& sig,
                       int n_parts,
                       vec new_inf,
                       vec pars,
                       vec q
                       ){
  
  for(int i = 0; i< n_parts; i++){
    // weights(i) = R::dbinom(obs, parts(1,i), pars(2), 1);
    weights(i) = R::dpois(obs, new_inf(i)*q(i), 1) + d_truncnorm(q(i), pars(2), pars(3),  0, 1,1) - d_truncnorm(q(i), mu(i), sig(i),  0, 1,1);
  }
}

void normaliseWeights(vec& tWeights,
                      vec& nWeights) {
  nWeights = exp(tWeights) / sum(exp(tWeights));
}

double computeLikelihood(vec& tWeights,
                         double lWeightsMax,
                         int n_particles){
  double tLikelihood = lWeightsMax + log(sum(exp(tWeights))) - log(n_particles);
  return tLikelihood;
}

void resampleParticles(mat& rParticles,
                       mat& particles,
                       vec& weights,
                       int n_particles){
  Rcpp::IntegerVector indices = Rcpp::sample(n_particles,n_particles, true, as<NumericVector>(wrap(weights)), false);
  for(int i = 0; i < n_particles ; i++) {
    rParticles.col(i) = particles.col(indices(i));
  }
}



// [[Rcpp::export(name=Particle_likelihood_SIR)]]
double Particle_likelihood_SIR( vec y,
                                vec init_pop,
                               arma::vec params,
                               int n_particles){
  
  double logLikelihood = 0;
  int n = accu(init_pop);
  int t = y.size();
  mat particles(3,n_particles);
  mat resampled_particles(3,n_particles);
  vec new_inf(n_particles);
  vec q(n_particles);
  vec mu(n_particles);
  vec sig(n_particles);
  for(int i = 0; i<n_particles; i++){
    resampled_particles.col(i) = init_pop;
  }
  vec  weights(n_particles);
  vec tWeights;
  vec nWeights(n_particles);
  double lWeightsMax;
  
  
  for(int i = 0; i < t; i++){
     // cout << "hi1" << endl;
    propagate_particles(resampled_particles, particles, new_inf, q, mu, sig, params, y(i), n_particles, n);
     // cout << "hi2" << endl;
     // cout << mu << endl;
    weight_particles(weights, y(i), mu, sig, n_particles, new_inf,params, q);
    // cout << weights << endl;
     // cout << "hi3" << endl;
    lWeightsMax = max(weights);
     // cout << "hi4" << endl;
    tWeights = weights - lWeightsMax;
    //  cout << "hi5" << endl;
    normaliseWeights(tWeights, nWeights);
     // cout << "hi6" << endl;
    double ll = computeLikelihood(tWeights, lWeightsMax, n_particles);
     // cout << "hi7" << endl;
    logLikelihood += ll;
    //  cout << "hi8" << endl;
    // cout << nWeights << endl;
    resampleParticles(resampled_particles, particles, nWeights,n_particles);
    // cout << "hi9" << endl;
  }
  return(logLikelihood);
}