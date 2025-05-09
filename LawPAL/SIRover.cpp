// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;



// params is of the form (beta, rho, gamma,h)
// [[Rcpp::export(name=SIRsim_over)]]
List SIRsim_over(int T,
                  int n,
                  arma::vec init_dist,
                  arma::vec xi,
                  arma::vec params,
                  vec q,
                  int tau,
                  double h) {
  
  // define initial conditions
  int total_pop = n;
  mat state(3,T+1);
  vec init_pop = n*init_dist;
  vec eta = init_pop/total_pop;
  vec new_exp(T);
  vec new_inf(T);
  vec new_rem(T);
  vec obs(T);
  
  state.col(0) = init_pop;
  
  for(int i = 0; i < T; i++) {
    
    // simulate new movements between compartments
    new_inf(i) = R::rbinom(state(0,i),1-exp(-h*xi(i)*params(0)*eta(1)));
    new_rem(i) = R::rbinom(state(1,i),1-exp(-h*params(1)));
    
    // update states at next timestep

    state(0,i+1) = state(0,i) - new_inf(i);;
    state(1,i+1) = state(1,i) + new_inf(i) - new_rem(i);
    state(2,i+1) = state(2,i) +  new_rem(i);

    
    
    // update eta

    eta = state.col(i+1)/total_pop;
  
  }
  
  for(int t = 0; t<T; t++){
    for(int j = 0; j < 4 ; j++){
      obs(t) = R::rbinom(new_inf(t), q(t));
    }
  }
  
    List L = List::create(Named("state") = state, Named("obs") = obs);
  
  return(L);
}


