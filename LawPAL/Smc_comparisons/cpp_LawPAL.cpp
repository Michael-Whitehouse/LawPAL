// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export(name=SIR_approx_lik_over)]]
double SIR_approx_lik(vec y,
                      vec init_dist,
                      vec params){
  
  int t = y.size();
  int n = sum(init_dist);
  double l;
  mat lambda_(3,t);
  mat lambda(3,t+1);
  mat K(3,3, fill::zeros);
  mat transitions(3,3, fill::zeros);
  vec w(t, fill::zeros);
  lambda.col(0) = init_dist;
  vec q_est(t, fill::zeros);
  vec q_var(t, fill::zeros);
  colvec ones(3, fill::ones);
  K(1,2) = 1 - exp(-params[1]);
  K(1,1) = 1 - K(1,2);
  K(2,2) = 1;
  
  for(int i =0; i < t; i++){
    K(0,0) = exp(-params(0)*lambda(1,i)/n);
    K(0,1) = 1 - K(0,0);

    transitions = ((ones*lambda.col(i).t()).t())%K;
    // update step

    l = transitions(0,1);
    q_est(i) = (0.5)*((params(2) -   l*params(3)*params(3)) + sqrt((l*params(3)*params(3) - params(2))*(l*params(3)*params(3) - params(2)) + 4*y(i)*params(3)*params(3)));

    q_var(i) = sqrt(1/(y(i)/(q_est(i)*q_est(i)) + 1/(params(3)*params(3))));
    
    transitions(0,1) = y(i) + (1-q_est(i))*l;

    lambda.col(i+1) = sum(transitions,0).t();
    
    w(i) = R::dpois(y(i), l*q_est(i) ,1) + d_truncnorm(q_est(i),params(2),params(3),0,1,1) + log(q_var(i)*sqrt(2*3.14159));
  }
  double lik = accu(w);
  
  return(lik);
}