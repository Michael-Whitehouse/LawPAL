# LawPAL

Implementations associated with the paper "Accelerated inference for stochastic compartmental models with over-dispersed partial observations"


An assumed density approximate likelihood is derived for a class of partially observed
stochastic compartmental models which permit observational over-dispersion. This is
achieved by treating time-varying reporting probabilities as latent variables and integrating them out using Laplace approximations within Poisson Approximate Likelihoods (LawPAL), resulting in a fast deterministic approximation to the marginal likelihood and filtering
distributions. We derive an asymptotically exact filtering result in the large population
regime, demonstrating the approximation’s ability to recover latent disease states and reporting probabilities. Through simulations we: 1) demonstrate favorable behavior of the
maximum approximate likelihood estimator in the large population and time horizon regime
in terms of ground truth recovery; 2) demonstrate order of magnitude computational speed
gains over a sequential Monte Carlo likelihood based approach, and explore the statistical
compromises our approximation implicitly makes. We conclude by embedding our methodol-
ogy within the probabilistic programming language Stan for automated Bayesian inference to
develop a model of practical interest using data from the Covid-19 outbreak in Switzerland.
