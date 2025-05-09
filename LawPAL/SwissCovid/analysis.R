setwd("~/Documents/R/disease_transmission_workflow-main")
library(tidybayes)
library(gridExtra)
library(tidyverse)
library(rstan)
library(bayesplot)
df_swiss <- read_csv("data/swiss_agg_data.csv")


df_swiss %>% 
  ggplot() + 
  geom_bar(mapping = aes(x = date, y = report_dt),  stat = "identity") +
  labs(y="Number of reported cases")

cases_raw <- df_swiss$report_dt

df_swiss %>% 
  ggplot() + 
  geom_bar(mapping = aes(x = date, y = report_dt), fill = c_mid, color = c_dark, stat = "identity") +
  labs(y="Number of reported cases")

# start from index case
cases <- cases_raw[24:132]
ts.plot(cases)
data <- list(y = cases)

fit_PAL <- sampling(model_PAL, 
                    data, 
                    iter=50000,
                    seed = 5,
                    chains = 4)

traceplot(fit_PAL, pars = c("gamma", "beta","rho","d","alpha","b","q","i0","e0"), chains = 1)
mean(samples_eq$beta/samples_eq$gamma)
samples_eq <- extract(fit_PAL)


summary <- summary(fit_PAL)$summary
means <- summary[, "mean"]
sds <- summary[, "sd"]




model_LawPAL <- stan_model("~/Documents/R/LawPAL/mcmc/Stan/SwissCovidLawPAL.stan")

fit_LawPAL <- sampling(model_LawPAL, 
                       data, 
                       iter=50000,
                       seed = 4,
                       chains = 4
                    )

traceplot(fit_LawPAL, pars = c("gamma", "beta","rho","alpha","d","b","q","sigma_q","i0","e0"), n_cols = 2)

summary <- summary(fit_LawPAL)$summary
means <- summary[, "mean"]
sds <- summary[, "sd"]


samples <- extract(fit_LawPAL)


