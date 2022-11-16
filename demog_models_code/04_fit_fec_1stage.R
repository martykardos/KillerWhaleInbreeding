
# use other subset of data
prior_pop = dplyr::filter(demog, sexF1M2==1, includeFec==1,
                          !is.na(Froh_min10Mb), !is.na(alive), 
                          !is.na(gave_birth),year>1975,
                          birth > 1960,
                          age >= 10, age <= 42)
# 43
mod = cmdstan_model("fecundity/2stage_fec_prior_bernoulli.stan")

for(ii in 2:4) {
  if(ii==1) prior_pop$froh = rep(0, nrow(X))
  if(ii==2) prior_pop$froh = prior_pop$Froh_min1Mb
  if(ii==3) prior_pop$froh = prior_pop$Froh_min5Mb
  if(ii==4) prior_pop$froh = prior_pop$Froh_min10Mb
  
  formula = gave_birth ~ s(year,k=7) + age + age2 + froh
  
  prior_pop$age2 = prior_pop$age^2
  stan_dat = make_standata(formula, data = prior_pop, family="bernoulli")
  
  stan_dat$prior_mean = rep(0,4)
  stan_dat$prior_cov = diag(4)
  stan_dat$Zs_mean = rep(0,5)
  stan_dat$Zs_cov = diag(5)
  stan_dat$bs_mean = 0
  stan_dat$bs_cov = 1
  
  fit = mod$sample(seed=123, chains=4, data=stan_dat, 
                   adapt_delta=.999,
                   iter_sampling = 5000,
                   iter_warmup = 5000,
                   parallel_chains=4,
                   max_treedepth=20)
  draws <- as.matrix(as_draws_df(fit$draws()))
  saveRDS(draws, file=paste0("fecundity/1stage_fec_",ii,"_bernoulli.rds"))
}