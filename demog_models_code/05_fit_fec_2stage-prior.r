set.seed(2021)
demog$age2 = demog$age^2
base_pop = dplyr::filter(demog, sexF1M2==1, includeFec==1,
                         is.na(Froh_min10Mb), !is.na(alive), 
                         !is.na(gave_birth), year>1975, 
                         birth > 1960,
                         age >= 10, age <= 42)
formula = gave_birth ~ s(year,k=7) + age + age2
# 19
mod = cmdstan_model("fecundity/2stage_fec_prior_bernoulli2.stan")
stan_dat = make_standata(formula, data = base_pop, family="bernoulli")
stan_dat$prior_mean = rep(0, 3)
stan_dat$prior_cov = diag(3)
stan_dat$Zs_mean = rep(0, ncol(stan_dat$Zs_1_1))
stan_dat$Zs_cov = cov(stan_dat$Zs_1_1)
stan_dat$bs_mean = 0
stan_dat$bs_cov = 1

fit = mod$sample(seed=123, chains=4, data=stan_dat, 
                 adapt_delta=.999,
                 iter_sampling = 5000,
                 iter_warmup = 5000,
                 parallel_chains=4,
                 max_treedepth=20)

draws <- as.matrix(as_draws_df(fit$draws()))

saveRDS(draws, "fecundity/fit_2stage_prior_bernoulli.rds")