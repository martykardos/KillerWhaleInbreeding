base_fit <- readRDS("fecundity/fit_2stage_prior_bernoulli.rds")
# summarize posteriors for new prior
pars = base_fit[,2:4]
mu = c(apply(pars, 2, mean),0)
Sigma = rbind(cbind(cov(pars),0),0)
Sigma[4,4] = 1
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
  # 47
  #fit <- brm(formula, family="bernoulli", 
  #           data=prior_pop, chains=1, iter=10)
  prior_pop$age2 = prior_pop$age^2
  stan_dat = make_standata(formula, data = prior_pop, family="bernoulli")
  
  # fixed effects 
  indx = grep("b\\[", colnames(base_fit))
  stan_dat$prior_mean = c(apply(base_fit[,indx],2,mean),0)
  stan_dat$prior_cov = cov(base_fit[,indx])
  stan_dat$prior_cov = cbind(rbind(stan_dat$prior_cov,0),0)
  stan_dat$prior_cov[4,4]=1
  indx = grep("zs_1_1", colnames(base_fit))
  stan_dat$Zs_mean = apply(base_fit[,indx], 2,mean)
  stan_dat$Zs_cov = cov(base_fit[,indx])
  indx = grep("bs", colnames(base_fit))
  stan_dat$bs_mean = mean(base_fit[,indx])
  stan_dat$bs_cov = sd(base_fit[,indx])
  
  fit = mod$sample(seed=123, chains=4, data=stan_dat, 
                   adapt_delta=.999,
                   iter_sampling = 5000,
                   iter_warmup = 5000,
                   parallel_chains=4,
                   max_treedepth=20)
  draws <- as.matrix(as_draws_df(fit$draws()))
  saveRDS(draws, file=paste0("fecundity/2stage_fec_",ii,"_bernoulli.rds"))
}