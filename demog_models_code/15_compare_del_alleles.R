dat$sexF1M2 = as.factor(dat$sexF1M2)

mod = cmdstan_model("survival/age_based_survival_bernoulli.stan")
fits = list()
for(ii in 2:4) {
  
  if(ii==1) dat$froh = dat$homs#rep(0, nrow(dat))
  if(ii==2) dat$froh = dat$Froh_min1Mb
  if(ii==3) dat$froh = dat$Froh_min5Mb
  if(ii==4) dat$froh = dat$Froh_min10Mb
  dat$age2 = dat$age^2
  #formula = alive ~ s(year,k=7)+s(age,k=5)+sexF1M2 + froh
  formula = alive ~ s(year,k=7)+age + age2 +sexF1M2 + froh
  #formula = alive ~ s(year,k=7)+age*sexF1M2 + age2*sexF1M2 + froh
  
  fit = brm(formula, family="bernoulli", data=dat,chains=3,iter=3000)
  fits[ii] = loo(fit)
  # newdf = data.frame(age=1:60,sexF1M2=1,froh=0,year=2000)
  # newdf$age2 = newdf$age^2
  # pred = predict(fit,newdf)
  #stan_dat = make_standata(formula, data = dat, family="bernoulli")
  
  #fit[[i]] = mod$sample(seed=123, chains=4, data=stan_dat, 
  #                 adapt_delta=.999,
  #                 iter_sampling = 5000,
  #                 iter_warmup = 5000,
  #                 parallel_chains=4,
  #                 max_treedepth=20)
  
  #draws <- as.matrix(as_draws_df(fit$draws()))
  
  #saveRDS(draws, file=paste0("survival/1stage_surv_",ii,"_bernoulli.rds"))
}
