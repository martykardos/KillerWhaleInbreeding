
mod = cmdstan_model("survival/age_based_survival_bernoulli.stan")
mod2 = cmdstan_model("survival/age_based_survival_bernoulli2.stan")

formula = alive ~ s(year,k=7)+s(age,k=5)+sexF1M2
stan_dat = make_standata(formula, data = prior_dat, family="bernoulli")
fit = mod$sample(seed=123, chains=4, data=stan_dat, 
                 adapt_delta=.999,
                 iter_sampling = 5000,
                 iter_warmup = 5000,
                 parallel_chains=4,
                 max_treedepth=20)  
prior_draws <- as.matrix(as_draws_df(fit$draws()))
saveRDS(prior_draws, file=paste0("survival/surv_priors.rds"))

for(ii in 2:4) {
  
  if(ii==1) dat$froh = rep(0, nrow(dat))
  if(ii==2) dat$froh = dat$Froh_min1Mb
  if(ii==3) dat$froh = dat$Froh_min5Mb
  if(ii==4) dat$froh = dat$Froh_min10Mb
  formula = alive ~ s(year,k=7)+s(age,k=5) + sexF1M2 + froh
  stan_dat = make_standata(formula, data = dat, family="bernoulli")
  # add prior info
  indx = grep("b\\[", colnames(prior_draws))
  stan_dat$b_mean = c(apply(prior_draws[,indx],2,mean),0)
  stan_dat$b_cov = cov(prior_draws[,indx]) 
  stan_dat$b_cov = rbind(cbind(stan_dat$b_cov,0),0)
  stan_dat$b_cov[ncol(stan_dat$b_cov),ncol(stan_dat$b_cov)]=1
  
  # year effects
  stan_dat$Zs1_mean = rep(0,5)#apply(prior_draws[,8:12],2,mean)
  stan_dat$Zs1_cov = diag(5)#cov(prior_draws[,8:12])
  
  indx = grep("bs", colnames(prior_draws))
  stan_dat$bs_mean = apply(prior_draws[,indx],2,mean)
  stan_dat$bs_cov = cov(prior_draws[,indx])  
  # age effect
  indx = grep("zs_2_1", colnames(prior_draws))
  stan_dat$Zs2_mean = apply(prior_draws[,indx],2,mean)
  stan_dat$Zs2_cov = cov(prior_draws[,indx])
  #stan_dat$sd1_mean = mean(log(prior_draws[,13]))
  #stan_dat$sd1_cov = sd(log(prior_draws[,13]))
  
  fit = mod2$sample(seed=123, chains=4, data=stan_dat, 
                    adapt_delta=.999,
                    iter_sampling = 5000,
                    iter_warmup = 5000,
                    parallel_chains=4,
                    max_treedepth=20)
  
  draws <- as.matrix(as_draws_df(fit$draws()))
  
  saveRDS(draws, file=paste0("survival/2stage_surv_",ii,"_bernoulli.rds"))
}
