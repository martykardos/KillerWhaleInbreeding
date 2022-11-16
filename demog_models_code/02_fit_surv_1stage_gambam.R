dat$sexF1M2 = as.factor(dat$sexF1M2)

#mod = cmdstan_model("survival/age_based_survival_bernoulli.stan")
gam_est = matrix(0, 3, 4)
bam_est = matrix(0, 3, 4)

for(ii in c(2,4)) {
  
  if(ii==1) dat$froh = rep(0, nrow(dat))
  if(ii==2) dat$froh = dat$Froh_min1Mb
  if(ii==3) dat$froh = dat$Froh_min5Mb
  if(ii==4) dat$froh = dat$Froh_min10Mb
  dat$age2 = dat$age^2
  formula = alive ~ s(year,k=7)+s(age,k=5)+sexF1M2 + froh
  #formula = alive ~ s(year,k=7)+age + age2 +sexF1M2 + froh
  g <- gam(formula, data=dat)
  #formula = alive ~ s(year,k=7)+age*sexF1M2 + age2*sexF1M2 + froh
  gam_est[ii-1,] = c(summary(g)$p.coeff[3], summary(g)$se[3], summary(g)$p.t[3], summary(g)$p.pv[3])
  
  fit = brm(formula, family="bernoulli", data=dat,chains=3,iter=2000)
  bfroh = rstan::extract(fit$fit, "b_froh")$b_froh
  bam_est[ii-1,] = c(mean(bfroh), quantile(bfroh, 0.025), quantile(bfroh, 0.975), length(which(bfroh < 0))/length(bfroh))
}

colnames(gam_est) = c("Estimate","Std. Error","t value","Pr(>|t|)")
colnames(bam_est) = c("mean","lower95","upper95","Pr_less_0")

saveRDS(gam_est,file="table_S9.rds")
saveRDS(bam_est,file="table_S11.rds")