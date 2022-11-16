
# use other subset of data
prior_pop = dplyr::filter(demog, sexF1M2==1, includeFec==1,
                          !is.na(Froh_min10Mb), !is.na(alive), 
                          !is.na(gave_birth),year>1975,
                          birth > 1960,
                          age >= 10, age <= 42)

for(ii in 2:4) {
  if(ii==1) prior_pop$froh = rep(0, nrow(X))
  if(ii==2) prior_pop$froh = prior_pop$Froh_min1Mb
  if(ii==3) prior_pop$froh = prior_pop$Froh_min5Mb
  if(ii==4) prior_pop$froh = prior_pop$Froh_min10Mb
  
  formula = gave_birth ~ s(year,k=7) + age + age2 + froh
  prior_pop$age2 = prior_pop$age^2
 
  g <- gam(formula, data=prior_pop)
  #formula = alive ~ s(year,k=7)+age*sexF1M2 + age2*sexF1M2 + froh
  gam_est[ii-1,] = c(summary(g)$p.coeff[4], summary(g)$se[4], summary(g)$p.t[4], summary(g)$p.pv[4])
  
  fit = brm(formula, family="bernoulli", data=prior_pop,chains=3,iter=2000)
  bfroh = rstan::extract(fit$fit, "b_froh")$b_froh
  bam_est[ii-1,] = c(mean(bfroh), quantile(bfroh, 0.025), quantile(bfroh, 0.975), length(which(bfroh < 0))/length(bfroh))
}


# colnames(gam_est) = c("Estimate","Std. Error","t value","Pr(>|t|)")
# colnames(bam_est) = c("mean","lower95","upper95","Pr_less_0")
# 
# saveRDS(gam_est,file="table_S9.rds")
# saveRDS(bam_est,file="table_S11.rds")