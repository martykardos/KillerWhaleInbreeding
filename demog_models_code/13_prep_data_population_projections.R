# use other subset of data
prior_pop = dplyr::filter(demog, sexF1M2==1, includeFec==1,
                          !is.na(Froh_min10Mb), !is.na(alive), 
                          !is.na(gave_birth),year>1975,
                          birth > 1960,
                          age >= 10, age <= 42)
prior_pop$froh = prior_pop$Froh_min10Mb

# use brms for this so it's straightforward to do predictions in 
# also incorporates error via posterior predictive dist
fit = brms::brm(gave_birth ~ s(year,k=7) + age + I(age^2)+froh,
                family="bernoulli",data =prior_pop,chains=4)

saveRDS(fit, "data/fec_mod.rds")
# prep survival data similarly
dat = dplyr::filter(demog, includeSurv==1,
                    !is.na(alive), birth>1960,
                    !is.na(Froh_min10Mb))
dat$froh = dat$Froh_min10Mb
formula = alive ~ s(year,k=7)+s(age,k=5)+sexF1M2 + froh

fit = brm(formula, family="bernoulli", data=dat,chains=4)
saveRDS(fit, "data/surv_mod.rds")

# bundle things into 1 workspace
surv_mod = readRDS("data/surv_mod.rds")
fec_mod = readRDS("data/fec_mod.rds")
save(surv_mod,fec_mod,file="demog_models_projection.Rdata")
