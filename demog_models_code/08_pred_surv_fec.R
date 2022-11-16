for(ii in 2:4) {
  n_draws = 10000
  df = expand.grid(alive = 1, 
                   year = seq(min(demog$year,na.rm=T),max(demog$year,na.rm=T)),
                   age = seq(0,50,by=5),
                   sexF1M2 = 1:2,
                   froh = seq(-4,4,length.out=21),
                   est = NA,
                   lo=NA, 
                   hi=NA)
  df = rbind(df, expand.grid(alive = 1, 
                             year = 2000,
                             age = 20,
                             sexF1M2 = 1:2,
                             froh = c(-1.262,-1,0),
                             est = NA,
                             lo=NA, 
                             hi=NA))
  
  df$sexF1M2 = as.factor(df$sexF1M2)
  df$age2 = df$age^2
  formula = alive ~ s(year,k=7)+s(age,k=5) + sexF1M2+froh
  stan_dat = make_standata(formula, data = df, family="bernoulli")
  
  draws <- readRDS(paste0("survival/1stage_surv_",ii,"_bernoulli.rds"))
  
  b_indx = grep("b\\[", colnames(draws))
  bs_indx = grep("bs", colnames(draws))
  zs_1_indx = grep("zs_1_1", colnames(draws))
  zs_2_indx = grep("zs_2_1", colnames(draws))
  
  df$row <- 1:nrow(df)
  # save 2000 mcmc draws 
  all_draws <- data.frame(row=NA, pred_logit = NA, 
                          pred_invlogit = NA, draws=NA)
  all_draws_20 <- data.frame(row=NA, pred_logit = NA, 
                             pred_invlogit = NA, draws=NA)
  for(j in 1:nrow(df)) {
    print(j)
    # mu = Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1;

    logit = draws[,b_indx] %*% stan_dat$X[j,] + 
      draws[,bs_indx] %*% stan_dat$Xs[j,] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[j,] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[j,]
    
    invlogit = plogis(as.numeric(unlist(logit)))
    df$est[j] = mean(invlogit)
    df$lo[j] = quantile(invlogit,0.025)
    df$hi[j] = quantile(invlogit,0.975)
    
    # save draws
    save_draws <- data.frame(draws=1:n_draws,row=j, pred_logit = logit[1:n_draws,1], 
                              pred_invlogit = invlogit[1:n_draws])
    if(df$froh[j]==0){
      all_draws <- rbind(all_draws, save_draws)
    }
    if(df$age[j] == 20) {
      all_draws_20 <- rbind(all_draws_20, save_draws)
    }
  }
  saveRDS(df, paste0("survival_output_",ii,".rds"))
  saveRDS(all_draws, paste0("survival_output_draws_",ii,".rds"))
  saveRDS(all_draws_20, paste0("survival_output_draws_20_",ii,".rds"))
  
  # do same thing with fecundity models
  # Make plots of survival and fecundity @ age
  df = expand.grid(gave_birth = 1, 
                   year = seq(min(demog$year,na.rm=T),max(demog$year,na.rm=T)),
                   age = 10:42,
                   sexF1M2 = 1,
                   froh = seq(-4,4,length.out=21),
                   est = NA,
                   lo=NA, 
                   hi=NA)
  
  formula = gave_birth ~ s(year,k=7) + age + age2 + froh
  df$age2 = df$age^2
  stan_dat = make_standata(formula, data = df, family="bernoulli")
  
  draws <- readRDS(paste0("fecundity/1stage_fec_",ii,"_bernoulli.rds"))
  
  b_indx = grep("b\\[", colnames(draws))
  bs_indx = grep("bs", colnames(draws))
  zs_1_indx = grep("zs_1_1", colnames(draws))
  
  all_draws <- data.frame(row=NA, pred_logit = NA, 
                           pred_invlogit = NA, draws=NA)
  all_draws_20 <- data.frame(row=NA, pred_logit = NA, 
                          pred_invlogit = NA, draws=NA)
  for(j in 1:nrow(df)) {
    # mu = Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1;
    
    logit = draws[,b_indx] %*% stan_dat$X[j,] + 
      draws[,bs_indx] * stan_dat$Xs[j,] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[j,]
    
    invlogit = plogis(as.numeric(unlist(logit)))
    df$est[j] = mean(invlogit)
    df$lo[j] = quantile(invlogit,0.025)
    df$hi[j] = quantile(invlogit,0.975)
    
    # save draws
    save_draws <- data.frame(draws=1:n_draws,row=j, pred_logit = logit[1:n_draws,1], 
                             pred_invlogit = invlogit[1:n_draws])
    if(df$froh[j]==0){
      all_draws <- rbind(all_draws, save_draws)
    }
    if(df$age[j] == 20) {
      all_draws_20 <- rbind(all_draws_20, save_draws)
    }
  }
  saveRDS(df, paste0("fecundity_output_",ii,".rds"))
  saveRDS(all_draws, paste0("fecundity_output_draws_",ii,".rds"))
  saveRDS(all_draws_20, paste0("fecundity_output_draws_20_",ii,".rds"))
}