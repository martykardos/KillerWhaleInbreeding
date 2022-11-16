df = expand.grid(alive = 1, 
                 year = seq(min(demog$year,na.rm=T),max(demog$year,na.rm=T)),
                 age = 1:40,
                 sexF1M2 = 1:2,
                 froh = 0)
df$sexF1M2 = as.factor(df$sexF1M2)
formula = alive ~ s(year,k=7)+s(age,k=5)+sexF1M2+froh

stan_dat = make_standata(formula, data = df, family="bernoulli")

means = c(0.2748935,0.09359155,0.03234788)#c(-1.307479,-2.481795,0.03234788)
sds = c(0.05170584,0.04568827,0.02563815)# c(0.1772409,0.4826107,0.02563815)

d = read.csv("data/FswithLHandPed.csv") 
d = dplyr::select(d, 
                  id, Froh_min1Mb,
                  Froh_min5Mb,
                  Froh_min10Mb)

for(ii in 1:3) {
  
  df_out = data.frame("cov"=paste0(c(1,5,10)[ii],"Mb"), 
                      "Froh" = seq(-4,4,length.out=50),
                      "est"=NA, "lo"=NA, "hi"=NA, "med"=NA,
                      "med_age"=NA,"med_lo"=NA,"med_hi"=NA)
  df_out$Froh_raw = sds[ii] * df_out$Froh + means[ii]
  #if(ii < 3) df_out$Froh_raw = exp(df_out$Froh_raw)
  
  df_out_fem = df_out  
  draws <- readRDS(paste0("survival/1stage_surv_",ii+1,"_bernoulli.rds"))
  
  b_indx = grep("b\\[", colnames(draws))
  bs_indx = grep("bs", colnames(draws))
  zs_1_indx = grep("zs_1_1", colnames(draws))
  zs_2_indx = grep("zs_2_1", colnames(draws))
  
  m = matrix(0,nrow(draws), max(df$age))
  fem = matrix(0,nrow(draws), max(df$age))
  for(j in 1:nrow(df_out)) {
    
    for(a in 1:ncol(m)) {
      # mu = Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 *   s_3_1;
      indx = which(df$sexF1M2==2 & df$age==a & df$year == 2000)
      logit = draws[,b_indx[1:(length(b_indx)-1)]] %*% stan_dat$X[indx,1:2] + 
        draws[,bs_indx] %*% stan_dat$Xs[indx,] + 
        draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx,] + 
        draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx,] +
        draws[,b_indx[length(b_indx)]] * df_out$Froh[j]
      
      m[,a] = plogis(as.numeric(unlist(logit)))
      
      indx = which(df$sexF1M2==1 & df$age==a & df$year == 2000)
      logit = draws[,b_indx[1:(length(b_indx)-1)]] %*% stan_dat$X[indx,1:2] + 
        draws[,bs_indx] %*% stan_dat$Xs[indx,] + 
        draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx,] + 
        draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx,] +
        draws[,b_indx[length(b_indx)]] * df_out$Froh[j]
      
      fem[,a] = plogis(as.numeric(unlist(logit)))
    }
    
    invlogit_m <- apply(m[,1:40],1,prod)
    df_out$est[j] = mean(invlogit_m)
    df_out$lo[j] = quantile(invlogit_m,0.025)
    df_out$hi[j] = quantile(invlogit_m,0.975)
    df_out$med[j] = quantile(invlogit_m,0.5)
    
    invlogit_f <- apply(fem[,1:40],1,prod)
    df_out_fem$est[j] = mean(invlogit_f)
    df_out_fem$lo[j] = quantile(invlogit_f,0.025)
    df_out_fem$hi[j] = quantile(invlogit_f,0.975)
    df_out_fem$med[j] = quantile(invlogit_f,0.5)  
    
    draws_m = data.frame(pred_invlogit = invlogit_m, row = j, ii = ii)
    draws_m$draw = seq(1,nrow(draws_m))
    draws_f = data.frame(pred_invlogit = invlogit_f, row = j, ii=ii)
    draws_f$draw = seq(1,nrow(draws_f))   
    if(j==1) {
      all_draws_m = draws_m
      all_draws_f = draws_f
    } else {
      all_draws_m = rbind(all_draws_m, draws_m)
      all_draws_f = rbind(all_draws_f, draws_f)
    }
    
  }
  
  if(ii==1) {
    df_all = df_out
    df_all_fem = df_out_fem
    df_all_draws_m = all_draws_m
    df_all_draws_f = all_draws_f
  } else {
    df_all = rbind(df_all,df_out)
    df_all_fem = rbind(df_all_fem, df_out_fem)
    df_all_draws_m = rbind(df_all_draws_m, all_draws_m)
    df_all_draws_f = rbind(df_all_draws_f, all_draws_f)
  }
}

saveRDS(df_all,"cumulative_survival_males.rds")
saveRDS(df_all_fem,"cumulative_survival_females.rds")

saveRDS(df_all_draws_m,"cumulative_survival_alldraws_males.rds")
saveRDS(df_all_draws_f,"cumulative_survival_alldraws_females.rds")