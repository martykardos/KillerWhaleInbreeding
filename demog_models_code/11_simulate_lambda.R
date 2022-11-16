df = expand.grid(alive = 1, 
                 year = seq(min(demog$year,na.rm=T),max(demog$year,na.rm=T)),
                 age = 1:100,
                 sexF1M2 = 1:2,
                 froh = 0,
                 gave_birth =1)
df$age2 = df$age^2
df$sexF1M2 = as.factor(df$sexF1M2)
# make surv dat
formula = alive ~ s(year,k=7)+s(age,k=5)+sexF1M2+froh
stan_dat = make_standata(formula, data = df, family="bernoulli")
# make fecundity dat
formula = gave_birth ~ s(year,k=7) + age + age2 + froh
birth_dat = make_standata(formula, data = df, family="bernoulli")

means = c(-1.307479,-2.481795,0.03234788)
sds = c(0.1772409,0.4826107,0.02563815)

d = read.csv("data/FswithLHandPed.csv") 
d = dplyr::select(d, 
                  id, Froh_min1Mb,
                  Froh_min5Mb,
                  Froh_min10Mb)

for(ii in 1:3) {
  
  df_out = data.frame("cov"=paste0(c(1,5,10)[ii],"Mb"), 
                      "Froh" = seq(-4,4,length.out=100),
                      "est"=NA, "lo"=NA, "hi"=NA)
  df_out$Froh_raw = sds[ii] * df_out$Froh + means[ii]
  if(ii < 3) df_out$Froh_raw = exp(df_out$Froh_raw)
  
  draws <- readRDS(paste0("survival/1stage_surv_",ii+1,"_bernoulli.rds"))
  draws2 <- readRDS(paste0("fecundity/1stage_fec_",ii+1,"_bernoulli.rds"))

  mat = matrix(0, 100, 100)
  lambda = 0
  for(j in 1:nrow(df_out)) {
    for(dd in 1:nrow(draws)) {
      
      for(a in 1:ncol(mat)) {
        # mu = Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 *   s_3_1;
        indx = which(df$sexF1M2=="1" & df$age==a & df$year == 2000)
        logit = draws[dd,2:3] %*% stan_dat$X[indx,1:2] + 
          draws[dd,5:6] %*% stan_dat$Xs[indx,] + 
          draws[dd,7:11] %*% stan_dat$Zs_1_1[indx,] + 
          draws[dd,13:15] %*% stan_dat$Zs_2_1[indx,] +
          draws[dd,4] * df_out$Froh[j]
        mat[a,a-1] = plogis(as.numeric(unlist(logit)))
        
        indx = which(df$sexF1M2==1 & df$age==a & df$year == 2000)
        logit = draws2[dd,2:5] %*% birth_dat$X[indx,] + 
          draws2[dd,6] * birth_dat$Xs[indx,] + 
          draws2[dd,7:11] %*% birth_dat$Zs_1_1[indx,]
        
        # do same thing with fecundity
        if(a %in% 10:42) mat[1,a] = plogis(as.numeric(unlist(logit)))    
      }
      lambda[dd] = eigen(mat)$values[1]
    }
    df_out$est[j] = mean(lambda)
    df_out$lo[j] = quantile(lambda,0.025)
    df_out$hi[j] = quantile(lambda,0.975)
  }
  
  if(ii==1) {
    df_all = df_out
  } else {
    df_all = rbind(df_all,df_out)
  }
}
saveRDS(df_all,"lambda_output.rds")
