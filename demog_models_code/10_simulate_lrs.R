df = expand.grid(alive = 1, 
                 year = seq(min(demog$year,na.rm=T),max(demog$year,na.rm=T)),
                 age = 1:45,
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

means = c(0.2748935,0.09359155,0.03234788)#c(-1.307479,-2.481795,0.03234788)
sds = c(0.05170584,0.04568827,0.02563815)# c(0.1772409,0.4826107,0.02563815)

d = read.csv("data/FswithLHandPed.csv") 
d = dplyr::select(d, 
                  id, Froh_min1Mb,
                  Froh_min5Mb,
                  Froh_min10Mb)

for(ii in 1:3) {
  
  df_out = data.frame("cov"=paste0(c(1,5,10)[ii],"Mb"), 
                      "Froh" = seq(-4,4,length.out=100),
                      "lrs"=NA, "lrs_lo"=NA, "lrs_hi"=NA)
  df_out$Froh_raw = sds[ii] * df_out$Froh + means[ii]
  #if(ii < 3) df_out$Froh_raw = exp(df_out$Froh_raw)
  
  draws <- readRDS(paste0("survival/1stage_surv_",ii+1,"_bernoulli.rds"))
  draws2 <- readRDS(paste0("fecundity/1stage_fec_",ii+1,"_bernoulli.rds"))
  m = matrix(0,nrow(draws), max(df$age))
  f = matrix(0,nrow(draws), max(df$age))
  
  b_indx = grep("b\\[", colnames(draws))
  bs_indx = grep("bs", colnames(draws))
  zs_1_indx = grep("zs_1_1", colnames(draws))
  zs_2_indx = grep("zs_2_1", colnames(draws))
  b_indx2 = grep("b\\[", colnames(draws2))
  bs_indx2 = grep("bs", colnames(draws2))
  zs_1_indx2 = grep("zs_1_1", colnames(draws2))

  for(j in 1:nrow(df_out)) {
    for(a in 1:ncol(m)) {
      # mu = Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 *   s_3_1;
      indx = which(df$sexF1M2==1 & df$age==a & df$year == 2000)
      logit = draws[,b_indx[1:(length(b_indx)-1)]] %*% stan_dat$X[indx,1:2] + 
        draws[,bs_indx] %*% stan_dat$Xs[indx,] + 
        draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx,] + 
        draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx,] +
        draws[,b_indx[length(b_indx)]] * df_out$Froh[j]
      m[,a] = plogis(as.numeric(unlist(logit)))
      
      indx = which(df$sexF1M2==1 & df$age==a & df$year == 2000)
      logit = draws2[,b_indx2] %*% birth_dat$X[indx,] + 
        draws2[,bs_indx2] * birth_dat$Xs[indx,] + 
        draws2[,zs_1_indx2] %*% birth_dat$Zs_1_1[indx,]
      
      # do same thing with fecundity
      f[,a] = plogis(as.numeric(unlist(logit)))    
    }
    f[,1:9] = 0 # youngest animals can't reproduce
    f[,43:45] = 0
    # simulate LRS
    ncolm = ncol(m)
    lrs = matrix(0, nrow(m), ncol(m))
    for(i in 1:nrow(lrs)) {
      # get indices for animals that die
      lrs[i,] = m[cbind(sample(1:nrow(m),size=ncolm,replace=T),1:ncolm)]
      indx = which(runif(ncolm) < 1-lrs[i,])
      if(length(indx) != 0) {
        lrs[i,1:(indx[1]-1)] = 1
        lrs[i,indx[1]:ncolm] = 0
      } else {
        lrs[i,] = 1
      }
      fs = f[cbind(sample(1:nrow(m),size=ncolm,replace=T),1:ncolm)]
      fs[1] = ifelse(runif(1) < fs[1], 1, 0)
      for(iii in 2:ncolm) {
        if(fs[iii-1]==0) {
          fs[iii] = ifelse(runif(1) < fs[iii], 1, 0)
        } else {
          fs[iii] = 0
        }
      }
      lrs[i,] = lrs[i,] * fs
    }
    
    df_out$lrs[j] = mean(apply(lrs,1,sum))
    df_out$lrs_lo[j] = quantile(apply(lrs,1,sum),0.025)
    df_out$lrs_hi[j] = quantile(apply(lrs,1,sum),0.975)
  }
  
  if(ii==1) {
    df_all = df_out
  } else {
    df_all = rbind(df_all,df_out)
  }
}
saveRDS(df_all,"sim_output.rds")
