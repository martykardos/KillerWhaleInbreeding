for(ii in 2:4) {
  df = expand.grid(alive = 1, 
                   year = seq(min(demog$year,na.rm=T),max(demog$year,na.rm=T)),
                   age = 11:40,
                   sexF1M2 = 1:2,
                   froh = c(-1.262,-1,0),
                   est = NA,
                   lo=NA, 
                   hi=NA)
  
  df$sexF1M2 = as.factor(df$sexF1M2)
  df$age2 = df$age^2
  formula = alive ~ s(year,k=7)+s(age,k=5) + sexF1M2+froh
  stan_dat = make_standata(formula, data = df, family="bernoulli")
  
  draws <- readRDS(paste0("survival/1stage_surv_",ii,"_bernoulli.rds"))
  
  b_indx = grep("b\\[", colnames(draws))
  bs_indx = grep("bs", colnames(draws))
  zs_1_indx = grep("zs_1_1", colnames(draws))
  zs_2_indx = grep("zs_2_1", colnames(draws))
  
  logit_f1 = matrix(0, nrow(draws), 1)
  logit_f2 = matrix(0, nrow(draws), 1)
  logit_f3 = matrix(0, nrow(draws), 1)
  logit_m1 = matrix(0, nrow(draws), 1)
  logit_m2 = matrix(0, nrow(draws), 1)
  logit_m3 = matrix(0, nrow(draws), 1)
  for(j in 11:40) {
    # mu = Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1;
    indx = which(df$year==2000 & df$sexF1M2==1 & df$age == j)
    
    logit_f1 = logit_f1 + draws[,b_indx] %*% stan_dat$X[indx[1],] + 
      draws[,bs_indx] %*% stan_dat$Xs[indx[1],] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx[1],] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx[1],]
    logit_f2 = logit_f2 + draws[,b_indx] %*% stan_dat$X[indx[2],] + 
      draws[,bs_indx] %*% stan_dat$Xs[indx[2],] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx[2],] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx[2],]
    logit_f3 = logit_f3 + draws[,b_indx] %*% stan_dat$X[indx[3],] + 
      draws[,bs_indx] %*% stan_dat$Xs[indx[3],] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx[3],] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx[3],]
    
    indx = which(df$year==2000 & df$sexF1M2==2 & df$age == j)
    
    logit_m1 = logit_m1 + draws[,b_indx] %*% stan_dat$X[indx[1],] + 
      draws[,bs_indx] %*% stan_dat$Xs[indx[1],] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx[1],] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx[1],]
    logit_m2 = logit_m2 + draws[,b_indx] %*% stan_dat$X[indx[2],] + 
      draws[,bs_indx] %*% stan_dat$Xs[indx[2],] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx[2],] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx[2],]
    logit_m3 = logit_m3 + draws[,b_indx] %*% stan_dat$X[indx[3],] + 
      draws[,bs_indx] %*% stan_dat$Xs[indx[3],] + 
      draws[,zs_1_indx] %*% stan_dat$Zs_1_1[indx[3],] + 
      draws[,zs_2_indx] %*% stan_dat$Zs_2_1[indx[3],]    

  }
  
  logit_f1 = logit_f1 / length(11:40)
  logit_f2 = logit_f2 / length(11:40)
  logit_f3 = logit_f3 / length(11:40)
  logit_m1 = logit_m1 / length(11:40)
  logit_m2 = logit_m2 / length(11:40)
  logit_m3 = logit_m3 / length(11:40)
  
  invlogit_f1 = plogis(as.numeric(unlist(logit_f1)))
  invlogit_f2 = plogis(as.numeric(unlist(logit_f2)))
  invlogit_f3 = plogis(as.numeric(unlist(logit_f3)))
  invlogit_m1 = plogis(as.numeric(unlist(logit_m1)))
  invlogit_m2 = plogis(as.numeric(unlist(logit_m2)))
  invlogit_m3 = plogis(as.numeric(unlist(logit_m3)))
  
  df1 = expand.grid("Sex" = "Female", 
                    "FROH" = c("0","ARKW","SRKW"),
                    "Mb" = c(NA,1,5,10)[ii],
                    "Mean"=0,"Low95"=0,"Hi95"=0)
  df1$Mean[1] = mean(invlogit_f1)
  df1$Low95[1] = quantile(invlogit_f1,0.025)
  df1$Hi95[1] = quantile(invlogit_f1,0.975)
  df1$Mean[2] = mean(invlogit_f2)
  df1$Low95[2] = quantile(invlogit_f2,0.025)
  df1$Hi95[2] = quantile(invlogit_f2,0.975) 
  df1$Mean[3] = mean(invlogit_f3)
  df1$Low95[3] = quantile(invlogit_f3,0.025)
  df1$Hi95[3] = quantile(invlogit_f3,0.975)    

  dfm = expand.grid("Sex" = "Male", 
                    "FROH" = c("0","ARKW","SRKW"),
                    "Mb" = c(NA,1,5,10)[ii],
                    "Mean"=0,"Low95"=0,"Hi95"=0)
  dfm$Mean[1] = mean(invlogit_m1)
  dfm$Low95[1] = quantile(invlogit_m1,0.025)
  dfm$Hi95[1] = quantile(invlogit_m1,0.975)
  dfm$Mean[2] = mean(invlogit_m2)
  dfm$Low95[2] = quantile(invlogit_m2,0.025)
  dfm$Hi95[2] = quantile(invlogit_m2,0.975) 
  dfm$Mean[3] = mean(invlogit_m3)
  dfm$Low95[3] = quantile(invlogit_m3,0.025)
  dfm$Hi95[3] = quantile(invlogit_m3,0.975)     
  
  if(ii== 2) df_all = rbind(df1,dfm)
  if(ii > 2) df_all = rbind(df_all, rbind(df1,dfm))

}

write.csv(df_all,paste0("average_survival",".csv"))