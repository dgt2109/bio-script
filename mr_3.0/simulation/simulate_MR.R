
num_indiv = 10000
num_instr = 500
sd_u = 0.01 
sd_Ex = 0.01 

rate = 75

if(TEST_SIM == T) {
  n_sim = 20 # usu 100
  beta_X = 0.1
  beta_U = 0.1
  sd_Ey = 0.6
}
#OVERLAP = 0 # set a % overlap by mixing u1 and u2

# record results: obs, ivw_bonf, egger_bonf, ivw_fdr, egger_fdr
positives <- rep(0,5)
coverage <- rep(0,5)
regression_type <- c("obs","ivw_bonf","ivw_fdr","egger_bonf","egger_fdr")
names(positives) <- regression_type
names(coverage) <- regression_type
betas = matrix(0,n_sim,5)
colnames(betas) <- regression_type

print(paste("sd_Ey:",sd_Ey,"beta_U:",beta_U,"beta_X:",beta_X,"x ",num_indiv,",",num_instr,sep=' '))
for (i in 1:n_sim) {
#  print(i)
  maf = runif(num_instr, 0.05, 0.5) #allele freq
  a = rexp(num_instr,rate) # fit to HDL distribution
  
  g=matrix(nrow=num_instr,ncol=num_indiv)
  for (j in 1:num_instr) {
    g[j,] = rbinom(num_indiv, 2, maf[j])
  }
  
  u1 = rnorm(num_indiv,0,sd_u)
  u2 = rnorm(num_indiv,0,sd_u)
  #    u = c(u1[0:round(OVERLAP*num_indiv)], # mixing code
  #                 u2[0:round((1-OVERLAP)*num_indiv)])
  
  # calculate exposure and outcome
  x = vector(mode="double",length=num_indiv)
  x = u1 + rnorm(num_indiv,0,sd_Ex) + colSums(a * g)
  
  y_obs = beta_X*x + beta_U*u1 + rnorm(num_indiv,0,sd_Ey)
  y = beta_X*x + beta_U*u2 + rnorm(num_indiv,0,sd_Ey)
  
  # generate regression results per instrument
  rnames = c("rf","p_rf","rf_se","do","do_se")
  for (var in rnames) {
    eval(parse(text=paste0(var,"<-vector(mode=\"double\",length=num_instr)")))
  }
  
  for (j in 1:num_instr) {
    res_rf = summary(lm(x~g[j,]))
    rf[j] = res_rf$coefficients[2] 
    rf_se[j] = res_rf$coefficients[4]
    p_rf[j] = res_rf$coefficients[8]
    res_do = summary(lm(y~g[j,]))
    do[j] = res_do$coefficients[2]
    do_se[j] = res_do$coefficients[4]
  }
  
  # compute exposure:outcome regressions
  res_obs = summary(lm(y_obs ~ x+0))
  betas[i,1] = res_obs$coef[1]
  
  if(res_obs$coef[4] < 0.05 && betas[i,1] > 0) { positives[1]=positives[1]+1 }
  if(beta_X > betas[i,1] - qt(0.975,df=length(x)-1)*res_obs$coef[1,2] &&
     beta_X < betas[i,1] + qt(0.975,df=length(x)-1)*res_obs$coef[1,2]) {
    coverage[1]=coverage[1]+1
  }
  
  bonf_sig = which(p_rf < 0.05/num_instr)
  sr_bonf = rf[bonf_sig]
  so_bonf = do[bonf_sig]
  err_bonf = do_se[bonf_sig]
  bonf_ivw = summary(lm(so_bonf ~ sr_bonf+0,weights=err_bonf^-2))
  bonf_egger = summary(lm(so_bonf ~ sr_bonf,weights=err_bonf^-2))
  betas[i,2] = bonf_ivw$coef[1]
  betas[i,4] = bonf_egger$coef[2,1]
  
  if(bonf_ivw$coef[4] < 0.05 && betas[i,2] > 0) { positives[2]=positives[2]+1 }
  if(bonf_egger$coef[2,4] < 0.05 && betas[i,4] > 0) { positives[4]=positives[4]+1 }
  if(beta_X > betas[i,2] - qt(0.975,df=length(bonf_sig)-1)*bonf_ivw$coef[1,2] &&
     beta_X < betas[i,2] + qt(0.975,df=length(bonf_sig)-1)*bonf_ivw$coef[1,2]) {
    coverage[2]=coverage[2]+1
  }
  if(beta_X > betas[i,4] - qt(0.975,df=length(bonf_sig)-1)*bonf_egger$coef[2,2] &&
     beta_X < betas[i,4] + qt(0.975,df=length(bonf_sig)-1)*bonf_egger$coef[2,2]) {
    coverage[4]=coverage[4]+1
  }
  
  q_rf = p.adjust(p_rf,method="fdr")
  fdr_sig = which(q_rf < 0.05)
  sr_fdr = rf[fdr_sig]
  so_fdr = do[fdr_sig]
  err_fdr = do_se[fdr_sig]
  fdr_ivw = summary(lm(so_fdr ~ sr_fdr+0,weights=err_fdr^-2))
  fdr_egger = summary(lm(so_fdr ~ sr_fdr,weights=err_fdr^-2))
  betas[i,3] = fdr_ivw$coef[1]
  betas[i,5] = fdr_egger$coef[2,1]
  
  if(fdr_ivw$coef[4] < 0.05 && betas[i,3] > 0) { positives[3]=positives[3]+1 }
  if(fdr_egger$coef[2,4] < 0.05 && betas[i,5] > 0) { positives[5]=positives[5]+1 }
  if(beta_X > betas[i,3] - qt(0.975,df=length(fdr_sig)-1)*fdr_ivw$coef[1,2] &&
     beta_X < betas[i,3] + qt(0.975,df=length(fdr_sig)-1)*fdr_ivw$coef[1,2]) {
    coverage[3]=coverage[3]+1
  }
  if(beta_X > betas[i,5] - qt(0.975,df=length(fdr_sig)-1)*fdr_egger$coef[2,2] &&
     beta_X < betas[i,5] + qt(0.975,df=length(fdr_sig)-1)*fdr_egger$coef[2,2]) {
    coverage[5]=coverage[5]+1
  }

}

posrate = positives/n_sim
covrate = coverage/n_sim
print(posrate)


