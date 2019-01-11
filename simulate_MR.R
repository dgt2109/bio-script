
# simulate mendelian randomization

num_indiv = 20000
n_sim = 100
#OVERLAP = 0 # set a % overlap by mixing u1 and u2

# record results: obs, ivw_bonf, egger_bonf, ivw_fdr, egger_fdr
pos <- vector(mode="double",length=5)
betas = matrix(0,n_sim,5)

print(paste0("PARAMETER: ",param))
for (i in 1:n_sim) {
  if(i%%10 == 0) { print(i) }
  
  # set parameters
  num_instr = 200
  
  sd_u = 1
  sd_Ex = 1
  sd_Ey = param
  
  str_min=0
  rate = 20
  
  beta_X = 1
  beta_U = 0.5
  
  maf = runif(num_instr, 0.05, 0.5) #allele freq
  a = str_min + rexp(num_instr,rate) # instrument strength
  
  # set genotypes
  g=matrix(nrow=num_indiv,ncol=num_instr)
  for (j in 1:num_indiv) {
    for (k in 1:num_instr) {
      g[j,k] = rbinom(1, 2, maf[k])
    }
  }
  
  u1 = rnorm(num_indiv,0,sd_u)
  u2 = rnorm(num_indiv,0,sd_u)
  #    u = c(u1[0:round(OVERLAP*num_indiv)], # mixing code
  #                 u2[0:round((1-OVERLAP)*num_indiv)])
  
  # calculate exposure and outcome
  x = vector(mode="double",length=num_indiv)
  for (m in 1:num_indiv) {
    x[m] = u1[m] + rnorm(1,0,sd_Ex) + sum(a*g[m,])
  }
  
  y_obs = beta_X*x + beta_U*u1 + rnorm(num_indiv,0,sd_Ey)
  y = beta_X*x + beta_U*u2 + rnorm(num_indiv,0,sd_Ey)
  
  # generate regression results per instrument
  rnames = c("rf","p_rf","rf_se","do","do_se")
  for (var in rnames) {
    eval(parse(text=paste0(var,
         "<-vector(mode=\"double\",length=num_instr)")))
  }
  
  for (k in 1:num_instr) {
    res_rf = summary(lm(x~g[,k]))
    rf[k] = res_rf$coefficients[2] 
    rf_se[k] = res_rf$coefficients[4]
    p_rf[k] = res_rf$coefficients[8]
    res_do = summary(lm(y~g[,k]))
    do[k] = res_do$coefficients[2]
    do_se[k] = res_do$coefficients[4]
  }
  
  # compute exposure:outcome regressions
  res_obs = summary(lm(y_obs ~ x+0))
  betas[i,1] = res_obs$coef[1]
  
  if(res_obs$coef[4] < 0.05 && betas[i,1] > 0) { pos[1]=pos[1]+1 }
  
  bonf_sig = which(p_rf < 0.05/num_instr)
  sr_bonf = rf[bonf_sig]
  so_bonf = do[bonf_sig]
  err_bonf = do_se[bonf_sig]
  bonf_ivw = summary(lm(so_bonf ~ sr_bonf+0,weights=err_bonf^-2))
  bonf_egger = summary(lm(so_bonf ~ sr_bonf,weights=err_bonf^-2))
  betas[i,2] = bonf_ivw$coef[1]
  betas[i,4] = bonf_egger$coef[2,1]
  
  if(bonf_ivw$coef[4] < 0.05 && betas[i,2] > 0) { pos[2]=pos[2]+1 }
  if(bonf_egger$coef[2,4] < 0.05 && betas[i,4] > 0) { pos[4]=pos[4]+1 }
  
  q_rf = p.adjust(p_rf,method="fdr")
  fdr_sig = which(q_rf < 0.05)
  sr_fdr = rf[fdr_sig]
  so_fdr = do[fdr_sig]
  err_fdr = do_se[fdr_sig]
  fdr_ivw = summary(lm(so_fdr ~ sr_fdr+0,weights=err_fdr^-2))
  fdr_egger = summary(lm(so_fdr ~ sr_fdr,weights=err_fdr^-2))
  betas[i,3] = fdr_ivw$coef[1]
  betas[i,5] = fdr_egger$coef[2,1]
    
  if(fdr_ivw$coef[4] < 0.05 && betas[i,3] > 0) { pos[3]=pos[3]+1 }
  if(fdr_egger$coef[2,4] < 0.05 && betas[i,5] > 0) { pos[5]=pos[5]+1 }
}

power = pos/n_sim
print(power)
