TEST_SIM = F
n_sim = 250

param <- matrix(c(0.5,0.0,0.6,0.1,0.7,0.2),nrow=2) 
TPR <- matrix(nrow=6,ncol=5)
FPR <- matrix(nrow=6,ncol=5)
CRCI <- matrix(nrow=6,ncol=5)
sim_betas <- list()

# 6 values/lists are saved:
# for increasing meas. error and then increasing confounding
for (k in 1:2) { # iterate over oY or U
  for (m in 1:3) { # iterate over parameter range
    if (k == 1) { sd_Ey = param[1,m]; beta_U = param[2,2] }
    if (k == 2) { sd_Ey = param[1,2]; beta_U = param[2,m] }
    
    beta_X = 0.1
    source("simulation/simulate_MR.R")
    TPR[m+(k-1)*3,] <- posrate
    CRCI[m+(k-1)*3,] <- covrate
    sim_betas[[m+(k-1)*3]] <- betas
    
    beta_X = 0
    source("simulation/simulate_MR.R")
    FPR[m+(k-1)*3,] <- posrate
  }
}
