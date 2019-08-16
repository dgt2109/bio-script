
byvar_results <- list()
for(k in 1:3) {
  # store matrix of results w/beta se heterog. intercept for ivw & egger
  per_lipid_res <- matrix(nrow=611,ncol=7)
  for(m in 5:615) { # replace with the exact # after presso variant removal
    PRESSO = T
    PARTIAL = T
    part <- 1:m
    source("mvmr_3.0.R")
    ivw_res <- summary(ivw)
    egger_res <- summary(eggers[[k]])
    per_lipid_res[m-4,] <- c(ivw_res$coef[k],ivw_res$coef[k,2],I[1,k],egger_res$coef[k+1],
                             egger_res$coef[k+1,2],I[2,k],egger_res$coef[1])
    if(m%%100 == 0) { print(m) }
  }
  byvar_results[[k]] <- per_lipid_res
}
