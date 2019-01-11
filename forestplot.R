
# vectors contain est for weak, moderate, strong in order

# CI weak < 10 -0.2126808, 0.01790345, 0.2484877
# CI moderate 10-100 -0.2603336, -0.1617508, -0.06316794
# CI strong > 100 -0.1369203, -0.01182121, 0.1132778

beta_LDL <- c(0.1413372,0.2158929,0.4578809,0.3113675)
lower_LDL <- c(-0.08202917,0.1354332,0.3426157,0.2552795)
upper_LDL <- c(0.3647036,0.2963525,0.5731462,0.3674555)

beta_HDL <- c(0.01790345, -0.1656635, -0.01182121, -0.1097544)
lower_HDL <- c(-0.2126808, -0.2663317, -0.1369203, -0.1695207)
upper_HDL <- c(0.2484877,-0.06499535,0.1132778, -0.04998817)

beta_TG <- c(0.04946334,0.1475094,-0.04221028,0.1178565)
lower_TG <- c(-0.2198454,0.04737182,-0.285871,0.04402893)
upper_TG <- c(0.3187721,0.2476469,0.2014504,0.191684)

OR_mean <- c(NA,exp(beta_LDL),NA,exp(beta_HDL),NA,exp(beta_TG))
OR_lower <- c(NA,exp(lower_LDL),NA,exp(lower_HDL),NA,exp(lower_TG))
OR_upper <- c(NA,exp(upper_LDL),NA,exp(upper_HDL),NA,exp(upper_TG))

ni <- c(396,217,18,631,318,296,17,631,367,252,12,631)

labels <- cbind(
  c("",rep("LDL-C",4),"",rep("HDL-C",4),"",rep("TG",4),""),
  c("F-statistic",rep(c("F<10","10<F<100","F>100","All",""),3)),
  c("  # Instr.",ni[1:4],"",ni[5:8],"",ni[9:12],""),
  c("OR",round(exp(beta_LDL),digits=2),"",round(exp(beta_HDL),digits=2),
    "",round(exp(beta_TG),digits=2),"")
  )
base_dir="C:/Users/David/Desktop/genetics/revision/snp680/strength/"
svg(paste0(base_dir,"byStr_fp.svg"),width=5,height=4)
forestplot(labels,OR_mean, OR_lower,OR_upper,
           is.summary=c(TRUE,rep(c(FALSE,FALSE,FALSE,TRUE,TRUE),3)),
           zero=1)
dev.off()