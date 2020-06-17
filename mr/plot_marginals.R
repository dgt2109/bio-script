
library(ggplot2)

# read ivw and egger models and data, perform plotting
# adjust CAD according to ivw model offsets

data <- read.table("C:/Genetics/SNP_Table/table_sixtrait.txt", header=TRUE, stringsAsFactors=FALSE)

outlier_data <- read.table("C:/Genetics/MRPRESSO/sixtrait_MRPRESSO.txt", header = TRUE, stringsAsFactors=FALSE)
rsid_MRPRESSO <- outlier_data[outlier_data$Pvalue < 0.05, 1]
data <- data[!(data$rsid %in% rsid_MRPRESSO), ]

param <- read.table("C:/Genetics/Results/sixtrait_MRPRESSO-results.txt", header=TRUE, stringsAsFactors=FALSE,
                    row.names=NULL)
rownames(param) <- make.unique(param[,1])
param <- param[,-1]

beta_indexes <- (1:7*3)-1
plots <- list()
exposures <- c("LDL-C Effect Size","HDL-C Effect Size","TG Effect Size",
               "BMI Effect Size","SBP Effect Size","T2D Effect Size")
legend_ypos <- c(0.15, 0.85, 0.15, 0.15, 0.15, 0.15)
for (ref_col_index in beta_indexes[1:6]) {
  # reorient to trait of interest
  for (i in beta_indexes) {
    if (i == ref_col_index) { next }
    data[, i] <- ifelse(data[, ref_col_index]>0,
                        data[, i],
                        data[, i]*-1)
  }
  data[, ref_col_index] <- abs(data[, ref_col_index])
  
  # adjust CAD estimate for covariates
  ref_index <- (ref_col_index+1)/3
  nonref_indexes <- beta_indexes[-ref_index]
  coef_indexes <- (1:6)[-ref_index]
  tmp <- colnames(data)[ref_col_index]
  colnames(data)[ref_col_index] <- "x"
  data$y <- data$CAD_Beta
  for (j in 1:5) { data$y <- data$y - data[ , nonref_indexes[j]] * log(param[coef_indexes[j], 1]) }
  data$y_ivw <- data[ , ref_col_index] * log(param[ref_index, 1])
  data$y_egger <- data[ , ref_col_index] * log(param[ref_index+6, 1]) + param[ref_index+6, 7]
  
  # plotting
  commonTheme <- list(labs(x=exposures[ref_index],y="Adjusted CAD Effect Size"), 
                      theme_bw()+theme(legend.position=c(0.7,legend_ypos[ref_index]),
                                       legend.title=element_blank(),
                                       legend.spacing.y = unit(0,"mm"),
                                       legend.background=element_blank(),
                                       legend.box.background=element_rect(colour="black")))
  plots[[ref_index]] <- ggplot(data=data, aes(x,y)) + geom_point() + 
    geom_line(data=data, aes(x, y_ivw, linetype="Multivariable MR-IVW")) + 
    geom_line(data=data, aes(x, y_egger, linetype="Multivariable MR-Egger")) +
    scale_linetype_manual(name="",values=c("Multivariable MR-IVW"=1,"Multivariable MR-Egger"=2),
                          guide=guide_legend(reverse=TRUE)) +
    commonTheme
  
  colnames(data)[ref_col_index] <- tmp
}

tiff("C:/Genetics/Figures/f3a.tif",width=4,height=3.5,units="in",res=500)
plots[[1]]
dev.off()

tiff("C:/Genetics/Figures/f3b.tif",width=4,height=3.5,units="in",res=500)
plots[[2]]
dev.off()

tiff("C:/Genetics/Figures/f3c.tif",width=4,height=3.5,units="in",res=500)
plots[[3]]
dev.off()

tiff("C:/Genetics/Figures/f3d.tif",width=4,height=3.5,units="in",res=500)
plots[[4]]
dev.off()

tiff("C:/Genetics/Figures/f3e.tif",width=4,height=3.5,units="in",res=500)
plots[[5]]
dev.off()

tiff("C:/Genetics/Figures/f3f.tif",width=4,height=3.5,units="in",res=500)
plots[[6]]
dev.off()