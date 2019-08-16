library(MendelianRandomization)

base_dir="C:/Users/David/Desktop/genetics/data/"
table_id="simple/HDL_simple"

data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
data=data[data$Common=="TRUE", ]

h <- data$lipid_beta
herr <- data$lipid_se
y <- data$CAD_beta
yerr <- data$CAD_se

MRInputObject <- mr_input(bx = h,
                          bxse = herr,
                          by = y,
                          byse = yerr)

IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)

print(IVWObject)

WeightedMedianObject <- mr_median(MRInputObject,
                                  weighting = "weighted",
                                  distribution = "normal",
                                  alpha = 0.05,
                                  iterations = 10000,
                                  seed = 314159265)

print(WeightedMedianObject)

EggerObject <- mr_egger(MRInputObject,
                        robust = FALSE,
                        penalized = FALSE,
                        correl = FALSE,
                        distribution = "normal",
                        alpha = 0.05)

print(EggerObject)

MaxLikObject <- mr_maxlik(MRInputObject,
                          model = "default",
                          correl = FALSE,
                          psi = 0,
                          distribution = "normal",
                          alpha = 0.05)

print(MaxLikObject)
