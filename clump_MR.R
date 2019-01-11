# use MRbase to extract SNPs with min p/q & low LD r2

library(TwoSampleMR)

base_dir = "C:/Users/david/Desktop/genetics/revision/clump/"

data <- read.delim(paste(base_dir, "lipid_5col.txt", sep=""), header=TRUE)
data$rsid <- make.names(data$rsid, unique=TRUE)
data2 <- format_data(data, type="exposure", 
                     snp_col = "rsid", 
                     effect_allele_col = "A1",
                     pval_col = "FDR_q")
clump <- clump_data(data2, clump_kb = 1000, clump_r2 = 0.001)

#write.csv(clump, paste(base_dir, "lipid_5col.csv", sep=""), row.names=FALSE)