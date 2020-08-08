library(stargazer)

# import SNP lists -> make table with all other associations
P_THRESHOLD <- 2e-5
exposures <- c("LDL", "HDL", "TG", "BMI", "T2D", "SBP")

SNP_tables <- list()
# for (i in 1:length(exposures)) {
#   SNP_list <- read.table(paste0("C:/Genetics/Variants/Significant/",exposures[i],"_SNPs.txt"),
#                          header=T, stringsAsFactors=F)
#   SNP_tables[[i]] <- SNP_list
# }

# process 185 historical SNPs
data <- read.table("C:/Genetics/Variants/do_suppl.txt", header=TRUE, stringsAsFactors=FALSE)
for (i in 1:3) {
  lipid_SNP <- data[data[, 5+i*2] < 1e-8,]
  SNP_list <- lipid_SNP[, c(1,5+i*2)]
  colnames(SNP_list) <- c("SNP","pval.exposure")
  SNP_tables[[i]] <- SNP_list
}
for (i in 4:6) {
  SNP_list <- read.table(paste0("C:/Genetics/Variants/Significant/",exposures[i],"_SNPs.txt"),
                         header=T, stringsAsFactors=F)
  SNP_tables[[i]] <- SNP_list
}

for (i in 1:length(exposures)) {
  print(paste0("Processing ",exposures[i],"..."))
  dataset <- read.table(paste0("C:/Genetics/Datasets/",exposures[i],"_data.txt"),
                        header=T, stringsAsFactors=F)
  # remove duplicate entry for SBP SNP rs4965529 in BMI and T2D datasets
  if (i == 4 | i == 5) {
    dataset <- dataset[-which(dataset$rsid == "rs4965529" & tolower(dataset$EA) == "t"),]
  }
  minitable <- dataset[, c("rsid","Beta","SE","P")]
  colnames(minitable)[2:4] <- paste0(exposures[i],"_",colnames(minitable)[2:4])
  for (j in 1:length(exposures)) {
    SNP_tables[[j]] <- merge(SNP_tables[[j]], minitable, by.x="SNP", by.y="rsid", all.x=TRUE)
  }
  rm(dataset)
}

overlap_num <- matrix(nrow=length(exposures),ncol=length(exposures)+1)
overlap_denom <- matrix(nrow=length(exposures),ncol=length(exposures)+1)
overlap_pct <- matrix(nrow=length(exposures),ncol=length(exposures)+1)

for (i in 1:length(exposures)) {
  SNP_table <- SNP_tables[[i]]
  for (j in 1:length(exposures)) {
    cross_ref <- SNP_table[, c(1:2,0:2+j*3)]
    cross_ref <- cross_ref[complete.cases(cross_ref), ]
    overlap_num[i,j] <- sum(cross_ref[, 5] < P_THRESHOLD)
    overlap_denom[i,j] <- nrow(cross_ref)
    overlap_pct[i,j] <- sum(cross_ref[, 5] < P_THRESHOLD) / nrow(cross_ref)
  }
  cross_ref <- SNP_table[, c(1:2,12:20)]
  cross_ref <- cross_ref[complete.cases(cross_ref), ]
  overlap_num[i,7] <- sum(cross_ref[, 5] < P_THRESHOLD | cross_ref[, 8] < P_THRESHOLD | cross_ref[, 11] < P_THRESHOLD)
  overlap_denom[i,7] <- nrow(cross_ref)
  overlap_pct[i,7] <- sum(cross_ref[, 5] < P_THRESHOLD | cross_ref[, 8] < P_THRESHOLD | cross_ref[, 11] < P_THRESHOLD) / nrow(cross_ref)
}

# chord plot
library(chorddiag)
m <- overlap_pct[, 1:6]*100
for (i in 1:length(exposures)) {
  SNP_table <- SNP_tables[[i]]
  cross_ref <- SNP_table[,-c(0:2+i*3)]
  cross_ref <- cross_ref[complete.cases(cross_ref), ]
  nonoverlap_num <- sum(cross_ref[, 5] > P_THRESHOLD & cross_ref[, 8] > P_THRESHOLD & cross_ref[, 11] > P_THRESHOLD
                        & cross_ref[, 14] > P_THRESHOLD & cross_ref[, 17] > P_THRESHOLD)
  nonoverlap_denom <- nrow(cross_ref)
  nonoverlap_pct <- sum(cross_ref[, 5] > P_THRESHOLD & cross_ref[, 8] > P_THRESHOLD & cross_ref[, 11] > P_THRESHOLD
                        & cross_ref[, 14] > P_THRESHOLD & cross_ref[, 17] > P_THRESHOLD) / nrow(cross_ref)
  m[i,i] <- nonoverlap_pct*100
}
dimnames(m) <- list(exposure = exposures, pleiotropic = exposures)

p <- chorddiag(m, groupnamePadding = 40)

library(htmlwidgets)
saveWidget(p, file="C:/Genetics/chord_interactive.html")
library(webshot)
webshot("C:/Genetics/chord_interactive.html","C:/Genetics/chord2.png", zoom=2)
