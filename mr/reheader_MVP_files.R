
traits <- c("LDL","HDL","TG")

for (trait in traits) {
  print(trait)
  data1 <- read.table(paste0("C:/genetics/datasets/lipid/MVP_",trait,"_data.txt"), 
                      header=TRUE, stringsAsFactors=FALSE)
  data1 <- data1[,-6]
  colnames(data1) <- c("rsid","Chr","Pos","EA","NEA","Freq","N","Beta","SE","P","Ethnic.Concordance")
  
  write.table(data1, paste0("C:/genetics/datasets/",trait,"_data.txt"), 
              sep='\t', row.names = FALSE)
}