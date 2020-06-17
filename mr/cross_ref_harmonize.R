
# clumping
exposures <- c("LDL","HDL","TG","BMI","SBP","T2D")
source("clump_with_MRBase.R") # set instruments

# cross-referencing
traits <- c("LDL","HDL","TG","BMI","SBP","T2D","CAD")
tables <- list()
for (i in 1:length(traits)) {
  print(paste0("Processing ", traits[i], " variants..."))
  database <- read.table(paste0("C:/Genetics/Datasets/", traits[i], "_data.txt"),
                         header=TRUE, stringsAsFactors = FALSE)
  
  tables[[i]] <- merge(database, instruments, by = "rsid")
  print(paste0(nrow(instruments)-nrow(tables[[i]]), " SNPs have missing data in ", traits[i]))
  rm(database)
}

# harmonization
harmonized <- list()
ref_EA_table <- tables[[2]][, c("rsid","EA")] # harmonize to HDL-C
colnames(ref_EA_table)[2] <- "Ref_EA"
for (i in 1:length(tables)) {
  this_table <- merge(tables[[i]], ref_EA_table, by = "rsid")
  this_table$Beta <- ifelse(tolower(this_table$EA) == this_table$Ref_EA,
                            this_table$Beta,
                            this_table$Beta * -1)
  harmonized[[i]] <- this_table
}

# merge harmonized tables
for (i in 1:length(harmonized)) {
  this_table <- harmonized[[i]][, c("rsid","Beta","SE","P")]
  colnames(this_table)[2:4] <- paste0(traits[i], "_", colnames(this_table)[2:4])
  if (i == 1) { table_data <- this_table }
  else { table_data <- merge(table_data, this_table, by = "rsid") }
}

table_data <- table_data[-c(959:961), ] # remove duplicate entries for rs4965529; there are 2 SNPs with this rsid
write.table(table_data, "C:/Genetics/SNP_Table/table_sixtrait.txt", sep='\t',
            row.names = FALSE)