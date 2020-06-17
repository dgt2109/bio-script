
#clumping
exposures <- c("LDL","HDL","TG")
source("clump_with_MRBase.R") # set instruments

# cross-referencing
tables <- list()
traits <- c("LDL","HDL","TG","CAD")
for (i in 1:length(traits)) {
  print(paste0("Processing ", traits[i], "data..."))
  database <- read.table(paste0("C:/Genetics/Datasets/", traits[i], "_data.txt"),
                         header=TRUE, stringsAsFactors = FALSE)
  
  tables[[i]] <- merge(database, instruments, by = "rsid")
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

data <- table_data
# write.table(data, "C:/Genetics/SNP_Table/table_three-trait.txt", sep='\t', 
#             row.names = FALSE)