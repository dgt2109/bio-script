
# specific data formatting tasks for dataset processing

# rename T2D headers / use notepad/text editor to rename other files to standard headers
# standard headers: Markername, Chr, Pos, rsid, EA, NEA, Freq, Beta, SE, P, N
t2d <- read.table("C:/Genetics/Datasets/T2D_data.txt", header = TRUE,
                  stringsAsFactors = FALSE) # path to original T2D_data file
colnames(t2d) <- c("rsid", "Chr", "Pos", "EA", "NEA", "Freq", "Beta", "SE", "P", "N_eff")
write.table(t2d,"C:/Genetics/Datasets/T2D_data_headerCorrected.txt", sep='\t', row.names = FALSE)

# remove _allele1_allele2 from "MarkerName" column in CAD_data
CAD_data <- read.table("C:/Genetics/Datasets/CAD_data.txt", header = TRUE,
                       stringsAsFactors = FALSE) # path to original CAD_data file
CAD_data$MarkerName <- substring(CAD_data$MarkerName, 1, regexpr("_", CAD_data$MarkerName)-1)
write.table(CAD_data, "C:/Genetics/Datasets/CAD_data_MarkerNameCorrected.txt", sep='\t', row.names = FALSE)

# for SBP & T2D, use the CAD table as a MarkerName:rsid reference to add an rsid column
rsid_reference <- CAD_data[, c("MarkerName","rsid")]
SBP_data <- read.table("C:/Genetics/Datasets/SBP_data.txt", header = TRUE,
                       stringsAsFactors = FALSE) # path to original SBP_data file
SBP_data <- merge(SBP_data, rsid_reference, all.x = TRUE)
write.table(SBP_data, "C:/Genetics/Datasets/SBP_data_with_rsid.txt", sep='\t', row.names = FALSE)

T2D_data <- read.table("C:/Genetics/Datasets/T2D_data.txt", header = TRUE,
                       stringsAsFactors = FALSE) # path to original T2D_data file
T2D_data <- merge(T2D_data, rsid_reference, all.x = TRUE)
write.table(T2D_data, "C:/Genetics/Datasets/T2D_data_with_rsid.txt", sep='\t', row.names = FALSE)
