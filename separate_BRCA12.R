##------ Wed Feb 27 11:15:34 2019 ------##
# originally jiachen's script, modified by Sakura
# Separates the BRCA expression matrix by BRCA1/2 mutation status

# set working directory
setwd('G:/My Drive/Jiachen Files/GDAC/')

# functions:
fix_names <- function(listofsampleid) {
  # set names so it would match normal and tumor names
  new_names <- substr(listofsampleid, 1, 12)
  new_names <- gsub('-', '.', new_names)
  return(new_names)
}

sep_matrices <- function(RNAmatrix, patientlist, filename) {
  # separate matrices by mutation and tissue subtypes
  new_colnames <- intersect(colnames(RNAmatrix), patientlist)
  tissuesubtype_df <- data.frame(RNAmatrix[,intersect(colnames(RNAmatrix),
                                                      patientlist)])
  colnames(tissuesubtype_df) <- new_colnames
  
  if (ncol(RNAmatrix) < 120) {
    write.table(tissuesubtype_df, file = paste(filename, '_n.txt', sep = ""),
                sep = '\t')
  } else {
    write.table(tissuesubtype_df, file = paste(filename, '_t.txt', sep = ""),
                sep = '\t')
  }
  print("file created")
}
# Main script

# Import normal
input <- read.table('BRCA_n.txt', sep = '\t', header = T, row.names = NULL)
rn <- input[,1]
rn <- make.unique(rn, sep = '.')
output <- input[,-1]
rownames(output) <- rn
names <- colnames(output)
for (i in 1:length(names)) {
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
}
colnames(output) <- names
exp_n <- output

# Import tumour
input <- read.table('BRCA_t.txt', sep = '\t', header = T, row.names = NULL)
rn <- input[,1]
rn <- make.unique(rn, sep = '.')
output <- input[,-1]
rownames(output) <- rn
names <- colnames(output)
for (i in 1:length(names)) {
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
}
colnames(output) <- names
exp_t <- output

# Import BRCA1/2 mutation data
# # only had to do once (minimize to needed info only)
# dir_sub <- "G:/My Drive/Jiachen Files/data/alterations_across_samples_BRCA_12.txt"
# # change the word "no alteration" to NA by using na.strings
# gen_info_BRCA1or2 <- read.table(dir_sub, sep = '\t', header = T,
#                                 row.names = NULL, na.strings = "no alteration")
# mutations_only <- gen_info_BRCA1or2[, c(1,2,7,11)]
# colnames(mutations_only)[c(3,4)] <- c("BRCA1", "BRCA2")
# # set row.names to false to not include row names in the final txt
# write.table(mutations_only,
#             file = 'G:/My Drive/Jiachen Files/data/BRCA_12.txt',
#             sep = '\t', row.names = FALSE)

dir_sub <- "G:/My Drive/Jiachen Files/data/BRCA_12.txt"
BRCA1or2 <- read.table(dir_sub, sep = '\t', header = T, row.names = NULL)
BRCA1_mut <- as.character(BRCA1or2$SAMPLE_ID[!is.na(BRCA1or2$BRCA1)])
BRCA2_mut <- as.character(BRCA1or2$SAMPLE_ID[!is.na(BRCA1or2$BRCA2)])
# samples with either BRCA1 and 2 mutation or both
mut <- c(BRCA1_mut, BRCA2_mut)
BRCA1_wt <- as.character(BRCA1or2$SAMPLE_ID[is.na(BRCA1or2$BRCA1)])
BRCA2_wt <- as.character(BRCA1or2$SAMPLE_ID[is.na(BRCA1or2$BRCA2)])
# samples with no mutations of either BRCA1 and 2 mutation
wt <- as.character(BRCA1or2$SAMPLE_ID[is.na(BRCA1or2$BRCA1) & is.na(BRCA1or2$BRCA2)])

# Fix names

BRCA1_mut <- fix_names(BRCA1_mut)
BRCA2_mut <- fix_names(BRCA2_mut)
mut <- fix_names(mut)
BRCA1_wt <- fix_names(BRCA1_wt)
BRCA2_wt <- fix_names(BRCA2_wt)
wt <- fix_names(wt)


# Separate matrices
# NOTE: <2 patients in BRCA1or2, BRCA2, and BRCA1 mut in normal tissues

sep_matrices(exp_n, BRCA1_mut, "BRCA1_mut")
sep_matrices(exp_n, BRCA1_wt, "BRCA1_wt")
sep_matrices(exp_t, BRCA1_mut, "BRCA1_mut")
sep_matrices(exp_t, BRCA1_wt, "BRCA1_wt")

sep_matrices(exp_n, BRCA2_mut, "BRCA2_mut")
sep_matrices(exp_n, BRCA2_wt, "BRCA2_wt")
sep_matrices(exp_t, BRCA2_mut, "BRCA2_mut")
sep_matrices(exp_t, BRCA2_wt, "BRCA2_wt")

sep_matrices(exp_n, mut, "BRCA1or2_mut")
sep_matrices(exp_n, wt, "BRCA1or2_wt")
sep_matrices(exp_t, mut, "BRCA1or2_mut")
sep_matrices(exp_t, wt, "BRCA1or2_wt")


setwd("G:/My Drive/Jiachen Files/codes_supp")
# save details of this R script
scriptname <- "separate_BRCA12"
savehistory(file = paste(scriptname, ".Rhistory", sep = ""))
sink(file = paste(Sys.Date(), scriptname, ".txt", sep = ""),
     type = c("output", "message"))
# details of the packages and other useful info.
print(sessionInfo()) 
sink(NULL)


