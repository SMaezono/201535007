# separate_subtypes.R
# modified by Sakura Maezono ##------ Wed Mar 13 13:20:34 2019 ------##
# Separates the KIRC expression matrix by VHL mutation status

root <- 'G:/My Drive/Jiachen Files/GDAC/'
setwd(root)

# Import normal
input <- read.table('KIRC_n.txt', sep = '\t', header = T, row.names = NULL)
rn <- input[,1]
rn <- make.unique(rn, sep = '.')
output <- input[,-1]
rownames(output) <- rn
names <- colnames(output)
for(i in 1:length(names)){
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
}
colnames(output) <- names
exp_n <- output

# Import tumour
input <- read.table('KIRC_t.txt', sep = '\t', header = T, row.names = NULL)
rn <- input[,1]
rn <- make.unique(rn, sep = '.')
output <- input[,-1]
rownames(output) <- rn
names <- colnames(output)
for(i in 1:length(names)){
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
}
colnames(output) <- names
exp_t <- output

# Import subtype data
dir_sub <- 'G:/My Drive/Jiachen Files/data/KIRC_vhl.txt'
subtypes <- read.table(dir_sub, sep = '\t', header = T, row.names = NULL)
mut <- as.character(subtypes$SAMPLE_ID[!is.na(subtypes$VHL)])
wt <- as.character(subtypes$SAMPLE_ID[is.na(subtypes$VHL)])

# Fix names
mut <- substr(mut, 1, 12)
wt <- substr(wt, 1, 12)
mut <- gsub('-', '.', mut)
wt <- gsub('-', '.', wt)

# Separate matrices
mut_n <- exp_n[,intersect(colnames(exp_n), mut)]
wt_n <- exp_n[,intersect(colnames(exp_n), wt)]
mut_t <- exp_t[,intersect(colnames(exp_t), mut)]
wt_t <- exp_t[,intersect(colnames(exp_t), wt)]

# Export
write.table(mut_n, file = 'KIRC vhl_mut_n.txt', sep = '\t')
write.table(wt_n, file = 'KIRC vhl_wt_n.txt', sep = '\t')
write.table(mut_t, file = 'KIRC vhl_mut_t.txt', sep = '\t')
write.table(wt_t, file = 'KIRC vhl_wt_t.txt', sep = '\t')
