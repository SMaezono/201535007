# separate_subtypes.R
# Jiachen Liang
# modified by Sakura Maezono ##------ Mon Mar 11 23:01:18 2019 ------##

# Separates the BRCA expression matrix into 4 subtypes
#   Basal-like
#   HER2-enriched
#   Luminal A
#   Luminal B

root <- 'G:/My Drive/Jiachen Files/GDAC/'
setwd(root)

# Import normal
input <- read.table('BRCA_n.txt', sep = '\t', header=T, row.names=NULL)
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
input <- read.table('BRCA_t.txt', sep = '\t', header=T, row.names=NULL)
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
dir_sub <- 'G:/My Drive/Jiachen Files/data/BRCA_subtypes.txt'
subtypes <- read.table(dir_sub, sep = '\t', header=T, row.names=1)
bas <- rownames(subtypes)[subtypes$Subtype == 'Basal-like']
her <- rownames(subtypes)[subtypes$Subtype == 'HER2-enriched']
luA <- rownames(subtypes)[subtypes$Subtype == 'Luminal A']
luB <- rownames(subtypes)[subtypes$Subtype == 'Luminal B']

# Separate matrices
bas_n <- exp_n[,intersect(colnames(exp_n), bas)]
her_n <- exp_n[,intersect(colnames(exp_n), her)]
luA_n <- exp_n[,intersect(colnames(exp_n), luA)]
luB_n <- exp_n[,intersect(colnames(exp_n), luB)]
bas_t <- exp_t[,intersect(colnames(exp_t), bas)]
her_t <- exp_t[,intersect(colnames(exp_t), her)]
luA_t <- exp_t[,intersect(colnames(exp_t), luA)]
luB_t <- exp_t[,intersect(colnames(exp_t), luB)]

# Export
write.table(bas_n, file = 'BRCA Basal-Like_n.txt', sep = '\t')
write.table(her_n, file = 'BRCA HER2-enriched_n.txt', sep = '\t')
write.table(luA_n, file = 'BRCA Luminal A_n.txt', sep = '\t')
write.table(luB_n, file = 'BRCA Luminal B_n.txt', sep = '\t')
write.table(bas_t, file = 'BRCA Basal-Like_t.txt', sep = '\t')
write.table(her_t, file = 'BRCA HER2-enriched_t.txt', sep = '\t')
write.table(luA_t, file = 'BRCA Luminal A_t.txt', sep = '\t')
write.table(luB_t, file = 'BRCA Luminal B_t.txt', sep = '\t')

