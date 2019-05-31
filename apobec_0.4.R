# apobec_0.4.R
# 1 thread

# Jiachen Liang

# Select for APOBEC genes that are:
#   Highly expressed
#   Differentially expressed in tumour vs normal
#   Correlated with hypoxia signature

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
#root <- "C:\\Users\\Jiachen Liang\\Documents\\GDAC\\"
setwd(root)

# Cancer list
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

#################################################################

# APOBEC gene list
genes_dir <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/"
apo <- as.character(read.table(
  file = paste(genes_dir, "tr_a.txt", sep=''),
  sep = '\t',
  header = TRUE)[,2])

# Incoming data matrices
exp_all <- matrix(nrow=length(apo), ncol=0)
exp_counts <- matrix(nrow=length(apo), ncol=0)
exp_de <- matrix(nrow=length(apo), ncol=0)

# Retrieve data
for(cancer in cancers){
  # Stats directory
  dir_stats <- paste(root, 'data/stats/', sep='')
  
  # Read expression values (tumour, apobec)
  exp_all_entry <- read.table(paste(dir_stats, cancer, '_t_a.txt', sep=''), sep='\t', header=TRUE, row.names=1)
  exp_all <- cbind(exp_all, exp_all_entry$Mean)
  
  # Read counts of expressed samples (tumour, apobec)
  exp_counts_entry <- read.table(paste(dir_stats, cancer, '_count_t_a.txt', sep=''), sep='\t', header=TRUE, row.names=1)
  exp_counts <- cbind(exp_counts, exp_counts_entry$Proportion)
  
  # Read differential expression values (tumour vs normal, apobec)
  exp_de_entry <- read.table(paste(dir_stats, cancer, '_de_a.txt', sep=''), sep='\t', header=TRUE, row.names=1)
  exp_de <- cbind(exp_de, exp_de_entry$DE.using.Mean)
}
rownames(exp_all) <- apo
colnames(exp_all) <- cancers
rownames(exp_counts) <- apo
colnames(exp_counts) <- cancers
rownames(exp_de) <- apo
colnames(exp_de) <- cancers

# Read correlation data
cor <- read.table(paste(root, 'output/2016/tumour/spearman/sig_mean_cor.txt', sep=''), sep='\t', header=TRUE, row.names=1)
cor <- cor[apo,]

# Select for gene/cancer couples that >=10% of samples are expressing
expressed <- exp_counts>0.1

# Select for gene/cancer couples that have >=0.2 DE
diff <- exp_de > 0.2

# Select for correlated gene/cancer couples
correlated <- (cor > 0.3 & !is.na(cor))

# Generate matrices
exp_sig <- exp_all
exp_sig[!expressed] <- NA

exp_sig_de <- exp_all
exp_sig_de[!(expressed&diff)] <- NA

exp_sig_cor <- exp_all
exp_sig_cor[!(expressed&correlated)] <- NA

exp_sig_de_cor <- exp_all
exp_sig_de_cor[!(expressed&diff&correlated)] <- NA

cor_sig <- cor
cor_sig[!expressed] <- NA

cor_sig_de <- cor
cor_sig_de[!(expressed&diff)] <- NA

cor_sig_cor <- cor
cor_sig_cor[!(expressed&correlated)] <- NA

cor_sig_de_cor <- cor
cor_sig_de_cor[!(expressed&diff&correlated)] <- NA

# Write files
out.dir <- paste(root, 'output/0.4/', sep='')
write.table(exp_sig, paste(out.dir, 'exp_sig.txt', sep=''), sep='\t')
write.table(exp_sig_de, paste(out.dir, 'exp_sig_de.txt', sep=''), sep='\t')
write.table(exp_sig_cor, paste(out.dir, 'exp_sig_cor.txt', sep=''), sep='\t')
write.table(exp_sig_de_cor, paste(out.dir, 'exp_sig_de_cor.txt', sep=''), sep='\t')
write.table(cor_sig, paste(out.dir, 'cor_sig.txt', sep=''), sep='\t')
write.table(cor_sig_de, paste(out.dir, 'cor_sig_de.txt', sep=''), sep='\t')
write.table(cor_sig_cor, paste(out.dir, 'cor_sig_cor.txt', sep=''), sep='\t')
write.table(cor_sig_de_cor, paste(out.dir, 'cor_sig_de_cor.txt', sep=''), sep='\t')

#################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
