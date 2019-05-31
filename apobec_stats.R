# apobec_stats.R
# n threads in MacOS and Linux using FORK
# 8 threads in Windows

# Jiachen Liang

# Get statistics computed by gdac_preprocess.R for APOBEC and hypoxia

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/"
#root <- "C:\\Users\\Jiachen Liang\\Documents\\GDAC\\"
setwd(root)

# Destination directory
dest <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/stats/"
#dest <- 'C:\\Users\\Jiachen Liang\\OneDrive - OnTheHub - The University of Oxford\\Buffa Lab\\apobec\\data\\stats\\'

# Using Windows?
win <- FALSE

# Cancer list
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

# Genes
genes_dir <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/"
#genes_dir <- "C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/"
apo <- as.character(read.table(
  file = paste(genes_dir, "tr_a.txt", sep=''),
  sep = '\t',
  header = TRUE)[,2])
hyp <- as.character(read.table(
  file = paste(genes_dir, "tr_h.txt", sep=''),
  sep = '\t',
  header = TRUE)[,2])

# Import packages
if(win){
  library(doParallel)
}else{
  library(parallel)
}

stats <- function(import_index){
  cancer <- cancers[import_index]
  
  # Import
  stats_a <- read.table(file=paste(root, 'stats/', cancer, '_g_a.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  stats_t <- read.table(file=paste(root, 'stats/', cancer, '_g_t.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  stats_n <- read.table(file=paste(root, 'stats/', cancer, '_g_n.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  stats_de <- read.table(file=paste(root, 'stats/', cancer, '_de.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  count_a <- read.table(file=paste(root, 'stats/', cancer, '_count_a.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  count_t <- read.table(file=paste(root, 'stats/', cancer, '_count_t.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  count_n <- read.table(file=paste(root, 'stats/', cancer, '_count_n.txt', sep=''), sep='\t', header=TRUE, row.names=NULL)
  
  # Fix names
  names <- stats_a[,1]
  stats_a <- as.matrix(stats_a[,2:4])
  stats_t <- as.matrix(stats_t[,2:4])
  stats_n <- as.matrix(stats_n[,2:4])
  stats_de <- as.matrix(stats_de[,2:3])
  count_a <- as.matrix(count_a[,2:3])
  count_t <- as.matrix(count_t[,2:3])
  count_n <- as.matrix(count_n[,2:3])
  colnames(stats_a) <- c('Mean', 'Median', 'SD')
  colnames(stats_t) <- c('Mean', 'Median', 'SD')
  colnames(stats_n) <- c('Mean', 'Median', 'SD')
  colnames(stats_de) <- c('DE using Mean', 'DE using Median')
  rownames(stats_a) <- names
  rownames(stats_t) <- names
  rownames(stats_n) <- names
  rownames(stats_de) <- names
  rownames(count_a) <- names
  rownames(count_t) <- names
  rownames(count_n) <- names
  colnames(count_a) <- c("Count", "Proportion")
  colnames(count_t) <- c("Count", "Proportion")
  colnames(count_n) <- c("Count", "Proportion")
  
  # Select
  stats_a_a <- stats_a[apo, ]
  stats_a_h <- stats_a[hyp, ]
  stats_t_a <- stats_t[apo, ]
  stats_t_h <- stats_t[hyp, ]
  stats_n_a <- stats_n[apo, ]
  stats_n_h <- stats_n[hyp, ]
  stats_de_a <- stats_de[apo, ]
  stats_de_h <- stats_de[hyp, ]
  count_a_a <- count_a[apo, ]
  count_a_h <- count_a[hyp, ]
  count_t_a <- count_t[apo, ]
  count_t_h <- count_t[hyp, ]
  count_n_a <- count_n[apo, ]
  count_n_h <- count_n[hyp, ]
  
  # Write
  write.table(stats_a_a, paste(dest, cancer, '_a_a.txt', sep=''), sep='\t')
  write.table(stats_a_h, paste(dest, cancer, '_a_h.txt', sep=''), sep='\t')
  write.table(stats_t_a, paste(dest, cancer, '_t_a.txt', sep=''), sep='\t')
  write.table(stats_t_h, paste(dest, cancer, '_t_h.txt', sep=''), sep='\t')
  write.table(stats_n_a, paste(dest, cancer, '_n_a.txt', sep=''), sep='\t')
  write.table(stats_n_h, paste(dest, cancer, '_n_h.txt', sep=''), sep='\t')
  write.table(stats_de_a, paste(dest, cancer, '_de_a.txt', sep=''), sep='\t')
  write.table(stats_de_h, paste(dest, cancer, '_de_h.txt', sep=''), sep='\t')
  write.table(count_a_a, paste(dest, cancer, '_count_a_a.txt', sep=''), sep='\t')
  write.table(count_a_h, paste(dest, cancer, '_count_a_h.txt', sep=''), sep='\t')
  write.table(count_t_a, paste(dest, cancer, '_count_t_a.txt', sep=''), sep='\t')
  write.table(count_t_h, paste(dest, cancer, '_count_t_h.txt', sep=''), sep='\t')
  write.table(count_n_a, paste(dest, cancer, '_count_n_a.txt', sep=''), sep='\t')
  write.table(count_n_h, paste(dest, cancer, '_count_n_h.txt', sep=''), sep='\t')
  
  return(cancer)
}

# Parallelization
if(win){
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:length(cancers)) %dopar% stats(i)
  stopCluster(cl)
}else{
  f <- vector("list", length=length(cancers))
  for(i in 1:length(cancers)){
    f[[i]] <- mcparallel(stats(i))
  }
  mccollect(f)
}

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)