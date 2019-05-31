# mut_sig_roberts.R

# Jiachen Liang
# modified by Sakura Maezono
# Compare results from Roberts et al. 2013 and hypoxia.

t_start <- Sys.time()
library(matrixStats)
library(ggplot2)

# Set working directory
root <- "G:/My Drive/Jiachen Files/mutsig/"
setwd(root)

# Import Roberts data
print('Importing Roberts et al. data')
roberts <- read.table(file='mut_sig_context_roberts.txt', sep='\t',
                      header=TRUE, row.names=NULL)
names <- as.character(roberts[,1])
for(i in 1:length(names)){
  split <- strsplit(names[i], '-')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
}
rownames(roberts) <- names

# Import expression data
cancers <- as.character(read.table(file='cancers_list_roberts.txt',
                                   sep='\t', header=F)[,1])
exp_a_t <- vector('list', length=length(cancers))
exp_h_t <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  print(paste('Importing', cancers[i], 'data'))
  target_file <- paste('G:/My Drive/Jiachen Files/data/2016/tumour/apobec/',
                       cancers[i], '.txt', sep='')
  if (i %in% c(3,7)) {
    exp_a_t[[i]] <- read.table(file=target_file, sep=' ', header=T, row.names=1)
  } else {
    exp_a_t[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  }
 
  target_file <- paste('G:/My Drive/Jiachen Files/data/2016/tumour/hypoxia/',
                       cancers[i], '.txt', sep='')
  if (i %in% c(3,7)) {
    exp_h_t[[i]] <- read.table(file=target_file, sep=' ', header=T, row.names=1)
  } else {
    exp_h_t[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  }
  
  # Fix names
  names <- colnames(exp_a_t[[i]])
  for(j in 1:length(names)){
    split <- strsplit(names[j], '\\.')
    names[j] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(exp_a_t[[i]]) <- names
  names <- colnames(exp_h_t[[i]])
  for(j in 1:length(names)){
    split <- strsplit(names[j], '\\.')
    names[j] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(exp_h_t[[i]]) <- names
}

root <- "G:/My Drive/Jiachen Files/mutsig/output_enrichment/"
setwd(root)

# Get common samples
common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(roberts), 
                           intersect(colnames(exp_a_t[[i]]), 
                                     colnames(exp_h_t[[i]])))
}

# Build matrices
matrices <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  mat <- roberts[common[[i]],'APOBEC_enrich40']
  mat <- cbind(mat, t(exp_a_t[[i]]['APOBEC3A',common[[i]]]),
               t(exp_a_t[[i]]['APOBEC3B',common[[i]]]))
  mat <- cbind(mat, colMeans(exp_h_t[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c('APOBEC signature enrichment', 'APOBEC3A expression', 
                     'APOBEC3B expression', 'Hypoxia signature expression')
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers[i], '_enrichment.txt', sep=''), sep='\t')
}

# Plot
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`APOBEC signature enrichment`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`APOBEC signature enrichment`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`APOBEC signature enrichment`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  pdf(file=paste(cancers[i], '_A3A-enr.pdf', sep=''), width=8, height=8)
  print(sp1)
  dev.off()
  pdf(file=paste(cancers[i], '_A3B-enr.pdf', sep=''), width=8, height=8)
  print(sp2)
  dev.off()
  pdf(file=paste(cancers[i], '_hyp-enr.pdf', sep=''), width=8, height=8)
  print(sp3)
  dev.off()
}

# Correlation
cor_mat <- matrix(nrow=length(cancers), ncol=6)
rownames(cor_mat) <- cancers
colnames(cor_mat) <- c('A3A', 'A3A pv', 'A3B', 'A3B pv', 'hyp', 'hyp pv')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  test <- cor.test(df$`APOBEC signature enrichment`, 
                   df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  test <- cor.test(df$`APOBEC signature enrichment`, 
                   df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  test <- cor.test(df$`APOBEC signature enrichment`, 
                   df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
}
write.table(cor_mat, file='correlation.txt', sep='\t')

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
