# mut_sig_counts.R

# Jiachen Liang
# modified by Sakura
# Compare the raw counts of TCW -> TTW/TGW and hypoxia.

t_start <- Sys.time()
library(matrixStats)
library(ggplot2)
library(deconstructSigs)

# Set working directory
root <- "G:/My Drive/Jiachen Files/Alexandrov_mut_sig/"
setwd(root)

# Import mutation data
print('Importing Alexandrov data')
cancers_alex <- as.character(read.table(file='cancer_list_alex.txt',
                                        sep='\t', header=F)[,1])
cancers <- as.character(read.table(file='cancer_list.txt', sep='\t',
                                   header=F)[,1])
mut <- data.frame(nrow=96, ncol=0)
for(c in cancers_alex){
  path <- paste('mutational_catalogs/exomes/', c, '/', c, '_exomes_mutational_catalog_96_subs.txt', sep='')
  this_mut <- read.table(file=path, sep='\t', header=TRUE, row.names=1)
  mut <- cbind(mut, this_mut)
}
names <- colnames(mut)
for(i in 1:length(names)){
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
}
colnames(mut) <- names

# Subsize APOBEC mutational signature as defined by Roberts et al. 2013
apo_sign <- c('T[C>T]A', 'T[C>G]A', 'T[C>T]T', 'T[C>G]T')
mut_apo <- mut[apo_sign,]

# Create list of sum of all APOBEC mutations
mut_sums_apo <- matrix(nrow=ncol(mut_apo), ncol=1)
mut_sums_apo[,1] <- colSums(mut_apo)
rownames(mut_sums_apo) <- colnames(mut_apo)

# Log scale list
mut_sums_apo_log <- log10(mut_sums_apo)
mut_sums_apo_log[is.infinite(mut_sums_apo_log)] <- NA

# Import expression data
exp_a_t <- vector('list', length=length(cancers))
exp_h_t <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  print(paste('Importing', cancers[i], 'data'))
  target_file <- paste('G:/My Drive/Jiachen Files/data/2016/tumour/apobec/', cancers[i], '.txt', sep='')
  exp_a_t[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  target_file <- paste('G:/My Drive/Jiachen Files/data/2016/tumour/hypoxia/', cancers[i], '.txt', sep='')
  exp_h_t[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
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

# Counts analysis ####################################################################################

root <- "G:/My Drive/Jiachen Files/Alexandrov_mut_sig/output_counts/"
setwd(root)

print('Analyzing mutational signatures using counts')

# Get common samples
common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(mut_sums_apo), intersect(colnames(exp_a_t[[i]]),
                                                             colnames(exp_h_t[[i]])))
}

# Build matrices
matrices <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  mat <- mut_sums_apo[common[[i]],1]
  mat <- cbind(mat, t(exp_a_t[[i]]['APOBEC3A',common[[i]]]),
               t(exp_a_t[[i]]['APOBEC3B',common[[i]]]))
  mat <- cbind(mat, colMeans(exp_h_t[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c('APOBEC mutation count', 'APOBEC3A expression',
                     'APOBEC3B expression', 'Hypoxia signature expression')
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers[i], '_mut_summary.txt', sep=''), sep='\t')
}

# Plot
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`APOBEC mutation count`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`APOBEC mutation count`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`APOBEC mutation count`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  pdf(file=paste(cancers[i], '_A3A-mut_count.pdf', sep=''), width=8, height=8)
  print(sp1)
  dev.off()
  pdf(file=paste(cancers[i], '_A3B-mut_count.pdf', sep=''), width=8, height=8)
  print(sp2)
  dev.off()
  pdf(file=paste(cancers[i], '_hyp-mut_count.pdf', sep=''), width=8, height=8)
  print(sp3)
  dev.off()
}

# Correlation
cor_mat <- matrix(nrow=length(cancers), ncol=6)
rownames(cor_mat) <- cancers
colnames(cor_mat) <- c('A3A', 'A3A pv', 'A3B', 'A3B pv', 'hyp', 'hyp pv')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  if (length(df$`APOBEC3A expression`[is.na(df$`APOBEC3A expression`)]) >= (length(df$`APOBEC3A expression`)-1)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`APOBEC mutation count`, df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  }
  if (length(df$`APOBEC3B expression`[is.na(df$`APOBEC3B expression`)]) >= (length(df$`APOBEC3B expression`)-1)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`APOBEC mutation count`, df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  }
  if (length(df$`Hypoxia signature expression`[is.na(df$`Hypoxia signature expression`)]) >= (length(df$`Hypoxia signature expression`)-1)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`APOBEC mutation count`, df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
  }
}
write.table(cor_mat, file='correlation.txt', sep='\t')

# Log scale analysis ###################################################################################

root <- "G:/My Drive/Jiachen Files/Alexandrov_mut_sig/output_counts/log/"
setwd(root)

print('Analyzing mutational signatures using counts (log)')

# Build matrices
matrices <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  mat <- mut_sums_apo_log[common[[i]],1]
  mat <- cbind(mat, t(exp_a_t[[i]]['APOBEC3A',common[[i]]]), t(exp_a_t[[i]]['APOBEC3B',common[[i]]]))
  mat <- cbind(mat, colMeans(exp_h_t[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c('APOBEC mutation count', 'APOBEC3A expression', 
                     'APOBEC3B expression', 'Hypoxia signature expression')
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers[i], '_mut_summary.txt', sep=''), sep='\t')
}

# Plot
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`APOBEC mutation count`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`APOBEC mutation count`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`APOBEC mutation count`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  pdf(file=paste(cancers[i], '_A3A-mut_count.pdf', sep=''), width=8, height=8)
  print(sp1)
  dev.off()
  pdf(file=paste(cancers[i], '_A3B-mut_count.pdf', sep=''), width=8, height=8)
  print(sp2)
  dev.off()
  pdf(file=paste(cancers[i], '_hyp-mut_count.pdf', sep=''), width=8, height=8)
  print(sp3)
  dev.off()
}

# Correlation
cor_mat <- matrix(nrow=length(cancers), ncol=6)
rownames(cor_mat) <- cancers
colnames(cor_mat) <- c('A3A', 'A3A pv', 'A3B', 'A3B pv', 'hyp', 'hyp pv')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  if (length(df$`APOBEC3A expression`[is.na(df$`APOBEC3A expression`)]) >= (length(df$`APOBEC3A expression`)-1)) {
    print("skipping, not enough samples")
  } else {
    test <- cor.test(df$`APOBEC mutation count`, 
                     df$`APOBEC3A expression`, 
                     method = 'spearman')
    cor_mat[i,1] <- as.numeric(test$estimate)
    cor_mat[i,2] <- as.numeric(test$p.value)
  }
  if (length(df$`APOBEC3B expression`[is.na(df$`APOBEC3B expression`)]) >= (length(df$`APOBEC3B expression`)-1)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`APOBEC mutation count`, 
                   df$`APOBEC3B expression`,
                   method = 'spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  }
  if (length(df$`Hypoxia signature expression`[is.na(df$`Hypoxia signature expression`)]) >= (length(df$`Hypoxia signature expression`)-1)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`APOBEC mutation count`,
                   df$`Hypoxia signature expression`, 
                   method = 'spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
  }
}
write.table(cor_mat, file='correlation.txt', sep='\t')

# Prepare for deconstructSigs (fractionize) ############################################################

root <- "G:/My Drive/Jiachen Files/Alexandrov_mut_sig/"
setwd(root)

# mut_total <- colSums(mut)
# mut_frac <- mut
# for(i in 1:nrow(mut_frac)){
#   for(j in 1:ncol(mut_frac)){
#     mut_frac[i,j] <- mut_frac[i,j]/as.numeric(mut_total[j])
#   }
# }
# mut_frac <- mut_frac[,-c(1,2)]
# mut_frac <- as.data.frame(t(mut_frac))
mut_counts <- mut[,-c(1,2)]
mut_counts <- as.data.frame(t(mut_counts))

# deconstructSigs ######################################################################################

# dcs_frac <- whichSignatures(mut_frac, rownames(mut_frac))
# dcs_counts <- whichSignatures(mut_counts, rownames(mut_counts), contexts.needed=T, tri.counts.method='exome')

if(file.exists('dcs_all.txt')){
  print('Reading from pre-analyzed deconstructSigs')
  dcs_counts <- read.table(file='dcs_all.txt', header=T, row.names=1, sep='\t')
}else{
  print('Analyzing mutational signatures using deconstructSigs')
  dcs_counts <- as.data.frame(matrix(nrow=0, ncol=27))
  for(i in 1:nrow(mut_counts)){
    print(paste('Analyzing sample', i, 'of', nrow(mut_counts)))
    dcs_counts <- rbind(
      dcs_counts,
      whichSignatures(mut_counts,
                      rownames(mut_counts)[i],
                      contexts.needed=T,
                      tri.counts.method='exome'
                      )$weights
      )
  }
  # Output matrix
  write.table(dcs_counts, file='dcs_all.txt', sep='\t')
}

# Subset for signatures 2 and 13
dcs_apo <- dcs_counts[,c(3,14)]

# Get common samples
common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(dcs_apo), intersect(colnames(exp_a_t[[i]]), colnames(exp_h_t[[i]])))
}

root <- "G:/My Drive/Jiachen Files/Alexandrov_mut_sig/output_dcs/"
setwd(root)

# Build matrices
matrices <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  mat <- dcs_apo[common[[i]],1]
  mat <- cbind(mat, dcs_apo[common[[i]],2])
  mat <- cbind(mat, (mat[,1] + mat[,2])/2, pmax(mat[,1], mat[,2]))
  mat <- cbind(mat, t(exp_a_t[[i]]['APOBEC3A',common[[i]]]), t(exp_a_t[[i]]['APOBEC3B',common[[i]]]))
  mat <- cbind(mat, pmax(mat[,5], mat[,6]))
  mat <- cbind(mat, colMeans(exp_h_t[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c(
    'Signature 2', 
    'Signature 13', 
    'Average of 2 and 13', 
    'Max of 2 and 13', 
    'APOBEC3A expression', 
    'APOBEC3B expression', 
    'Max of 3A and 3B exp',
    'Hypoxia signature expression'
    )
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers[i], '_signature_summary.txt', sep=''), sep='\t')
}

# Plot
print('Plotting')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp4 <- ggplot(df, aes(x=`Max of 3A and 3B exp`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  pdf(file=paste(cancers[i], '_A3A-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp1)
  dev.off()
  pdf(file=paste(cancers[i], '_A3B-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp2)
  dev.off()
  pdf(file=paste(cancers[i], '_hyp-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp3)
  dev.off()
  pdf(file=paste(cancers[i], '_max3A3B-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp4)
  dev.off()
}

# Correlation
print('Analyzing correlation')
cor_mat <- matrix(nrow=length(cancers), ncol=8)
rownames(cor_mat) <- cancers
colnames(cor_mat) <- c('A3A', 'A3A pv', 'A3B', 'A3B pv', 'hyp', 'hyp pv', 'max3A3B', 'max3A3B pv')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  if (length(df$`APOBEC3A expression`[is.na(df$`APOBEC3A expression`)]) >= (length(df$`APOBEC3A expression`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  }
  if (length(df$`APOBEC3B expression`[is.na(df$`APOBEC3B expression`)]) >= (length(df$`APOBEC3B expression`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2) | i == 5) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  }
  if (length(df$`Hypoxia signature expression`[is.na(df$`Hypoxia signature expression`)]) >= (length(df$`Hypoxia signature expression`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2)) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`Max of 2 and 13`, df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
  }
  if (length(df$`Max of 3A and 3B exp`[is.na(df$`Max of 3A and 3B exp`)]) >= (length(df$`Max of 3A and 3B exp`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2) | i == 5) {
    print("skipping, not enough samples")
  } else {
  test <- cor.test(df$`Max of 2 and 13`, df$`Max of 3A and 3B exp`, method='spearman')
  cor_mat[i,7] <- as.numeric(test$estimate)
  cor_mat[i,8] <- as.numeric(test$p.value)
  }
}
write.table(cor_mat, file='correlation.txt', sep='\t')

#################################
# Filtering

root <- "G:/My Drive/Jiachen Files/Alexandrov_mut_sig/output_dcs/filtered/"
setwd(root)

# Remove zeroes
dcs_apo_rm0 <- dcs_apo
dcs_apo_rm0[dcs_apo_rm0==0] <- NA

# Get expression thresholds
thr_mat <- read.table("G:/My Drive/Jiachen Files/data/exp_med.txt", 
                      sep='\t', header=T, row.names=1)

# Build matrices
matrices <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  mat <- dcs_apo_rm0[common[[i]],1]
  mat <- cbind(mat, dcs_apo_rm0[common[[i]],2])
  mat <- cbind(mat, (mat[,1] + mat[,2])/2, pmax(mat[,1], mat[,2], na.rm = T))
  mat <- cbind(mat, t(exp_a_t[[i]]['APOBEC3A',common[[i]]]),
               t(exp_a_t[[i]]['APOBEC3B',common[[i]]]))
  
  # Remove unexpressed
  thr <- as.numeric(thr_mat[cancers[i], 2])
  for(j in 1:length(mat[,5])){
    if(mat[j,5]<thr | is.na(mat[j,5])){
      mat[j,5] <- NA
    }
  }
  for(j in 1:length(mat[,6])){
    if(mat[j,6]<thr | is.na(mat[j,6])){
      mat[j,6] <- NA
    }
  }
  
  mat <- cbind(mat, pmax(mat[,5], mat[,6], na.rm = T))
  mat <- cbind(mat, colMeans(exp_h_t[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c(
    'Signature 2', 
    'Signature 13', 
    'Average of 2 and 13', 
    'Max of 2 and 13', 
    'APOBEC3A expression', 
    'APOBEC3B expression', 
    'Max of 3A and 3B exp',
    'Hypoxia signature expression'
  )
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers[i], '_signature_summary.txt',
                              sep=''), sep='\t')
}

# Plot
print('Plotting')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  sp4 <- ggplot(df, aes(x=`Max of 3A and 3B exp`, y=`Max of 2 and 13`)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  pdf(file=paste(cancers[i], '_A3A-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp1)
  dev.off()
  pdf(file=paste(cancers[i], '_A3B-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp2)
  dev.off()
  pdf(file=paste(cancers[i], '_hyp-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp3)
  dev.off()
  pdf(file=paste(cancers[i], '_max3A3B-mut_dcs.pdf', sep=''), width=8, height=8)
  print(sp4)
  dev.off()
}

# Correlation
print('Analyzing correlation')
cor_mat <- matrix(nrow=length(cancers), ncol=8)
rownames(cor_mat) <- cancers
colnames(cor_mat) <- c('A3A', 'A3A pv', 'A3B', 'A3B pv', 'hyp',
                       'hyp pv', 'max3A3B', 'max3A3B pv')
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  if (length(df$`APOBEC3A expression`[is.na(df$`APOBEC3A expression`)]) >= (length(df$`APOBEC3A expression`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2)) {
    print("skipping, not enough samples")
  } else {
    test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3A expression`, method='spearman')
    cor_mat[i,1] <- as.numeric(test$estimate)
    cor_mat[i,2] <- as.numeric(test$p.value)
  }
  if (length(df$`APOBEC3B expression`[is.na(df$`APOBEC3B expression`)]) >= (length(df$`APOBEC3B expression`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2) | i == 5) {
    print("skipping, not enough samples")
  } else {
    test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3B expression`, method='spearman')
    cor_mat[i,3] <- as.numeric(test$estimate)
    cor_mat[i,4] <- as.numeric(test$p.value)
  }
  if (length(df$`Hypoxia signature expression`[is.na(df$`Hypoxia signature expression`)]) >= (length(df$`Hypoxia signature expression`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2)) {
    print("skipping, not enough samples")
  } else {
    test <- cor.test(df$`Max of 2 and 13`, df$`Hypoxia signature expression`, method='spearman')
    cor_mat[i,5] <- as.numeric(test$estimate)
    cor_mat[i,6] <- as.numeric(test$p.value)
  }
  if (length(df$`Max of 3A and 3B exp`[is.na(df$`Max of 3A and 3B exp`)]) >= (length(df$`Max of 3A and 3B exp`)-2) | length(df$`Max of 2 and 13`[is.na(df$`Max of 2 and 13`)]) >= (length(df$`Max of 2 and 13`)-2) | i == 5) {
    print("skipping, not enough samples")
  } else {
    test <- cor.test(df$`Max of 2 and 13`, df$`Max of 3A and 3B exp`, method='spearman')
    cor_mat[i,7] <- as.numeric(test$estimate)
    cor_mat[i,8] <- as.numeric(test$p.value)
  }
}
write.table(cor_mat, file='correlation.txt', sep='\t')

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)