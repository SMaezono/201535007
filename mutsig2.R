# mutsig2.R

# Calculate mutational signature using 4 methods

# Windows directory: /Users/jiachen
# MacOS directory: /Users/jiachen

#####################################################################################

t_start <- Sys.time()
t_last <- Sys.time()

library(matrixStats)
library(ggplot2)
library(deconstructSigs)

# Number of controls
ctrlcount <- 100

# Set working directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/"
setwd(root)

# Import mutation data
print('Importing Alexandrov data')
cancers_alex <- as.character(read.table(file='cancer_list_alex.txt', sep='\t', header=F)[,1])
cancers <- as.character(read.table(file='cancer_list.txt', sep='\t', header=F)[,1])
############################################################# EXCLUDE COAD AND UCEC (too few samples)
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
exp_all <- vector('list', length=length(cancers))
med <- c()

for(i in 1:length(cancers)){
  
  print(paste('Importing', cancers[i], 'data'))
  target_file <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/apobec/', cancers[i], '.txt', sep='')
  exp_a_t[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  target_file <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/hypoxia/', cancers[i], '.txt', sep='')
  exp_h_t[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  target_file <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/', cancers[i], '_t.txt', sep='')
  input <- read.table(file=target_file, sep='\t', header=T, row.names=NULL)
  
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
  
  # Fix names (exp)
  rn <- input[,1]
  rn <- make.unique(rn, sep='.')
  
  exp <- input[,-1]
  rownames(exp) <- rn
  
  names <- colnames(exp)
  for(j in 1:length(names)){
    split <- strsplit(names[j], '\\.')
    names[j] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(exp) <- names
  
  # Median
  exp_med <- as.numeric(as.matrix(exp))
  exp_med[is.na(exp_med)] <- 0
  med <- c(med, median(exp_med, na.rm=T))
  # Filter
  exp_out <- exp[(rowCounts(!is.na(exp))/ncol(exp))>0.5,]
  exp_out <- exp_out[rowMeans(exp, na.rm=T)>=med[i],]
  exp_all[[i]] <- exp_out
}

# Matrices with correlation scores
all_matrices <- vector('list', length=4)

# Counts analysis ####################################################################################

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/output_counts/"
setwd(root)

print(Sys.time()-t_last)
t_last <- Sys.time()
print('Analyzing mutational signatures using counts (log)')

# Get common samples
common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(mut_sums_apo), intersect(colnames(exp_a_t[[i]]), colnames(exp_h_t[[i]])))
}

# Build matrices
matrices <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  mat <- mut_sums_apo_log[common[[i]],1]
  mat <- cbind(mat, t(exp_a_t[[i]]['APOBEC3A',common[[i]]]), t(exp_a_t[[i]]['APOBEC3B',common[[i]]]))
  mat <- cbind(mat, colMeans(exp_h_t[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c('APOBEC mutation count', 'APOBEC3A expression', 'APOBEC3B expression', 'Hypoxia signature expression')
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers[i], '_mut_summary.txt', sep=''), sep='\t')
}

# Plot
for(i in 1:length(cancers)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`APOBEC mutation count`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`APOBEC mutation count`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`APOBEC mutation count`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
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
  test <- cor.test(df$`APOBEC mutation count`, df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  test <- cor.test(df$`APOBEC mutation count`, df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  test <- cor.test(df$`APOBEC mutation count`, df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
}
write.table(cor_mat, file='correlation.txt', sep='\t')
all_matrices[[1]] <- cor_mat

# Control
control_mat <- matrix(nrow=ctrlcount, ncol=length(cancers))
control.gene <- c()
control.cancer <- c()
control.cor <- c()
control.pv <- c()
control.cor.f <- c()
control.method <- c()

common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(mut_sums_apo), colnames(exp_all[[i]]))
}
for(i in 1:length(cancers)){
  this.exp <- exp_all[[i]]
  this.exp <- this.exp[,common[[i]]]
  for(j in 1:ctrlcount){
    this.gene <- NULL
    ok <- F
    while(!ok){
      ok <- F
      this.gene <- rownames(this.exp)[sample(1:nrow(this.exp), 1)]
      if(!is.na(this.gene)){ # filter out garbage genes
        if(substr(this.gene, 1, 2)!='NA'){
          if(!(this.gene %in% control_mat[,i])){ # check if unique
            ok <- T
          }
        }
      }
    }
    control_mat[j,i] <- this.gene
    control.gene <- c(control.gene, 'Control')
    control.cancer <- c(control.cancer, cancers[i])
    control.method <- c(control.method, 'Counts')
    this.test <- cor.test(as.numeric(this.exp[this.gene,]), as.numeric(mut_sums_apo_log[common[[i]],1]))
    this.cor <- this.test$estimate
    this.pv <- this.test$p.value
    this.cor.f <- this.cor
    if(this.pv>=0.05){
      this.cor.f <- NA
    }
    control.cor <- c(control.cor, this.cor)
    control.pv <- c(control.pv, this.pv)
    control.cor.f <- c(control.cor.f, this.cor.f)
  }
}
df_counts <- data.frame(control.gene, control.cancer, control.cor, control.pv, control.cor.f, control.method)

for(i in 1:length(cancers)){
  actual <- as.data.frame(t(c('Hypoxia', cancers[i], cor_mat[i,5], cor_mat[i,6], cor_mat[i,5], 'Counts')))
  colnames(actual) <- colnames(df_counts)
  a3a <- as.data.frame(t(c('A3A', cancers[i], cor_mat[i,1], cor_mat[i,2], cor_mat[i,1], 'Counts')))
  colnames(a3a) <- colnames(df_counts)
  a3b <- as.data.frame(t(c('A3B', cancers[i], cor_mat[i,3], cor_mat[i,4], cor_mat[i,3], 'Counts')))
  colnames(a3b) <- colnames(df_counts)
  df_counts <- rbind(df_counts, actual, a3a, a3b)
}

# Enrichment analysis ################################################################################

# Set working directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/"
setwd(root)

# Import Roberts data
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Importing Roberts et al. data')
cancers_c <- as.character(read.table(file='cancers_list_roberts.txt', sep='\t', header=F)[,1])
roberts <- read.table(file='mut_sig_context_roberts.txt', sep='\t', header=TRUE, row.names=NULL)
names <- as.character(roberts[,1])
for(i in 1:length(names)){
  split <- strsplit(names[i], '-')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
}
rownames(roberts) <- names

# Import expression data
exp_a_t_r <- vector('list', length=length(cancers_c))
exp_h_t_r <- vector('list', length=length(cancers_c))
exp_all_r <- vector('list', length=length(cancers_c))
for(i in 1:length(cancers_c)){
  print(paste('Importing', cancers_c[i], 'data'))
  target_file <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/apobec/', cancers_c[i], '.txt', sep='')
  exp_a_t_r[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  target_file <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/hypoxia/', cancers_c[i], '.txt', sep='')
  exp_h_t_r[[i]] <- read.table(file=target_file, sep='\t', header=T, row.names=1)
  # Fix names
  names <- colnames(exp_a_t_r[[i]])
  for(j in 1:length(names)){
    split <- strsplit(names[j], '\\.')
    names[j] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(exp_a_t_r[[i]]) <- names
  names <- colnames(exp_h_t_r[[i]])
  for(j in 1:length(names)){
    split <- strsplit(names[j], '\\.')
    names[j] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(exp_h_t_r[[i]]) <- names
  
  target_file <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/', cancers_c[i], '_t.txt', sep='')
  input <- read.table(file=target_file, sep='\t', header=T, row.names=NULL)
  
  # Fix names (exp)
  rn <- input[,1]
  rn <- make.unique(rn, sep='.')
  
  exp <- input[,-1]
  rownames(exp) <- rn
  
  names <- colnames(exp)
  for(j in 1:length(names)){
    split <- strsplit(names[j], '\\.')
    names[j] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(exp) <- names
  
  # Median
  exp_med <- as.numeric(as.matrix(exp))
  exp_med[is.na(exp_med)] <- 0
  med <- c(med, median(exp_med, na.rm=T))
  # Filter
  exp_out <- exp[(rowCounts(!is.na(exp))/ncol(exp))>0.5,]
  exp_out <- exp_out[rowMeans(exp, na.rm=T)>=med[i],]
  exp_all_r[[i]] <- exp_out
}

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/output_enrichment/"
setwd(root)

# Get common samples
common <- vector('list', length=length(cancers_c))
for(i in 1:length(cancers_c)){
  common[[i]] <- intersect(rownames(roberts), intersect(colnames(exp_a_t_r[[i]]), colnames(exp_h_t_r[[i]])))
}

# Build matrices
matrices <- vector('list', length=length(cancers_c))
for(i in 1:length(cancers_c)){
  mat <- roberts[common[[i]],'APOBEC_enrich40']
  mat <- cbind(mat, t(exp_a_t_r[[i]]['APOBEC3A',common[[i]]]), t(exp_a_t_r[[i]]['APOBEC3B',common[[i]]]))
  mat <- cbind(mat, colMeans(exp_h_t_r[[i]][,common[[i]]], na.rm=T))
  colnames(mat) <- c('APOBEC signature enrichment', 'APOBEC3A expression', 'APOBEC3B expression', 'Hypoxia signature expression')
  matrices[[i]] <- mat
  write.table(mat, file=paste(cancers_c[i], '_enrichment.txt', sep=''), sep='\t')
}

# Plot
for(i in 1:length(cancers_c)){
  df <- as.data.frame(matrices[[i]])
  sp1 <- ggplot(df, aes(x=`APOBEC3A expression`, y=`APOBEC signature enrichment`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`APOBEC signature enrichment`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`APOBEC signature enrichment`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  pdf(file=paste(cancers_c[i], '_A3A-enr.pdf', sep=''), width=8, height=8)
  print(sp1)
  dev.off()
  pdf(file=paste(cancers_c[i], '_A3B-enr.pdf', sep=''), width=8, height=8)
  print(sp2)
  dev.off()
  pdf(file=paste(cancers_c[i], '_hyp-enr.pdf', sep=''), width=8, height=8)
  print(sp3)
  dev.off()
}

# Correlation
cor_mat <- matrix(nrow=length(cancers_c), ncol=6)
rownames(cor_mat) <- cancers_c
colnames(cor_mat) <- c('A3A', 'A3A pv', 'A3B', 'A3B pv', 'hyp', 'hyp pv')
for(i in 1:length(cancers_c)){
  df <- as.data.frame(matrices[[i]])
  test <- cor.test(df$`APOBEC signature enrichment`, df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  test <- cor.test(df$`APOBEC signature enrichment`, df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  test <- cor.test(df$`APOBEC signature enrichment`, df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
}
write.table(cor_mat, file='correlation.txt', sep='\t')
all_matrices[[2]] <- cor_mat

# Control
control_mat <- matrix(nrow=ctrlcount, ncol=length(cancers_c))
control.gene <- c()
control.cancer <- c()
control.cor <- c()
control.pv <- c()
control.cor.f <- c()
control.method <- c()

common <- vector('list', length=length(cancers))
for(i in 1:length(cancers_c)){
  common[[i]] <- intersect(rownames(roberts), colnames(exp_all_r[[i]]))
}
for(i in 1:length(cancers_c)){
  this.exp <- exp_all_r[[i]]
  this.exp <- this.exp[,common[[i]]]
  for(j in 1:ctrlcount){
    this.gene <- NULL
    ok <- F
    while(!ok){
      ok <- F
      this.gene <- rownames(this.exp)[sample(1:nrow(this.exp), 1)]
      if(!is.na(this.gene)){ # filter out garbage genes
        if(substr(this.gene, 1, 2)!='NA'){
          if(!(this.gene %in% control_mat[,i])){ # check if unique
            ok <- T
          }
        }
      }
    }
    control_mat[j,i] <- this.gene
    control.gene <- c(control.gene, 'Control')
    control.cancer <- c(control.cancer, cancers_c[i])
    control.method <- c(control.method, 'Enrichment')
    this.test <- cor.test(as.numeric(this.exp[this.gene,]), as.numeric(roberts[common[[i]],1]))
    this.cor <- this.test$estimate
    this.pv <- this.test$p.value
    this.cor.f <- this.cor
    if(this.pv>=0.05){
      this.cor.f <- NA
    }
    control.cor <- c(control.cor, this.cor)
    control.pv <- c(control.pv, this.pv)
    control.cor.f <- c(control.cor.f, this.cor.f)
  }
}
df_enrichment <- data.frame(control.gene, control.cancer, control.cor, control.pv, control.cor.f, control.method)

for(i in 1:length(cancers_c)){
  actual <- as.data.frame(t(c('Hypoxia', cancers_c[i], cor_mat[i,5], cor_mat[i,6], cor_mat[i,5], 'Enrichment')))
  colnames(actual) <- colnames(df_enrichment)
  a3a <- as.data.frame(t(c('A3A', cancers_c[i], cor_mat[i,1], cor_mat[i,2], cor_mat[i,1], 'Enrichment')))
  colnames(a3a) <- colnames(df_enrichment)
  a3b <- as.data.frame(t(c('A3B', cancers_c[i], cor_mat[i,3], cor_mat[i,4], cor_mat[i,3], 'Enrichment')))
  colnames(a3b) <- colnames(df_enrichment)
  df_enrichment <- rbind(df_enrichment, actual, a3a, a3b)
}

# Prepare for deconstructSigs (fractionize) ############################################################

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/"
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

# Filter for only TCGA data
mut_counts <- mut_counts[substr(rownames(mut_counts),1,4)=='TCGA',]

# DCS all signatures #################################################################################

# dcs_frac <- whichSignatures(mut_frac, rownames(mut_frac))
# dcs_counts <- whichSignatures(mut_counts, rownames(mut_counts), contexts.needed=T, tri.counts.method='exome')

print(Sys.time()-t_last)
t_last <- Sys.time()

if(file.exists('dcs_all.txt')){
  print('Reading from pre-analyzed deconstructSigs')
  dcs_counts <- read.table(file='dcs_all.txt', header=T, row.names=1, sep='\t')
}else{
  print('Analyzing mutational signatures using deconstructSigs')
  dcs_counts <- as.data.frame(matrix(nrow=0, ncol=30))
  
  # sigpath <- 'signatures_apobec.txt'
  # sigpath <- read.table('signatures_probabilities.txt', header=T, row.names=NULL, sep='\t')
  # sigpath <- t(sigpath)
  # colnames(sigpath) <- as.character(sigpath['Somatic.Mutation.Type',])
  # sigpath <- sigpath[4:33,]
  
  sigpath <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/signatures alt.txt"
  
  for(i in 1:nrow(mut_counts)){
    print(paste('Analyzing sample', i, 'of', nrow(mut_counts)))
    dcs_counts <- rbind(
      dcs_counts,
      whichSignatures(mut_counts,
                      rownames(mut_counts)[i],
                      signatures.ref = sigpath,
                      contexts.needed=T,
                      tri.counts.method='exome'
      )$weights
    )
  }
  # Output matrix
  write.table(dcs_counts, file='dcs_all.txt', sep='\t')
}

# Subset for signatures 2 and 13
dcs_apo <- dcs_counts[,c(2,13)]

# Get common samples
common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(dcs_apo), intersect(colnames(exp_a_t[[i]]), colnames(exp_h_t[[i]])))
}

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/output_dcs/"
setwd(root)

# Plot out DCS
plot.cancer <- c()
plot.value <- c()
plot.sig <- c()
plot.matrices <- vector('list', length=length(common))
for(c in 1:length(cancers)){
  cur.mat <- dcs_counts[common[[c]],]
  for(i in 1:nrow(cur.mat)){
    for(j in 1:ncol(cur.mat)){
      plot.cancer <- c(plot.cancer, cancers[c])
      plot.value <- c(plot.value, cur.mat[i,j])
      if(nchar(colnames(cur.mat)[j])==11){
        split <- strsplit(colnames(cur.mat)[j], '\\.')
        sig_name <- paste(split[[1]][1], '.0', split[[1]][2], sep='')
        cur.col.names <- colnames(cur.mat)
        cur.col.names[j] <- sig_name
        colnames(cur.mat) <- cur.col.names
      }else{
        sig_name <- colnames(cur.mat)[j]
      }
      plot.sig <- c(plot.sig, sig_name)
    }
  }
  plot.matrices[[c]] <- cur.mat
}
plot.df <- data.frame(plot.cancer, plot.value, plot.sig)
for(c in 1:length(plot.matrices)){
  plot_out <- ggplot(plot.df[plot.df$plot.cancer==cancers[c],], aes(x=plot.sig, y=plot.value)) +
    geom_jitter(alpha=0.1) +
    ggtitle(paste('DCS (all) -', cancers[c])) +
    xlab("Signature") + ylab("DCS score") +
    ylim(0,1) +
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=45, hjust=1))
  pdf(file=paste(cancers[c], '_DCS_all.pdf', sep=''), width=10, height=4)
  print(plot_out)
  dev.off()
  
  # percentage of non-null values
  nnp <- matrix(nrow=ncol(plot.matrices[[c]]), ncol=1)
  rownames(nnp) <- colnames(plot.matrices[[c]])
  colnames(nnp) <- '% non-null'
  for(j in 1:ncol(plot.matrices[[c]])){
    nnp[j,1] <- 
      count(plot.df[plot.df$plot.cancer==cancers[c] & plot.df$plot.sig==colnames(plot.matrices[[c]])[j],]$plot.value>0, na.rm=T)/
      length(plot.df[plot.df$plot.cancer==cancers[c] & plot.df$plot.sig==colnames(plot.matrices[[c]])[j],]$plot.value)
    nnp[j,1] <- nnp[j,1]*100
  }
  write.table(nnp, file=paste(cancers[c], '_nnp.txt', sep=''), sep='\t')
}

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
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`Max of 2 and 13`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`Max of 2 and 13`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp4 <- ggplot(df, aes(x=`Max of 3A and 3B exp`, y=`Max of 2 and 13`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
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
  test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  test <- cor.test(df$`Max of 2 and 13`, df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
  test <- cor.test(df$`Max of 2 and 13`, df$`Max of 3A and 3B exp`, method='spearman')
  cor_mat[i,7] <- as.numeric(test$estimate)
  cor_mat[i,8] <- as.numeric(test$p.value)
}
write.table(cor_mat, file='correlation.txt', sep='\t')
all_matrices[[3]] <- cor_mat

# Control
control_mat <- matrix(nrow=ctrlcount, ncol=length(cancers))
control.gene <- c()
control.cancer <- c()
control.cor <- c()
control.pv <- c()
control.cor.f <- c()
control.method <- c()

common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(dcs_apo), colnames(exp_all[[i]]))
}
for(i in 1:length(cancers)){
  this.exp <- exp_all[[i]]
  this.exp <- this.exp[,common[[i]]]
  for(j in 1:ctrlcount){
    this.gene <- NULL
    ok <- F
    while(!ok){
      ok <- F
      this.gene <- rownames(this.exp)[sample(1:nrow(this.exp), 1)]
      if(!is.na(this.gene)){ # filter out garbage genes
        if(substr(this.gene, 1, 2)!='NA'){
          if(!(this.gene %in% control_mat[,i])){ # check if unique
            ok <- T
          }
        }
      }
    }
    control_mat[j,i] <- this.gene
    control.gene <- c(control.gene, 'Control')
    control.cancer <- c(control.cancer, cancers[i])
    control.method <- c(control.method, 'DCS all')
    this.test <- cor.test(as.numeric(this.exp[this.gene,]), as.numeric(dcs_apo[common[[i]],1]))
    this.cor <- this.test$estimate
    this.pv <- this.test$p.value
    this.cor.f <- this.cor
    if(this.pv>=0.05){
      this.cor.f <- NA
    }
    control.cor <- c(control.cor, this.cor)
    control.pv <- c(control.pv, this.pv)
    control.cor.f <- c(control.cor.f, this.cor.f)
  }
}
df_dcs_all <- data.frame(control.gene, control.cancer, control.cor, control.pv, control.cor.f, control.method)

for(i in 1:length(cancers)){
  actual <- as.data.frame(t(c('Hypoxia', cancers[i], cor_mat[i,5], cor_mat[i,6], cor_mat[i,5], 'DCS all')))
  colnames(actual) <- colnames(df_dcs_all)
  a3a <- as.data.frame(t(c('A3A', cancers[i], cor_mat[i,1], cor_mat[i,2], cor_mat[i,1], 'DCS all')))
  colnames(a3a) <- colnames(df_dcs_all)
  a3b <- as.data.frame(t(c('A3B', cancers[i], cor_mat[i,3], cor_mat[i,4], cor_mat[i,3], 'DCS all')))
  colnames(a3b) <- colnames(df_dcs_all)
  df_dcs_all <- rbind(df_dcs_all, actual, a3a, a3b)
}

# DCS curated #################################################################################

# dcs_frac <- whichSignatures(mut_frac, rownames(mut_frac))
# dcs_counts <- whichSignatures(mut_counts, rownames(mut_counts), contexts.needed=T, tri.counts.method='exome')

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/"
setwd(root)

print(Sys.time()-t_last)
t_last <- Sys.time()

print('Analyzing mutational signatures using deconstructSigs with curated signatures')

dcs_counts_all <- vector('list', length=length(cancers))
dcs_counts <- as.data.frame(matrix(nrow=0, ncol=2))

# IF WE NEED TO INPUT OUR OWN MATRIX
# sigpath <- 'signatures_apobec.txt'
# sigpath <- read.table('signatures_probabilities.txt', header=T, row.names=NULL, sep='\t')
# sigpath <- t(sigpath)
# colnames(sigpath) <- as.character(sigpath['Somatic.Mutation.Type',])
# sigpath <- sigpath[4:33,]

# sigpath <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/signatures alt.txt"

# Subset matrix into different cancer types
mut_counts_pc <- vector('list', length=length(cancers))
for(j in 1:length(cancers)){
  mut_counts_pc[[j]] <- mut_counts[intersect(rownames(mut_counts), intersect(colnames(exp_a_t[[j]]), colnames(exp_h_t[[j]]))),]
}

for(j in 1:length(cancers)){
  sigpath <- paste('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/signatures_', cancers[j], '.txt', sep='')
  dcs_counts_pc <- as.data.frame(matrix(nrow=0, ncol=0))
  if(file.exists(paste('dcs_curated_', cancers[j], '.txt', sep=''))){
    print(paste(cancers[j], 'Reading from pre-analyzed deconstructSigs'))
    dcs_counts_pc <- read.table(file=paste('dcs_curated_', cancers[j], '.txt', sep=''), header=T, row.names=1, sep='\t')
    dcs_counts_all[[j]] <- dcs_counts_pc
  }else{
    print(paste(cancers[j], 'Analyzing using deconstructSigs'))
    for(i in 1:nrow(mut_counts_pc[[j]])){
      print(paste('Analyzing', cancers[j], 'sample', i, 'of', nrow(mut_counts_pc[[j]])))
      this.dcs <- whichSignatures(mut_counts_pc[[j]],
                                  rownames(mut_counts_pc[[j]])[i],
                                  signatures.ref = sigpath,
                                  contexts.needed=T,
                                  tri.counts.method='exome'
                                  )$weights
      dcs_counts_pc <- rbind(dcs_counts_pc,this.dcs)
    }
    
    dcs_counts_all[[j]] <- dcs_counts_pc
    
    # Output matrix
    write.table(dcs_counts_pc, file=paste('dcs_curated_', cancers[j], '.txt', sep=''), sep='\t')
  }
}
print('Read complete')

# Plot out signatures
plot.cancer <- c()
plot.value <- c()
plot.sig <- c()
for(c in 1:length(dcs_counts_all)){
  for(i in 1:nrow(dcs_counts_all[[c]])){
    for(j in 1:ncol(dcs_counts_all[[c]])){
      plot.cancer <- c(plot.cancer, cancers[c])
      plot.value <- c(plot.value, dcs_counts_all[[c]][i,j])
      if(nchar(colnames(dcs_counts_all[[c]])[j])==11){
        split <- strsplit(colnames(dcs_counts_all[[c]])[j], '\\.')
        sig_name <- paste(split[[1]][1], '.0', split[[1]][2], sep='')
        cur.col.names <- colnames(dcs_counts_all[[c]])
        cur.col.names[j] <- sig_name
        colnames(dcs_counts_all[[c]]) <- cur.col.names
      }else{
        sig_name <- colnames(dcs_counts_all[[c]])[j]
      }
      plot.sig <- c(plot.sig, sig_name)
    }
  }
}
plot.df <- data.frame(plot.cancer, plot.value, plot.sig)
for(c in 1:length(dcs_counts_all)){
  plot_out <- ggplot(plot.df[plot.df$plot.cancer==cancers[c],], aes(x=plot.sig, y=plot.value)) +
    geom_jitter(alpha=0.1) +
    ggtitle(paste('DCS (curated) -', cancers[c])) +
    xlab("Signature") + ylab("DCS score") +
    ylim(0,1) +
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=45, hjust=1))
  pdf(file=paste('output_dcs_curated/', cancers[c], '_DCS_cur.pdf', sep=''), width=6, height=4)
  print(plot_out)
  dev.off()
  
  # percentage of non-null values
  nnp <- matrix(nrow=ncol(dcs_counts_all[[c]]), ncol=1)
  rownames(nnp) <- colnames(dcs_counts_all[[c]])
  colnames(nnp) <- '% non-null'
  for(j in 1:ncol(dcs_counts_all[[c]])){
    nnp[j,1] <- 
      count(plot.df[plot.df$plot.cancer==cancers[c] & plot.df$plot.sig==colnames(dcs_counts_all[[c]])[j],]$plot.value>0, na.rm=T)/
      length(plot.df[plot.df$plot.cancer==cancers[c] & plot.df$plot.sig==colnames(dcs_counts_all[[c]])[j],]$plot.value)
    nnp[j,1] <- nnp[j,1]*100
  }
  write.table(nnp, file=paste('output_dcs_curated/', cancers[c], '_nnp.txt', sep=''), sep='\t')
}

# Subset for signatures 2 and 13 and combine
for(i in 1:length(cancers)){
  this.dcs.counts <- dcs_counts_all[[i]][,c('Signature.02', 'Signature.13')]
  dcs_counts <- rbind(dcs_counts, this.dcs.counts)
}
dcs_apo <- dcs_counts

# Get common samples
common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(dcs_apo), intersect(colnames(exp_a_t[[i]]), colnames(exp_h_t[[i]])))
}

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/output_dcs_curated/"
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
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp2 <- ggplot(df, aes(x=`APOBEC3B expression`, y=`Max of 2 and 13`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp3 <- ggplot(df, aes(x=`Hypoxia signature expression`, y=`Max of 2 and 13`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
  sp4 <- ggplot(df, aes(x=`Max of 3A and 3B exp`, y=`Max of 2 and 13`)) +
    geom_point(alpha=0.6, pch=16) + theme(text = element_text(size=24))
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
  test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3A expression`, method='spearman')
  cor_mat[i,1] <- as.numeric(test$estimate)
  cor_mat[i,2] <- as.numeric(test$p.value)
  test <- cor.test(df$`Max of 2 and 13`, df$`APOBEC3B expression`, method='spearman')
  cor_mat[i,3] <- as.numeric(test$estimate)
  cor_mat[i,4] <- as.numeric(test$p.value)
  test <- cor.test(df$`Max of 2 and 13`, df$`Hypoxia signature expression`, method='spearman')
  cor_mat[i,5] <- as.numeric(test$estimate)
  cor_mat[i,6] <- as.numeric(test$p.value)
  test <- cor.test(df$`Max of 2 and 13`, df$`Max of 3A and 3B exp`, method='spearman')
  cor_mat[i,7] <- as.numeric(test$estimate)
  cor_mat[i,8] <- as.numeric(test$p.value)
}
write.table(cor_mat, file='correlation.txt', sep='\t')
all_matrices[[4]] <- cor_mat

# Control
control_mat <- matrix(nrow=ctrlcount, ncol=length(cancers))
control.gene <- c()
control.cancer <- c()
control.cor <- c()
control.pv <- c()
control.cor.f <- c()
control.method <- c()

common <- vector('list', length=length(cancers))
for(i in 1:length(cancers)){
  common[[i]] <- intersect(rownames(dcs_apo), colnames(exp_all[[i]]))
}
for(i in 1:length(cancers)){
  this.exp <- exp_all[[i]]
  this.exp <- this.exp[,common[[i]]]
  for(j in 1:ctrlcount){
    this.gene <- NULL
    ok <- F
    while(!ok){
      ok <- F
      this.gene <- rownames(this.exp)[sample(1:nrow(this.exp), 1)]
      if(!is.na(this.gene)){ # filter out garbage genes
        if(substr(this.gene, 1, 2)!='NA'){
          if(!(this.gene %in% control_mat[,i])){ # check if unique
            ok <- T
          }
        }
      }
    }
    control_mat[j,i] <- this.gene
    control.gene <- c(control.gene, 'Control')
    control.cancer <- c(control.cancer, cancers[i])
    control.method <- c(control.method, 'DCS curated')
    this.test <- cor.test(as.numeric(this.exp[this.gene,]), as.numeric(dcs_apo[common[[i]],1]))
    this.cor <- this.test$estimate
    this.pv <- this.test$p.value
    this.cor.f <- this.cor
    if(this.pv>=0.05){
      this.cor.f <- NA
    }
    control.cor <- c(control.cor, this.cor)
    control.pv <- c(control.pv, this.pv)
    control.cor.f <- c(control.cor.f, this.cor.f)
  }
}
df_dcs_cur <- data.frame(control.gene, control.cancer, control.cor, control.pv, control.cor.f, control.method)

for(i in 1:length(cancers)){
  actual <- as.data.frame(t(c('Hypoxia', cancers[i], cor_mat[i,5], cor_mat[i,6], cor_mat[i,5], 'DCS curated')))
  colnames(actual) <- colnames(df_dcs_cur)
  a3a <- as.data.frame(t(c('A3A', cancers[i], cor_mat[i,1], cor_mat[i,2], cor_mat[i,1], 'DCS curated')))
  colnames(a3a) <- colnames(df_dcs_cur)
  a3b <- as.data.frame(t(c('A3B', cancers[i], cor_mat[i,3], cor_mat[i,4], cor_mat[i,3], 'DCS curated')))
  colnames(a3b) <- colnames(df_dcs_cur)
  df_dcs_cur <- rbind(df_dcs_cur, actual, a3a, a3b)
}

# Controls #################################################################################

print(Sys.time()-t_last)
t_last <- Sys.time()

print('Plotting controls')

df_counts$control.cor <- as.numeric(df_counts$control.cor)
df_counts$control.pv <- as.numeric(df_counts$control.pv)
df_counts$control.cor.f <- as.numeric(df_counts$control.cor.f)

df_enrichment$control.cor <- as.numeric(df_enrichment$control.cor)
df_enrichment$control.pv <- as.numeric(df_enrichment$control.pv)
df_enrichment$control.cor.f <- as.numeric(df_enrichment$control.cor.f)

df_dcs_all$control.cor <- as.numeric(df_dcs_all$control.cor)
df_dcs_all$control.pv <- as.numeric(df_dcs_all$control.pv)
df_dcs_all$control.cor.f <- as.numeric(df_dcs_all$control.cor.f)

df_dcs_cur$control.cor <- as.numeric(df_dcs_cur$control.cor)
df_dcs_cur$control.pv <- as.numeric(df_dcs_cur$control.pv)
df_dcs_cur$control.cor.f <- as.numeric(df_dcs_cur$control.cor.f)

df_all <- as.data.frame(rbind(df_counts, df_enrichment, df_dcs_all, df_dcs_cur))

root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/mutsig/output_controls/"
setwd(root)

control_plot_counts <- ggplot() +
  geom_boxplot(data=df_counts[df_counts$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_counts[df_counts$control.gene=='Hypoxia',], aes(x=control.cancer, y=control.cor), color='red') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - Counts') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
control_plot_enrichment <- ggplot() +
  geom_boxplot(data=df_enrichment[df_enrichment$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_enrichment[df_enrichment$control.gene=='Hypoxia',], aes(x=control.cancer, y=control.cor), color='red') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - Enrichment') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
control_plot_dcs_all <- ggplot() +
  geom_boxplot(data=df_dcs_all[df_dcs_all$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_dcs_all[df_dcs_all$control.gene=='Hypoxia',], aes(x=control.cancer, y=control.cor), color='red') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - DCS all') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
control_plot_dcs_cur <- ggplot() +
  geom_boxplot(data=df_dcs_cur[df_dcs_cur$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_dcs_cur[df_dcs_cur$control.gene=='Hypoxia',], aes(x=control.cancer, y=control.cor), color='red') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - DCS curated') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))

pdf(file=paste('grouped by method/counts.pdf', sep=''), width=4, height=4)
print(control_plot_counts)
dev.off()
pdf(file=paste('grouped by method/enrichment.pdf', sep=''), width=4, height=4)
print(control_plot_enrichment)
dev.off()
pdf(file=paste('grouped by method/dcs_all.pdf', sep=''), width=4, height=4)
print(control_plot_dcs_all)
dev.off()
pdf(file=paste('grouped by method/dcs_cur.pdf', sep=''), width=4, height=4)
print(control_plot_dcs_cur)
dev.off()

control_apo_counts <- ggplot() +
  geom_boxplot(data=df_counts[df_counts$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_counts[df_counts$control.gene=='A3A',], aes(x=control.cancer, y=control.cor), color='green') +
  geom_point(data=df_counts[df_counts$control.gene=='A3B',], aes(x=control.cancer, y=control.cor), color='blue') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - Counts') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
control_apo_enrichment <- ggplot() +
  geom_boxplot(data=df_enrichment[df_enrichment$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_enrichment[df_enrichment$control.gene=='A3A',], aes(x=control.cancer, y=control.cor), color='green') +
  geom_point(data=df_enrichment[df_enrichment$control.gene=='A3B',], aes(x=control.cancer, y=control.cor), color='blue') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - Enrichment') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
control_apo_dcs_all <- ggplot() +
  geom_boxplot(data=df_dcs_all[df_dcs_all$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_dcs_all[df_dcs_all$control.gene=='A3A',], aes(x=control.cancer, y=control.cor), color='green') +
  geom_point(data=df_dcs_all[df_dcs_all$control.gene=='A3B',], aes(x=control.cancer, y=control.cor), color='blue') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - DCS all') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
control_apo_dcs_cur <- ggplot() +
  geom_boxplot(data=df_dcs_cur[df_dcs_cur$control.gene=='Control',], aes(x=control.cancer, y=control.cor)) +
  geom_point(data=df_dcs_cur[df_dcs_cur$control.gene=='A3A',], aes(x=control.cancer, y=control.cor), color='green') +
  geom_point(data=df_dcs_cur[df_dcs_cur$control.gene=='A3B',], aes(x=control.cancer, y=control.cor), color='blue') +
  ylim(-0.5,0.5) +
  ggtitle('Control plot - DCS curated') +
  xlab("Cancer") + ylab("Spearman's rho") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))

pdf(file=paste('apobec/counts.pdf', sep=''), width=4, height=4)
print(control_apo_counts)
dev.off()
pdf(file=paste('apobec/enrichment.pdf', sep=''), width=4, height=4)
print(control_apo_enrichment)
dev.off()
pdf(file=paste('apobec/dcs_all.pdf', sep=''), width=4, height=4)
print(control_apo_dcs_all)
dev.off()
pdf(file=paste('apobec/dcs_cur.pdf', sep=''), width=4, height=4)
print(control_apo_dcs_cur)
dev.off()

for(i in 1:length(cancers_c)){
  plot_control <- ggplot() +
    geom_boxplot(data=df_all[(df_all$control.cancer==cancers_c[i] & df_all$control.gene!='Hypoxia'),], aes(x=control.method, y=control.cor)) +
    geom_point(data=df_all[(df_all$control.cancer==cancers_c[i] & df_all$control.gene=='Hypoxia'),], aes(x=control.method, y=control.cor), color='red') +
    ylim(-0.5,0.5) +
    ggtitle(paste('Control plot -', cancers_c[i])) +
    xlab("Method") + ylab("Spearman's rho") +
    theme(text = element_text(size=16), axis.text.x = element_text(angle=45, hjust=1))
  pdf(file=paste('grouped by cancer/', cancers_c[i], '.pdf', sep=''), width=4, height=4)
  print(plot_control)
  dev.off()
}

z_counts <- matrix(nrow=length(cancers), ncol=3)
rownames(z_counts) <- cancers
colnames(z_counts) <- c('A3A', 'A3B', 'Hypoxia')
for(i in 1:length(cancers)){
  dist <- as.numeric(df_counts[(df_counts$control.gene=='Control' & df_counts$control.cancer==cancers[i]),]$control.cor)
  val_a <- as.numeric(df_counts[(df_counts$control.gene=='A3A' & df_counts$control.cancer==cancers[i]),]$control.cor)
  val_b <- as.numeric(df_counts[(df_counts$control.gene=='A3B' & df_counts$control.cancer==cancers[i]),]$control.cor)
  val_h <- as.numeric(df_counts[(df_counts$control.gene=='Hypoxia' & df_counts$control.cancer==cancers[i]),]$control.cor)
  z_counts[i,1] <- (val_a-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_counts[i,2] <- (val_b-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_counts[i,3] <- (val_h-mean(dist, na.rm=T))/sd(dist, na.rm=T)
}
write.table(z_counts, file='z_counts.txt', sep='\t')

z_enrichment <- matrix(nrow=length(cancers_c), ncol=3)
rownames(z_enrichment) <- cancers_c
colnames(z_enrichment) <- c('A3A', 'A3B', 'Hypoxia')
for(i in 1:length(cancers_c)){
  dist <- as.numeric(df_enrichment[(df_enrichment$control.gene=='Control' & df_enrichment$control.cancer==cancers_c[i]),]$control.cor)
  val_a <- as.numeric(df_enrichment[(df_enrichment$control.gene=='A3A' & df_enrichment$control.cancer==cancers_c[i]),]$control.cor)
  val_b <- as.numeric(df_enrichment[(df_enrichment$control.gene=='A3B' & df_enrichment$control.cancer==cancers_c[i]),]$control.cor)
  val_h <- as.numeric(df_enrichment[(df_enrichment$control.gene=='Hypoxia' & df_enrichment$control.cancer==cancers_c[i]),]$control.cor)
  z_enrichment[i,1] <- (val_a-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_enrichment[i,2] <- (val_b-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_enrichment[i,3] <- (val_h-mean(dist, na.rm=T))/sd(dist, na.rm=T)
}
write.table(z_enrichment, file='z_enrichment.txt', sep='\t')

z_dcs_all <- matrix(nrow=length(cancers), ncol=3)
rownames(z_dcs_all) <- cancers
colnames(z_dcs_all) <- c('A3A', 'A3B', 'Hypoxia')
for(i in 1:length(cancers)){
  dist <- as.numeric(df_dcs_all[(df_dcs_all$control.gene=='Control' & df_dcs_all$control.cancer==cancers[i]),]$control.cor)
  val_a <- as.numeric(df_dcs_all[(df_dcs_all$control.gene=='A3A' & df_dcs_all$control.cancer==cancers[i]),]$control.cor)
  val_b <- as.numeric(df_dcs_all[(df_dcs_all$control.gene=='A3B' & df_dcs_all$control.cancer==cancers[i]),]$control.cor)
  val_h <- as.numeric(df_dcs_all[(df_dcs_all$control.gene=='Hypoxia' & df_dcs_all$control.cancer==cancers[i]),]$control.cor)
  z_dcs_all[i,1] <- (val_a-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_dcs_all[i,2] <- (val_b-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_dcs_all[i,3] <- (val_h-mean(dist, na.rm=T))/sd(dist, na.rm=T)
}
write.table(z_dcs_all, file='z_dcs_all.txt', sep='\t')

z_dcs_cur <- matrix(nrow=length(cancers), ncol=3)
rownames(z_dcs_cur) <- cancers
colnames(z_dcs_cur) <- c('A3A', 'A3B', 'Hypoxia')
for(i in 1:length(cancers)){
  dist <- as.numeric(df_dcs_cur[(df_dcs_cur$control.gene=='Control' & df_dcs_cur$control.cancer==cancers[i]),]$control.cor)
  val_a <- as.numeric(df_dcs_cur[(df_dcs_cur$control.gene=='A3A' & df_dcs_cur$control.cancer==cancers[i]),]$control.cor)
  val_b <- as.numeric(df_dcs_cur[(df_dcs_cur$control.gene=='A3B' & df_dcs_cur$control.cancer==cancers[i]),]$control.cor)
  val_h <- as.numeric(df_dcs_cur[(df_dcs_cur$control.gene=='Hypoxia' & df_dcs_cur$control.cancer==cancers[i]),]$control.cor)
  z_dcs_cur[i,1] <- (val_a-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_dcs_cur[i,2] <- (val_b-mean(dist, na.rm=T))/sd(dist, na.rm=T)
  z_dcs_cur[i,3] <- (val_h-mean(dist, na.rm=T))/sd(dist, na.rm=T)
}
write.table(z_dcs_cur, file='z_dcs_cur.txt', sep='\t')

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)