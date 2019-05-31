# master.R


# Jiachen Liang
# 2018
# modified by Sakura Maezono for gdrive ##------ Wed Feb 27 16:40:48 2019 ------##

# Buffa Laboratory
# Computational Biology and Integrative Genomics Group
# Department of Oncology
# University of Oxford

###############################################################################################

# DESCRIPTION
# Similar to signatures.R 
# master.R incorporates everything signatures does,
# minus I think the plotting of the huge xCell heatmaps
# and a few other not-so-useful graphs, 
# and adds differential expression, plus a bit of clean up on global variables
###############################################################################################
# Initiation
###############################################################################################

t_start <- Sys.time()
t_last <- Sys.time()

library(parallel)
library(matrixStats)
library(gplots)
library(ggplot2)

###############################################################################################
# File locations
###############################################################################################

# Set working directory
root <- "G:/My Drive/Jiachen Files/"
setwd(root)

# Directory: output
dir.create('G:/My Drive/Jiachen Files/Sakura/output/')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/master')
dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/master/'

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'

# create directories [IMPORTANT, files cannot be saved to non-existing
# directory]
dir.create(paste(dir_o, 'exp boxplots per signature/', sep = ""))
dir.create(paste(dir_o, 'exp boxplots across cancers/', sep = ""))
dir.create(paste(dir_o, 'cor_xc/', sep = ""))
dir.create(paste(dir_o, 'cor single gene/', sep = ""))
dir.create(paste(dir_o, 'cor mean master/', sep = ""))
dir.create(paste(dir_o, 'cor mean master scatterplots/', sep = ""))
dir.create(paste(dir_o, 'cor mean master scatterplots/_NT', sep = ""))
dir.create(paste(dir_o, 'cor mean master/clustered/', sep = ""))
dir.create(paste(dir_o, 'xCell heatmaps/', sep = ""))
dir.create(paste(dir_o, 'xCell cor/', sep = ""))
dir.create(paste(dir_o, 'xCell cor/unclustered/', sep = ""))
dir.create(paste(dir_o, 'xCell cor all/', sep = ""))
dir.create(paste(dir_o, 'xCell cor all/unclustered/', sep = ""))
dir.create(paste(dir_o, 'heatmaps master individual genes/', sep = ""))
dir.create(paste(dir_o, 
                 'heatmaps master individual genes/signature clustered/',
                 sep = ""))
dir.create(paste(dir_o, 
                 'heatmaps master individual genes/full clustered/',
                 sep = ""))
###############################################################################################
# Parameters
###############################################################################################

# Number of tested signatures
nts <- 4

# Boxplot dimensions
w <- 8 # width
h <- 6 # height
bpts <- 16 # text size

# Heatmap parameters
margins <- c(12,8) # margins for output image
key <- FALSE # display legend
skey <- TRUE # symmetric key
keysize <- 0.1
fontsize <- 1.2
hm_w <- 12 # width
hm_h <- 12 # height

# Significance threshold
thr <- 0.05

# Using means (T) vs medians(F)
um <- T

###############################################################################################
# Import expression matrix from file system
###############################################################################################

exp_import <- function(cancer, state) {
  
  dir <- paste(dir_e, cancer, '_', state, '.txt', sep='')
  input <- read.table(dir, sep='\t', header=T, row.names=NULL)
  
  rn <- input[,1]
  rn <- make.unique(rn, sep='.')
  
  output <- as.data.frame(input[,-1])
  rownames(output) <- rn
  
  # Fix col names
  names <- colnames(output)
  for(i in 1:length(names)){
    split <- strsplit(names[i], '\\.')
    names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(output) <- names
  
  return(output)
  
}

###############################################################################################
# Expression box plots for every gene in a given signature
###############################################################################################

ebp_sign <- function(cancer_index){
  
  # get data
  cancer <- cancers[cancer_index]
  exp_n <- exp_n_all[[cancer_index]]
  exp_t <- exp_t_all[[cancer_index]]
  
  # mean expression per sample (for next step in analysis)
  meps_n <- matrix(nrow=length(signatures), ncol=ncol(exp_n))
  meps_t <- matrix(nrow=length(signatures), ncol=ncol(exp_t))
  rn <- c()
  
  # for each signature
  for(i in 1:length(signatures)){
    
    # signature name
    name <- colnames(signatures[[i]])
    
    # subsize the expression matrices
    exp_n_sub <- exp_n[as.character(signatures[[i]][,1]),]
    exp_t_sub <- exp_t[as.character(signatures[[i]][,1]),]
    
    # build data frame
    values <- c()
    genes <- c()
    status <- c()
    
    for(j in 1:nrow(exp_n_sub)){
      for(k in 1:ncol(exp_n_sub)){
        values <- c(values, exp_n_sub[j,k])
        genes <- c(genes, rownames(exp_n_sub)[j])
        status <- c(status, 'Normal')
      }
    }
    for(j in 1:nrow(exp_t_sub)){
      for(k in 1:ncol(exp_t_sub)){
        values <- c(values, exp_t_sub[j,k])
        genes <- c(genes, rownames(exp_t_sub)[j])
        status <- c(status, 'Tumour')
      }
    }
    df <- data.frame(values, genes, status)
    
    # output box plot
    library(ggplot2)
    plot_out <- ggplot(df, aes(x=genes, y=values, fill=status)) +
      geom_boxplot() +
      ggtitle(paste('Expression of', name, 'genes in', cancer)) +
      xlab("Gene") + ylab("mRNA count (log2, normalized)") +
      theme(text = element_text(size=bpts),
            axis.text.x = element_text(angle=45, hjust=1))
    pdf(file=paste(dir_o, 'exp boxplots per signature/', cancer, '_', name, '.pdf', sep=''), width=w, height=h)
    print(plot_out)
    dev.off()
    
    # perform t-test for each gene
    pv <- matrix(nrow=nrow(exp_n_sub), ncol=1)
    for(j in 1:nrow(exp_n_sub)){
      if((length(na.omit(as.numeric(exp_n_sub[j,]))) < 3) | (length(na.omit(as.numeric(exp_t_sub[j,]))) < 3)){ # filter out genes that have less than 3 values
        pv[j,1] <- NA
      }else{
        pv[j,1] <- t.test(exp_n_sub[j,], exp_t_sub[j,])$p.value
      }
    }
    rownames(pv) <- rownames(exp_n_sub)
    colnames(pv) <- 'p-value (t test)'
    write.table(pv, file=paste(dir_o, 'exp boxplots per signature/', cancer, '_', name, '_', 'pv.txt', sep=''), sep='\t')
    
    # get mean expression per sample (for next step in analysis)
    if(um){
      meps_n[i,] <- colMeans(exp_n_sub, na.rm = T)
      meps_t[i,] <- colMeans(exp_t_sub, na.rm = T)
      rn <- c(rn, name)
    }else{
      meps_n[i,] <- colMedians(exp_n_sub, na.rm = T)
      meps_t[i,] <- colMedians(exp_t_sub, na.rm = T)
      rn <- c(rn, name)
    }
    
  }
  
  # output mean expression per sample for all signatures (for next step in analysis)
  rownames(meps_n) <- rn
  colnames(meps_n) <- colnames(exp_n)
  rownames(meps_t) <- rn
  colnames(meps_t) <- colnames(exp_t)
  meps <- list(meps_n, meps_t)
  
  return(meps)
  
}

###############################################################################################
# Expression box plots across cancers
###############################################################################################

ebp_cancers <- function(s_i){
  
  # signature name
  name <- colnames(signatures[[s_i]])
  
  # build data frame
  v <- c()
  c <- c()
  state <- c()
  
  # for each cancer
  for(i in 1:length(meps_all)){
    
    # for each sample in signature
    for(j in 1:ncol(meps_all[[i]][[1]])){
      v <- c(v, meps_all[[i]][[1]][s_i, j])
      c <- c(c, cancers[i])
      state <- c(state, 'Normal')
    }
    for(j in 1:ncol(meps_all[[i]][[2]])){
      v <- c(v, meps_all[[i]][[2]][s_i, j])
      c <- c(c, cancers[i])
      state <- c(state, 'Tumour')
    }
  }
  df <- data.frame(v, c, state)
  
  # output box plot
  if(um){
    plot_out <- ggplot(df, aes(x=c, y=v, fill=state)) +
      geom_boxplot() +
      ggtitle(paste('Expression of', name, 'genes across cancers')) +
      xlab("Cancer") + ylab("mean mRNA count (log2, normalized)") +
      theme(text = element_text(size=bpts),
            axis.text.x = element_text(angle=45, hjust=1))
  }else{
    plot_out <- ggplot(df, aes(x=c, y=v, fill=state)) +
      geom_boxplot() +
      ggtitle(paste('Expression of', name, 'genes across cancers')) +
      xlab("Cancer") + ylab("median mRNA count (log2, normalized)") +
      theme(text = element_text(size=bpts),
            axis.text.x = element_text(angle=45, hjust=1))
  }
  pdf(file=paste(dir_o, 'exp boxplots across cancers/', name, '.pdf', sep=''), width=w, height=h)
  print(plot_out)
  dev.off()
  
  # perform t-test for each gene
  pv <- matrix(nrow=length(meps_all), ncol=1)
  for(i in 1:length(meps_all)){
    if((length(na.omit(as.numeric(meps_all[[i]][[1]][s_i,]))) < 3) | (length(na.omit(as.numeric(meps_all[[i]][[2]][s_i,]))) < 3)){ # filter out < 3 values
      pv[i,1] <- NA
    }else{
      pv[i,1] <- t.test(meps_all[[i]][[1]][s_i,], meps_all[[i]][[2]][s_i,])$p.value
    }
  }
  rownames(pv) <- cancers
  colnames(pv) <- 'p-value (t test)'
  
  write.table(pv, file=paste(dir_o, 'exp boxplots across cancers/', name, 'pv.txt', sep=''), sep='\t')
  
  model <- aov(formula = v ~ c + state, data = df)
  tukey <- anova(model)
  
  write.table(tukey, file=paste(dir_o, 'exp boxplots across cancers/', name, 'anova.txt', sep=''), sep='\t')
  
}

###############################################################################################
# Correlation heatmaps (single genes)
###############################################################################################

cor_hm <- function(cancer_index){
  
  # combine all tested signatures into one list
  combined_sign <- c()
  # for(i in 2:length(signatures)){
  #   for(j in 1:nrow(signatures[[i]])){
  #     combined_sign <- c(combined_sign, as.character(signatures[[i]][j,1]))
  #   }
  # }
  combined_sign <- as.character(signatures[[2]][,1]) # OVERRIDE: DO NOT INCLUDE OTHER SIGNATURES
  master_sign <- as.character(signatures[[1]][,1]) # master signature
  
  # subset data
  cancer <- cancers[cancer_index]
  exp_t <- exp_t_all[[cancer_index]]
  exp_t_sub <- exp_t[combined_sign,] # tested signatures
  exp_t_master <- exp_t[master_sign,] # master signature
  exp_n <- exp_n_all[[cancer_index]]
  exp_n_sub <- exp_n[combined_sign,] # tested signatures
  exp_n_master <- exp_n[master_sign,] # master signature
  
  # correlation analysis (spearman)
  table_cor <- matrix(nrow=nrow(exp_t_master), ncol=nrow(exp_t_sub))
  table_pv <- matrix(nrow=nrow(exp_t_master), ncol=nrow(exp_t_sub))
  for(i in 1:nrow(exp_t_master)){
    for(j in 1:nrow(exp_t_sub)){
      test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(exp_t_sub[j,]), method='spearman')
      table_cor[i,j] <- test$estimate
      table_pv[i,j] <- test$p.value
    }
  }
  
  # remove insignificant values
  table_cor[table_pv >= 0.05] <- 0
  
  # naming
  rownames(table_cor) <- master_sign
  colnames(table_cor) <- combined_sign
  rownames(table_pv) <- master_sign
  colnames(table_pv) <- combined_sign
  
  # plot heatmap
  pdf(file=paste(dir_o, 'cor single gene/', cancer, '.pdf', sep=''), width=12, height=4)
  hm <- heatmap.2(
    data.matrix(table_cor),
    #main=paste("1. Mean -", m),
    Rowv=NA,
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=margins,
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  # write tables
  write.table(table_cor, file=paste(dir_o, 'cor single gene/', cancer, '_cor.txt', sep=''), sep='\t')
  write.table(table_pv, file=paste(dir_o, 'cor single gene/', cancer, '_pv.txt', sep=''), sep='\t')
  
  # get mean of signatures (for downstream analysis)
  table_means <- matrix(nrow=length(signatures)-1, ncol=ncol(exp_t_sub))
  for(i in 2:length(signatures)){
    if(um){
      table_means[i-1,] <- colMeans(exp_t[as.character(signatures[[i]][,1]),])
    }else{
      table_means[i-1,] <- colMedians(exp_t[as.character(signatures[[i]][,1]),])
    }
  }
  table_means_n <- matrix(nrow=length(signatures)-1, ncol=ncol(exp_n_sub))
  for(i in 2:length(signatures)){
    if(um){
      table_means_n[i-1,] <- colMeans(exp_n[as.character(signatures[[i]][,1]),])
    }else{
      table_means_n[i-1,] <- colMedians(exp_n[as.character(signatures[[i]][,1]),])
    }
  }
  
  # correlate means with master signature
  table_cor_means <- matrix(nrow=nrow(exp_t_master), ncol=nrow(table_means))
  table_pv_means <- matrix(nrow=nrow(exp_t_master), ncol=nrow(table_means))
  for(i in 1:nrow(exp_t_master)){
    for(j in 1:nrow(table_means)){
      test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(table_means[j,]), method='spearman')
      table_cor_means[i,j] <- test$estimate
      table_pv_means[i,j] <- test$p.value
    }
  }
  
  # Generate scatter plots
  for(i in 1:nrow(exp_t_master)){
    out_dir <- paste(dir_o, 'cor mean master scatterplots/', cancer, '-', as.character(signatures[[1]][,1])[i], '_cor.pdf', sep='')
    a <- as.numeric(exp_t_master[i,])
    h <- as.numeric(table_means[1,])
    df <- data.frame(a, h)
    plot <- ggplot(df, aes(x=a, y=h)) +
      geom_point(alpha=0.6) +
      ggtitle(paste(as.character(signatures[[1]][,1])[i], 'vs Hypoxia -', cancer, '(Tumour)')) +
      xlab(as.character(signatures[[1]][,1])[i]) +
      ylab('mean Hypoxia signature') +
      theme(text = element_text(size=12))
    pdf(file=out_dir, width=6, height=4)
    print(plot)
    dev.off()
  }
  for(i in 1:nrow(exp_t_master)){
    out_dir <- paste(dir_o, 'cor mean master scatterplots/_NT/', cancer, '-', as.character(signatures[[1]][,1])[i], '_cor_NT.pdf', sep='')
    a <- c(as.numeric(exp_t_master[i,]), as.numeric(exp_n_master[i,]))
    h <- c(as.numeric(table_means[1,]), as.numeric(table_means_n[1,]))
    c <- c(rep('tumour', length(as.numeric(exp_t_master[i,]))), rep('normal', length(as.numeric(exp_n_master[i,]))))
    df <- data.frame(a, h, c)
    plot <- ggplot(df, aes(x=a, y=h, color=c)) +
      geom_point(alpha=0.6) +
      ggtitle(paste(as.character(signatures[[1]][,1])[i], 'vs Hypoxia -', cancer, '(Tumour)')) +
      xlab(as.character(signatures[[1]][,1])[i]) +
      ylab('mean Hypoxia signature') +
      theme(text = element_text(size=12))
    pdf(file=out_dir, width=6, height=4)
    print(plot)
    dev.off()
  }
  
  # remove insignificant values
  table_cor_means[table_pv_means >= 0.05] <- 0
  
  # get signature names
  sign_names <- c()
  for(i in 2:length(signatures)){
    sign_names <- c(sign_names, colnames(signatures[[i]]))
  }
  
  # naming
  rownames(table_cor_means) <- master_sign
  colnames(table_cor_means) <- sign_names
  rownames(table_pv_means) <- master_sign
  colnames(table_pv_means) <- sign_names
  
  # output
  out <- list(table_cor_means, table_pv_means)
  return(out)
  
}

###############################################################################################
# Correlation heatmaps (signature means)
###############################################################################################

cor_mean_hm <- function(mean_cor){
  
  for(i in 2:length(signatures)){
    
    sign_name <- colnames(signatures[[i]])
    
    # assemble matrix
    mat_cor <- matrix(nrow=nrow(signatures[[1]]), ncol=length(cancers))
    mat_pv <- matrix(nrow=nrow(signatures[[1]]), ncol=length(cancers))
    for(j in 1:length(cancers)){
      mc <- mean_cor[[j]]
      mc_cor <- mc[[1]]
      mc_pv <- mc[[2]]
      mat_cor[,j] <- mc_cor[,i-1]
      mat_pv[,j] <- mc_pv[,i-1]
    }
    
    # names
    rownames(mat_cor) <- as.character(signatures[[1]][,1])
    colnames(mat_cor) <- cancers
    rownames(mat_pv) <- as.character(signatures[[1]][,1])
    colnames(mat_pv) <- cancers
    
    # plot heatmap
    pdf(file=paste(dir_o, 'cor mean master/', sign_name, '.pdf', sep=''), width=9, height=9)
    hm <- heatmap.2(
      data.matrix(mat_cor),
      #main=paste("1. Mean -", m),
      Rowv=NA, Colv=NA,
      dendrogram="none",
      col=greenred,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # plot heatmap (clustered)
    pdf(file=paste(dir_o, 'cor mean master/clustered/', sign_name, '.pdf', sep=''), width=9, height=9)
    hm <- heatmap.2(
      data.matrix(mat_cor),
      #main=paste("1. Mean -", m),
      col=greenred,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # write files
    write.table(mat_cor, file=paste(dir_o, 'cor mean master/', sign_name, '_cor.txt', sep=''), sep='\t')
    write.table(mat_pv, file=paste(dir_o, 'cor mean master/', sign_name, '_pv.txt', sep=''), sep='\t')
    
  }
  
}

###############################################################################################
# xCell
###############################################################################################

import_xcell <- function(){
  xc_raw <- read.table(file=paste(dir_p, 'xCell_TCGA_RSEM.txt', sep=''), sep='\t', header=T, row.names=1)
  names <- colnames(xc_raw)
  for(i in 1:length(names)){
    split <- strsplit(names[i], '\\.')
    names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
  }
  colnames(xc_raw) <- names
  xc <- list('vector', length=length(cancers))
  for(i in 1:length(cancers)){
    cn <- colnames(exp_t_all[[i]])
    common <- intersect(cn, colnames(xc_raw))
    entry <- xc_raw[,common]
    xc[[i]] <- entry
  }
  return(xc)
}

xcell <- function(cancer_index){
  
  cancer <- cancers[cancer_index]
  
  # plot heatmap
  pdf(file=paste(dir_o, 'xCell heatmaps/', cancer, '.pdf', sep=''), width=12, height=12)
  hm <- heatmap.2(
    data.matrix(xc[[cancer_index]]),
    #main=paste("1. Mean -", m),
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=margins,
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  # row means for further analysis
  if(um){
    xc_rm <- rowMeans(xc[[cancer_index]])
  }else{
    xc_rm <- rowMedians(xc[[cancer_index]])
  }
  
  # Prepare for correlation analysis
  cancer <- cancers[cancer_index]
  exp_t <- exp_t_all[[cancer_index]]
  common <- intersect(colnames(xc[[cancer_index]]), colnames(exp_t))
  exp_t <- exp_t[,common]
  
  # APOBEC vs xCell
  exp_t_sub <- exp_t[as.character(signatures[[1]][,1]),]
  exp_t_master <- xc[[cancer_index]]
  
  # correlation analysis (spearman)
  table_cor <- matrix(nrow=nrow(exp_t_master), ncol=nrow(exp_t_sub))
  table_pv <- matrix(nrow=nrow(exp_t_master), ncol=nrow(exp_t_sub))
  for(i in 1:nrow(exp_t_master)){
    for(j in 1:nrow(exp_t_sub)){
      test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(exp_t_sub[j,]), method='spearman')
      table_cor[i,j] <- test$estimate
      table_pv[i,j] <- test$p.value
    }
  }
  
  # remove insignificant values
  table_cor[table_pv >= 0.05] <- 0
  
  # naming
  rownames(table_cor) <- rownames(xc[[cancer_index]])
  colnames(table_cor) <- rownames(exp_t_sub)
  rownames(table_pv) <- rownames(xc[[cancer_index]])
  colnames(table_pv) <- rownames(exp_t_sub)
  
  # plot heatmap
  pdf(file=paste(dir_o, 'xCell cor/', cancer, '_', sig_names[1], '.pdf', sep=''), width=6, height=16)
  hm <- heatmap.2(
    data.matrix(table_cor),
    #main=paste("1. Mean -", m),
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  pdf(file=paste(dir_o, 'xCell cor/unclustered/', cancer, '_', sig_names[1], '.pdf', sep=''), width=6, height=16)
  hm <- heatmap.2(
    data.matrix(table_cor),
    #main=paste("1. Mean -", m),
    Rowv=F, Colv=F,
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  # write tables
  write.table(table_cor, file=paste(dir_o, 'xCell cor/', cancer, '_', sig_names[1], '_cor.txt', sep=''), sep='\t')
  write.table(table_pv, file=paste(dir_o, 'xCell cor/', cancer, '_', sig_names[1], '_pv.txt', sep=''), sep='\t')
  
  # get mean of signatures (for downstream analysis)
  table_means <- matrix(nrow=2, ncol=ncol(exp_t_sub))
  if(um){
    table_means[1,] <- colMeans(exp_t_sub, na.rm=T)
  }else{
    table_means[1,] <- colMedians(exp_t_sub, na.rm=T)
  }
  
  # Hypoxia vs xCell
  exp_t_sub <- exp_t[as.character(signatures[[2]][,1]),]
  exp_t_master <- xc[[cancer_index]]
  
  # correlation analysis (spearman)
  table_cor <- matrix(nrow=nrow(exp_t_master), ncol=nrow(exp_t_sub))
  table_pv <- matrix(nrow=nrow(exp_t_master), ncol=nrow(exp_t_sub))
  for(i in 1:nrow(exp_t_master)){
    for(j in 1:nrow(exp_t_sub)){
      test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(exp_t_sub[j,]), method='spearman')
      table_cor[i,j] <- test$estimate
      table_pv[i,j] <- test$p.value
    }
  }
  
  # remove insignificant values
  table_cor[table_pv >= 0.05] <- 0
  
  # naming
  rownames(table_cor) <- rownames(xc[[cancer_index]])
  colnames(table_cor) <- rownames(exp_t_sub)
  rownames(table_pv) <- rownames(xc[[cancer_index]])
  colnames(table_pv) <- rownames(exp_t_sub)
  
  # plot heatmap
  pdf(file=paste(dir_o, 'xCell cor/', cancer, '_', sig_names[2], '.pdf', sep=''), width=12, height=16)
  hm <- heatmap.2(
    data.matrix(table_cor),
    #main=paste("1. Mean -", m),
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  pdf(file=paste(dir_o, 'xCell cor/unclustered/', cancer, '_', sig_names[2], '.pdf', sep=''), width=12, height=16)
  hm <- heatmap.2(
    data.matrix(table_cor),
    #main=paste("1. Mean -", m),
    Rowv=F, Colv=F,
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  # write tables
  write.table(table_cor, file=paste(dir_o, 'xCell cor/', cancer, '_', sig_names[2], '_cor.txt', sep=''), sep='\t')
  write.table(table_pv, file=paste(dir_o, 'xCell cor/', cancer, '_', sig_names[2], '_pv.txt', sep=''), sep='\t')
  
  # get mean of signatures (for downstream analysis)
  if(um){
    table_means[2,] <- colMeans(exp_t_sub, na.rm=T)
  }else{
    table_means[2,] <- colMedians(exp_t_sub, na.rm=T)
  }
  
  # correlate means with master signature
  table_cor_means <- matrix(nrow=nrow(exp_t_master), ncol=nrow(table_means))
  table_pv_means <- matrix(nrow=nrow(exp_t_master), ncol=nrow(table_means))
  for(i in 1:nrow(exp_t_master)){
    for(j in 1:nrow(table_means)){
      test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(table_means[j,]), method='spearman')
      table_cor_means[i,j] <- test$estimate
      table_pv_means[i,j] <- test$p.value
    }
  }
  
  # remove insignificant values from means
  table_cor_means[table_pv_means >= 0.05] <- NA
  
  # naming means
  rownames(table_cor_means) <- rownames(xc[[cancer_index]])
  rownames(table_pv_means) <- rownames(xc[[cancer_index]])
  colnames(table_cor_means) <- c(sig_names[1], sig_names[2])
  colnames(table_pv_means) <- c(sig_names[1], sig_names[2])
  
  out <- list(table_cor_means, table_pv_means, xc_rm)
  return(out)
  
}

xcell_all <- function(xc_out){
  
  # get means data from previous analysis
  score_means <- matrix(nrow=length(xc_out[[1]][[3]]), ncol=length(xc_out))
  for(i in 1:length(xc_out)){
    score_means[,i] <- xc_out[[i]][[3]]
  }
  rownames(score_means) <- rownames(xc[[1]])
  colnames(score_means) <- cancers
  
  # remove insignificant values
  score_means[is.na(score_means)] <- 0
  
  # file naming
  file.name <- c()
  if(um){
    file.name <- 'xCell cor all/mean_xcell_scores.pdf'
  }else{
    file.name <- 'xCell cor all/median_xcell_scores.pdf'
  }
  
  # graph
  pdf(file=paste(dir_o, file.name, sep=''), width=6, height=18)
  hm <- heatmap.2(
    score_means,
    #main=paste("1. Mean -", m),
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  score_means <- matrix(nrow=length(as.numeric(xc_out[[1]][[1]][,1])), ncol=length(xc_out))
  for(i in 1:length(xc_out)){
    score_means[,i] <- as.numeric(xc_out[[i]][[1]][,1])
  }
  rownames(score_means) <- rownames(xc[[i]])
  colnames(score_means) <- cancers
  score_means[is.na(score_means)] <- 0
  pdf(file=paste(dir_o, 'xCell cor all/cor_', sig_names[1], '.pdf', sep=''), width=6, height=18)
  hm <- heatmap.2(
    score_means,
    #main=paste("1. Mean -", m),
    col=redgreen,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  pdf(file=paste(dir_o, 'xCell cor all/unclustered/cor_', sig_names[1], '.pdf', sep=''), width=6, height=18)
  hm <- heatmap.2(
    score_means,
    #main=paste("1. Mean -", m),
    Rowv=F, Colv=F,
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  score_means <- matrix(nrow=length(as.numeric(xc_out[[1]][[1]][,2])), ncol=length(xc_out))
  for(i in 1:length(xc_out)){
    score_means[,i] <- as.numeric(xc_out[[i]][[1]][,2])
  }
  rownames(score_means) <- rownames(xc[[i]])
  colnames(score_means) <- cancers
  score_means[is.na(score_means)] <- 0
  pdf(file=paste(dir_o, 'xCell cor all/cor_', sig_names[2], '.pdf', sep=''), width=6, height=18)
  hm <- heatmap.2(
    score_means,
    #main=paste("1. Mean -", m),
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  pdf(file=paste(dir_o, 'xCell cor all/unclustered/cor_', sig_names[2], '.pdf', sep=''), width=6, height=18)
  hm <- heatmap.2(
    score_means,
    #main=paste("1. Mean -", m),
    Rowv=F, Colv=F,
    col=greenred,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=c(10,20),
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
}

###############################################################################################
# Heatmaps of correlation of every gene in master signature with others + xCell
###############################################################################################

hm_master <- function(master_index){
  
  this.gene <- as.character(signatures[[1]][master_index,1])
  
  # for each signature
  for(i in 2:length(signatures)){
    
    # matrix to store correlation values
    mat.cor <- matrix(nrow=nrow(signatures[[i]])+1, ncol=length(cancers)+1)
    mat.pv <- matrix(nrow=nrow(signatures[[i]])+1, ncol=length(cancers)+1)
    if(um){
      rownames(mat.cor) <- c(as.character(signatures[[i]][,1]), 'Signature mean')
      rownames(mat.pv) <- c(as.character(signatures[[i]][,1]), 'Signature mean')
      colnames(mat.cor) <- c(cancers, 'Pan-cancer mean')
      colnames(mat.pv) <- c(cancers, 'Pan-cancer mean')
    }else{
      rownames(mat.cor) <- c(as.character(signatures[[i]][,1]), 'Signature median')
      rownames(mat.pv) <- c(as.character(signatures[[i]][,1]), 'Signature median')
      colnames(mat.cor) <- c(cancers, 'Pan-cancer median')
      colnames(mat.pv) <- c(cancers, 'Pan-cancer median')
    }
    
    # for each cancer
    for(cancer_index in 1:length(cancers)+1){
      
      # get the expression for this master gene in one cancer
      # exception: the last index refers to across cancers
      exp_this.gene <- NULL
      if(cancer_index<=length(cancers)){
        exp_this.gene <- as.numeric(exp_t_all[[cancer_index]][this.gene,])
      }else{
        exp_this.gene <- as.numeric(exp_t_full[this.gene,])
      }
      
      # for each gene in the signature
      for(j in 1:nrow(signatures[[i]])){
        
        # get the expression for the gene
        sig.gene <- as.character(signatures[[i]][j,1])
        if(cancer_index<=length(cancers)){
          exp_sig.gene <- as.numeric(exp_t_all[[cancer_index]][sig.gene,])
        }else{
          exp_sig.gene <- as.numeric(exp_t_full[sig.gene,])
        }
        
        # correlate the expression of master gene with signature gene
        test <- cor.test(exp_this.gene, exp_sig.gene, method='spearman')
        mat.cor[j,cancer_index] <- test$estimate
        mat.pv[j,cancer_index] <- test$p.value
        
      }
      
      # mean of the signature
      sigmean <- c()
      if(um){
        if(cancer_index<=length(cancers)){ # individual cancers
          sigmean <- colMeans(exp_t_all[[cancer_index]][as.character(signatures[[i]][,1]),])
        }else{ # across cancers
          sigmean <- colMeans(exp_t_full[as.character(signatures[[i]][,1]),])
        }
      }else{
        if(cancer_index<=length(cancers)){
          sigmean <- colMedians(exp_t_all[[cancer_index]][as.character(signatures[[i]][,1]),])
        }else{
          sigmean <- colMedians(exp_t_full[as.character(signatures[[i]][,1]),])
        }
      }
      sigmean <- as.numeric(sigmean)
      
      # correlate the expressions using the mean of the signature
      test <- cor.test(exp_this.gene, sigmean, method='spearman')
      mat.cor[j+1,cancer_index] <- test$estimate
      mat.pv[j+1,cancer_index] <- test$p.value
      
    }
    
    # filter out p>0.05
    mat.cor.filtered <- mat.cor
    mat.cor.filtered[mat.pv>thr] <- 0
    mat.cor.filtered[is.na(mat.cor.filtered)] <- 0
    
    # plot heatmaps (unclustered)
    pdf(file=paste(dir_o, 'heatmaps master individual genes/',
                   this.gene, '_', sig_names[i], '.pdf', sep=''), width=6, height=18)
    hm <- heatmap.2(
      mat.cor.filtered,
      Rowv=F, Colv=F,
      col=greenred,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(10,10),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # plot heatmaps (signature clustered)
    pdf(file=paste(dir_o, 'heatmaps master individual genes/signature clustered/',
                   this.gene, '_', sig_names[i], '.pdf', sep=''), width=6, height=18)
    hm <- heatmap.2(
      mat.cor.filtered,
      Colv=F,
      col=greenred,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(10,10),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # plot heatmaps (fully clustered)
    pdf(file=paste(dir_o, 'heatmaps master individual genes/full clustered/', this.gene, '_', sig_names[i], '.pdf', sep=''), width=6, height=18)
    hm <- heatmap.2(
      mat.cor.filtered,
      col=greenred,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(10,10),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # output tables
    write.table(mat.cor, file=paste(dir_o, 'heatmaps master individual genes/', this.gene, '_', sig_names[i], '_cor.txt', sep=''))
    write.table(mat.pv, file=paste(dir_o, 'heatmaps master individual genes/', this.gene, '_', sig_names[i], '_pv.txt', sep=''))
    
  }
  
}

###############################################################################################
# Main function
###############################################################################################

# Import parameters
print('Importing parameters')
cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep=''), header=T, sep='\t')[,1])
signatures <- vector('list', nts+1)
signatures[[1]] <- read.table(paste(dir_p, 'sign_master.txt', sep=''), header=T, sep='\t')
sig_names <- as.character(colnames(signatures[[1]]))
for(i in 1:nts){
  signatures[[i+1]] <- read.table(paste(dir_p, 'sign_', i, '.txt', sep=''), header=T, sep='\t')
  sig_names <- c(sig_names, as.character(colnames(signatures[[i+1]])))
}

# Import expression matrices (parallel)
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Importing expression matrices')
f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(exp_import(cancers[i], 'n'))
}
exp_n_all <- mccollect(f) # all normal exp matrices
f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(exp_import(cancers[i], 't'))
}
exp_t_all <- mccollect(f) # tumour exp matrices

# Import xCell data
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Importing xCell data')
xc <- import_xcell()

# Rearrange all samples across cancers in one list
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Pre-processing data')
exp_t_full <- exp_t_all[[1]]
for(i in 2:length(exp_t_all)){
  exp_t_full <- cbind(exp_t_full, exp_t_all[[i]])
}

# Median expression of all cancers
exp_med <- matrix(nrow=length(cancers), ncol=2)
for(i in 1:length(cancers)){
  median_n <- exp_n_all[[i]]
  median_n[is.na(median_n)] <- 0
  median_t <- exp_t_all[[i]]
  median_t[is.na(median_t)] <- 0
  exp_med[i,1] <- median(as.matrix(median_t), na.rm = T)
  exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)
}
rownames(exp_med) <- cancers
colnames(exp_med) <- c('Normal', 'Tumour')
write.table(exp_med, file='data/exp_med.txt', sep='\t')

# Generate expression box plots for every gene in all signatures, in all cancers
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Generating expression box plots for genes in every signature')
f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(ebp_sign(i))
}
meps_all <- mccollect(f)

# Generate expression boxplots across cancers
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Generating expression box plots across cancers')
f <- vector("list", length=length(signatures))
for(i in 1:length(signatures)){
  f[[i]] <- mcparallel(ebp_cancers(i))
}
mccollect(f)

# Generate single gene correlation heatmaps
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Generating single gene correlation heatmaps')
f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(cor_hm(i))
}
mean_cor <- mccollect(f)

# Generate signature mean correlation heatmaps
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Generating signature mean correlation heatmaps')
cor_mean_hm(mean_cor)

# xCell
print(Sys.time()-t_last)
t_last <- Sys.time()
print('xCell')
f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(xcell(i))
}
xc_out <- mccollect(f)
xcell_all(xc_out)

# Generate heatmaps of correlation of every gene in master signature with others
print(Sys.time()-t_last)
t_last <- Sys.time()
print('Generating heatmaps of correlation of every gene in master signature with others')
f <- vector("list", length=nrow(signatures[[1]]))
for(i in 1:nrow(signatures[[1]])){
  f[[i]] <- mcparallel(hm_master(i))
}
mccollect(f)

###############################################################################################
# Correlation partial to xCell
###############################################################################################

print(Sys.time()-t_last)
t_last <- Sys.time()
print('Calculating correlation partial to xCell')

# Directory: output
dir_o <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/master/cor_xc/'

library(matrixStats)
library(ggplot2)
library(DAAG)
library(gplots)
library(ppcor)

# Significance threshold
thr_sig <- 0.05

# Threshold for xCell inclusion (xCell enrichment)
thr_xc <- 0.2

# Threshold for xCell inclusion (% of samples enriched)
thr_xc_s <- 0.2

# Minimum proportion of non-NA samples
thr_na <- 0.8

# APOBEC and hypoxia signatures
apo.genes <- as.character(read.table(file="data/tr_a.txt", sep="\t", header=TRUE, row.names=1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file="data/tr_h.txt", sep="\t", header=TRUE, row.names=1)$Gene.Symbol)

# Signatures to exclude
exc <- c(
  'Adipocytes',
  'Epithelial cells',
  'Smooth muscle',
  'Sebocytes',
  'Hepatocytes',
  'Keratinocytes',
  'Endothelial cells'
)

# Import xCell data
exp_all <- exp_t_all
xc_raw <- read.table(file=paste(dir_p, 'xCell_TCGA_RSEM.txt', sep=''), sep='\t', header=T, row.names=1)
names <- colnames(xc_raw)
for(i in 1:length(names)){
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
}
colnames(xc_raw) <- names
xc <- list('vector', length=length(cancers))
for(i in 1:length(cancers)){
  cn <- colnames(exp_all[[i]])
  common <- intersect(cn, colnames(xc_raw))
  entry <- xc_raw[,common]
  xc[[i]] <- entry
  exp_all[[i]] <- exp_all[[i]][,common]
}

# Main function

# Matrices to store correlation values
cor_p <- matrix(nrow=length(apo.genes), ncol=length(cancers)) # rho (partial)
rownames(cor_p) <- apo.genes
colnames(cor_p) <- cancers
pv_p <- cor_p # pv (partial)
cor_r <- cor_p # rho (regular)
pv_r <- cor_p # pv (regular)
diff <- cor_p # differential between partial and normal

for(c in 1:length(cancers)){
  
  cancer <- cancers[c]
  
  #################################################################################################
  # Determine which xCell signatures to keep (using thresholds)
  
  xc_mat_full <- xc[[c]]
  to_keep <- c()
  for(i in 1:nrow(xc_mat_full)){
    if(
      count(as.numeric(xc_mat_full[i,])>=thr_xc, na.rm=T) # Number of samples with significant signature enrichment
      >=
      (thr_xc_s*ncol(xc_mat_full)) # thr_xc_s% of samples in the dataset
    ){
      to_keep <- c(to_keep, i) # keep this signature
    }
  }
  xc_mat <- xc_mat_full[to_keep,]
  xc_mat <- xc_mat[setdiff(rownames(xc_mat), exc),]
  
  #################################################################################################
  # Generate mean hypoxia signature
  
  exp_hyp <- exp_all[[c]][hyp.genes,]
  hyp_mean <- colMeans(exp_hyp, na.rm=T)
  
  #################################################################################################
  # Subset APOBECs
  
  exp_apo <- exp_all[[c]][apo.genes,]
  
  #################################################################################################
  # Correlation
  
  for(i in 1:length(apo.genes)){
    gene <- apo.genes[i]
    # only use rows with no missing values
    temp_mat <- cbind(as.numeric(exp_all[[c]][gene,]), as.numeric(hyp_mean), t(xc_mat))
    not_missing <- !is.na(rowMeans(temp_mat))
    # only test if there are > 3 non-null rows
    if(count(not_missing)>=3){
      x <- as.numeric(exp_all[[c]][gene,])[not_missing]
      y <- as.numeric(hyp_mean)[not_missing]
      z <- t(xc_mat)[not_missing,]
      
      # Partial
      test <- pcor.test(x, y, z, method="spearman")
      cor_p[i,c] <- test$estimate
      pv_p[i,c] <- test$p.value
      
      # Regular
      test <- cor.test(x, y, method="spearman")
      cor_r[i,c] <- test$estimate
      pv_r[i,c] <- test$p.value
      
    }else{
      cor_p[i,c] <- NA
      pv_p[i,c] <- NA
    }
  }
  pv_p[,c] <- as.numeric(p.adjust(pv_p[,c], method = 'fdr'))
  pv_r[,c] <- as.numeric(p.adjust(pv_r[,c], method = 'fdr'))
  
}

# Unfiltered rho values
cor_p_unfiltered <- cor_p
cor_r_unfiltered <- cor_r

# Filter rho values for p values
cor_p[pv_p>=thr_sig] <- NA
cor_r[pv_r>=thr_sig] <- NA

# For plotting, replace NAs with zeroes
cor_p_plot <- cor_p
cor_r_plot <- cor_r
cor_p_plot[is.na(cor_p_plot)] <- 0
cor_r_plot[is.na(cor_r_plot)] <- 0

# Differential (log10)
diff <- cor_p/cor_r
diff[is.nan(diff)] <- NA
diff[is.infinite(diff)] <- NA
diff <- log10(diff)
diff[is.na(diff)] <- 0

write.table(cor_p, file=paste(dir_o, 'cor_p.txt', sep=''), sep='\t')
write.table(cor_r, file=paste(dir_o, 'cor_r.txt', sep=''), sep='\t')
write.table(cor_p_unfiltered, file=paste(dir_o, 'cor_p_unfiltered.txt', sep=''), sep='\t')
write.table(cor_r_unfiltered, file=paste(dir_o, 'cor_r_unfiltered.txt', sep=''), sep='\t')
write.table(pv_p, file=paste(dir_o, 'pv_p.txt', sep=''), sep='\t')
write.table(pv_r, file=paste(dir_o, 'pv_r.txt', sep=''), sep='\t')
write.table(diff, file=paste(dir_o, 'diff.txt', sep=''), sep='\t')

pdf(file=paste(dir_o, 'cor_p.pdf', sep=''), width=8, height=8)
hm <- heatmap.2(
  cor_p_plot,
  Rowv=NA, Colv=NA,
  xlab = 'Cancers',
  main = 'Spearman Correlation (partial to xCell)',
  dendrogram="none", density.info="none", trace="none",
  col=greenred,
  breaks=seq(-1, 1, length.out=101),
  margins=c(7,7)
)
dev.off()
pdf(file=paste(dir_o, 'cor_r.pdf', sep=''), width=8, height=8)
hm <- heatmap.2(
  cor_r_plot,
  Rowv=NA, Colv=NA,
  xlab = 'Cancers',
  main = 'Spearman Correlation',
  dendrogram="none", density.info="none", trace="none",
  col=greenred,
  breaks=seq(-1, 1, length.out=101),
  margins=c(7,7)
)
dev.off()
pdf(file=paste(dir_o, 'diff.pdf', sep=''), width=8, height=8)
hm <- heatmap.2(
  diff,
  Rowv=NA, Colv=NA,
  xlab = 'Cancers',
  main = 'Differential (log10)',
  dendrogram="none", density.info="none", trace="none",
  col=greenred,
  breaks=seq(-1, 1, length.out=101),
  margins=c(7,7)
)
dev.off()

###############################################################################################
# Self-correlation
###############################################################################################

print(Sys.time()-t_last)
t_last <- Sys.time()
print('Calculating self-correlation')

# Directory: expression matrices
dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/apobec/'

# Directory: output
dir_o <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/master/autocorr/'

# Import parameters
cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep=''), header=T, sep='\t')[,1])

for(c in cancers){
  exp <- read.table(paste(dir_e, c, '.txt', sep=''), sep='\t', header=T, row.names=1)
  exp <- exp[order(rownames(exp)),]
  mat <- matrix(nrow=nrow(exp), ncol=nrow(exp))
  rownames(mat) <- rownames(exp)
  colnames(mat) <- rownames(exp)
  pv <- mat
  for(i in 1:nrow(exp)){
    for(j in 1:nrow(exp)){
      # t <- cor.test(as.numeric(exp[i,]), as.numeric(exp[j,]), method = 'spearman')
      # mat[i,j] <- t$estimate
      # pv[i,j] <- t$p.value
      result = tryCatch({
        t <- cor.test(as.numeric(exp[i,]), as.numeric(exp[j,]), method = 'spearman')
        mat[i,j] <- t$estimate
        pv[i,j] <- t$p.value
      }, error = function(e) {
        mat[i,j] <- NA
        pv[i,j] <- NA
      }, finally = {
        
      })
    }
  }
  write.table(mat, file=paste(dir_o, c, '_cor.txt', sep=''), sep='\t')
  write.table(pv, file=paste(dir_o, c, '_pv.txt', sep=''), sep='\t')
  mat[pv>=0.05] <- NA
  mat[is.na(mat)] <- 0
  pdf(file=paste(dir_o, c, '_autocorr.pdf', sep=''), width=5, height=5)
  hm <- heatmap.2(
    mat,
    main=c,
    Rowv=NA, Colv=NA,
    dendrogram="none",
    col=greenred,
    margins=c(8,8),
    density.info="none", trace="none",
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
}

###############################################################################################
# Termination
###############################################################################################

print(Sys.time()-t_last)

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
