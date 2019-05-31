# signatures.R

# Jiachen Liang
# 2018

# Buffa Laboratory
# Computational Biology and Integrative Genomics Group
# Department of Oncology
# University of Oxford

###############################################################################################

# Master function to compare a signature with one or many others

# Inputs:
#   List of cancers
#   Master signature (sign. 1)
#   Tested signatures (sign. 2 to n-1)
#   Housekeeping genes (sign. n)
#   Expression matrix (normal)
#   Expression matrix (tumour)

# Outputs:
#   In each cancer
#     Boxplot of master signature expression (n/t)
#       + t-test results
#     Boxplots of tested signatures expression (n/t)
#       + t-test results
#     Heatmap of correlation between all genes in master signature vs all genes in all tested signatures
#       + associated rho and p-value tables
#   Boxplot of mean master signature expression across cancers (n/t)
#     + t-test results
#   Boxplots of mean tested signatures expression across cancers (n/t)
#     t-test results
#   Heatmap of expression of all genes across all signautres, in all samples
#   Heatmap of expression of all genes across all signatures, in all samples, normalized to housekeeping genes
#   Heatmap of correlation between the mean expression of master signautre in all samples vs each gene in all tested signatures
#     + associated rho and p-value tables

###############################################################################################
# Initiation
###############################################################################################

t_start <- Sys.time()
t_last <- Sys.time()

# Using Windows?
win <- F

###############################################################################################
# Associated library functions
###############################################################################################

if(win){
  library(doParallel)
}else{
  library(parallel)
}

library(matrixStats)
library(gplots)
library(ggplot2)

###############################################################################################
# File locations
###############################################################################################

if(win){
  # Set working directory
  root <- "C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
  setwd(root)
  
  # Directory: output
  dir_o <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/signatures_output/'
  
  # Directory: parameters
  dir_p <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'
}else{
  # Set working directory
  root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
  setwd(root)
  
  # Directory: output
  dir_o <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/signatures_output/'
  
  # Directory: parameters
  dir_p <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'
}

###############################################################################################
# Parameters
###############################################################################################

# Number of tested signatures
nts <- 4

# Boxplot dimensions
w <- 8 # width
h <- 6 # height

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

exp_import <- function(cancer, state){
  
  dir <- paste(dir_e, cancer, '_', state, '.txt', sep='')
  input <- read.table(dir, sep='\t', header=T, row.names=NULL)
  
  rn <- input[,1]
  rn <- make.unique(rn, sep='.')
  
  output <- input[,-1]
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
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
      theme(text = element_text(size=16),
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
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
      theme(text = element_text(size=16),
            axis.text.x = element_text(angle=45, hjust=1))
  }else{
    plot_out <- ggplot(df, aes(x=c, y=v, fill=state)) +
      geom_boxplot() +
      ggtitle(paste('Expression of', name, 'genes across cancers')) +
      xlab("Cancer") + ylab("median mRNA count (log2, normalized)") +
      theme(text = element_text(size=16),
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
# Generate expression heatmap
###############################################################################################

exp_hm <- function(cancer_index){
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
  # combine all signatures into one list
  combined_sign <- c()
  for(i in 1:length(signatures)){
    for(j in 1:nrow(signatures[[i]])){
      combined_sign <- c(combined_sign, as.character(signatures[[i]][j,1]))
    }
  }
  
  # subset data
  cancer <- cancers[cancer_index]
  exp_t <- exp_t_all[[cancer_index]]
  exp_t_sub <- exp_t[combined_sign,]
  
  # remove NAs
  exp_t_sub[is.na(exp_t_sub)] <- 0
  
  # heatmap limits
  med <- exp_med[cancer_index, 2]
  ul <- 2*med
  ll <- 0
  
  # plot heatmap
  pdf(file=paste(dir_o, 'exp heatmaps/', cancer, '.pdf', sep=''), width=hm_w, height=hm_h)
  hm <- heatmap.2(
    data.matrix(exp_t_sub),
    #main=paste("1. Mean -", m),
    Rowv=NA,
    dendrogram="none",
    col=redgreen,
    key=key, symkey=skey, keysize=keysize,
    labCol=NA,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=margins,
    breaks=seq(ll, ul, length.out=101)
  )
  dev.off()
  
  # normalize for housekeeping genes
  hk_mean <- t(colMeans(exp_t_sub[as.character(signatures[[length(signatures)]][,1]),], na.rm = T))
  exp_t_sub_norm <- exp_t_sub
  for(i in 1:ncol(exp_t_sub_norm)){
    exp_t_sub_norm[,i] <- exp_t_sub_norm[,i] - hk_mean[i]
  }
  
  # heatmap limits
  ul <- -10
  ll <- 10
  
  # plot heatmap
  pdf(file=paste(dir_o, 'exp heatmaps/norm_hk/', cancer, '.pdf', sep=''), width=hm_w, height=hm_h)
  hm <- heatmap.2(
    data.matrix(exp_t_sub_norm),
    #main=paste("1. Mean -", m),
    Rowv=NA,
    dendrogram="none",
    col=redgreen,
    key=key, symkey=skey, keysize=keysize,
    labCol=NA,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=margins,
    breaks=seq(ll, ul, length.out=101)
  )
  dev.off()
  
}

###############################################################################################
# Correlation heatmaps (single genes)
###############################################################################################

cor_hm <- function(cancer_index){
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
    col=redgreen,
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
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
      col=redgreen,
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
      col=redgreen,
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
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
  cancer <- cancers[cancer_index]
  
  # plot heatmap
  pdf(file=paste(dir_o, 'xCell heatmaps/', cancer, '.pdf', sep=''), width=12, height=12)
  hm <- heatmap.2(
    data.matrix(xc[[cancer_index]]),
    #main=paste("1. Mean -", m),
    col=redgreen,
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
    col=redgreen,
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
    col=redgreen,
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
    col=redgreen,
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
    col=redgreen,
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
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
    col=redgreen,
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
    col=redgreen,
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
    col=redgreen,
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
    col=redgreen,
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
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
    pdf(file=paste(dir_o, 'heatmaps master individual genes/', this.gene, '_', sig_names[i], '.pdf', sep=''), width=6, height=18)
    hm <- heatmap.2(
      mat.cor.filtered,
      Rowv=F, Colv=F,
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(10,10),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # plot heatmaps (signature clustered)
    pdf(file=paste(dir_o, 'heatmaps master individual genes/signature clustered/', this.gene, '_', sig_names[i], '.pdf', sep=''), width=6, height=18)
    hm <- heatmap.2(
      mat.cor.filtered,
      Colv=F,
      col=redgreen,
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
      col=redgreen,
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

if(win){ ### WINDOWS ##########################################################################
  
  # Import parameters
  print('Importing parameters')
  t_last <- Sys.time()
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
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  exp_n_all <- foreach(i=1:length(cancers)) %dopar% exp_import(cancers[i], 'n')
  stopCluster(cl)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  exp_t_all <- foreach(i=1:length(cancers)) %dopar% exp_import(cancers[i], 't')
  stopCluster(cl)
  
  # Import xCell data
  xc <- import_xcell()
  
  # Rearrange all samples across cancers in one matrix
  exp_t_full <- exp_t_all[[1]]
  for (i in 2:length(exp_t_all)) {
    exp_t_full <- cbind(exp_t_full, exp_t_all[[i]])
  }
  
  # Median expression of all cancers
  exp_med <- matrix(nrow = length(cancers), ncol = 2)
  for (i in 1:length(cancers)) {
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
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  meps_all <- foreach(i=1:length(cancers)) %dopar% ebp_sign(i)
  stopCluster(cl)
  
  # Generate expression boxplots across cancers
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Generating expression box plots across cancers')
  cores <- min(detectCores(), length(signatures))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:length(signatures)) %dopar% ebp_cancers(i)
  stopCluster(cl)
  
  # Generate expression heatmaps
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Generating expression heat maps')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:length(cancers)) %dopar% exp_hm(i)
  stopCluster(cl)
  
  # Generate single gene correlation heatmaps
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Generating single gene correlation heatmaps')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  mean_cor <- foreach(i=1:length(cancers)) %dopar% cor_hm(i)
  stopCluster(cl)
  
  # Generate signature mean correlation heatmaps
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Generating signature mean correlation heatmaps')
  cor_mean_hm(mean_cor)
  
  # xCell
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('xCell')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  xc_out <- foreach(i=1:length(cancers)) %dopar% xcell(i)
  stopCluster(cl)
  xcell_all(xc_out)
  
  # Generate heatmaps of correlation of every gene in master signature with others
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Generating heatmaps of correlation of every gene in master signature with others')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:nrow(signatures[[1]])) %dopar% hm_master(i)
  stopCluster(cl)
  
  # Terminate
  print(Sys.time()-t_last)
  
}else{ ##### Mac OS ###########################################################################
  
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
  
  # Rearrange all samples across cancers in one matrix
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
  
  # Generate expression heatmaps
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Generating expression heat maps')
  f <- vector("list", length=length(cancers))
  for(i in 1:length(cancers)){
    f[[i]] <- mcparallel(exp_hm(i))
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
  
  # Terminate
  print(Sys.time()-t_last)
}

###############################################################################################
# Termination
###############################################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
