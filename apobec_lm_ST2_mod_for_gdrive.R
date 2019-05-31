# apobec_lm.R

# Jiachen Liang
# Modified by Sakura Maezono

# added ifelse to make sure OV doesn't crash (considering it has
# no normal tissues)

# Compared to apobec_lm_ST_mod_for_gdrive.R,
# this uses GDC and LOGGED


# Generate linear models for APOBEC

# SINGLE THREADED

# Inputs:
#   List of cancers
#   Expression data (tumour and normal)

# Outputs:
#   For each cancer type
#     For tumour and normal
#       PCA rotation table
#       PCA rotated data
#       PCA sdev bar plot
#       Scatter plots for original data vs PCs (with sdev > 1)
#       Coefficients and residuals of simple LM with PCA (filtered for PC's sdev > 1) with k-fold cross-validation
#       Work in progress: Coefficients and residuals of penalized LM with k-fold cross-validation

###############################################################################################
# Initiation
###############################################################################################

t_start <- Sys.time()
t_last <- Sys.time()

# Using Windows?
win <- T

###############################################################################################
# Associated library functions
###############################################################################################

library(matrixStats)
library(ggplot2)
library(DAAG)

if(win) {
  library(doParallel)
} else {
  library(parallel)
}

###############################################################################################
# File locations
###############################################################################################

if(win) {
  # Set working directory
  root <- "G:/My Drive/Jiachen Files/"
  setwd(root)
  
  # Directory: output
  dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2')
  dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/'
  
  # Directory: parameters
  dir_p <- 'G:/My Drive/Jiachen Files/data/'
  
  # Directory: expression matrices
  dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'
} else {
  # Set working directory
  root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
  setwd(root)
  
  # Directory: output
  dir_o <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/lm/'
  
  # Directory: parameters
  dir_p <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'
}


# create directories
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/normal')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/normal/PCA')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/tumour')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/tumour/PCA')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/normal/scatterplots')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/tumour/scatterplots')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/normal/LM')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/tumour/LM')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/normal/hyp vs pc')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/lm_ST2/tumour/hyp vs pc')
###############################################################################################
# Parameters
###############################################################################################

# Cancer list
cancers <- as.character(read.table(paste(dir_p, 
                                         'cancer_list.txt', sep = ''),
                                   header = T, sep = '\t')[,1])

# APOBEC and hypoxia genes
apo.genes <- as.character(read.table(file = "data/tr_a.txt", sep = "\t", 
                                     header = TRUE, row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file = "data/tr_h.txt", sep = "\t", 
                                     header = TRUE, row.names = 1)$Gene.Symbol)

# Minimum proportion of non-NA samples
thr_na <- 0.8

# Folds for cross-validation
k_folds <- 5

# Restrict to A3 family?
restrict <- F

###############################################################################################
# Import expression matrix from file system and preprocess
###############################################################################################

exp_import <- function(cancer) {
  
  library(matrixStats)
  if (cancer == "OV") {
    # no normal tissue samples
    dir_t <- paste(dir_e, cancer, '_t.txt', sep = '')
    input <- as.matrix(read.table(dir_t, sep = '\t', header = T, row.names = NULL))
    rn <- input[,1]
    rn <- make.unique(rn, sep = '.')
    output <- as.data.frame(input[,-1])
    rownames(output) <- rn
    names <- colnames(output)
    for(i in 1:length(names)) {
      split <- strsplit(names[i], '\\.')
      names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
    }
    colnames(output) <- names
    output <- log2(data.matrix(output))
    # according to Francesca, log the data
    input_t <- output
    
    hyp_t <- input_t[hyp.genes,]
    
    input_t <- input_t[apo.genes,]
    
    input_t <- input_t[order(rownames(input_t)),]
    
    if(restrict) {
      input_t <- input_t[c('APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 
                           'APOBEC3F', 'APOBEC3G', 'APOBEC3H'),]
    }
    
    # Replace NA by median of gene
    med_t <- rowMedians(input_t, na.rm = T)
    for(i in 1:nrow(input_t)) {
      for(j in 1:ncol(input_t)) {
        if(is.na(input_t[i,j])) {
          input_t[i,j] <- med_t[i]
        }
      }
    }
    
    output <- list(input_t, hyp_t)
  } else {
    dir_n <- paste(dir_e, cancer, '_n.txt', sep = '')
    input <- as.matrix(read.table(dir_n, sep = '\t', 
                                  header = T, row.names = NULL))
    rn <- input[,1]
    rn <- make.unique(rn, sep = '.')
    output <- as.data.frame(input[,-1])
    rownames(output) <- rn
    names <- colnames(output)
    for(i in 1:length(names)) {
      split <- strsplit(names[i], '\\.')
      names[i] <- paste(split[[1]][1], split[[1]][2], 
                        split[[1]][3], sep = '.')
    }
    colnames(output) <- names
    output <- log2(data.matrix(output))
    input_n <- output
    
    dir_t <- paste(dir_e, cancer, '_t.txt', sep = '')
    input <- as.matrix(read.table(dir_t, sep = '\t', 
                                  header = T, row.names = NULL))
    rn <- input[,1]
    rn <- make.unique(rn, sep = '.')
    output <- as.data.frame(input[,-1])
    rownames(output) <- rn
    names <- colnames(output)
    for(i in 1:length(names)) {
      split <- strsplit(names[i], '\\.')
      names[i] <- paste(split[[1]][1], split[[1]][2], 
                        split[[1]][3], sep = '.')
    }
    colnames(output) <- names
    output <- log2(data.matrix(output))
    input_t <- output
    
    hyp_n <- input_n[hyp.genes,]
    hyp_t <- input_t[hyp.genes,]
    
    input_n <- input_n[apo.genes,]
    input_t <- input_t[apo.genes,]
    
    input_n <- input_n[order(rownames(input_n)),]
    input_t <- input_t[order(rownames(input_t)),]
    
    if(restrict) {
      input_n <- input_n[c('APOBEC3A', 'APOBEC3B', 'APOBEC3C', 
                           'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H'),]
      input_t <- input_t[c('APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H'),]
    }
    # Filter out genes non-NA in less than x proportion of samples
    na_filter <-  
      (rowCounts(!is.na(input_n))/ncol(input_n)>= thr_na) &
      (rowCounts(!is.na(input_t))/ncol(input_t)>= thr_na)
    input_n <- input_n[na_filter,]
    input_t <- input_t[na_filter,]
    
    # Replace NA by median of gene
    med_n <- rowMedians(input_n, na.rm = T)
    med_t <- rowMedians(input_t, na.rm = T)
    for(i in 1:nrow(input_n)) {
      for(j in 1:ncol(input_n)) {
        if(is.na(input_n[i,j])) {
          input_n[i,j] <- med_n[i]
        }
      }
    }
    for(i in 1:nrow(input_t)) {
      for(j in 1:ncol(input_t)) {
        if(is.na(input_t[i,j])) {
          input_t[i,j] <- med_t[i]
        }
      }
    }
    
    output <- list(input_n, input_t, hyp_n, hyp_t)
  }
  
  return(output)
  
}

if (win) { ### WINDOWS ##########################################################################
  
  # Import expression matrices (parallel)
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Importing expression matrices')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  exp_all <- foreach(i = 1:length(cancers)) %dopar% exp_import(cancers[i])
  stopCluster(cl)
  
} else {##### Mac OS ###########################################################################
  
  # Import expression matrices (parallel)
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Importing expression matrices')
  f <- vector("list", length = length(cancers))
  for(i in 1:length(cancers)) {
    f[[i]] <- mcparallel(exp_import(cancers[i]))
  }
  exp_all <- mccollect(f)
  
}

###############################################################################################
# PCA and linear model
###############################################################################################

for(cancer_index in 1:length(cancers)) {
  
  cancer <- cancers[cancer_index]
  if (cancer == "OV") {
    # Get expression matrices
    mat_t <- exp_all[[cancer_index]][[1]]
    hyp_t <- exp_all[[cancer_index]][[2]]
    
    
    # PCA
    pca_t <- prcomp(mat_t)
    
    # Output rotation data and rotated values
    write.table(pca_t$rotation, file=paste(dir_o, 'tumour/PCA/', cancer, '_rotation.txt', sep=''), sep='\t')
    write.table(pca_t$x, file=paste(dir_o, 'tumour/PCA/', cancer, '_x.txt', sep=''), sep='\t')
    
    # Output bar plot of sdev -------------------------------------------
    
    sdev_df_t <- as.data.frame(pca_t$sdev)
    colnames(sdev_df_t) <- 'Standard Deviation'
    plot_t <- ggplot(data=sdev_df_t, aes(x=rownames(sdev_df_t), y=`Standard Deviation`)) +
      geom_bar(stat='identity') +
      xlab('Principal Component') +
      ggtitle(paste('Standard deviation of principal components -', cancer, '(Tumour)')) +
      theme(text = element_text(size=12))
    
    pdf(file=paste(dir_o, 'tumour/PCA/', cancer, '_sdev.pdf', sep=''), width=6, height=4)
    print(plot_t)
    dev.off()
    
    # Scatter plots for original data vs PCs ----------------------------
    
    # filter for PCs with sdev > 1
    good_pc_t <- pca_t$rotation[,pca_t$sdev>1]
    
    # correlation score and p-value matrices
    cor_t <- matrix(nrow=nrow(mat_t), ncol=ncol(good_pc_t))
    pv_t <- matrix(nrow=nrow(mat_t), ncol=ncol(good_pc_t))
    rownames(cor_t) <- rownames(mat_t)
    rownames(pv_t) <- rownames(mat_t)
    colnames(cor_t) <- colnames(good_pc_t)
    colnames(pv_t) <- colnames(good_pc_t)
    
    
    # TUMOUR
    for(i in 1:nrow(mat_t)){ # for every APOBEC
      for(j in 1:ncol(good_pc_t)){ # for every PC
        name <- paste(rownames(mat_t)[i], 'vs', colnames(good_pc_t)[j])
        dir <- paste(dir_o, 'tumour/scatterplots/', cancer, '_',
                     colnames(good_pc_t)[j], '_', 
                     rownames(mat_t)[i], '.pdf', sep='')
        Original <- c() # original data
        PC <- c() # principal component
        for(k in 1:ncol(mat_t)){ # for every sample
          Original <- c(Original, mat_t[i,k])
          PC <- c(PC, good_pc_t[k,j])
        }
        df <- data.frame(Original, PC)
        # plot
        plot <- ggplot(df, aes(x=Original, y=PC)) +
          geom_point(shape=1) +
          geom_smooth(method=lm , color="red", se=TRUE) +
          ggtitle(paste(rownames(mat_t)[i], 'vs', 
                        colnames(good_pc_t)[j], '-',
                        cancer, '(Tumour)')) +
          xlab(rownames(mat_t)[i]) +
          ylab(colnames(good_pc_t)[j]) +
          theme(text = element_text(size=12))
        pdf(file=dir, width=6, height=4)
        print(plot)
        dev.off()
        # correlation
        cor_t[i,j] <- as.numeric(cor.test(as.numeric(mat_t[i,]), 
                                          as.numeric(good_pc_t[,j]), 
                                          method = 'spearman')$estimate)
        pv_t[i,j] <- as.numeric(cor.test(as.numeric(mat_t[i,]),
                                         as.numeric(good_pc_t[,j]), 
                                         method = 'spearman')$p.value)
      }
      # FDR
      pv_t[i,] <- p.adjust(as.numeric(pv_t[i,]), method='fdr')
    }
    
    # Output rho and pv
    write.table(cor_t, file=paste(dir_o, 'tumour/scatterplots/',
                                  cancer, '_cor.txt', sep=''), sep='\t')
    write.table(pv_t, file=paste(dir_o, 'tumour/scatterplots/', 
                                 cancer, '_pv.txt', sep=''), sep='\t')
    
    # Filter rho with FDR-adjusted pv < 0.05
    cor_t_f <- cor_t
    cor_t_f[pv_t>=0.05] <- NA
    
    # Output filtered rho
    write.table(cor_t_f, file=paste(dir_o, 'tumour/scatterplots/', cancer, '_cor_filtered.txt', sep=''), sep='\t')
    
    # Correlation between PCA and hypoxia -------------------------------
    
    # filter for PCs with sdev > 1
    good_pc_t <- pca_t$rotation[,pca_t$sdev>1]
    
    # correlation score and p-value matrices
    cor_t <- matrix(nrow=ncol(good_pc_t), ncol=1)
    rownames(cor_t) <- colnames(good_pc_t)
    pv_t <- cor_t
    
    # matrices of mean hypoxia score
    mat_t <- t(as.matrix(colMeans(hyp_t, na.rm=T)))
    
    # TUMOUR
    for(j in 1:ncol(good_pc_t)){ # for every PC
      name <- paste('Hypoxia', 'vs', colnames(good_pc_t)[j])
      dir <- paste(dir_o, 'tumour/hyp vs pc/', cancer, '_', colnames(good_pc_t)[j],'.pdf', sep='')
      Original <- c() # original data
      PC <- c() # principal component
      for(k in 1:ncol(mat_t)){ # for every sample
        Original <- c(Original, mat_t[1,k])
        PC <- c(PC, good_pc_t[k,j])
      }
      df <- data.frame(Original, PC)
      # plot
      plot <- ggplot(df, aes(x=Original, y=PC)) +
        geom_point(shape=1) +
        geom_smooth(method=lm , color="red", se=TRUE) +
        ggtitle(paste('Hypoxia', 'vs', colnames(good_pc_t)[j], '-', cancer, '(Tumour)')) +
        xlab('Mean Hypoxia Score') +
        ylab(colnames(good_pc_t)[j]) +
        theme(text = element_text(size=12))
      pdf(file=dir, width=6, height=4)
      print(plot)
      dev.off()
      # correlation
      cor_t[j,1] <- as.numeric(cor.test(as.numeric(mat_t[1,]), 
                                        as.numeric(good_pc_t[,j]),
                                        method = 'spearman')$estimate)
      pv_t[j,1] <- as.numeric(cor.test(as.numeric(mat_t[1,]), 
                                       as.numeric(good_pc_t[,j]), 
                                       method = 'spearman')$p.value)
    }
    # FDR
    pv_t[,1] <- p.adjust(as.numeric(pv_t[,1]), method='fdr')
    
    # Output rho and pv
    write.table(cor_t, file=paste(dir_o, 'tumour/hyp vs pc/', 
                                  cancer, '_cor.txt', sep=''), sep='\t')
    write.table(pv_t, file=paste(dir_o, 'tumour/hyp vs pc/', 
                                 cancer, '_pv.txt', sep=''), sep='\t')
    
    # Filter rho with FDR-adjusted pv < 0.05
    cor_t_f <- cor_t
    cor_t_f[pv_t>=0.05] <- NA
    
    # Output filtered rho
    write.table(cor_t_f, file=paste(dir_o, 'tumour/hyp vs pc/', 
                                    cancer, '_cor_filtered.txt', sep=''), sep='\t')
    
    # LM with PCA -------------------------------------------------------
    
    # Data
    data_t <- as.data.frame(cbind(good_pc_t, rowMeans(t(hyp_t), na.rm=T)))
    colnames(data_t) <- c(colnames(good_pc_t), 'Hypoxia')
    
    # Formulae
    f_t <- paste('Hypoxia ~', paste(colnames(good_pc_t), collapse=' + '))
    
    # LM
    lm_t <- lm(f_t, data=data_t)
    
    # Directories
    dir_t <- paste(dir_o, 'tumour/LM/', cancer, '_', sep='')
    
    # Coefficients
    lm_t_c <- as.matrix(lm_t$coefficients)
    colnames(lm_t_c) <- cancer
    write.table(lm_t_c, file=paste(dir_t, 'LM_coef.txt', sep=''), sep='\t')
    
    # Significance
    sink(paste(dir_t, 'LM_summary.txt', sep=''))
    print(summary(lm_t))
    sink()
    
    # Residuals
    lm_t_r <- as.matrix(lm_t$residuals)
    colnames(lm_t_r) <- cancer
    write.table(lm_t_r, file=paste(dir_t, 'LM_res.txt', sep=''), sep='\t')
    
    # Plot residuals
    pdf(file=paste(dir_t, 'LM_res_hist.pdf', sep=''), width=6, height=4)
    print(
      hist(lm_t$residuals, main=paste("Residuals Histogram -", cancer, '(Tumour)'),
           ylab="Frequency", xlab="Value")
    )
    print(
      qqnorm(lm_t$residuals, main=paste("Residuals Q-Q plot -", cancer, '(Tumour)'))
    )
    print(
      qqline(lm_t$residuals)
    )
    dev.off()
    
    # What about analysis of error and variance?
    # http://www.learnbymarketing.com/tutorials/linear-regression-in-r/
    
    # k-fold cross-validation
    pdf(file=paste(dir_t, 'LM_cv.pdf', sep=''), width=6, height=6)
    print(cv.lm(data=data_t, lm_t, m=k_folds))
    dev.off()
    
    cv_t <- cv.lm(data=data_t, lm_t, m=k_folds)
    
    write.table(cv_t, file=paste(dir_t, 'LM_cv.txt', sep=''), sep='\t')
    
    # Penalized LM ------------------------------------------------------
    
  } else {
    # Get expression matrices
    mat_n <- exp_all[[cancer_index]][[1]]
    mat_t <- exp_all[[cancer_index]][[2]]
    hyp_n <- exp_all[[cancer_index]][[3]]
    hyp_t <- exp_all[[cancer_index]][[4]]
    
    # PCA
    pca_n <- prcomp(mat_n)
    pca_t <- prcomp(mat_t)
    
    # Output rotation data and rotated values
    write.table(pca_n$rotation, file=paste(dir_o, 'normal/PCA/', cancer, '_rotation.txt', sep=''), sep='\t')
    write.table(pca_t$rotation, file=paste(dir_o, 'tumour/PCA/', cancer, '_rotation.txt', sep=''), sep='\t')
    write.table(pca_n$x, file=paste(dir_o, 'normal/PCA/', cancer, '_x.txt', sep=''), sep='\t')
    write.table(pca_t$x, file=paste(dir_o, 'tumour/PCA/', cancer, '_x.txt', sep=''), sep='\t')
    
    # Output bar plot of sdev -------------------------------------------
    
    sdev_df_n <- as.data.frame(pca_n$sdev)
    sdev_df_t <- as.data.frame(pca_t$sdev)
    colnames(sdev_df_n) <- 'Standard Deviation'
    colnames(sdev_df_t) <- 'Standard Deviation'
    
    plot_n <- ggplot(data=sdev_df_n, aes(x=rownames(sdev_df_n), y=`Standard Deviation`)) +
      geom_bar(stat='identity') +
      xlab('Principal Component') +
      ggtitle(paste('Standard deviation of principal components -', cancer, '(Normal)')) +
      theme(text = element_text(size=12))
    plot_t <- ggplot(data=sdev_df_t, aes(x=rownames(sdev_df_t), y=`Standard Deviation`)) +
      geom_bar(stat='identity') +
      xlab('Principal Component') +
      ggtitle(paste('Standard deviation of principal components -', cancer, '(Tumour)')) +
      theme(text = element_text(size=12))
    
    pdf(file=paste(dir_o, 'normal/PCA/', cancer, '_sdev.pdf', sep=''), width=6, height=4)
    print(plot_n)
    dev.off()
    pdf(file=paste(dir_o, 'tumour/PCA/', cancer, '_sdev.pdf', sep=''), width=6, height=4)
    print(plot_t)
    dev.off()
    
    # Scatter plots for original data vs PCs ----------------------------
    
    # filter for PCs with sdev > 1
    good_pc_n <- as.matrix(pca_n$rotation[,pca_n$sdev>1])
    good_pc_t <- as.matrix(pca_t$rotation[,pca_t$sdev>1])
    
    # correlation score and p-value matrices
    cor_n <- matrix(nrow=nrow(mat_n), ncol=ncol(good_pc_n))
    pv_n <- matrix(nrow=nrow(mat_n), ncol=ncol(good_pc_n))
    cor_t <- matrix(nrow=nrow(mat_t), ncol=ncol(good_pc_t))
    pv_t <- matrix(nrow=nrow(mat_t), ncol=ncol(good_pc_t))
    rownames(cor_n) <- rownames(mat_n)
    rownames(pv_n) <- rownames(mat_n)
    colnames(cor_n) <- colnames(good_pc_n)
    colnames(pv_n) <- colnames(good_pc_n)
    rownames(cor_t) <- rownames(mat_t)
    rownames(pv_t) <- rownames(mat_t)
    colnames(cor_t) <- colnames(good_pc_t)
    colnames(pv_t) <- colnames(good_pc_t)
    
    # NORMAL
    for(i in 1:nrow(mat_n)){ # for every APOBEC
      for(j in 1:ncol(good_pc_n)){ # for every PC
        name <- paste(rownames(mat_n)[i], 'vs', colnames(good_pc_n)[j])
        dir <- paste(dir_o, 'normal/scatterplots/', cancer, '_', 
                     colnames(good_pc_n)[j], '_', rownames(mat_n)[i],
                     '.pdf', sep='')
        Original <- c() # original data
        PC <- c() # principal component
        for(k in 1:ncol(mat_n)){ # for every sample
          Original <- c(Original, mat_n[i,k])
          PC <- c(PC, good_pc_n[k,j])
        }
        df <- data.frame(Original, PC)
        # plot
        plot <- ggplot(df, aes(x=Original, y=PC)) +
          geom_point(shape=1) +
          geom_smooth(method=lm , color="red", se=TRUE) +
          ggtitle(paste(rownames(mat_n)[i], 'vs', colnames(good_pc_n)[j], '-', cancer, '(Normal)')) +
          xlab(rownames(mat_n)[i]) +
          ylab(colnames(good_pc_n)[j]) +
          theme(text = element_text(size=12))
        pdf(file=dir, width=6, height=4)
        print(plot)
        dev.off()
        # correlation
        cor_n[i,j] <- as.numeric(cor.test(as.numeric(mat_n[i,]), as.numeric(good_pc_n[,j]), method = 'spearman')$estimate)
        pv_n[i,j] <- as.numeric(cor.test(as.numeric(mat_n[i,]), as.numeric(good_pc_n[,j]), method = 'spearman')$p.value)
      }
      # FDR
      pv_n[i,] <- p.adjust(as.numeric(pv_n[i,]), method='fdr')
    }
    # TUMOUR
    for(i in 1:nrow(mat_t)){ # for every APOBEC
      for(j in 1:ncol(good_pc_t)){ # for every PC
        name <- paste(rownames(mat_t)[i], 'vs', colnames(good_pc_t)[j])
        dir <- paste(dir_o, 'tumour/scatterplots/', cancer, '_', colnames(good_pc_t)[j], '_', rownames(mat_t)[i], '.pdf', sep='')
        Original <- c() # original data
        PC <- c() # principal component
        for(k in 1:ncol(mat_t)){ # for every sample
          Original <- c(Original, mat_t[i,k])
          PC <- c(PC, good_pc_t[k,j])
        }
        df <- data.frame(Original, PC)
        # plot
        plot <- ggplot(df, aes(x=Original, y=PC)) +
          geom_point(shape=1) +
          geom_smooth(method=lm , color="red", se=TRUE) +
          ggtitle(paste(rownames(mat_t)[i], 'vs', colnames(good_pc_t)[j], '-', cancer, '(Tumour)')) +
          xlab(rownames(mat_t)[i]) +
          ylab(colnames(good_pc_t)[j]) +
          theme(text = element_text(size=12))
        pdf(file=dir, width=6, height=4)
        print(plot)
        dev.off()
        # correlation
        cor_t[i,j] <- as.numeric(cor.test(as.numeric(mat_t[i,]), as.numeric(good_pc_t[,j]), method = 'spearman')$estimate)
        pv_t[i,j] <- as.numeric(cor.test(as.numeric(mat_t[i,]), as.numeric(good_pc_t[,j]), method = 'spearman')$p.value)
      }
      # FDR
      pv_t[i,] <- p.adjust(as.numeric(pv_t[i,]), method='fdr')
    }
    
    # Output rho and pv
    write.table(cor_n, file=paste(dir_o, 'normal/scatterplots/', cancer, '_cor.txt', sep=''), sep='\t')
    write.table(pv_n, file=paste(dir_o, 'normal/scatterplots/', cancer, '_pv.txt', sep=''), sep='\t')
    write.table(cor_t, file=paste(dir_o, 'tumour/scatterplots/', cancer, '_cor.txt', sep=''), sep='\t')
    write.table(pv_t, file=paste(dir_o, 'tumour/scatterplots/', cancer, '_pv.txt', sep=''), sep='\t')
    
    # Filter rho with FDR-adjusted pv < 0.05
    cor_n_f <- cor_n
    cor_t_f <- cor_t
    cor_n_f[pv_n>=0.05] <- NA
    cor_t_f[pv_t>=0.05] <- NA
    
    # Output filtered rho
    write.table(cor_n_f, file=paste(dir_o, 'normal/scatterplots/', cancer, '_cor_filtered.txt', sep=''), sep='\t')
    write.table(cor_t_f, file=paste(dir_o, 'tumour/scatterplots/', cancer, '_cor_filtered.txt', sep=''), sep='\t')
    
    # Correlation between PCA and hypoxia -------------------------------
    
    # filter for PCs with sdev > 1
    good_pc_n <- as.matrix(pca_n$rotation[,pca_n$sdev>1])
    good_pc_t <- as.matrix(pca_t$rotation[,pca_t$sdev>1])
    
    # correlation score and p-value matrices
    cor_n <- matrix(nrow=ncol(good_pc_n), ncol=1)
    rownames(cor_n) <- colnames(good_pc_n)
    pv_n <- cor_n
    cor_t <- matrix(nrow=ncol(good_pc_t), ncol=1)
    rownames(cor_t) <- colnames(good_pc_t)
    pv_t <- cor_t
    
    # matrices of mean hypoxia score
    mat_n <- t(as.matrix(colMeans(hyp_n, na.rm=T)))
    mat_t <- t(as.matrix(colMeans(hyp_t, na.rm=T)))
    
    # NORMAL
    for(j in 1:ncol(good_pc_n)){ # for every PC
      name <- paste('Hypoxia', 'vs', colnames(good_pc_n)[j])
      dir <- paste(dir_o, 'normal/hyp vs pc/', cancer, '_', colnames(good_pc_n)[j],'.pdf', sep='')
      Original <- c() # original data
      PC <- c() # principal component
      for(k in 1:ncol(mat_n)){ # for every sample
        Original <- c(Original, mat_n[1,k])
        PC <- c(PC, good_pc_n[k,j])
      }
      df <- data.frame(Original, PC)
      # plot
      plot <- ggplot(df, aes(x=Original, y=PC)) +
        geom_point(shape=1) +
        geom_smooth(method=lm , color="red", se=TRUE) +
        ggtitle(paste('Hypoxia', 'vs', colnames(good_pc_n)[j], '-', cancer, '(Normal)')) +
        xlab('Mean Hypoxia Score') +
        ylab(colnames(good_pc_n)[j]) +
        theme(text = element_text(size=12))
      pdf(file=dir, width=6, height=4)
      print(plot)
      dev.off()
      # correlation
      cor_n[j,1] <- as.numeric(cor.test(as.numeric(mat_n[1,]), as.numeric(good_pc_n[,j]), method = 'spearman')$estimate)
      pv_n[j,1] <- as.numeric(cor.test(as.numeric(mat_n[1,]), as.numeric(good_pc_n[,j]), method = 'spearman')$p.value)
    }
    # FDR
    pv_n[,1] <- p.adjust(as.numeric(pv_n[,1]), method='fdr')
    # TUMOUR
    for(j in 1:ncol(good_pc_t)){ # for every PC
      name <- paste('Hypoxia', 'vs', colnames(good_pc_t)[j])
      dir <- paste(dir_o, 'tumour/hyp vs pc/', cancer, '_', colnames(good_pc_t)[j],'.pdf', sep='')
      Original <- c() # original data
      PC <- c() # principal component
      for(k in 1:ncol(mat_t)){ # for every sample
        Original <- c(Original, mat_t[1,k])
        PC <- c(PC, good_pc_t[k,j])
      }
      df <- data.frame(Original, PC)
      # plot
      plot <- ggplot(df, aes(x=Original, y=PC)) +
        geom_point(shape=1) +
        geom_smooth(method=lm , color="red", se=TRUE) +
        ggtitle(paste('Hypoxia', 'vs', colnames(good_pc_t)[j], '-', cancer, '(Tumour)')) +
        xlab('Mean Hypoxia Score') +
        ylab(colnames(good_pc_t)[j]) +
        theme(text = element_text(size=12))
      pdf(file=dir, width=6, height=4)
      print(plot)
      dev.off()
      # correlation
      cor_t[j,1] <- as.numeric(cor.test(as.numeric(mat_t[1,]), as.numeric(good_pc_t[,j]), method = 'spearman')$estimate)
      pv_t[j,1] <- as.numeric(cor.test(as.numeric(mat_t[1,]), as.numeric(good_pc_t[,j]), method = 'spearman')$p.value)
    }
    # FDR
    pv_t[,1] <- p.adjust(as.numeric(pv_t[,1]), method='fdr')
    
    # Output rho and pv
    write.table(cor_n, file=paste(dir_o, 'normal/hyp vs pc/', cancer, '_cor.txt', sep=''), sep='\t')
    write.table(pv_n, file=paste(dir_o, 'normal/hyp vs pc/', cancer, '_pv.txt', sep=''), sep='\t')
    write.table(cor_t, file=paste(dir_o, 'tumour/hyp vs pc/', cancer, '_cor.txt', sep=''), sep='\t')
    write.table(pv_t, file=paste(dir_o, 'tumour/hyp vs pc/', cancer, '_pv.txt', sep=''), sep='\t')
    
    # Filter rho with FDR-adjusted pv < 0.05
    cor_n_f <- cor_n
    cor_t_f <- cor_t
    cor_n_f[pv_n>=0.05] <- NA
    cor_t_f[pv_t>=0.05] <- NA
    
    # Output filtered rho
    write.table(cor_n_f, file=paste(dir_o, 'normal/hyp vs pc/', cancer, '_cor_filtered.txt', sep=''), sep='\t')
    write.table(cor_t_f, file=paste(dir_o, 'tumour/hyp vs pc/', cancer, '_cor_filtered.txt', sep=''), sep='\t')
    
    # LM with PCA -------------------------------------------------------
    
    # Data
    data_n <- as.data.frame(cbind(good_pc_n, rowMeans(t(hyp_n), na.rm=T)))
    data_t <- as.data.frame(cbind(good_pc_t, rowMeans(t(hyp_t), na.rm=T)))
    
    if (ncol(good_pc_n) == 1) {
      colnames(data_n) <- c("PC1", 'Hypoxia')
    } else {
      colnames(data_n) <- c(colnames(good_pc_n), 'Hypoxia')
    }
      colnames(data_t) <- c(colnames(good_pc_t), 'Hypoxia')
    # Formulae
    if (ncol(good_pc_n) == 1) {
      f_n <- paste('Hypoxia ~', paste("PC1", collapse=' + '))
    } else {
      f_n <- paste('Hypoxia ~', paste(colnames(good_pc_n), collapse=' + '))
    }
    f_t <- paste('Hypoxia ~', paste(colnames(good_pc_t), collapse=' + '))
    
    # LM
    lm_n <- lm(f_n, data=data_n)
    lm_t <- lm(f_t, data=data_t)
    
    # Directories
    dir_n <- paste(dir_o, 'normal/LM/', cancer, '_', sep='')
    dir_t <- paste(dir_o, 'tumour/LM/', cancer, '_', sep='')
    
    # Coefficients
    lm_n_c <- as.matrix(lm_n$coefficients)
    colnames(lm_n_c) <- cancer
    lm_t_c <- as.matrix(lm_t$coefficients)
    colnames(lm_t_c) <- cancer
    write.table(lm_n_c, file=paste(dir_n, 'LM_coef.txt', sep=''), sep='\t')
    write.table(lm_t_c, file=paste(dir_t, 'LM_coef.txt', sep=''), sep='\t')
    
    # Significance
    sink(paste(dir_n, 'LM_summary.txt', sep=''))
    print(summary(lm_n))
    sink()
    sink(paste(dir_t, 'LM_summary.txt', sep=''))
    print(summary(lm_t))
    sink()
    
    # Residuals
    lm_n_r <- as.matrix(lm_n$residuals)
    colnames(lm_n_r) <- cancer
    lm_t_r <- as.matrix(lm_t$residuals)
    colnames(lm_t_r) <- cancer
    write.table(lm_n_r, file=paste(dir_n, 'LM_res.txt', sep=''), sep='\t')
    write.table(lm_t_r, file=paste(dir_t, 'LM_res.txt', sep=''), sep='\t')
    
    # Plot residuals
    pdf(file=paste(dir_n, 'LM_res_hist.pdf', sep=''), width=6, height=4)
    print(
      hist(lm_n$residuals, main=paste("Residuals Histogram -", cancer, '(Normal)'),
           ylab="Frequency", xlab="Value")
    )
    print(
      qqnorm(lm_n$residuals, main=paste("Residuals Q-Q plot -", cancer, '(Normal)'))
    )
    print(
      qqline(lm_n$residuals)
    )
    dev.off()
    pdf(file=paste(dir_t, 'LM_res_hist.pdf', sep=''), width=6, height=4)
    print(
      hist(lm_t$residuals, main=paste("Residuals Histogram -", cancer, '(Tumour)'),
           ylab="Frequency", xlab="Value")
    )
    print(
      qqnorm(lm_t$residuals, main=paste("Residuals Q-Q plot -", cancer, '(Tumour)'))
    )
    print(
      qqline(lm_t$residuals)
    )
    dev.off()
    
    # What about analysis of error and variance?
    # http://www.learnbymarketing.com/tutorials/linear-regression-in-r/
    
    # k-fold cross-validation
    pdf(file=paste(dir_n, 'LM_cv.pdf', sep=''), width=6, height=6)
    print(cv.lm(data=data_n, lm_n, m=k_folds))
    dev.off()
    pdf(file=paste(dir_t, 'LM_cv.pdf', sep=''), width=6, height=6)
    print(cv.lm(data=data_t, lm_t, m=k_folds))
    dev.off()
    
    cv_n <- cv.lm(data=data_n, lm_n, m=k_folds)
    cv_t <- cv.lm(data=data_t, lm_t, m=k_folds)
    
    write.table(cv_n, file=paste(dir_n, 'LM_cv.txt', sep=''), sep='\t')
    write.table(cv_t, file=paste(dir_t, 'LM_cv.txt', sep=''), sep='\t')
    
    # Penalized LM ------------------------------------------------------
  }
}

###############################################################################################
# Termination
###############################################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
