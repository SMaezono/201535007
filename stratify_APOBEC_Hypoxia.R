# Sakura Maezono # Date created: ##------ Tue Apr 30 21:06:59 2019 ------##
# Purpose: 
# stratify samples into high/low status quadrants: 
# high/low hypoxia vs high/low APOBEC

# NOT FINISHED, new version is stratify_into_Hypoxia_quadrants_APOBEC


###############################################################################################
# Initiation
###############################################################################################

t_start <- Sys.time()
t_last <- Sys.time()

###############################################################################################
# Associated library functions
###############################################################################################
# Windows or Mac
win <- TRUE

library(matrixStats)
if (win) {
  library(doParallel)
}else{
  library(parallel)
}
library(gplots)


# create folder(s) if it doesn't/ they don't exist

dir.create("G:/My Drive/Jiachen Files/Sakura/Hypoxia_APOBEC_Stratification/")

# Set working directory
root <- "G:/My Drive/Jiachen Files/Sakura/Hypoxia_APOBEC_Stratification/"
setwd(root)

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'


# APOBEC and hypoxia genes
dir_ext <- "G:/My Drive/Jiachen Files/"
apo.genes <- as.character(read.table(file = paste(dir_ext, 
                                                  "data/tr_a.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file = paste(dir_ext,
                                                  "data/tr_h.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)
# Minimum proportion of non-NA samples
thr_na <- 0.8

# "t" for tumour data and "n" for normal data
state <- "t"
 if (state == "t") {
   dir.create("G:/My Drive/Jiachen Files/Sakura/Hypoxia_APOBEC_Stratification/Buffa 52/tumour/")
   dir.create("G:/My Drive/Jiachen Files/Sakura/Hypoxia_APOBEC_Stratification/APOBEC/tumour/")
 } else {
   dir.create("G:/My Drive/Jiachen Files/Sakura/Hypoxia_APOBEC_Stratification/Buffa 52/normal/")
   dir.create("G:/My Drive/Jiachen Files/Sakura/Hypoxia_APOBEC_Stratification/APOBEC/normal/")
 }

# Restrict to A3 family?
restrict <- FALSE

# by quartiles or top and bottom %?
quartiles <- TRUE

if (!quartiles) {
  # How much in % are top and bottom
  percentindecimal <- 0.10
}


# Cancer list
cancers <- as.character(read.table(paste(dir_p, 
                                         'cancer_list.txt',
                                         sep = ''),
                                   header = T, sep = '\t')[,1])


###############################################################################
# Import expression matrix from file system
###############################################################################

import_cancer_exp <- function(cancer, state) {
  # cancer = c("BRCA", etc.)
  # state = c("t", "n") 
  
  dir <- paste(dir_e, cancer, '_', state, '.txt', sep = '')
  input <- read.table(dir, sep = '\t', header = T, row.names = NULL)
  
  rn <- input[,1]
  rn <- make.unique(rn, sep = '.')
  
  output <- as.data.frame(input[,-1])
  rownames(output) <- rn
  
  exp <- input[,-1]
  rownames(exp) <- rn
  
  # Fix col names
  names <- colnames(exp)
  for (i in 1:length(names)) {
    split <- strsplit(names[i], '\\.')
    names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
  }
  colnames(exp) <- names
  # output <- data.frame(output)
  # hypoxia genes only in order
  hyp_t <- as.matrix(output[hyp.genes,])
  # # remove rows that only have NAs
  # hyp_t <- rm_onlyNArow(hyp_t)
  na_filter <-
    (rowCounts(!is.na(hyp_t))/ncol(hyp_t)>=thr_na)
  hyp_t <- hyp_t[na_filter,]

  # Replace NA by median of gene
  # med_t <- rowMedians(hyp_t, na.rm = T)
  # for (i in 1:nrow(hyp_t)) {
  #   for (j in 1:ncol(hyp_t)) {
  #     if (is.na(hyp_t[i,j])) {
  #       hyp_t[i,j] <- med_t[i]
  #     }
  #   }
  # }
  input_t <- as.matrix(output[apo.genes,])
  if (restrict) {
    input_t <- input_t[c('APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
                         'APOBEC3F', 'APOBEC3G', 'APOBEC3H'),]
  }
  # input_t <- data.frame(output[order(rownames(input_t)),])
  # remove rows that only have NAs
  # input_t <- rm_onlyNArow(input_t)
  # Filter out genes non-NA in less than x proportion of samples
  na_filter <-  
    (rowCounts(!is.na(input_t))/ncol(input_t)>=thr_na)
  
  input_t <- input_t[na_filter,]
  
  # # Replace NA by median of gene
  # med_t <- rowMedians(input_t, na.rm = T)
  # for(i in 1:nrow(input_t)) {
  #   for(j in 1:ncol(input_t)) {
  #     if(is.na(input_t[i,j])) {
  #       input_t[i,j] <- med_t[i]
  #     }
  #   }
  # }
  
  # matrices of mean & median scores for every patient
  # hypoxia gene sig
  hyp_t <- t(hyp_t)
  hyp_metrics <- matrix(nrow = length(hyp_t[,1]), ncol = 2 )
  rownames(hyp_metrics) <- rownames(hyp_t)
  colnames(hyp_metrics) <- c("Mean", "Median")
  hyp_metrics[,1] <- rowMeans(as.matrix(hyp_t), na.rm=T)
  hyp_metrics[,2] <- rowMedians(as.matrix(hyp_t), na.rm=T)
  # reorder by mean
  hyp_metrics <- hyp_metrics[order(hyp_metrics[,1]),]
  # APOBEC gene sig
  input_t <- t(input_t)
  apo_metrics <- matrix(nrow = length(input_t[,1]), ncol = 2 )
  rownames(apo_metrics) <- rownames(input_t)
  colnames(apo_metrics) <- c("Mean", "Median")
  apo_metrics[,1] <- rowMeans(as.matrix(input_t), na.rm=T)
  apo_metrics[,2] <- rowMedians(as.matrix(input_t), na.rm=T)
  # reorder by mean
  apo_metrics <- apo_metrics[order(apo_metrics[,1]),]
  
  # Output_all, mean and median of gene sig dor each patient
  write.table(hyp_metrics, 
              file = paste(cancer,
                         '_meanormedian_of_hypoxia_genes.txt', 
                         sep = ''), sep = '\t')
  
  # Output_all, mean and median of gene sig dor each patient
  write.table(apo_metrics, 
              file = paste(cancer,
                           '_meanormedian_of_APOBEC_genes.txt', 
                           sep = ''), sep = '\t')
  # reorder based on median lowest to highest
  hyp_metrics_med <- hyp_metrics[order(hyp_metrics[,2]),]
  apo_metrics_med <- apo_metrics[order(apo_metrics[,2]),]
  
  cancer_output <- list(data.frame(hyp_metrics[,1]), 
                        data.frame(hyp_metrics_med[,2]),
                        data.frame(apo_metrics[,1]),
                        data.frame(apo_metrics_med[,2]))
  return(cancer_output)
}

###############################################################################
# 
###############################################################################
Strat_APOBEC_Hypoxia <- function(cancer) {
  
  cancer_output <- import_cancer_exp(cancer, state)
  
  Hypoxia_mean <- cbind(rownames(cancer_output[[1]]), cancer_output[[1]])
  Hypoxia_median <- cbind(rownames(cancer_output[[2]]), cancer_output[[2]])
  apo_mean <- cbind(rownames(cancer_output[[3]]), cancer_output[[3]])
  apo_median <- cbind(rownames(cancer_output[[4]]), cancer_output[[4]])
  
  
  extract_Top_and_bottom <- function(which_sig) {
    # ceiling () rounds up the values
    if (which_sig == "hypoxia") {
      if (percentindecimal < 1) {
        desired_percent <- ceiling(length(Hypoxia_mean[,1])*percentindecimal) 
      } else {
        desired_percent <- length(Hypoxia_mean[,1])
      }
    } else {
      if (percentindecimal < 1) {
        desired_percent <- ceiling(length(apo_mean[,1])*percentindecimal) 
      } else {
        desired_percent <- length(apo_mean[,1])
      }
    }
    # prints out how many samples were considered
    print(paste(desired_percent, 
                "sample ids were used for the top",
                "or bottom group based on", which_sig,"gene sig score",
                sep = " "))
    
    if (which_sig == "hypoxia") {
      mean <- Hypoxia_mean
      median <- Hypoxia_median
    } else {
      mean <- apo_mean
      median <- apo_median
    }
    ### Hypoxia gene sig ###
    # lowest to highest (default) order used to get top
    # order the patients based on mean/median of the hypoxia signature
    # and take top ?? %
    top_mean <- data.frame(Hypoxia_mean[order(mean[,2], 
                                              decreasing = TRUE),])
    top_mean <- data.frame(top_mean[1:desired_percent,]) 
    top_mean_patients <- top_mean[,1]
    
    top_median <- data.frame(median[order(Hypoxia_median[,2], 
                                          decreasing = TRUE),])
    top_median <- data.frame(top_median[1:desired_percent,]) 
    top_median_patients <- top_median[,1]
    
    # lowest to highest (default) order used to get bottom
    bottom_mean <- data.frame(mean[1:desired_percent,]) 
    bottom_median <- data.frame(median[1:desired_percent,]) 
    bottom_mean_patients <- bottom_mean[,1]
    bottom_median_patients <- bottom_median[,1]
    
    output <- list(top_mean, bottom_mean, top_median, bottom_median,
                   top_mean_patients, bottom_mean_patients,
                   top_median_patients, bottom_median_patients)
    
    return(output)
  }
  
  Strat_Hypoxia <- extract_Top_and_bottom("hypoxia")
  Strat_Apobec <-  extract_Top_and_bottom("apobec")

 
  # Output rotation data and rotated values (Hypoxia)
  

  # per patient (metric scores of each patient including all genes)
  output <- as.matrix(cancer_output[[4]])
  out_mean <- data.frame(colMeans(output, na.rm = T))
  out_median <- data.frame(colMedians(output, na.rm = T))
  output_metric_score <- cbind(out_mean, out_median)
  colnames(output_metric_score) <- c("Mean", "Median")
  
  
  if (state == "t") {
    # APOBEC vs full gene list metrics (mean and median)
    APOBEC_vs_full_gene_list <- correlation_fxn(pca_APOBEC_t)
    write.table(APOBEC_vs_full_gene_list[[1]], 
                file = paste('APOBEC/tumour/',
                             cancer, '_correlation_spearman.txt',
                             sep = ''), sep = '\t')
    write.table(APOBEC_vs_full_gene_list[[2]], 
                file = paste('APOBEC/tumour/',
                             cancer, '_correlation_p_val.txt', 
                             sep = ''), sep = '\t')
    # Hypoxia genes vs full gene list metrics (mean and median)
    Hypoxia_vs_full_gene_list <- correlation_fxn(pca_Hypoxia_t)
    write.table(Hypoxia_vs_full_gene_list[[1]], 
                file = paste('Buffa 52/tumour/',
                             cancer, '_correlation_spearman.txt', sep = ''),
                sep = '\t')
    write.table(Hypoxia_vs_full_gene_list[[2]], 
                file = paste('Buffa 52/tumour/',
                             cancer, '_correlation_p_val.txt', sep = ''), 
                sep = '\t')
  } else {
    # APOBEC vs full gene list metrics (mean and median)
    APOBEC_vs_full_gene_list <- correlation_fxn(pca_APOBEC_t)
    write.table(APOBEC_vs_full_gene_list[[1]], 
                file = paste('APOBEC/normal/',
                             cancer, '_correlation_spearman.txt',
                             sep = ''), sep = '\t')
    write.table(APOBEC_vs_full_gene_list[[2]], 
                file = paste('APOBEC/normal/',
                             cancer, '_correlation_p_val.txt', 
                             sep = ''), sep = '\t')
    # Hypoxia genes vs full gene list metrics (mean and median)
    Hypoxia_vs_full_gene_list <- correlation_fxn(pca_Hypoxia_t)
    write.table(Hypoxia_vs_full_gene_list[[1]], 
                file = paste('Buffa 52/normal/',
                             cancer, '_correlation_spearman.txt', sep = ''),
                sep = '\t')
    write.table(Hypoxia_vs_full_gene_list[[2]], 
                file = paste('Buffa 52/normal/',
                             cancer, '_correlation_p_val.txt', sep = ''), 
                sep = '\t')
  }
  
  
  output <- list(APOBEC_vs_full_gene_list, Hypoxia_vs_full_gene_list)
  return(output)
}
