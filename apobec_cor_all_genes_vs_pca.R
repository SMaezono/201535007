# Sakura Maezono # date created: Sat Apr 20 18:12:38 2019 
# calculate the correlation between the principal components of
# APOBEC expression and the full gene list (all genes) 
# and all hypoxia in the different cancer types
# In performing PCA, all NAs were replaced with median
# extract list of genes that correlate with APOBEC PC1 and PC2


###############################################################################################
# Initiation
###############################################################################################
rm(list = ls())
t_start <- Sys.time()
t_last <- Sys.time()

###############################################################################################
# Associated library functions
###############################################################################################
# Windows or Mac
win <- TRUE

library(matrixStats)
library(gtools)
if (win) {
  library(doParallel)
}else{
  library(parallel)
}
library(gplots)
library(RDAVIDWebService)


###############################################################################################
# File locations
###############################################################################################

# create folder(s) if it doesn't/ they don't exist

dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/APOBEC/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/APOBEC/tumour/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/APOBEC/normal/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/Buffa 52/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/Buffa 52/tumour/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/Buffa 52/normal/")
# Set working directory
root <- "G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/"
setwd(root)

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'

# Directory: for output
dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/'


###############################################################################################
# Parameters
###############################################################################################

# Significance threshold
thr <- 0.05

# Correlation threshold (filtered for DAVID)
cor_thr <- 0.3

# Filter for expression?
do_filter <- T


# Minimum proportion of non-NA samples
thr_na <- 0.8

# Folds for cross-validation
k_folds <- 5

# Restrict to A3 family?
restrict <- FALSE

# Restrict to first 2 PCs?
restrict_to_PC1_PC2 <- TRUE


# Heatmap parameters (for correlation data)
key <- TRUE # display legend
skey <- TRUE # symmetric key
fontsize <- 0.8
hm_w <- 16 # width
hm_h <- 14 # height

# "t" for tumour data and "n" for normal data
state <- "t"

# Cancer list
cancers <- as.character(read.table(paste(dir_p, 
                                         'cancer_list.txt',
                                         sep = ''),
                                   header = T, sep = '\t')[,1])

if (state == "n") {
  # "OV" has no normal tissues
  cancers <- cancers[-length(cancers)]
}


### In house functions ####

# remove rows that are only NA
rm_onlyNArow <- function(df) df[rowSums(is.na(df)) != ncol(df),]
rm_onlyNAcol <- function(df) df[colSums(is.na(df)) != nrow(df),]
rm_onlyNANrow <- function(df) df[rowSums(is.nan(df)) != ncol(df),]
rm_row_without_sig <- function(df) df[rowSums(df > 0.05) != ncol(df),]


# put in decimal the desired percentage for top or bottom group
percent_in_decimal <- 1
###############################################################################
# Import expression matrix from file system
###############################################################################
# APOBEC and hypoxia genes
dir_ext <- "G:/My Drive/Jiachen Files/data/"
apo.genes <- as.character(read.table(file = paste(dir_ext, 
                                                  "tr_a.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file = paste(dir_ext,
                                                  "tr_h.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)


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
  # output <- as.matrix(output)
  
  
  # Filter out genes non-NA in less than x proportion of samples
  na_filter <-  
    (rowCounts(!is.na(output))/ncol(output)>=thr_na)
  output <- output[na_filter,]
  
  # Replace NA by median of gene
  output_mat <- as.matrix(output)
  med_t <- rowMedians(output_mat, na.rm = T)
  for (i in 1:nrow(output_mat)) {
    for (j in 1:ncol(output_mat)) {
      if (is.na(output_mat[i,j])) {
        output_mat[i,j] <- med_t[i]
      }
    }
  }
  
  # remove rows that only have NAs or NAN
  output_mat <- rm_onlyNArow(output_mat)
  output_mat <- rm_onlyNANrow(output_mat)
 
  
  hyp_t <- as.matrix(output[hyp.genes,])

  na_filter <- (rowCounts(!is.na(hyp_t))/ncol(hyp_t)>=thr_na)
  hyp_t <- hyp_t[na_filter,]

  # Replace NA by median of gene
  med_t <- rowMedians(hyp_t, na.rm = T)
  for (i in 1:nrow(hyp_t)) {
    for (j in 1:ncol(hyp_t)) {
      if (is.na(hyp_t[i,j])) {
        hyp_t[i,j] <- med_t[i]
      }
    }
  }
  # remove rows that only have NAs
  hyp_t <- rm_onlyNArow(hyp_t)
  hyp_t <- rm_onlyNANrow(hyp_t)
  
  input_t <- as.matrix(output[apo.genes,])
  if (restrict) {
    input_t <- input_t[c('APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
                         'APOBEC3F', 'APOBEC3G', 'APOBEC3H'),]
  }
  # input_t <- data.frame(output[order(rownames(input_t)),])
  # remove rows that only have NAs
  # input_t <- rm_onlyNArow(input_t)
  # Filter out genes non-NA in less than x proportion of samples
  na_filter <- (rowCounts(!is.na(input_t))/ncol(input_t)>=thr_na)
  input_t <- input_t[na_filter,]
  # 
  # Replace NA by median of gene
  med_t <- rowMedians(input_t, na.rm = T)
  for(i in 1:nrow(input_t)) {
    for(j in 1:ncol(input_t)) {
      if(is.na(input_t[i,j])) {
        input_t[i,j] <- med_t[i]
      }
    }
  }
  # remove rows that only have NAs
  input_t <- rm_onlyNArow(input_t)
  input_t <- rm_onlyNANrow(input_t)
  # matrices of mean & median scores for every gene
  output_all <- matrix(nrow = length(output[,1]), ncol = 2 )
  rownames(output_all) <- rownames(output)
  colnames(output_all) <- c("Mean", "Median")
  output_all[,1] <- rowMeans(output, na.rm=T)
  output_all[,2] <- rowMedians(as.matrix(output), na.rm=T)
  
  # Output_all mean and median of full list of genes 
  write.table(output_all, 
              file = paste(dir_o, cancer,"_", state,"_",
                         'meanormedian_of_full_list_genes.txt', 
                         sep = ""), sep = '\t')
  
  cancer_output <- list(output_all, data.frame(hyp_t), 
                        data.frame(input_t), output)
  return(cancer_output)
}

###############################################################################
# Correlate the Principal components of gene signatures with hypoxia genes or 
# full gene list 
###############################################################################
## correlation ###
correlation_plot <- function(testtable, filename) {
  # put all the correlation (non-param) significance 
  # between different PCs of gene sig 
  # vs. full gene list of the different cancer types
  if (length(rownames(testtable)) >= 100) {
    fontsize <- 0.3
  }
  testtable1 <- data.matrix(testtable)
  # testtable should be correlation table (spearman rho)
  pdf(file = paste(dir_o, filename, "_", state, '.pdf', sep = ""), 
      width = 12, height = 4)
  ht <- try(heatmap.2(testtable1,
                  Rowv = TRUE, Colv = TRUE,
                  #main=paste("1. Mean -", m)
                  col = redblue, key = key, symkey = skey,
                  cexRow = fontsize, cexCol = fontsize,
                  density.info = "none", trace = "none",
                  breaks = seq(-1, 1, length.out = 101),
                  na.color = "white"))
  if (class(ht) == "try-error") {    
    ht <- heatmap.2(testtable1,
                    Rowv = FALSE, Colv = FALSE,
                    #main=paste("1. Mean -", m)
                    col = redblue, key = key, symkey = skey,
                    cexRow = fontsize, cexCol = fontsize,
                    density.info = "none", trace = "none",
                    breaks = seq(-1, 1, length.out = 101),
                    na.color = "white")
    cat("caught an error in row and/or column dendogram") 
  }

  dev.off()
}
correlation_plot_2 <- function(testtable, filename) {
  # put all the correlation (non-param) significance 
  # between different PCs of gene sig 
  # vs. full gene list of the different cancer types
  
  if (length(rownames(testtable)) >= 100) {
    fontsize <- 0.3
  }
  testtable1 <- data.matrix(testtable)
  # testtable should be correlation table (spearman rho)
  pdf(file = paste(dir_o, filename, "_", state, '.pdf', sep = ""))
  ht <- try(heatmap.2(data.matrix(testtable1),
                      col = redblue, cexRow = fontsize, 
                      cexCol = fontsize, density.info = "none", 
                      trace = "none", 
                      breaks = seq(-1, 1, length.out = 101),
                      na.color = "white")) 
  if (class(ht) == "try-error") {    
    ht <- heatmap.2(testtable1,
                    Rowv = FALSE, Colv = FALSE,
                    #main=paste("1. Mean -", m)
                    col = redblue, key = key, symkey = skey,
                    cexRow = fontsize, cexCol = fontsize,
                    density.info = "none", trace = "none",
                    breaks = seq(-1, 1, length.out = 101),
                    na.color = "white")
    cat("caught an error in row and/or column dendogram") 
  }
  dev.off()

}
main_fxn_pca_vs_full_gene_list <- function(cancer) { 
  
  cancer_output <- import_cancer_exp(cancer, state)
  
  # PCA of Hypoxia genes
  only_Hypoxia_t <- cancer_output[[2]]
  pca_Hypoxia_t <- prcomp(only_Hypoxia_t,
                          center = TRUE, scale. = TRUE)
  
  # PCA of APOBEC genes
  only_APOBEC_t <- cancer_output[[3]]
  pca_APOBEC_t <- prcomp(only_APOBEC_t, 
                         center = TRUE, scale. = TRUE) 
  
  # Output rotation data and rotated values (Hypoxia)
  if (state == "t") {
    write.table(pca_Hypoxia_t$rotation, 
                file = paste(dir_o, 'Buffa 52/tumour/', cancer,
                             '_rotation.txt', sep = ''), sep = '\t')
    write.table(pca_Hypoxia_t$x, 
                file = paste(dir_o, 'Buffa 52/tumour/',
                             cancer, '_x.txt', sep = ''), sep = '\t')
    
    # Output rotation data and rotated values (APOBEC)
    write.table(pca_APOBEC_t$rotation, 
                file = paste(dir_o, 'APOBEC/tumour/',
                             cancer, '_rotation.txt', sep = ''), sep = '\t')
    write.table(pca_APOBEC_t$x, 
                file = paste(dir_o, 'APOBEC/tumour/', 
                             cancer, '_x.txt', sep = ''), sep = '\t')
  } else {
    write.table(pca_Hypoxia_t$rotation, 
                file = paste(dir_o, 'Buffa 52/normal/', cancer,
                             '_rotation.txt', sep = ''), sep = '\t')
    write.table(pca_Hypoxia_t$x, 
                file = paste(dir_o, 'Buffa 52/normal/',
                             cancer, '_x.txt', sep = ''), sep = '\t')
    
    # Output rotation data and rotated values (APOBEC)
    write.table(pca_APOBEC_t$rotation, 
                file = paste(dir_o, 'APOBEC/normal/',
                             cancer, '_rotation.txt', sep = ''), sep = '\t')
    write.table(pca_APOBEC_t$x, 
                file = paste(dir_o,'APOBEC/normal/', 
                             cancer, '_x.txt', sep = ''), sep = '\t')
  }
  
  # 
  # # per patient (metric scores of each genes)
  # output <- as.matrix(cancer_output[[4]])
  # out_mean <- data.frame(rowMeans(output, na.rm = T))
  # out_median <- data.frame(rowMedians(output, na.rm = T))
  # output_metric_score <- cbind(out_mean, out_median)
  # colnames(output_metric_score) <- c("Mean", "Median")
  # rownames(output_metric_score) <- rownames(output)
  
  # graph biplot
  require(factoextra)
  pca_APOBEC_transpose <- prcomp(t(only_APOBEC_t), 
                         center = TRUE, scale. = TRUE)
  pdf(file = paste(dir_o, cancer, "_APOBEC_biplot_", state, '.pdf', sep = ""))
  pc_APOBEC_x <- fviz_pca_biplot(pca_APOBEC_transpose, label = "var", col.var = "orange",
                                 habillage = "none",
                                 addEllipses = TRUE, ellipse.level = 0.95) +
    labs(title = "PCA", x = "PC1", y = "PC2")
  print(pc_APOBEC_x) 
  
  dev.off()
  
  # scree plot (which PC is important)
  pdf(file = paste(dir_o, cancer, "_APOBEC_screeplot_", state, '.pdf', sep = ""))
  pca_APOBEC_transpose_scree <- fviz_screeplot(pca_APOBEC_transpose, ncp = 9) + 
    labs(title = "PCA Scree plot", x = "Principal Components (PC)",
         y = "Percentage of explained variances (%)")
  print(pca_APOBEC_transpose_scree) 
  dev.off()
  
  require(corrplot)
  # looking at correlation between different APOBECs and PCs
  pdf(file = paste(dir_o,cancer, "_APOBEC_corrplot_pc1topc3_", state, '.pdf', sep = ""))
  corrplot(pca_APOBEC_transpose$rotation[,c(1:3)], "pie")
  dev.off()
  
  pc_APOBEC_t <- t(pca_APOBEC_t$rotation)
  if (restrict_to_PC1_PC2) {
    pc_APOBEC_t <- pc_APOBEC_t[c(1:2),]  
  }
  
  output_full_list <- cancer_output[[4]]
  
  ### PCA of gene sig vs. hypoxia or full gene list in cancer ###
  correlation_fxn <- function(gene_list) {
    # correlate PCs of APOBEC vs. hypoxia or full gene list
    # gene_list = c(only_Hypoxia_t, output_full_list)
    
    # correlation analysis (spearman) rho
    table_cor <- matrix(nrow = nrow(as.matrix(pc_APOBEC_t)), ncol = nrow(as.matrix(gene_list)))
    table_pv <- matrix(nrow = nrow(as.matrix(pc_APOBEC_t)), ncol = nrow(as.matrix(gene_list)))
    # naming
    rownames(table_cor) <- rownames(pc_APOBEC_t)
    colnames(table_cor) <- rownames(gene_list)
    rownames(table_pv) <- rownames(pc_APOBEC_t)
    colnames(table_pv) <- rownames(gene_list)
    # calculate median
    exp_med <- as.numeric(as.matrix(gene_list))
    exp_med[is.na(exp_med)] <- 0
    med <- median(exp_med)
    
    for (i in 1:nrow(as.matrix(pc_APOBEC_t))) {
      for (j in 1:nrow(as.matrix(gene_list))) { # for all genes
        if ((mean(as.matrix(gene_list)[j,], na.rm = T) < med) & do_filter) {
          table_cor[1,j] <- NA
          table_pv[1,j] <- NA
          } else { 
            test <- cor.test(as.numeric(pc_APOBEC_t[i,]), 
                             as.numeric(as.matrix(gene_list[j,])),
                             method = 'spearman')
            table_cor[i,j] <- test$estimate
            table_pv[i,j] <- test$p.value
          }
      }
    }
    
    
    # remove genes (columns) with NA in both PC1 and PC2
    table_cor <- t(rm_onlyNArow(t(table_cor)))
  
    # FDR
    table_pv[1,] <- p.adjust(table_pv[1,], method = 'fdr')
    table_pv[2,] <- p.adjust(table_pv[2,], method = 'fdr')
    table_pv <- as.matrix(t(rm_onlyNArow(t(table_pv))))
    
    # replace insignificant values with zero
    table_cor[table_pv >= thr] <- 0
    # #################################################################
    # # filter based on rho (corr np)
    # # Make binary form (correlated or not)
    # table_bin <- t(table_cor)
    # table_bin[table_bin>cor_thr] <- 1
    # table_bin[table_bin<1] <- 0
    # table_bin[is.na(table_bin)] <- 0
    # table_bin <- table_bin[order(rowSums(table_bin), decreasing=T),]
    # bin_sum <- tablerix(as.numeric(rowSums(table_bin)), nrow=nrow(table_bin), ncol=1)
    # colnames(bin_sum) <- 'Count'
    # table_bin <- cbind(table_bin, bin_sum)
    # table_bin <- as.data.frame(table_bin)
    # 
    # # Make summed form (sum of all rho values)
    # table_sum <- as.tablerix(rowSums(t(table_cor), na.rm=T))
    # rownames(table_sum) <- colnames(table_cor)
    # colnames(table_sum) <- 'Sum of Rho'
    # # filter for genes at least correlated with one gene in signature
    # table_sum <- table_sum[rownames(table_bin[table_bin$Count>0,]),] 
    # table_sum <- table_sum[order(-table_sum)]
    # ####################################################################################
    ### visualize correlation (heatmap) ###
    if (length(gene_list[,1]) <= length(only_Hypoxia_t[,1])) {
      correlation_plot(table_cor, 
                       paste("cor_PCs_of_APOBEC_vs_Hypoxia_genes", 
                             cancer, sep = "_"))
      gene_sig <- "vs_hypoxia_sig"
      
    } else {
      correlation_plot(table_cor, 
                       paste("cor_PCs_of_APOBEC_vs_all_genes",
                             cancer, sep = "_"))
      gene_sig <- "vs_all_genes"
      
    }
    
    if (state == "t") {
 
      # write tables
      write.table(table_cor, file = paste(dir_o, cancer, "_tumour_", gene_sig, "_",
                                          'cor.txt', sep = ""), sep = '\t')
      write.table(table_pv, file = paste(dir_o, cancer, "_tumour_", gene_sig, "_",
                                         'pv.txt', sep = ""), sep = '\t')
    } else {
      # write tables
      write.table(table_cor, file = paste(dir_o, cancer, "_normal_", gene_sig, "_",
                                          'cor.txt', sep = "_"), sep = '\t')
      write.table(table_pv, file = paste(dir_o, cancer, "_normal_", gene_sig, "_",
                                         'pv.txt', sep = "_"), sep = '\t')
    }
    output <- list(data.frame(table_cor), data.frame(table_pv))
    return(output)
  }
  output_all_genes <- correlation_fxn(output_full_list)
  output_hypoxia <- correlation_fxn(only_Hypoxia_t)
  output <- list(output_all_genes, output_hypoxia)
  return(output)
  
}

### to visualize the summary of the correlation in the different cancers ###

## based on p-val of correlation ###
complex_heatmap_sum <- function(testtable, heatmapfilename, columntitle) {
  # put all the correlation (non-param) significance 
  # between different PCs of gene sig 
  # vs. full gene list of the different cancer types
  # heatmap with a histogram annotation (non-param)
  
  # testtable should be pval table and  
  # heatmapfilename == name of the file 
  
  # columntitle should be the name of signature with number of genes
  #i.e. "Buffa Hypoxia genes (52)"
  
  ### heatmap annotation###
  test_table1 <- data.frame(testtable)
  # save rownames to be used later to maintain exact chars
  saverownames <-  as.character(rownames(test_table1)) 
  
  for_histogram <- lapply(c(1:ncol(test_table1)), function(x) {
    # get the length (# of cancer in which a PC of gene sig is significantly different 
    # from full gene list metrics (mean/median))
    
    len <- length(test_table1[,x][test_table1[,x] < thr]) 
    df <-  data.frame(len)
    colnames(df) <- colnames(test_table1)[x]
    
    #set rownames to have uniform rownames needed to bind them later
    rownames(df) <-  "frequency" 
    return(df)
  }) 
  
  # to combine the data with number of cancers in which hypoxia genes are 
  # significantly different when comparing tumor vs cancer for barplot (bp)
  require(dplyr)
  anno_bp <-  data.frame(bind_cols(for_histogram)) 
  anno_bp_reverse <-  t(anno_bp) # to transpose so it can be ordered based on frequency
  
  # to rearrange the gene datasets based on frequency
  
  # cbind does not match by rownames unless combined with match... 
  # merge(x,y, by= "row.names") is another method but it makes the rownames 
  # a column, all.x means keep all rownames
  test_table1_reverse <- merge(anno_bp_reverse, t(test_table1),
                               by = "row.names", all.x = TRUE, all.y = TRUE) 
  
  # put the rownames back 
  rownames(test_table1_reverse) <-  test_table1_reverse$Row.names 
  
  # remove rownames as a column as it has been put back as rownames
  test_table1_reverse <- test_table1_reverse[,-1] 
  
  # order by frequency of the cancers
  test_table1_reverse <- data.frame(
    test_table1_reverse[order(test_table1_reverse[,1]),]) 
  
  # make this the new test_table1 that has been transposed back and 
  # arrange the genes from low to high freq
  test_table1 <- data.frame(t(test_table1_reverse))
  
  # barplots only allow numeric vectors/matrices;
  # numbers match with columns of the heatmap
  anno_bpval <- as.numeric(test_table1[1,]) 
  
  # remove freq row
  test_table1 <- test_table1[-1,]
  rownames(test_table1) <- saverownames
  #test_table1$no_of_sig_genes <- length(test_table1[x,][test_table1[x,]<0.05]) 
  
  for_main <- lapply(c(1:length(test_table1[,1])), function(x) {
    len <- length(test_table1[x,][test_table1[x,] < thr]) 
    df <- data.frame(len)
    rownames(df) <- rownames(test_table1)[x]
    colnames(df) <- "frequency" 
    return(df)
  }) 
  
  maindf <-  data.frame(bind_rows(for_main)) 
  test_table1 <-  cbind(test_table1,maindf)
  test_table1 <-  test_table1[order(test_table1$frequency),]
  #test_table2 <<-  test_table1[-length(test_table1[1,])] # for it to go to the global environment
  # libraries needed for the next few ;ines of code
  require(ComplexHeatmap)
  require(circlize)
  #to create the histogram annotation (no of cancer in which a particular gene is significant)
  ha = HeatmapAnnotation(show_annotation_name = TRUE,  
                         annotation_name_offset = unit(10, "mm"), 
                         annotation_name_side = "left",
                         annotation_name_rot = 90, 
                         Number_of_Cancer_types = anno_barplot(anno_bpval,
                                                               which = "column", 
                                                               border = TRUE, 
                                                               bar_width = 0.6, 
                                                               gp = gpar(fill = "#CCCCCC"),
                                                               ylim = c(0, 23), 
                                                               axis = TRUE, 
                                                               axis_side = "left", 
                                                               axis_gp = gpar(fontsize = 9), 
                                                               axis_direction = "normal"))
  
  ### Main heatmap ###
  
  # cluster_columns set to false so histograms/barplots would match with the ht
  
  # the Euclidean distance is used by default to measure the dissimilarity 
  # between each pair of observations.  clustering_distance_rows = "euclidean",
  
  # top_annotation_height could be adjusted in the main heatmap
  
  
  ht = Heatmap(test_table1[-length(test_table1[1,])], top_annotation = ha, 
               top_annotation_height = unit(5, "cm"), 
               heatmap_legend_param = list(title = "p_value", legend_direction = "vertical",
                                           title_gp = gpar(fontsize = 14, fontface = "bold"),
                                           legend_height = unit(24,"cm"), 
                                           labels_gp = gpar(fontsize = 14),
                                           color_bar = "continuous", 
                                           labels = c("5e-02", "1e-02", "1e-03", "1e-04"), 
                                           at = c(0.05,0.01, 0.001,0.0001)), 
               na_col = "red", col = colorRamp2(c(.05, 0.01, 0.001, .0001), 
                                                c("white","lightblue", "blue", "blue4")), 
               cluster_rows = FALSE, cluster_columns = FALSE, column_title = columntitle, 
               column_title_side = "bottom", column_title_gp = gpar(fontsize = 20),
               column_title_rot = 0) 
  
  #To draw the heatmap and save as png
  png(paste(dir_o, heatmapfilename, "_", state, "_", "complexheatmap.png", sep = ""), width = 1000, 
      height = 1000, units = "px", pointsize = 12) 
  
  #draw is to output the heatmaps from ht_list
  #padding = unit(c(bottom,left,top,right), "cm") for the sides of the image 
  draw(ht, padding = unit(c(1,2,1,1), "cm")) 
  dev.off()
}



### MAIN ######################################################################
cancer_output <- lapply(cancers, main_fxn_pca_vs_full_gene_list)


per_pc_corr <- function(which_genelist, which_values, which_PC) {
  # create correlation graphs including all cancers
  # per principal component of APOBEC
  # for which_genelist, PCs vs. all genes would be 1, and 
  # PCs vs. hypoxia genes would be 2
  # for which_values, expression values would be 1, and 
  # significant difference values would be 2
  # for which_PC, PC1 would be 1, and 
  # PC2 would be 2, etc...
  
  output_PC <- lapply(c(1:length(cancers)), 
                      function(x) cancer_output[[x]][[which_genelist]][[which_values]][which_PC,])
  
  # smartbind from the package gtools does rbind for df with unequal col numbers
  if (state == "t") {
    output_PC <- smartbind(output_PC[[1]], output_PC[[2]], output_PC[[3]],
                           output_PC[[4]], output_PC[[5]], output_PC[[6]],
                           output_PC[[7]], output_PC[[8]], output_PC[[9]],
                           output_PC[[10]], output_PC[[11]], output_PC[[12]],
                           output_PC[[13]], output_PC[[14]], output_PC[[15]])
    
  } else {
    output_PC <- smartbind(output_PC[[1]], output_PC[[2]], output_PC[[3]],
                           output_PC[[4]], output_PC[[5]], output_PC[[6]],
                           output_PC[[7]], output_PC[[8]], output_PC[[9]],
                           output_PC[[10]], output_PC[[11]], output_PC[[12]],
                           output_PC[[13]], output_PC[[14]])
    
  }
 
  
  
  rownames(output_PC) <- cancers
  # Filter out genes non-NA in less than x proportion of samples
  # na_filter <-  
  #   (rowCounts(!is.na(output_PC))/ncol(output_PC)>=thr_na)
  # output_PC <- output_PC[na_filter,]
  # na_filter <-  
  #   (colCounts(!is.na(output_PC))/nrow(output_PC)>=thr_na)
  # output_PC <- output_PC[,na_filter]
  output_PC <- as.matrix(output_PC)
  
  sig_output_PC <- lapply(c(1:length(cancers)), 
                      function(x) cancer_output[[x]][[which_genelist]][[2]][which_PC,])
  if (state == "t") {
    sig_output_PC_0.05 <- smartbind(sig_output_PC[[1]], sig_output_PC[[2]], sig_output_PC[[3]],
                                    sig_output_PC[[4]], sig_output_PC[[5]], sig_output_PC[[6]],
                                    sig_output_PC[[7]], sig_output_PC[[8]], sig_output_PC[[9]],
                                    sig_output_PC[[10]], sig_output_PC[[11]], sig_output_PC[[12]],
                                    sig_output_PC[[13]], sig_output_PC[[14]], sig_output_PC[[15]])
    
  } else {
    sig_output_PC_0.05 <- smartbind(sig_output_PC[[1]], sig_output_PC[[2]], sig_output_PC[[3]],
                                    sig_output_PC[[4]], sig_output_PC[[5]], sig_output_PC[[6]],
                                    sig_output_PC[[7]], sig_output_PC[[8]], sig_output_PC[[9]],
                                    sig_output_PC[[10]], sig_output_PC[[11]], sig_output_PC[[12]],
                                    sig_output_PC[[13]], sig_output_PC[[14]])
    
  }
  rownames(sig_output_PC_0.05) <- cancers
  
  sig_output_PC_0.05 <- t(sig_output_PC_0.05)
  sig_output_PC_0.05 <- rm_onlyNArow(sig_output_PC_0.05)
  sig_output_PC_0.05 <- rm_onlyNANrow(sig_output_PC_0.05)
  sig_output_PC_0.05 <- rm_row_without_sig(sig_output_PC_0.05)
  sig_output_PC_0.05 <- data.frame(sig_output_PC_0.05)
  sig_output_PC_0.05 <- rm_onlyNArow(sig_output_PC_0.05)

  overlapping_genes <- intersect(colnames(output_PC), rownames(sig_output_PC_0.05))
  # genes that were significantly correlated in at least 1 cancer
  # for which_genelist, PCs vs. all genes would be 1, and 
  # PCs vs. hypoxia genes would be 2
  # for which_values, expression values would be 1, and 
  # significant difference values would be 2
  # for which_PC, PC1 would be 1, and 
  # PC2 would be 2, etc...
  which_genelist_filename <- ifelse(which_genelist == 1, "PCs_vs_all_genes",
                                    "PCs_vs_hypoxia genes")
  
  which_PC_filename <- ifelse(which_PC == 1, "PC1", "PC2")
  write.table(overlapping_genes, 
              file = paste(dir_o, state, "_", which_genelist_filename, "_",
                           which_PC_filename, "_",
                           "for_gene_enrichment_pathway.txt", 
                           sep = ""), sep = '\t')
  if (which_genelist == 1) {
    output_PC <- data.frame(output_PC[,overlapping_genes])
  } else {
    output_PC <- data.frame(output_PC)
  }
 
  
  
  return(output_PC)
}

# PC1 vs. all genes, expression values
PC1_exp_allgenes <- per_pc_corr(1,1,1)
# PC2 vs. all genes, expression values
PC2_exp_allgenes <- per_pc_corr(1,1,2)
# PC1 vs. hypoxia genes, expression values
PC1_exp_hypgenes <- per_pc_corr(2,1,1)
# PC2 vs. hypoxia genes, expression values
PC2_exp_hypgenes <- per_pc_corr(2,1,2)

### visualize correlation (heatmap) ###
correlation_plot_2(PC1_exp_allgenes, "APOBEC_PC1_vs_all_genes_rho") 
correlation_plot_2(PC2_exp_allgenes, "APOBEC_PC2_vs_all_genes_rho") 
correlation_plot_2(PC1_exp_hypgenes, "APOBEC_PC1_vs_hyp_genes_rho") 
correlation_plot_2(PC2_exp_hypgenes, "APOBEC_PC2_vs_hyp_genes_rho") 


# PC1 vs. all genes, sig diff values
PC1_signif_allgenes <- per_pc_corr(1,2,1)
# PC2 vs. all genes, sig diff values
PC2_signif_allgenes <- per_pc_corr(1,2,2)

# PC1 vs. hypoxia genes, sig diff values
PC1_signif_hypgenes <- per_pc_corr(2,2,1)
# PC2 vs. hypoxia genes, sig diff values
PC2_signif_hypgenes <- per_pc_corr(2,2,2)


### visualize summarized correlation significance (heatmap) ###

complex_heatmap_sum(PC1_signif_allgenes,"APOBEC_PC1_signif_allgenes", 
                    "gene expressions")

complex_heatmap_sum(PC2_signif_allgenes,"APOBEC_PC2_signif_allgenes",
                    "gene expressions")

complex_heatmap_sum(PC1_signif_hypgenes,"APOBEC_PC1_signif_hypoxiagenes", 
                    "Hypoxia gene expressions")

complex_heatmap_sum(PC2_signif_hypgenes,"APOBEC_PC2_signif_hypoxiagenes",
                    "Hypoxia gene expressions")


root <- "G:/My Drive/Jiachen Files/codes_supp/"
setwd(root)
gc()

### extract list of genes that correlate with APOBEC pc1 and pc2 (tumor)

extract_list_of_genes_corr_to_PCs <- function(cancer) {
  # extract list of genes that correlate with APOBEC pc1 and pc2 (tumor)
  extract_genes <- function(gene_sig) {
    # extract genes that are significantly different based on p-value alone 
    # and those that are both filtered by pval and correlation
    
    # extract genes based on p-value
    list_of_genes <- read.table(file = paste(dir_o, cancer, "_", "tumour", "_", gene_sig,
                                             "_", 'pv.txt', sep = ""), sep = '\t')
    list_of_genes <- data.frame(t(list_of_genes))
    
    # get list of genes significantly correlated with PC1 based on pval
    PC1_corr_genes <- list_of_genes[list_of_genes$PC1 < thr,]
    PC1_corr_genes <- PC1_corr_genes[order(PC1_corr_genes$PC1, decreasing = FALSE),]
    PC1_corr_genes <- rownames(PC1_corr_genes)
    write.table(PC1_corr_genes, file = paste(dir_o, cancer, "_", 
                                             "tumour", "_", gene_sig,
                                             "_",  "cor_to_PC1.txt",
                                             sep = ""), sep = '\t')
    
    
    # extract genes based on correlation rho 
    list_of_genes_cor <- read.table(file = paste(dir_o, cancer, "_", "tumour", "_", gene_sig,
                                            "_", 'cor.txt', sep = ""), sep = '\t')
    list_of_genes_cor <- data.frame(t(list_of_genes_cor))
    
    # get list of genes significantly correlated with PC1 based on rho
    PC1_corr_genes_rho <- list_of_genes_cor[abs(list_of_genes_cor$PC1) >= cor_thr,]
    PC1_corr_genes_rho <- PC1_corr_genes_rho[order(PC1_corr_genes_rho$PC1, 
                                                   decreasing = FALSE),]
    PC1_corr_genes_rho <- rownames(PC1_corr_genes_rho)
    if (length(PC1_corr_genes_rho) >= 1) {
      write.table(PC1_corr_genes_rho, file = paste(dir_o, cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_",  "cor_to_PC1_", cor_thr, ".txt",
                                                   sep = ""), sep = '\t')
    } else {
      write.table(PC1_corr_genes_rho, file = paste(dir_o, "EMPTY_", cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_",  "cor_to_PC1_", cor_thr, ".txt",
                                                   sep = ""), sep = '\t')
    }

    if (length(PC1_corr_genes_rho) >= 1 & length(PC1_corr_genes) >= 1) { 
    # get list of genes significantly correlated with PC1 filtered
    # by both pval AND rho
      PC1_corr_genes_filtered <- intersect(PC1_corr_genes, PC1_corr_genes_rho)
      write.table(PC1_corr_genes_rho, file = paste(dir_o, cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_",  "double_filter_cor_to_PC1_",
                                                   cor_thr, ".txt",
                                                   sep = ""), sep = '\t')
    } else {
      print(paste("Double filter could not be done in", cancer, gene_sig, 
                  "for APOBEC PC1"))
    }
    
    if (percent_in_decimal < 1) {
      desired_percent <- ceiling(length(PC1_corr_genes)*percent_in_decimal) 
    } else {
      desired_percent <- length(PC1_corr_genes)
    }
    
    # prints out how many samples were considered
    print(paste(desired_percent, 
                "gene(s) make up the top", percent_in_decimal*100, "%",
                "of the group correlated with APOBEC PC1 in",
                cancer, "based on p-value alone"))
    
    # and take top desired percent
    PC1_corr_genes <- PC1_corr_genes[1:desired_percent] 
    write.table(PC1_corr_genes, file = paste(dir_o, cancer, "_", 
                                             "tumour", "_", gene_sig,
                                             "_", 
                                             "cor_to_PC1_top_",
                                             percent_in_decimal*100,
                                             "%.txt",
                                             sep = ""), sep = '\t')
    
    
    
    # get list of genes correlated with PC2
    PC2_corr_genes <- list_of_genes[list_of_genes$PC2 < thr,]
    PC2_corr_genes <- PC2_corr_genes[order(PC2_corr_genes$PC2, decreasing = FALSE),]
    PC2_corr_genes <- rownames(PC2_corr_genes)
    write.table(PC2_corr_genes, file = paste(dir_o, cancer, "_",
                                             "tumour", "_", gene_sig,
                                             "_", "cor_to_PC2.txt", sep = ""), 
                sep = '\t')
    
    # extract genes based on correlation rho 
    list_of_genes <- read.table(file = paste(dir_o, cancer, "_", "tumour", "_", gene_sig,
                                             "_", 'cor.txt', sep = ""), sep = '\t')
    list_of_genes_cor <- data.frame(t(list_of_genes))
    
    # get list of genes significantly correlated with PC2 based on rho
    PC2_corr_genes_rho <- list_of_genes_cor[abs(list_of_genes_cor$PC2) >= cor_thr,]
    PC2_corr_genes_rho <- PC2_corr_genes_rho[order(PC2_corr_genes_rho$PC2,
                                                   decreasing = FALSE),]
    PC2_corr_genes_rho <- rownames(PC2_corr_genes_rho)
    if (length(PC2_corr_genes_rho) >= 1) {
      write.table(PC2_corr_genes_rho, file = paste(dir_o, cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_",  "cor_to_PC2_", cor_thr, ".txt",
                                                   sep = ""), sep = '\t')
    } else {
      write.table(PC2_corr_genes_rho, file = paste(dir_o, "EMPTY_", cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_",  "cor_to_PC2_", cor_thr, ".txt",
                                                   sep = ""), sep = '\t')
    }
    
    if (length(PC2_corr_genes_rho) >= 1 & length(PC2_corr_genes) >= 1) { 
      # get list of genes significantly correlated with PC2 filtered
      # by both pval AND rho
      PC2_corr_genes_filtered <- intersect(PC2_corr_genes, PC2_corr_genes_rho)
      write.table(PC2_corr_genes_rho, file = paste(dir_o, cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_",  "double_filter_cor_to_PC2_",
                                                   cor_thr, ".txt",
                                                   sep = ""), sep = '\t')
    } else {
      print(paste("Double filter could not be done in", cancer, gene_sig, 
                  "for APOBEC PC2"))
    }
    if (percent_in_decimal < 1) {
      desired_percent <- ceiling(length(PC2_corr_genes)*percent_in_decimal) 
    } else {
      desired_percent <- length(PC2_corr_genes)
    }
    
    # prints out how many samples were considered
    print(paste(desired_percent, 
                "gene(s) make up the top", percent_in_decimal*100, "%",
                "of the group correlated with APOBEC PC2",
                cancer))
    
    # and take top desired percentage
    PC2_corr_genes <- PC2_corr_genes[1:desired_percent] 
    write.table(PC2_corr_genes, file = paste(dir_o, cancer, "_", 
                                             "tumour", "_", gene_sig,
                                             "_", 
                                             "cor_to_PC2_top",
                                             percent_in_decimal*100, "%.txt",
                                             sep = ""), sep = '\t')
    
    
    
    # get list of genes correlated with BOTH PC1 and PC2
    PC1andPC2_corr_genes <- intersect(PC1_corr_genes,PC2_corr_genes)
    write.table(PC1andPC2_corr_genes, file = paste(dir_o, cancer, "_" ,"tumour", gene_sig,
                                                   "_", "cor_to_PC1andPC2.txt", sep = ""), 
                sep = '\t')
    print(paste("Saved list of genes correlated with PCs in", cancer, gene_sig))
    if (percent_in_decimal < 1) {
      desired_percent <- ceiling(length(PC1andPC2_corr_genes)*percent_in_decimal) 
    } else {
      desired_percent <- length(PC1andPC2_corr_genes)
    }
    
    # prints out how many samples were considered
    print(paste(desired_percent, 
                "gene(s) make up the top", percent_in_decimal*100, "%",
                "of the group correlated with APOBEC PC1andPC2",
                cancer))
    
    # and take top 10% 
    PC1andPC2_corr_genes <- PC1andPC2_corr_genes[1:desired_percent] 
    write.table(PC1andPC2_corr_genes, file = paste(dir_o, cancer, "_", 
                                                   "tumour", "_", gene_sig,
                                                   "_", 
                                                   "cor_to_PC1andPC2_top_",
                                                   percent_in_decimal*100, "%.txt",
                                                   sep = ""), sep = '\t')
  }
  lapply(c("vs_hypoxia_sig", "vs_all_genes"), 
         extract_genes)
}

lapply(cancers, extract_list_of_genes_corr_to_PCs)

# reannotate_hngc_to_entrez_fxn <- function(genelist, up_or_down) {
#   # reannotate the ensembl ids from WeiChen 
#   # (for BRCA, PRAD, PAAD, COAD) using biomaRt
#   
#   ensembl_entrez <- getBM_hsapiens_ensembl("ensembl_gene_id", genelist)
#   biomart_SF_eg <- na.omit(ensembl_entrez) 
#   
#   write.csv(biomart_SF_eg, paste("Annotations",up_or_down,
#                                  "SF_WeiChen.csv", sep = "_"))
#   
#   # match ensembl ids
#   matching_affyid <- intersect(genelist,
#                                unique(biomart_SF_eg$ensembl_gene_id)) 
#   print(paste(length(matching_affyid), 
#               " ensembl id from WeiChen have other latest biomaRt annotations"))
#   
#   # detemine the ensembl ids missing
#   missing_affy_id <- genelist[!genelist %in% biomart_SF_eg$ensembl_gene_id] 
#   
#   if (length(missing_affy_id) > 0) {
#     # reannotate to hgnc symbol
#     print(paste(missing_affy_id,
#                 "is/are not present in the latest biomaRt"))
#   } else {
#     print(paste("All ensembl ids from WeiChen have the latest biomaRt")) 
#   }
#   
#   # only keep the entrez and ensembl id columns
#   biomart_SF_eg <- biomart_SF_eg[c(2:4)] 
#   
#   # only keep the unique rows 
#   biomart_SF_eg <- unique(biomart_SF_eg) 
#   
#   # order based on entrez symbol
#   biomart_SF_eg <- biomart_SF_eg[order(biomart_SF_eg$entrezgene),] 
#   rownames(biomart_SF_eg) <- biomart_SF_eg$entrezgene
#   
#   # SF entrez gene version
#   SF_entrez_mRNA <- as.character(unique(biomart_SF_eg$entrezgene)) 
#   print(paste(length(SF_entrez_mRNA), "entrez ids"))
#   
#   # unique ensembl ids version
#   SF_en_mRNA <- as.character(unique(biomart_SF_eg$ensembl_gene_id.1)) 
#   print(paste(length(SF_en_mRNA), "ensembl ids")) 
#   
#   # unique gene symbol version
#   SF_symbol_mRNA <- as.character(unique(biomart_SF_eg$hgnc_symbol)) 
#   print(paste(length(SF_symbol_mRNA), "hgnc symbols"))
#   
#   return(biomart_SF_eg)
# }


# save details of this R script
scriptname <- "apobec_cor_all_genes_vs_pca"
savehistory(file = paste(scriptname, ".Rhistory", sep = ""))
sink(file = paste(Sys.Date(), scriptname, ".txt", sep = ""),
     type = c("output", "message"))
# details of the packages and other useful info.
print(sessionInfo()) 
sink(NULL)

# save globalenv/workspace
save.image("apobec_cor_all_genes_vs_pca.Rdata")  
