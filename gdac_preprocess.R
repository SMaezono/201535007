# gdac_preprocess.R
# n threads in MacOS and Linux using FORK
# 8 threads in Windows

# Jiachen Liang

# Pre-process GDAC data

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/"
#root <- "C:\\Users\\Jiachen Liang\\Documents\\GDAC\\"
setwd(root)

# Using Windows?
win <- FALSE

# Cancer list
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC",
             "LUAD", "LUSC", "OV", "PRAD", "STAD", "THCA", "UCEC")

# Parameters
gene_names <- TRUE # Set to TRUE to export using gene names; FALSE to export using gene IDs
write_all <- TRUE # Export a file with all samples
write_t <- TRUE # Export a file with only solid tumour samples
write_n <- TRUE # Export a file with only normal samples
write_stats <- TRUE # Export mean, median, sd of samples
q <- 0.5 # Threshold quantile for expression (0.5 for 50% quantile)

#################################################################

# Import packages
if(win){
  library(doParallel)
}else{
  library(parallel)
}
library(matrixStats)

# MAIN FUNCTION
main <- function(import_index){
  library(matrixStats)
  cancer <- cancers[import_index]
  print(cancer)
  print('   reading')
  if(win){
    dir <- paste(root, 'TCGA\\', cancer, ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", sep="")
  }else{
    dir <- paste(root, 'TCGA/', cancer, ".txt", sep="")
  }
  input <- read.table(dir, sep="\t", header = TRUE, row.names=1)
  input <- input[-1,]
  
  # Fix row names
  # Log transform
  name <- rownames(input)
  col <- colnames(input)
  input <- matrix(as.numeric(as.matrix(input)), nrow=nrow(input), ncol=ncol(input))
  input <- log2(input)
  input[is.infinite(input)] <- NA
  input[is.nan(input)] <- NA
  out.name <- input # output matrix using gene names
  out.id <- input # using gene IDs
  rn_name <- c()
  rn_id <- c()
  for(i in 1:length(name)){ # get names
    split <- strsplit(name[i], "\\|")
    if(split[[1]][1]=='?'){
      rn_name <- c(rn_name, paste('UNKNOWN_', as.character(i), sep=''))
    }else{
      rn_name <- c(rn_name, split[[1]][1])
    }
    rn_id <- c(rn_id, split[[1]][2])
  }
  rownames(out.name) <- rn_name
  rownames(out.id) <- rn_id
  colnames(out.name) <- col
  colnames(out.id) <- col
  
  # Write file with all samples
  if(write_all){
    print('   writing processed')
    if(gene_names){
      write.table(out.name, paste(root, cancer, '_all.txt', sep=''), sep='\t')
    }else{
      write.table(out.id, paste(root, cancer, '_all.txt', sep=''), sep='\t')
    }
  }
  
  # Separate tumour and normal
  all_col <- colnames(out.id)
  normal <- c()
  tumour <- c()
  for(i in all_col){
    char14 <- substr(i, 14, 15)
    if(char14=="01"){ # tumour
      tumour <- c(tumour, i)
    }
    else if(char14=="11"){ # normal
      normal <- c(normal, i)
    }
  }
  if(gene_names){
    out.t <- out.name[,tumour]
    out.n <- out.name[,normal]
  }else{
    out.t <- out.id[,tumour]
    out.n <- out.id[,normal]
  }
  
  # Fix col names
  colnames(out.t) <- substr(colnames(out.t), 1, 12)
  colnames(out.n) <- substr(colnames(out.n), 1, 12)
  
  if(write_t){
    print('   writing tumour')
    write.table(out.t, paste(root, cancer, '_t.txt', sep=''), sep='\t')
  }
  if(write_n){
    print('   writing normal')
    write.table(out.n, paste(root, cancer, '_n.txt', sep=''), sep='\t')
  }
  
  if(write_stats){
    
    # Stats for all
    all_mean <- mean(input, na.rm=TRUE)
    all_median <- median(input, na.rm=TRUE)
    all_sd <- sd(input, na.rm=TRUE)
    stats_all <- rbind(all_mean, all_median, all_sd)
    rownames(stats_all) <- c('Mean', 'Median', 'SD')
    colnames(stats_all) <- 'Whole matrix'
    
    # Threshold calculation
    thr_a <- quantile(input, prob=q, na.rm=TRUE)
    thr_t <- quantile(out.t, prob=q, na.rm=TRUE)
    thr_n <- quantile(out.n, prob=q, na.rm=TRUE)
    
    # Stats for each sample
    s_a_mean <- colMeans(input, na.rm=TRUE)
    s_a_median <- colMedians(input, na.rm=TRUE)
    s_t_mean <- colMeans(out.t, na.rm=TRUE)
    s_t_median <- colMedians(out.t, na.rm=TRUE)
    s_n_mean <- colMeans(out.n, na.rm=TRUE)
    s_n_median <- colMedians(out.n, na.rm=TRUE)
    stats_s_a <- cbind(s_a_mean, s_a_median)
    stats_s_t <- cbind(s_t_mean, s_t_median)
    stats_s_n <- cbind(s_n_mean, s_n_median)
    rownames(stats_s_a) <- col
    rownames(stats_s_t) <- tumour
    rownames(stats_s_n) <- normal
    colnames(stats_s_a) <- c('Mean', 'Median')
    colnames(stats_s_t) <- c('Mean', 'Median')
    colnames(stats_s_n) <- c('Mean', 'Median')
    
    # Stats for each gene
    g_a_mean <- rowMeans(input, na.rm=TRUE)
    g_a_median <- rowMedians(input, na.rm=TRUE)
    g_a_sd <- rowSds(input, na.rm=TRUE)
    g_t_mean <- rowMeans(out.t, na.rm=TRUE)
    g_t_median <- rowMedians(out.t, na.rm=TRUE)
    g_t_sd <- rowSds(out.t, na.rm=TRUE)
    g_n_mean <- rowMeans(out.n, na.rm=TRUE)
    g_n_median <- rowMedians(out.n, na.rm=TRUE)
    g_n_sd <- rowSds(out.n, na.rm=TRUE)
    stats_g_a <- cbind(g_a_mean, g_a_median, g_a_sd)
    stats_g_t <- cbind(g_t_mean, g_t_median, g_t_sd)
    stats_g_n <- cbind(g_n_mean, g_n_median, g_n_sd)
    rownames(stats_g_a) <- rownames(out.t)
    rownames(stats_g_t) <- rownames(out.t)
    rownames(stats_g_n) <- rownames(out.t)
    colnames(stats_g_a) <- c('Mean', 'Median', 'SD')
    colnames(stats_g_t) <- c('Mean', 'Median', 'SD')
    colnames(stats_g_n) <- c('Mean', 'Median', 'SD')
    
    # Differential expression
    common <- intersect(colnames(out.t), colnames(out.n))
    de_mean <- rowMeans(out.t[,common]-out.n[,common], na.rm=T)
    de_median <- g_t_mean-g_n_mean
    de <- cbind(de_mean, de_median)
    de[is.infinite(de)] <- NA
    de[is.nan(de)] <- NA
    rownames(de) <- rownames(out.t)
    colnames(de) <- c('DE using Mean', 'DE unpaired')
    
    # Filter out genes that have NA in over 75% of samples in either normal or tumour
    for(g in rownames(de)){
      prop_n <- count(is.na(out.n[g,]))/ncol(out.n)
      prop_t <- count(is.na(out.t[g,]))/ncol(out.t)
      prop_max <- max(prop_n, prop_t)
      if(prop_max > 0.75){
        de[g,] <- c(NA, NA)
      }
    }
    
    # Counts of appearances of genes with significant expression
    sel_a <- input
    sel_t <- out.t
    sel_n <- out.n
    sel_a[sel_a<thr_a] <- NA
    sel_t[sel_t<thr_t] <- NA
    sel_n[sel_n<thr_n] <- NA
    count_a <- matrix(nrow=nrow(out.t), ncol=2)
    count_t <- matrix(nrow=nrow(out.t), ncol=2)
    count_n <- matrix(nrow=nrow(out.t), ncol=2)
    rownames(count_a) <- rownames(out.t)
    rownames(count_t) <- rownames(out.t)
    rownames(count_n) <- rownames(out.t)
    colnames(count_a) <- c('Count', 'Proportion')
    colnames(count_t) <- c('Count', 'Proportion')
    colnames(count_n) <- c('Count', 'Proportion')
    for(i in 1:nrow(sel_a)){
      count <- 0
      for(j in 1:ncol(sel_a)){
        if(!is.na(sel_a[i,j])){
          count <- count + 1
        }
      }
      count_a[i,1] <- count
      count_a[i,2] <- count/ncol(sel_a)
    }
    for(i in 1:nrow(sel_t)){
      count <- 0
      for(j in 1:ncol(sel_t)){
        if(!is.na(sel_t[i,j])){
          count <- count + 1
        }
      }
      count_t[i,1] <- count
      count_t[i,2] <- count/ncol(sel_t)
    }
    for(i in 1:nrow(sel_n)){
      count <- 0
      for(j in 1:ncol(sel_n)){
        if(!is.na(sel_n[i,j])){
          count <- count + 1
        }
      }
      count_n[i,1] <- count
      count_n[i,2] <- count/ncol(sel_n)
    }
    
    # Write
    write.table(stats_all, paste(root, 'stats/', cancer, '_all.txt', sep=''), sep='\t')
    write.table(stats_s_a, paste(root, 'stats/', cancer, '_s_a.txt', sep=''), sep='\t')
    write.table(stats_s_t, paste(root, 'stats/', cancer, '_s_t.txt', sep=''), sep='\t')
    write.table(stats_s_n, paste(root, 'stats/', cancer, '_s_n.txt', sep=''), sep='\t')
    write.table(stats_g_a, paste(root, 'stats/', cancer, '_g_a.txt', sep=''), sep='\t')
    write.table(stats_g_t, paste(root, 'stats/', cancer, '_g_t.txt', sep=''), sep='\t')
    write.table(stats_g_n, paste(root, 'stats/', cancer, '_g_n.txt', sep=''), sep='\t')
    write.table(de, paste(root, 'stats/', cancer, '_de.txt', sep=''), sep='\t')
    write.table(count_a, paste(root, 'stats/', cancer, '_count_a.txt', sep=''), sep='\t')
    write.table(count_t, paste(root, 'stats/', cancer, '_count_t.txt', sep=''), sep='\t')
    write.table(count_n, paste(root, 'stats/', cancer, '_count_n.txt', sep=''), sep='\t')
  }
  
  return(cancer)
}

# Parallelization
if(win){
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:length(cancers)) %dopar% main(i)
  stopCluster(cl)
}else{
  f <- vector("list", length=length(cancers))
  for(i in 1:length(cancers)){
    f[[i]] <- mcparallel(main(i))
  }
  mccollect(f)
}

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)

