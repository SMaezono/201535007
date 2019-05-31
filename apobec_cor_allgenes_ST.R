# apobec_cor_allgenes_ST.R

# Jiachen Liang

# Generate correlation scores for the expression of each APOBEC gene with every gene in the TCGA expression matrix for all cancers

# !!!!! SINGLE THREADED VERSION !!!!!
print('!!!!! SINGLE THREADED VERSION !!!!!')
print('EXPECTED RUNTIME: 12 HOURS')
print('RUN apobec_cor_allgenes.R FOR PARALLEL VERSION')

# Output:
#   Heatmaps per cancer of APOBECs vs all genes
#     Heatmap and list of all genes
#     Heatmap and list of top 100 genes
#     Heatmap and list of top 50 genes
#     Heatmap and list of top 50 genes for A3B only
#   Heatmaps across cancers of each APOBEC vs all genes
#     Heatmap and list of all genes
#     Heatmap and list of top 100 genes
#     Heatmap and list of top 50 genes
#     Heatmap and list of top 50 genes for A3B only
#   List of genes correlated with all A3 family genes (all correlation scores and p-values)
#   List of genes correlated with more than 5 A3 family genes (all correlation scores and p-values)
#   List of genes correlated with A3B (all correlation scores and p-values)

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
library(ppcor)

###############################################################################################
# Parameters
###############################################################################################

# Heatmap parameters
margins <- c(8,8) # margins for output image
key <- FALSE # display legend
skey <- TRUE # symmetric key
keysize <- 0.1
fontsize <- 1.2
hm_w <- 12 # width
hm_h <- 12 # height

# Significance threshold
thr <- 0.05

# Correlation threshold
cor_thr <- 0.3

# Filter expressed genes?
filter <- F # set to F to analyze all genes regardless of expression level
filter_quantile <- 0.5 # by how much? 0.5 is median

# plot?
doplot <- F

# Partial correlation for hypoxia signature?
p_cor <- F

###############################################################################################
# File locations
###############################################################################################

if(win){
  # Set working directory
  root <- "C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/cor_all/"
  setwd(root)
  
  # Directory: parameters
  dir_p <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'
}else{
  # Set working directory
  root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/cor_all/"
  setwd(root)
  
  # Directory: parameters
  dir_p <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'
}

# List of cancers
cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep=''), header=T, sep='\t')[,1])
cancers <- 'BRCA'

# List of APOBECs
apobecs <- as.character(read.table(paste(dir_p, 'sign_master.txt', sep=''), header=T, sep='\t')[,1])
apobecs <- apobecs[order(apobecs)]

# Hypoxia signature (for partial correlation)
hyp_sign <- as.character(read.table(paste(dir_p, 'sign_1.txt', sep=''), header=T, sep='\t')[,1])

###############################################################################################
# Log system
###############################################################################################

# log <- matrix(nrow=1, ncol=2)
# log[1,1] <- '00:00:00'
# log[1,2] <- 'Progress timing log'
# 
# atl <- function(input){
#   entry <- matrix(nrow=1, ncol=2)
#   
#   time <- as.numeric(Sys.time())-as.numeric(t_start)
#   time_h <- as.character(floor(time/3600))
#   time_m <- as.character(floor((time%%3600)/60))
#   time_s <- as.character(floor((time%%3600)%%60))
#   if(nchar(time_h)==1){time_h<-paste('0',time_h,sep='')}
#   if(nchar(time_m)==1){time_m<-paste('0',time_m,sep='')}
#   if(nchar(time_m)==1){time_s<-paste('0',time_s,sep='')}
#   
#   time_entry <- paste(
#     time_h, ':',
#     time_m, ':',
#     time_s,
#     sep=''
#   )
#   
#   entry[1,1] <- time_entry
#   entry[1,2] <- input
#   log <- rbind(log, entry)
# }

to_time <- function(){
  time <- as.numeric(Sys.time())-as.numeric(t_start)
  time_h <- as.character(floor(time/3600))
  time_m <- as.character(floor((time%%3600)/60))
  time_s <- as.character(floor((time%%3600)%%60))
  if(nchar(time_h)==1){time_h<-paste('0',time_h,sep='')}
  if(nchar(time_m)==1){time_m<-paste('0',time_m,sep='')}
  if(nchar(time_s)==1){time_s<-paste('0',time_s,sep='')}

  time_entry <- paste(
    time_h, ':',
    time_m, ':',
    time_s,
    sep=''
  )
  
  return(time_entry)
}

###############################################################################################
# Import expression matrix from file system
###############################################################################################

exp_import <- function(cancer, state){
  
  message(paste(to_time(), '       ', cancer, '- Importing'))
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  
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
  
  message(paste(to_time(), '       ', cancer, '- Importing done'))
  
  return(output)
  
}

###############################################################################################
# Parallel function calls
###############################################################################################

if(win){ ### WINDOWS ##########################################################################
  
  # Import expression matrices (parallel)
  print(Sys.time()-t_last)
  t_last <- Sys.time()
  print('Importing expression matrices')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores, outfile='log_1.txt')
  registerDoParallel(cl)
  exp_t_all <- foreach(n=1:length(cancers)) %dopar% exp_import(cancers[n], 't')
  stopCluster(cl)
  
}else{ ##### Mac OS ###########################################################################
  
  # Import expression matrices (parallel)
  print('Importing expression matrices')
  f <- vector("list", length=length(cancers))
  for(n in 1:length(cancers)){
    f[[n]] <- mcparallel(exp_import(cancers[n], 't'))
  }
  exp_t_all <- mccollect(f)
  
}

###############################################################################################
# Heatmaps per cancer
###############################################################################################

for(cancer_index in 1:length(cancers)){ # Main heatmap function
  
  library(matrixStats)
  library(gplots)
  library(ggplot2)
  library(ppcor)
  
  print(paste(to_time(), '       ', 'Starting', cancers[cancer_index]))
  
  # filter out genes with no data
  # to.del <- c()
  # for(i in 1:nrow(exp)){
  #   if(sum(!is.na(exp[i,]))<5){
  #     to.del <- c(to.del, i)
  #   }
  # }
  # exp <- exp[-(to.del),]
  
  # get expression matrix for specific cancer
  exp <- exp_t_all[[cancer_index]]
  
  # get median of expression matrix
  exp_na <- as.numeric(as.matrix(exp))
  exp_na[is.na(exp_na)] <- 0
  med <- as.numeric(quantile(exp_na, filter_quantile, na.rm=T))
  
  # do we filter?
  if(!filter){
    med <- 0
  }
  
  # hypoxia signature matrix (for partial correlation)
  hyp_mat <- t(exp[hyp_sign,])
  
  # output matrix
  mat_cor <- matrix(nrow=length(apobecs), ncol=nrow(exp))
  mat_pv <- matrix(nrow=length(apobecs), ncol=nrow(exp))
  rownames(mat_cor) <- apobecs
  rownames(mat_pv) <- apobecs
  colnames(mat_cor) <- rownames(exp)
  colnames(mat_pv) <- rownames(exp)
  
  # calculate correlation (parallel)
  
  if(file.exists(paste('per cancer/', cancers[cancer_index], '_cor.txt', sep=''))){ # Read from file system if available
    print(paste(to_time(), '       ', cancers[cancer_index], '- Reading from file system'))
    mat_cor <- as.matrix(read.table(paste('per cancer/', cancers[cancer_index], '_cor.txt', sep=''), sep='\t', header=T, row.names=1))
    mat_pv <- as.matrix(read.table(paste('per cancer/', cancers[cancer_index], '_pv.txt', sep=''), sep='\t', header=T, row.names=1))
    print(paste(to_time(), '       ', cancers[cancer_index], '- Reading done'))
  }else{ # Compute if not
    
    print(paste(to_time(), '       ', cancers[cancer_index], '- Calculating correlation'))
    
    for(i in 1:length(apobecs)){
      
      print(paste(to_time(), '       ', cancers[cancer_index], '- Calculating correlation -', apobecs[i]))
      
      # list of cor and pv for current APOBEC vs all genes
      cor <- matrix(nrow=1, ncol=nrow(exp))
      pv <- matrix(nrow=1, ncol=nrow(exp))
      
      for(j in 1:nrow(exp)){ # for all genes
        result = tryCatch({
          if(mean(as.numeric(exp[j,]), na.rm=T)>med){ # filter for expression
            if(p_cor){ # partial correlation for hypoxia signature
              # only use rows with no missing values
              temp_mat <- cbind(as.numeric(exp[apobecs[i],]), as.numeric(exp[j,]), hyp_mat)
              not_missing <- !is.na(rowMeans(temp_mat))
              # only test if there are > 3 non-null rows
              if(count(not_missing)>=3){
                x <- as.numeric(exp[apobecs[i],])[not_missing]
                y <- as.numeric(exp[j,])[not_missing]
                z <- hyp_mat[not_missing,]
                test <- pcor.test(x, y, z, method="spearman")
                cor[1,j] <- test$estimate
                pv[1,j] <- test$p.value
              }else{
                cor[1,j] <- NA
                pv[1,j] <- NA
              }
            }else{
              test <- cor.test(as.numeric(exp[apobecs[i],]), as.numeric(exp[j,]))
              cor[1,j] <- test$estimate
              pv[1,j] <- test$p.value
            }
          }else{
            cor[1,j] <- NA
            pv[1,j] <- NA
          }
        }, warning = function(w) {
          
        }, error = function(e) {
          cor[1,j] <- NA
          pv[1,j] <- NA
        }, finally = {
          
        })
      }
      
      cor <- as.numeric(cor[1,])
      pv <- as.numeric(pv[1,])
      
      # FDR
      pv <- p.adjust(pv, method='fdr')
      
      mat_cor[i,] <- cor
      mat_pv[i,] <- pv
      
    }
    
    # cores <- min(detectCores(), length(apobecs))
    # cl2 <- makeCluster(cores, outfile='log_3.txt')
    # registerDoParallel(cl2)
    # cor_apo <- foreach(i=1:length(apobecs)) %dopar% cor_cal(i)
    # stopCluster(cl2)
    # 
    # for(i in 1:nrow(mat_cor)){
    #   mat_cor[i,] <- cor_apo[[i]][[1]]
    #   mat_pv[i,] <- cor_apo[[i]][[2]]
    # }
    # 
    # output
    write.table(mat_cor, file=paste('per cancer/', cancers[cancer_index], '_cor.txt', sep=''), sep='\t')
    write.table(mat_pv, file=paste('per cancer/', cancers[cancer_index], '_pv.txt', sep=''), sep='\t')
    
    print(paste(to_time(), '       ', cancers[cancer_index], '- Calculation done'))
  }
  
  # filter out insignificant values
  mat_cor[mat_pv>thr] <- NA
  
  # Version with NAs set to 0
  mat_cor_NA <- mat_cor
  mat_cor_NA[is.na(mat_cor_NA)] <- 0
  mat_cor_100_NA <- mat_cor_NA[,order(colMeans(mat_cor_NA), decreasing = TRUE)]
  mat_cor_100_NA <- mat_cor_100_NA[,1:100]
  
  # truncate correlation matrix
  mat_cor_ranked <- mat_cor[,order(colMeans(mat_cor), decreasing = TRUE)]
  mat_cor_ranked_a3b <- mat_cor[,order(mat_cor['APOBEC3B',], decreasing = TRUE)]
  mat_cor_100 <- mat_cor_ranked[,1:100]
  mat_cor_100_a3b <- mat_cor_ranked_a3b[,1:100]
  mat_cor_auto <- mat_cor[,rownames(mat_cor)]
  
  if(doplot){
    print(paste(to_time(), '       ', 'Plotting', cancers[cancer_index]))
    # generate heatmap (all)
    # pdf(file=paste('per cancer/', cancers[cancer_index], '_all.pdf', sep=''), width=50, height=10)
    # hm <- heatmap.2(
    #   mat_cor_NA,
    #   Rowv=F, 
    #   col=redgreen,
    #   key=T, symkey=skey, keysize=1,
    #   cexRow=fontsize, cexCol=fontsize,
    #   density.info="none", trace="none",
    #   margins=c(8,8),
    #   breaks=seq(-1, 1, length.out=101)
    # )
    # dev.off()
    
    # Get hclust output
    # hm_out <- hm$colInd
    
    # Get top 1000 genes from hclust
    # hc1000 <- hm_out[(length(hm_out)-999):(length(hm_out))]
    # mat_cor_hc_1000 <- mat_cor[,hc1000]
    
    # generate heatmap (apobecs)
    pdf(file=paste('per cancer/', cancers[cancer_index], '_apo.pdf', sep=''), width=8, height=8)
    hm <- heatmap.2(
      mat_cor_auto,
      Rowv=F, Colv=F, dendrogram = 'none',
      col=redgreen,
      key=T, symkey=skey, keysize=1,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # generate heatmap (top 100)
    pdf(file=paste('per cancer/', cancers[cancer_index], '_100.pdf', sep=''), width=25, height=6)
    hm <- heatmap.2(
      mat_cor_100,
      Rowv=F, Colv=F, dendrogram = 'none',
      col=redgreen,
      key=T, symkey=skey, keysize=1,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(8,8),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # generate heatmap (top 100, NAs included)
    pdf(file=paste('per cancer/', cancers[cancer_index], '_100_NA.pdf', sep=''), width=25, height=6)
    hm <- heatmap.2(
      mat_cor_100_NA,
      Rowv=F, Colv=F, dendrogram = 'none',
      col=redgreen,
      key=T, symkey=skey, keysize=1,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(8,8),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
    
    # generate heatmap (top 100 A3B)
    pdf(file=paste('per cancer/', cancers[cancer_index], '_100_A3B.pdf', sep=''), width=25, height=6)
    hm <- heatmap.2(
      mat_cor_100_a3b,
      Rowv=F, Colv=F, dendrogram = 'none',
      col=redgreen,
      key=T, symkey=skey, keysize=1,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=c(8,8),
      breaks=seq(-1, 1, length.out=101)
    )
    dev.off()
  }
  
  # Make binary form (correlated or not)
  mat_bin <- t(mat_cor)
  mat_bin[mat_bin>cor_thr] <- 1
  mat_bin[mat_bin<1] <- 0
  mat_bin[is.na(mat_bin)] <- 0
  mat_bin <- mat_bin[order(rowSums(mat_bin), decreasing=T),]
  bin_sum <- matrix(as.numeric(rowSums(mat_bin)), nrow=nrow(mat_bin), ncol=1)
  colnames(bin_sum) <- 'Count'
  mat_bin <- cbind(mat_bin, bin_sum)
  mat_bin <- as.data.frame(mat_bin)
  write.table(mat_bin, file=paste('per cancer/', cancers[cancer_index], '_bin.txt', sep=''), sep='\t')
  write.table(mat_bin[mat_bin$Count>0,], file=paste('per cancer/', cancers[cancer_index], '_list.txt', sep=''), sep='\t')
  write.table(mat_bin[mat_bin$APOBEC3B>0,], file=paste('per cancer/', cancers[cancer_index], '_list_A3B.txt', sep=''), sep='\t')
  
  print(paste(to_time(), '       ', cancers[cancer_index], '- FINISHED'))
  
  return(mat_cor)
  
}

###############################################################################################
# Termination
###############################################################################################

print(Sys.time()-t_last)
print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)