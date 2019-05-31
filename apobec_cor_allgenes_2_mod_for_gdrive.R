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
if(win) {
  library(doParallel)
}else{
  library(parallel)
}

to_time <- function() {
  time <- as.numeric(Sys.time())-as.numeric(t_start)
  time_h <- as.character(floor(time/3600))
  time_m <- as.character(floor((time%%3600)/60))
  time_s <- as.character(floor((time%%3600)%%60))
  if(nchar(time_h) == 1) {time_h<-paste('0',time_h,sep = '')}
  if(nchar(time_m) == 1) {time_m<-paste('0',time_m,sep = '')}
  if(nchar(time_s) == 1) {time_s<-paste('0',time_s,sep = '')}
  
  time_entry <- paste(
    time_h, ':',
    time_m, ':',
    time_s,
    sep = ''
  )
  
  return(time_entry)
}

###############################################################################################
# Parameters
###############################################################################################

# Significance threshold
thr <- 0.05

# Correlation threshold
cor_thr <- 0.3

# Filter for expression?
do_filter <- T

###############################################################################################
# File locations
###############################################################################################
# create folder(s) if it doesn't/ they don't exist
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all/")
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cor_all/Buffa 52/")

# Set working directory
root <- "G:/My Drive/Jiachen Files/Sakura/output/cor_all/"
setwd(root)

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'

# List of tested genes
sign <- as.character(read.table(paste(dir_p, 'sign_1.txt', sep = ''), header=T,
                                sep = '\t')[,1])
sign <- sign[order(sign)]



# APOBEC and hypoxia genes
dir_ext <- "G:/My Drive/Jiachen Files/"
apo.genes <- as.character(read.table(file = paste(dir_ext, "data/tr_a.txt", sep=""),
                                     sep = "\t", header = TRUE, row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file = paste(dir_ext,"data/tr_h.txt", sep=""),
                                     sep = "\t", header = TRUE, row.names = 1)$Gene.Symbol)
###############################################################################################
# Import expression matrix from file system
###############################################################################################

cancer <- 'BRCA'
  
state <- 't'

dir <- paste(dir_e, cancer, '_', state, '.txt', sep = '')
input <- read.table(dir, sep = '\t', header=T, row.names=NULL)

rn <- input[,1]
rn <- make.unique(rn, sep = '.')

exp <- input[,-1]
rownames(exp) <- rn

# Fix col names
names <- colnames(exp)
for(i in 1:length(names)) {
  split <- strsplit(names[i], '\\.')
  names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
}
colnames(exp) <- names

# calculate median
exp_med <- as.numeric(as.matrix(exp))
exp_med[is.na(exp_med)] <- 0
med <- median(exp_med)

###############################################################################################
# Main function
###############################################################################################

# output matrix
mat_cor <- matrix(nrow=length(sign), ncol=nrow(exp))
mat_pv <- matrix(nrow=length(sign), ncol=nrow(exp))
rownames(mat_cor) <- sign
rownames(mat_pv) <- sign
colnames(mat_cor) <- rownames(exp)
colnames(mat_pv) <- rownames(exp)

# calculate correlation
calculate <- function(i) {
  # list of cor and pv for current gene vs all genes
  cor <- matrix(nrow=1, ncol=nrow(exp))
  pv <- matrix(nrow=1, ncol=nrow(exp))
  
  for(j in 1:nrow(exp)) { # for all genes
    # ERROR: mean(exp[j,], na.rm=T)<med results to NA instead of TRUE or false
    if((mean(exp[j,], na.rm=T)<med) & do_filter) {
      cor[1,j] <- NA
      pv[1,j] <- NA
    }else{
      result = tryCatch({
        test <- cor.test(as.numeric(exp[sign[i],]), as.numeric(exp[j,]))
        cor[1,j] <- test$estimate
        pv[1,j] <- test$p.value
      }, warning = function(w) {
        test <- cor.test(as.numeric(exp[sign[i],]), as.numeric(exp[j,]))
        cor[1,j] <- test$estimate
        pv[1,j] <- test$p.value
      }, error = function(e) {
        cor[1,j] <- NA
        pv[1,j] <- NA
      }, finally = {
        
      })
    }
  }
  
  cor <- as.numeric(cor[1,])
  pv <- as.numeric(pv[1,])
  
  # FDR
  pv <- p.adjust(pv, method='fdr')
  
  return(list(cor, pv))
}

# calculation call
cores <- min(detectCores(), length(sign))
cl <- makeCluster(cores)
registerDoParallel(cl)
result <- foreach(n=1:length(sign)) %dopar% calculate(n)
stopCluster(cl)

# put in matrix
for(i in 1:length(sign)) {
  mat_cor[i,] <- result[[i]][[1]]
  mat_pv[i,] <- result[[i]][[2]]
}

write.table(mat_cor, file = paste('Buffa52/', cancer, '_cor.txt',
                                  sep = ''), sep = '\t')
write.table(mat_pv, file = paste('Buffa52/', cancer, '_pv.txt', 
                                 sep = ''), sep = '\t')

print(paste(to_time(), '       ', cancer, '- Calculation done'))

# filter out insignificant values
mat_cor[mat_pv>thr] <- NA

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
write.table(mat_bin, file = paste('Buffa52/', cancer, '_bin.txt',
                                  sep = ''), sep = '\t')
write.table(mat_bin[mat_bin$Count>0,], file = paste('Buffa52/',
                                                    cancer, '_list.txt',
                                                    sep = ''), sep = '\t')

# Make summed form (sum of all rho values)
mat_sum <- as.matrix(rowSums(t(mat_cor), na.rm=T))
rownames(mat_sum) <- colnames(mat_cor)
colnames(mat_sum) <- 'Sum of Rho'
# filter for genes at least correlated with one gene in signature
mat_sum <- mat_sum[rownames(mat_bin[mat_bin$Count>0,]),] 
mat_sum <- mat_sum[order(-mat_sum)]
write.table(mat_sum, file = paste('Buffa52/', cancer, '_sum.txt', 
                                  sep = ''), sep = '\t')
write.table(mat_sum[1:100], file = paste('Buffa52/', cancer, 
                                         '_sum_top_100.txt', 
                                         sep = ''), sep = '\t')
write.table(mat_sum[1:500], file = paste('Buffa52/', cancer, 
                                         '_sum_top_500.txt', 
                                         sep = ''), sep = '\t')
write.table(mat_sum[1:1000], file = paste('Buffa52/', cancer,
                                          '_sum_top_1000.txt',
                                          sep = ''), sep = '\t')
write.table(mat_sum[mat_sum>=10], file = paste('Buffa52/', cancer, 
                                               '_sum_above10.txt', 
                                               sep = ''), sep = '\t')
write.table(mat_sum[mat_sum>=1], file = paste('Buffa52/', cancer,
                                              '_sum_above1.txt', 
                                              sep = ''), sep = '\t')
write(names(mat_sum[1:100]), file = paste('Buffa52/', cancer, 
                                          '_genes_top_100.txt', sep = ''))
write(names(mat_sum[1:500]), file = paste('Buffa52/', cancer, 
                                          '_genes_top_500.txt', sep = ''))
write(names(mat_sum[1:1000]), file = paste('Buffa52/', cancer,
                                           '_genes_top_1000.txt', sep = ''))
write(names(mat_sum[mat_sum>=10]), file = paste('Buffa52/', cancer, 
                                                '_genes_above10.txt', sep = ''))
write(names(mat_sum[mat_sum>=1]), file = paste('Buffa52/', cancer, 
                                               '_genes_above1.txt', sep = ''))
