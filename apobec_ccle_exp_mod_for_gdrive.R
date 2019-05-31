# apobec_ccle_exp.R

# Jiachen Liang
# modified by Sakura Maezono

# Plot out expression of APOBEC genes in CCLE data

t_start <- Sys.time()

# Root directory
root <- "G:/My Drive/Jiachen Files/data/"
setwd(root)
# Output directory
dir.create("G:/My Drive/Jiachen Files/Sakura/output/ccle/")
out_dir <- "G:/My Drive/Jiachen Files/Sakura/output/ccle/"

library(matrixStats)
library(ggplot2)
library(EnvStats)

# Boxplot dimensions
w <- 15 # width
h <- 10 # height

# remove genes in cancer with less than 5 values?
rm_genes <- F
###############################################################################

# Import

# Read
raw <- read.csv("ccle.csv", header = TRUE, row.names = 1)

# Sanitize
raw[raw == Inf] <- NA
raw[raw == '#NAME?'] <- NA
raw[raw == 0] <- NA

###############################################################################


cancers <- c('Prostate', 'Breast', 'Pancreas', 'Colon')

# Separate baseline and fold change

raw_bl <- matrix(nrow = nrow(raw), ncol = 0)
raw_fc <- matrix(nrow = nrow(raw), ncol = 0)

rn_bl <- c()
rn_fc <- c()

for(i in 1:ncol(raw)){
  sample_name <- colnames(raw)[i]
  split <- strsplit(sample_name, '\\.')[[1]]
  if(length(split) == 2){
    # correct for some samples having a dot in their name
    type <- split[2]
    name_cut <- split[1]
  }else{
    type <- split[3]
    name_cut <- paste(split[1], split[2], sep = '_')
  }
  if(type  ==  'FoldChange'){
    raw_fc <- cbind(raw_fc, raw[,i])
    rn_fc <- c(rn_fc, name_cut)
  }else{
    raw_bl <- cbind(raw_bl, raw[,i])
    rn_bl <- c(rn_bl, name_cut)
  }
}

rownames(raw_bl) <- rownames(raw)
rownames(raw_fc) <- rownames(raw)
colnames(raw_bl) <- rn_bl
colnames(raw_fc) <- rn_fc

###############################################################################

# Take log 2 of baseline

raw_bl <- log2(as.matrix(raw_bl))

###############################################################################

# Separate cancer types

separate <- function(input){
  
  exp_prostate <- matrix(nrow = nrow(raw), ncol = 0)
  exp_breast <- matrix(nrow = nrow(raw), ncol = 0)
  exp_pancreas <- matrix(nrow = nrow(raw), ncol = 0)
  exp_colon <- matrix(nrow = nrow(raw), ncol = 0)
  
  for(i in 1:ncol(input)){
    sample_name <- colnames(input)[i]
    split <- strsplit(sample_name, '_')[[1]]
    type <- split[1]
    if(type == 'prostate'){
      exp_prostate <- cbind(exp_prostate, input[,i])
      exp_prostate <- data.frame(exp_prostate) 
    }else if(type == 'breast'){
      exp_breast <- cbind(exp_breast, input[,i])
      exp_breast <- data.frame(exp_breast) 
    }else if(type == 'pancreas'){
      exp_pancreas <- cbind(exp_pancreas, input[,i])
      exp_pancreas <- data.frame(exp_pancreas) 
    }else if(type == 'colon'){
      exp_colon <- cbind(exp_colon, input[,i])
      exp_colon <- data.frame(exp_colon) 
    }
  }
  
  output  =  list(exp_prostate, exp_breast, exp_pancreas, exp_colon)
  
  if (rm_genes) {
    # Remove genes that have less than 5 data points
    for(i in 1:length(output)){ # cycle through cancers
      for(j in 1:nrow(output[[i]])){ # cycle through genes in each cancer
        if(count(!is.na(output[[i]][j,])) < 5){ # if a gene has less than 5 data points that are not NA
          for(k in 1:length(output[[i]][j,])){ # replace all values in that row by NA
            output[[i]][j,k] <- NA
          }
        }
      }
    }
  }

  
  return(output)
}

sep_bl <- separate(raw_bl)
sep_fc <- separate(raw_fc)

###############################################################################

# Plot

plot <- function(input, group, cancer_index){

  c <- cancers[cancer_index]
  
  if(group  ==  'Baseline'){
    y_label <- 'mRNA count (log2 normalized)'
  }else{
    y_label <- 'Fold change in mRNA count (hypoxia vs normoxia)'
  }
  
  values <- c()
  genes <- c()
  
  for(i in 1:nrow(input)){
    this_gene <- rownames(input)[i]
    for(j in 1:ncol(input)){
      values <- c(values, input[i,j])
      genes <- c(genes, this_gene)
    }
  }
  
  df <- data.frame(values, genes)
  
  plot_out <- ggplot(df, aes(x = genes, y = values)) +
    geom_boxplot() +
    ggtitle(paste('Expression of APOBEC genes in', c, '-', group)) +
    xlab("Gene") + ylab(y_label)
  plot_out <- plot_out + stat_n_text() 
  
  pdf(file = paste(out_dir, 'plot_', group, '_', c, '.pdf', sep = ''), width = w, height = h)
  print(plot_out)
  dev.off()
  
}

for(i in 1:length(sep_bl)){
  plot(sep_bl[[i]], 'Baseline', i)
}
for(i in 1:length(sep_fc)){
  plot(sep_fc[[i]], 'Fold change', i)
}

###############################################################################

# Write file to disk with mean/median expression of APOBEC genes

mean_bl <- matrix(nrow = nrow(raw), ncol = 4) # create matrices
mean_fc <- matrix(nrow = nrow(raw), ncol = 4)

rownames(mean_bl) <- rownames(raw)
rownames(mean_fc) <- rownames(raw)

colnames(mean_bl) <- c('Prostate', 'Breast', 'Pancreas', 'Colon')
colnames(mean_fc) <- c('Prostate', 'Breast', 'Pancreas', 'Colon')

median_bl <- matrix(nrow = nrow(raw), ncol = 4)
median_fc <- matrix(nrow = nrow(raw), ncol = 4)

rownames(median_bl) <- rownames(raw)
rownames(median_fc) <- rownames(raw)

colnames(median_bl) <- c('Prostate', 'Breast', 'Pancreas', 'Colon')
colnames(median_fc) <- c('Prostate', 'Breast', 'Pancreas', 'Colon')

for(i in 1:length(sep_bl)){ # fill with mean and median values
  mean_bl[,i] <- rowMeans(sep_bl[[i]], na.rm  =  T)
  median_bl[,i] <- rowMedians(sep_bl[[i]], na.rm  =  T)
}
for(i in 1:length(sep_fc)){
  mean_fc[,i] <- rowMeans(sep_fc[[i]], na.rm  =  T)
  median_fc[,i] <- rowMedians(sep_fc[[i]], na.rm  =  T)
}

mean_bl[is.nan(mean_bl)] <- NA
mean_fc[is.nan(mean_fc)] <- NA
median_bl[is.nan(median_bl)] <- NA
median_fc[is.nan(median_fc)] <- NA

write.table(mean_bl, file = paste(out_dir, 'mean_bl.txt', sep = ''), sep = '\t')
write.table(mean_fc, file = paste(out_dir, 'mean_fc.txt', sep = ''), sep = '\t')
write.table(median_bl, file = paste(out_dir, 'median_bl.txt', sep = ''), sep = '\t')
write.table(median_fc, file = paste(out_dir, 'median_fc.txt', sep = ''), sep = '\t')

###############################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)