# apobec_single_cor.R
# n threads

# Jiachen Liang

# Calculate the correlation of APOBEC genes with single genes

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
setwd(root)

# Cancer types
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
#cancers <- 'BRCA'

# Genes to be tested
genes <- c('CA9', 'VEGFA', 'VEGFB', 'VEGFC', 'ENO1', 'SLC2A1')

# Seed genes to be tested
seeds <- c('AKL3', 'LDHA', 'BNIP3', 'SLC2A1', 'CA9', 'VEGF', 'PGK1', 'ENO1', 'HK2', 'ADM')

# APOBEC genes
apo <- read.table(file = "data/tr_a.txt", sep = '\t', header = TRUE)[,2]

# Expression data directory
exp_dir <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'

# Heatmap parameters
margins <- c(8,8) # margins for output image
key <- FALSE # display legend
skey <- TRUE # symmetric key
keysize <- 0.1
fontsize <- 1.2
w <- 5 # width
h <- 5 # height

library(parallel)
library(matrixStats)
library(ggplot2)
library(gplots)

###############################################################################

# Main function

main <- function(c_i){
  
  cancer <- cancers[c_i]
  
  # Import
  dir <- paste(exp_dir, cancer, '_t.txt', sep='')
  raw_exp <- read.table(dir, header=T, row.names=NULL, sep='\t')
  names <- raw_exp[,1]
  names <- make.unique(names, sep='_')
  raw_exp <- raw_exp[,-1]
  rownames(raw_exp) <- names
  
  # Singular analysis
  singular <- sing(cancer, raw_exp)
  
}

###############################################################################

# Singular analysis

sing <- function(cancer, raw_exp){
  
  # Matrix to store correlation values (rho)
  cor <- matrix(nrow=length(apo), ncol=length(genes))
  
  # Matrix to store p values
  pv <- matrix(nrow=length(apo), ncol=length(genes))
  
  for(i in 1:length(apo)){
    for(j in 1:length(genes)){
      gene_a <- as.character(apo[i])
      gene_s <- as.character(genes[j])
      if((rowSums(!is.na(raw_exp[gene_a,])) > 3) & (rowSums(!is.na(raw_exp[gene_s,])) > 3)){ # deal with empty rows
        test <- cor.test(as.numeric(raw_exp[gene_a,]), as.numeric(raw_exp[gene_s,]), method='spearman')
        cor[i,j] <- as.numeric(test$estimate)
        pv[i,j] <- as.numeric(test$p.value)
      }else{
        cor[i,j] <- 0
        pv[i,j] <- 1
      }
    }
  }
  
  rownames(cor) <- apo
  colnames(cor) <- genes
  rownames(pv) <- apo
  colnames(pv) <- genes
  
  # Remove low significance
  for(i in 1:nrow(cor)){
    for(j in 1:ncol(cor)){
      if(pv[i,j] >= 0.05){
        cor[i,j] <- 0
      }
    }
  }
  
  # Output vector
  v <- vector("list", length=2)
  v[[1]] <- cor
  v[[2]] <- pv
  
  # Plot
  filename <- paste('output/singular/', cancer, '.pdf', sep='')
  pdf(file=filename, width=w, height=h)
  hm <- heatmap.2(
    data.matrix(cor),
    #main=paste(cancer, 'single gene correlation (Spearman)'),
    col=redgreen,
    key=key, symkey=skey, keysize=keysize,
    cexRow=fontsize, cexCol=fontsize,
    density.info="none", trace="none",
    margins=margins,
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
  
  return(v)
  
}

###############################################################################

# Seeds analysis



###############################################################################

# Parallelization

f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(main(i))
}
mccollect(f)

###############################################################################


print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)

