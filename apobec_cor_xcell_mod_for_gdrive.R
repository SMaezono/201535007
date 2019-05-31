# apobec_cor_xcell.R

# Jiachen Liang
# Modified by Sakura Maezono

# Can't get the same results as Jiachen's

# Correlation between APOBEC and hypoxia partial to significant xCell signatures

###################################################################################################
# Directories
# Set working directory
root <- "G:/My Drive/Jiachen Files/"
setwd(root)

# Directory: output
dir.create('G:/My Drive/Jiachen Files/Sakura/output/')
dir.create('G:/My Drive/Jiachen Files/Sakura/output/cor_xc/')

dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/cor_xc/'

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'

###################################################################################################
# Parameters

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

# Cancer list
cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep=''),
                                   header=T, sep='\t')[,1])

# APOBEC and hypoxia signatures
apo.genes <- as.character(read.table(file="data/tr_a.txt", sep="\t", 
                                     header=TRUE, row.names=1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file="data/tr_h.txt", sep="\t", 
                                     header=TRUE, row.names=1)$Gene.Symbol)

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

###################################################################################################
# Import expression data

exp_all <- list(length=length(cancers))

for(c in 1:length(cancers)){
  
  cancer <- cancers[c]
  
  dir <- paste(dir_e, cancer, '_t.txt', sep='')
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
  
  exp_all[[c]] <- output

}

###################################################################################################
# Import xCell data

xc_raw <- read.table(file=paste(dir_p, 'xCell_TCGA_RSEM.txt',
                                sep=''), sep='\t', header=T, row.names=1)
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

###################################################################################################
# Main function

# Matrices to store correlation values
cor_p <- matrix(nrow=length(apo.genes), ncol=length(cancers)) # rho (partial)
rownames(cor_p) <- apo.genes
colnames(cor_p) <- cancers
pv_p <- cor_p # pv (partial)
cor_r <- cor_p # rho (regular)
pv_r <- cor_p # pv (regular)
diff <- cor_p # differential between partial and normal

for (c in 1:length(cancers)) {
  
  cancer <- cancers[c]
  
  #################################################################################################
  # Determine which xCell signatures to keep (using thresholds)
  
  xc_mat_full <- xc[[c]]
  to_keep <- c()
  for (i in 1:nrow(xc_mat_full)) {
    if (
      # removed the line na.rm = T
      count(as.numeric(xc_mat_full[i,]) >= thr_xc) # Number of samples with significant signature enrichment
      >=
      (thr_xc_s*ncol(xc_mat_full)) # thr_xc_s% of samples in the dataset
    ) {
      to_keep <- c(to_keep, i) # keep this signature
    }
  }
  xc_mat <- xc_mat_full[to_keep,]
  xc_mat <- xc_mat[setdiff(rownames(xc_mat), exc),]
  
  #################################################################################################
  # Generate mean hypoxia signature
  
  exp_hyp <- exp_all[[c]][hyp.genes,]
  hyp_mean <- colMeans(exp_hyp, na.rm = T)
  
  #################################################################################################
  # Subset APOBECs
  
  exp_apo <- exp_all[[c]][apo.genes,]
  
  #################################################################################################
  # Correlation
  
  for (i in 1:length(apo.genes)) {
    gene <- apo.genes[i]
    # only use rows with no missing values
    temp_mat <- cbind(as.numeric(exp_all[[c]][gene,]), as.numeric(hyp_mean), t(xc_mat))
    not_missing <- !is.na(rowMeans(temp_mat))
    # only test if there are > 3 non-null rows
    if (count(not_missing) >= 3) {
      x <- as.numeric(exp_all[[c]][gene,])[not_missing]
      y <- as.numeric(hyp_mean)[not_missing]
      z <- t(xc_mat)[not_missing,]
      
      # Partial
      test <- pcor.test(x, y, z, method = "spearman")
      cor_p[i,c] <- test$estimate
      pv_p[i,c] <- test$p.value
      
      # Regular
      test <- cor.test(x, y, method = "spearman")
      cor_r[i,c] <- test$estimate
      pv_r[i,c] <- test$p.value
      
    } else {
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
  col=redgreen,
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
  col=redgreen,
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
  col=redgreen,
  breaks=seq(-1, 1, length.out=101),
  margins=c(7,7)
)
dev.off()
