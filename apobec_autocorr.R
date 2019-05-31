library(matrixStats)
library(gplots)
library(ggplot2)

# Set working directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
setwd(root)

# Directory: parameters
dir_p <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'

# Directory: expression matrices
dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/apobec/'

# Directory: output
dir_o <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/autocorr/'

# Import parameters
cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep=''), header=T, sep='\t')[,1])

for(c in cancers){
  exp <- read.table(paste(dir_e, c, '.txt', sep=''), sep='\t', header=T, row.names=1)
  exp <- exp[order(rownames(exp)),]
  mat <- matrix(nrow=nrow(exp), ncol=nrow(exp))
  rownames(mat) <- rownames(exp)
  colnames(mat) <- rownames(exp)
  pv <- mat
  for(i in 1:nrow(exp)){
    for(j in 1:nrow(exp)){
      # t <- cor.test(as.numeric(exp[i,]), as.numeric(exp[j,]), method = 'spearman')
      # mat[i,j] <- t$estimate
      # pv[i,j] <- t$p.value
      result = tryCatch({
        t <- cor.test(as.numeric(exp[i,]), as.numeric(exp[j,]), method = 'spearman')
        mat[i,j] <- t$estimate
        pv[i,j] <- t$p.value
      }, error = function(e) {
        mat[i,j] <- NA
        pv[i,j] <- NA
      }, finally = {

      })
    }
  }
  write.table(mat, file=paste(dir_o, c, '_cor.txt', sep=''), sep='\t')
  write.table(pv, file=paste(dir_o, c, '_pv.txt', sep=''), sep='\t')
  mat[pv>=0.05] <- NA
  mat[is.na(mat)] <- 0
  pdf(file=paste(dir_o, c, '_autocorr.pdf', sep=''), width=5, height=5)
  hm <- heatmap.2(
    mat,
    main=c,
    Rowv=NA, Colv=NA,
    dendrogram="none",
    col=redgreen,
    margins=c(8,8),
    density.info="none", trace="none",
    breaks=seq(-1, 1, length.out=101)
  )
  dev.off()
}
