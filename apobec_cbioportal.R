library(ggplot2)
library(matrixStats)

setwd('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/')

mut <- read.table('BRCA_apo_mut.txt', sep='\t', header=T, row.names=1)
amp <- read.table('BRCA_apo_amp.txt', sep='\t', header=T, row.names=1)

mut <- mut[,order(colnames(mut))]
amp <- amp[,order(colnames(amp))]

perc_m <- matrix(nrow=ncol(mut), ncol=1)
rownames(perc_m) <- colnames(mut)
colnames(perc_m) <- 'Mutated'

perc_a <- matrix(nrow=ncol(mut), ncol=1)
rownames(perc_a) <- colnames(mut)
colnames(perc_a) <- 'Amplified'
perc_d <- matrix(nrow=ncol(mut), ncol=1)
rownames(perc_d) <- colnames(mut)
colnames(perc_d) <- 'Deleted'

for(i in 1:length(colnames(mut))){
  perc_m[i,1] <- count(!is.na(mut[,i]))
  perc_a[i,1] <- count(amp[,i]>0)
  perc_d[i,1] <- count(amp[,i]<0)
}

perc_m <- perc_m/nrow(mut)
perc_a <- perc_a/nrow(amp)
perc_d <- perc_d/nrow(amp)