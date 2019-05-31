# a3a_a3b.R

# Jiachen Liang

# Compare and correlate the expression of APOBEC 3A and 3B genes in all cancer types

# Import libraries
library(ggplot2)
library(gplots)

# Parameters
margins <- c(8,8) # margins for output image
key <- FALSE # display legend
skey <- TRUE # symmetric key
keysize <- 0.1
fontsize <- 1.2
hm_w <- 12 # width
hm_h <- 12 # height

win <- T

thr <- 0.05 # Significance threshold

if(win){
  # Set working directory
  root <- "C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/a3a_a3b/"
  setwd(root)
  
  # Directory: parameters
  dir_p <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- 'C:/Users/Jiachen Liang/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/apobec/'
}else{
  # Set working directory
  root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/a3a_a3b/"
  setwd(root)
  
  # Directory: parameters
  dir_p <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
  
  # Directory: expression matrices
  dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/apobec/data/2016/tumour/apobec/'
}

mat <- matrix(nrow=length(cancers), ncol=2)
rownames(mat) <- cancers
colnames(mat) <- c('Rho', 'p value')

cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep=''), header=T, sep='\t')[,1])
for(i in 1:length(cancers)){
  
  print(cancers[i])
  
  # Import data
  path <- paste(dir_e, cancers[i], '.txt', sep='')
  exp <- read.table(file=path, header=T, row.names=1, sep='\t')
  
  # Data frame
  df <- as.data.frame(cbind(t(exp['APOBEC3A',]), t(exp['APOBEC3B',])))
  colnames(df) <- c('APOBEC3A', 'APOBEC3B')
  
  # Scatter plot
  sc_plot <- ggplot(df, aes(x=APOBEC3A, y=APOBEC3B)) +
    geom_point(shape=1) +
    geom_smooth(method=lm , color="red", se=TRUE)
  pdf(file=paste(cancers[i], '.pdf', sep=''), width=8, height=8)
  print(sc_plot)
  dev.off()
  
  # Correlate
  test <- cor.test(as.numeric(exp['APOBEC3A',]), as.numeric(exp['APOBEC3B',]))
  mat[i,1] <- test$estimate
  mat[i,2] <- test$p.value
}

write.table(mat, file='correlation.txt', sep='\t')