# apobec_cna.R
# Jiachen Liang

# Calculates the correlation between copy number alterations in APOBEC genes and the hypoxia signature

library(miscTools)
setwd("/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/")

cancers <- c("BLCA", "BRCA", "COADREAD", "HNSC", "KIRC", "LUAD", "LUSC", "UCEC")
#cancers <- c("BLCA")

apo.genes <- read.table(file="data/tr_a.txt", sep="\t", header=TRUE, row.names=1)$Gene.Symbol
hyp.genes <- read.table(file="data/tr_h.txt", sep="\t", header=TRUE, row.names=1)$Gene.Symbol

m <- "pearson"

all.cor.mean <- matrix(nrow=12, ncol=0)
all.pv.mean <- matrix(nrow=12, ncol=0)
# all.cor.median <- matrix(nrow=12, ncol=0)
# all.pv.median <- matrix(nrow=12, ncol=0)
# all.cor.pc1 <- matrix(nrow=12, ncol=0)
# all.pv.pc1 <- matrix(nrow=12, ncol=0)

for(cancer in cancers){
  
  cna.dir <- paste("/Users/Jiachen/Downloads/AdditionalFile3/TCGA/2014_01_15/", cancer, 
                   "/tumour/cna_log2_hgncsymbol.txt.gz", sep="")
  all.cna <- read.table(file = cna.dir, sep = '\t', header = TRUE, row.names = 1)
  
  apo.cna <- data.frame(matrix(nrow=0, ncol=ncol(all.cna)))
  colnames(apo.cna) <- colnames(all.cna)
  
  for(apo in apo.genes){
    entry <- data.frame(all.cna[apo,])
    rownames(entry) <- apo
    apo.cna <- rbind(apo.cna, entry)
  }
  
  hyp.cna <- data.frame(matrix(nrow=0, ncol=ncol(all.cna)))
  colnames(hyp.cna) <- colnames(hyp.cna)
  
  for(hyp in hyp.genes){
    entry <- data.frame(all.cna[hyp,])
    rownames(entry) <- hyp
    hyp.cna <- rbind(hyp.cna, entry)
  }
  
  # # Remove zeroes
  # for(i in 1:nrow(apo.cna)){
  #   for(j in 1:ncol(apo.cna)){
  #     if(!is.na(apo.cna[i,j])){
  #       if(apo.cna[i,j]==0){
  #         apo.cna[i,j] <- NA
  #       }
  #     }
  #   }
  # }
  # for(i in 1:nrow(hyp.cna)){
  #   for(j in 1:ncol(hyp.cna)){
  #     if(!is.na(hyp.cna[i,j])){
  #       if(hyp.cna[i,j]==0){
  #         hyp.cna[i,j] <- NA
  #       }
  #     }
  #   }
  # }
  # to.del <- c()
  # for(i in 1:nrow(hyp.cna)){
  #   if(is.nan(rowMeans(hyp.cna[i,], na.rm=TRUE))){
  #     to.del <- c(to.del, i)
  #   }
  # }
  # hyp.cna <- hyp.cna[-to.del,]
  # 
  # # Calculate hypoxia PC1
  # hyp.cna.pc <- hyp.cna
  # for(i in 1:nrow(hyp.cna.pc)){
  #   for(j in 1:ncol(hyp.cna.pc)){
  #     if(is.na(hyp.cna.pc[i,j])){
  #       hyp.cna.pc[i,j] <- rowMeans(hyp.cna.pc[i,], na.rm=TRUE)
  #     }
  #   }
  # }
  # pc1 <- princomp(x = t(hyp.cna.pc), scores = TRUE)$scores[,1] ### Error in cov.wt(z) : 'x' must contain finite values only
  # pc1 <- t(pc1)
  
  # Calculate hypoxia medians
  # medians <- data.frame(colMedians(hyp.cna, na.rm=TRUE))
  # medians <- t(medians)
  
  # Calculate hypoxia means
  means <- data.frame(colMeans(hyp.cna, na.rm=TRUE))
  means <- t(means)
  
  # Compute correlation
  total.cor.mean <- matrix(nrow=0, ncol=1)
  total.pv.mean <- matrix(nrow=0, ncol=1)
  # total.cor.median <- matrix(nrow=0, ncol=1)
  # total.pv.median <- matrix(nrow=0, ncol=1)
  # total.cor.pc1 <- matrix(nrow=0, ncol=1)
  # total.pv.pc1 <- matrix(nrow=0, ncol=1)
  
  for(apo in rownames(apo.cna)){
    
    test <- cor.test(as.numeric(apo.cna[apo,]), as.numeric(means), method=m)
    test.cor.mean <- matrix(as.numeric(test$estimate))
    test.pv.mean <- matrix(test$p.value)
    apo.cor.mean <- matrix(nrow=0, ncol=1)
    apo.pv.mean <- matrix(nrow=0, ncol=1)
    apo.cor.mean <- rbind(apo.cor.mean, test.cor.mean)
    apo.pv.mean <- rbind(apo.pv.mean, test.pv.mean)
    
    # test <- cor.test(as.numeric(apo.cna[apo,]), as.numeric(medians), method=m)
    # test.cor.median <- matrix(as.numeric(test$estimate))
    # test.pv.median <- matrix(test$p.value)
    # apo.cor.median <- matrix(nrow=0, ncol=1)
    # apo.pv.median <- matrix(nrow=0, ncol=1)
    # apo.cor.median <- rbind(apo.cor.median, test.cor.median)
    # apo.pv.median <- rbind(apo.pv.median, test.pv.median)
    # 
    # test <- cor.test(as.numeric(apo.cna[apo,]), as.numeric(pc1), method=m)
    # test.cor.pc1 <- matrix(as.numeric(test$estimate))
    # test.pv.pc1 <- matrix(test$p.value)
    # apo.cor.pc1 <- matrix(nrow=0, ncol=1)
    # apo.pv.pc1 <- matrix(nrow=0, ncol=1)
    # apo.cor.pc1 <- rbind(apo.cor.pc1, test.cor.pc1)
    # apo.pv.pc1 <- rbind(apo.pv.pc1, test.pv.pc1)
    
    total.cor.mean <- rbind(total.cor.mean, apo.cor.mean)
    total.pv.mean <- rbind(total.pv.mean, apo.pv.mean)
    # total.cor.median <- rbind(total.cor.median, apo.cor.median)
    # total.pv.median <- rbind(total.pv.median, apo.pv.median)
    # total.cor.pc1 <- rbind(total.cor.pc1, apo.cor.pc1)
    # total.pv.pc1 <- rbind(total.pv.pc1, apo.pv.pc1)
    
  }
  
  rownames(total.cor.mean) <- rownames(apo.cna)
  rownames(total.pv.mean) <- rownames(apo.cna)
  # rownames(total.cor.median) <- rownames(apo.cna)
  # rownames(total.pv.median) <- rownames(apo.cna)
  # rownames(total.cor.pc1) <- rownames(apo.cna)
  # rownames(total.pv.pc1) <- rownames(apo.cna)
  
  all.cor.mean <- cbind(all.cor.mean, total.cor.mean)
  all.pv.mean <- cbind(all.pv.mean, total.pv.mean)
  # all.cor.median <- cbind(all.cor.median, total.cor.median)
  # all.pv.median <- cbind(all.pv.median, total.pv.median)
  # all.cor.pc1 <- cbind(all.cor.pc1, total.cor.pc1)
  # all.pv.pc1 <- cbind(all.pv.pc1, total.pv.pc1)
  
}

colnames(all.cor.mean) <- cancers
colnames(all.pv.mean) <- cancers
# colnames(all.cor.median) <- cancers
# colnames(all.pv.median) <- cancers
# colnames(all.cor.pc1) <- cancers
# colnames(all.pv.pc1) <- cancers
