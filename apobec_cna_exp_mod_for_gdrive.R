# apobec_cna_exp.R

# Jiachen Liang

# Calculates correlation between APOBEC's CNA and expression

library(gplots)

setwd("G:/My Drive/Jiachen Files/")

# create directories
dir.create("G:/My Drive/Jiachen Files/Sakura/output/cna_cor/")

cancers <- c("BLCA", "BRCA", "HNSC", "KIRC", "LUAD", "LUSC", "UCEC", "OV")
# NO CESC in this study
#cancers <- c("BLCA")

apo.genes <- as.character(read.table(file = "data/tr_a.txt", sep = "\t", header = TRUE, row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]

mat_cor <- matrix(nrow = length(cancers), ncol = 12)
rownames(mat_cor) <- cancers
colnames(mat_cor) <- apo.genes
mat_pv <- matrix(nrow = length(cancers), ncol = 12)
rownames(mat_pv) <- cancers
colnames(mat_pv) <- apo.genes

for(c in 1:length(cancers)){
  
  cancer <- cancers[c]
  
  cna.dir <- paste("G:/My Drive/Jiachen Files/AdditionalFile3/TCGA/2014_01_15/", cancer, 
                   "/tumour/cna_log2_hgncsymbol.txt.gz", sep = "")
  all.cna <- read.table(file  =  cna.dir, sep  =  '\t', header  =  TRUE, row.names  =  1)
  
  exp.dir <- paste('G:/My Drive/Jiachen Files/data/2016/tumour/apobec/', cancer, '.txt', sep = '')
  all.exp <- read.table(file  =  exp.dir, sep  =  ' ', header  =  TRUE, row.names  =  1)
  
  apo.cna <- data.frame(matrix(nrow = 0, ncol = ncol(all.cna)))
  apo.exp <- data.frame(matrix(nrow = 0, ncol = ncol(all.exp)))
  
  for(apo in apo.genes){
    entry <- data.frame(all.cna[apo,])
    rownames(entry) <- apo
    apo.cna <- rbind(apo.cna, entry)
    entry <- data.frame(all.exp[apo,])
    rownames(entry) <- apo
    apo.exp <- rbind(apo.exp, entry)
  }
  
  colnames(apo.cna) <- colnames(all.cna)
  colnames(apo.exp) <- substr(colnames(all.exp), 1, 12)
  
  common.samples <- intersect(colnames(apo.cna), colnames(apo.exp))
  apo.cna <- apo.cna[,common.samples]
  apo.exp <- apo.exp[,common.samples]
  
  cor <- c()
  pv <- c()
  for(apo in apo.genes){
    test <- cor.test(as.numeric(apo.cna[apo,]), as.numeric(apo.exp[apo,]), method = 'spearman')
    cor <- c(cor, as.numeric(test$estimate))
    pv <- c(pv, as.numeric(test$p.value))
  }
  pv <- p.adjust(pv, method = 'fdr')
  
  mat_cor[c,] <- cor
  mat_pv[c,] <- pv
  
  ac_exp_cor <- matrix(nrow = 12, ncol = 12)
  rownames(ac_exp_cor) <- apo.genes
  colnames(ac_exp_cor) <- apo.genes
  ac_exp_pv <- ac_exp_cor
  ac_cna_cor <- ac_exp_cor
  ac_cna_pv <- ac_exp_pv
  
  for(apo1 in 1:length(apo.genes)){
    for(apo2 in 1:length(apo.genes)){
      test <- cor.test(as.numeric(apo.cna[apo1,]), as.numeric(apo.cna[apo2,]), method = 'spearman')
      ac_cna_cor[apo1,apo2] <- as.numeric(test$estimate)
      ac_cna_pv[apo1,apo2] <- as.numeric(test$p.value)
      test <- cor.test(as.numeric(apo.exp[apo1,]), as.numeric(apo.exp[apo2,]), method = 'spearman')
      ac_exp_cor[apo1,apo2] <- as.numeric(test$estimate)
      ac_exp_pv[apo1,apo2] <- as.numeric(test$p.value)
    }
    ac_cna_pv[apo1,] <- p.adjust(ac_cna_pv[apo1,], method = 'fdr')
    ac_exp_pv[apo1,] <- p.adjust(ac_exp_pv[apo1,], method = 'fdr')
  }
  
  ac_cna_cor_unf <- ac_cna_cor
  ac_exp_cor_unf <- ac_exp_cor
  
  ac_cna_cor[ac_cna_pv> = 0.05] <- NA
  ac_exp_cor[ac_exp_pv> = 0.05] <- NA
  

  write.table(ac_cna_cor, file = paste('Sakura/output/cna_cor/ac_cna_', cancer, '_cor.txt', sep = ''), sep = '\t')
  write.table(ac_cna_cor_unf, file = paste('Sakura/output/cna_cor/ac_cna_', cancer, '_cor_unf.txt', sep = ''), sep = '\t')
  write.table(ac_cna_pv, file = paste('Sakura/output/cna_cor/ac_cna_', cancer, '_pv.txt', sep = ''), sep = '\t')
  
  write.table(ac_exp_cor, file = paste('Sakura/output/cna_cor/ac_exp_', cancer, '_cor.txt', sep = ''), sep = '\t')
  write.table(ac_exp_cor_unf, file = paste('Sakura/output/cna_cor/ac_exp_', cancer, '_cor_unf.txt', sep = ''), sep = '\t')
  write.table(ac_exp_pv, file = paste('Sakura/output/cna_cor/ac_exp_', cancer, '_pv.txt', sep = ''), sep = '\t')
  
  ac_cna_hm <- ac_cna_cor
  ac_cna_hm[is.na(ac_cna_hm)] <- 0
  ac_exp_hm <- ac_exp_cor
  ac_exp_hm[is.na(ac_exp_hm)] <- 0
  
  pdf(file = paste('Sakura/output/cna_cor/ac_cna_', cancer, '.pdf', sep = ''), width = 8, height = 8)
  hm <- heatmap.2(
    ac_cna_hm,
    Rowv = NA, Colv = NA,
    dendrogram = "none", density.info = "none", trace = "none",
    col = redgreen,
    breaks = seq(-1, 1, length.out = 101),
    margins = c(7,7)
  )
  dev.off()
  pdf(file = paste('Sakura/output/cna_cor/ac_exp_', cancer, '.pdf', sep = ''), width = 8, height = 8)
  hm <- heatmap.2(
    ac_exp_hm,
    Rowv = NA, Colv = NA,
    dendrogram = "none", density.info = "none", trace = "none",
    col = redgreen,
    breaks = seq(-1, 1, length.out = 101),
    margins = c(7,7)
  )
  dev.off()
  
}

mat_cor_unf <- mat_cor

mat_cor[mat_pv> = 0.05] <- NA

write.table(mat_cor, file = 'Sakura/output/cna_cor/exp_cna_cor.txt', sep = '\t')
write.table(mat_cor_unf, file = 'Sakura/output/cna_cor/exp_cna_cor_unf.txt', sep = '\t')
write.table(mat_pv, file = 'Sakura/output/cna_cor/exp_cna_pv.txt', sep = '\t')

mat_hm <- mat_cor
mat_hm[is.na(mat_hm)] <- 0

pdf(file = 'Sakura/output/cna_cor/exp_cna.pdf', width = 8, height = 8)
hm <- heatmap.2(
  mat_hm,
  Rowv = NA, Colv = NA,
  dendrogram = "none", density.info = "none", trace = "none",
  col = redgreen,
  breaks = seq(-1, 1, length.out = 101),
  margins = c(7,7)
)
dev.off()
