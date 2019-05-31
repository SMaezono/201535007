# apobec_boxplots.R
# 4 threads

# Jiachen Liang

# Generates boxplots with data from apobec_0.3.R and apobec_exp.R
# Compares tumour/normal across cancer types

setwd("~/OneDrive - OnTheHub - The University of Oxford/Buffa Lab")

library(parallel)
library(miscTools)
library(plyr)
library(R.matlab)
library(ggplot2)

###############################################################################################
# Parameters
###############################################################################################

# Cancer types
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
# Omit: KIRC, OV, PRAD

# Cancer types to compare expression between T/N
#sel <- c("KIRC", "LIHC", "THCA")
sel <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

# Year of snapshot
y <- c("2016", "2015")

# State of sample
s <- c("tumour", "normal")

# Gene set
g <- c("apobec", "hypoxia")

# Method for correlation
m <- c("pearson", "spearman")

# Kind of data used for correlation
k <- c("mean", "median", "pc1")

# Gene names
apo.genes <- read.table(file = "apobec/data/tr_a.txt", sep = '\t', header = TRUE)[,2]
hyp.genes <- read.table(file = "apobec/data/tr_h.txt", sep = '\t', header = TRUE)[,2]

# Boxplot colors
bpcol <- c("blue", "red")

# Boxplot dimensions
w <- 15 # width
h <- 10 # height

###############################################################################################
# Incoming data storage
###############################################################################################

exp.mean.apo.n <- matrix(nrow=length(apo.genes), ncol=0)
rownames(exp.mean.apo.n) <- apo.genes
exp.mean.apo.t <- matrix(nrow=length(apo.genes), ncol=0)
rownames(exp.mean.apo.t) <- apo.genes

exp.mean.hyp.n <- matrix(nrow=length(hyp.genes), ncol=0)
rownames(exp.mean.hyp.n) <- hyp.genes
exp.mean.hyp.t <- matrix(nrow=length(hyp.genes), ncol=0)
rownames(exp.mean.hyp.t) <- hyp.genes

exp.median.apo.n <- matrix(nrow=length(apo.genes), ncol=0)
rownames(exp.median.apo.n) <- apo.genes
exp.median.apo.t <- matrix(nrow=length(apo.genes), ncol=0)
rownames(exp.median.apo.t) <- apo.genes

exp.median.hyp.n <- matrix(nrow=length(hyp.genes), ncol=0)
rownames(exp.median.hyp.n) <- hyp.genes
exp.median.hyp.t <- matrix(nrow=length(hyp.genes), ncol=0)
rownames(exp.median.hyp.t) <- hyp.genes

cor.mean.p.n <- NULL
cor.mean.p.t <- NULL

cor.mean.s.n <- NULL
cor.mean.s.t <- NULL

cor.median.p.n <- NULL
cor.median.p.t <- NULL

cor.median.s.n <- NULL
cor.median.s.t <- NULL

cor.pc1.p.n <- NULL
cor.pc1.p.t <- NULL

cor.pc1.s.n <- NULL
cor.pc1.s.t <- NULL

###############################################################################################
# Import functions
###############################################################################################

import_exp <- function(year, state, gene_set, cancer){
  dir <- paste('apobec/data/', year, '/', state, '/', gene_set, '/', cancer, '.txt', sep='')
  raw <- read.table(file=dir, sep='\t', header=TRUE, row.names=1)
  return(raw)
}
import_cor <- function(year, state, method, kind){
  dir <- paste('apobec/output/', year, '/', state, '/', method, '/sig_', kind, '_cor.txt', sep='')
  raw <- read.table(file=dir, sep='\t', header=TRUE, row.names=1)
  
  cancer.names <- colnames(raw)
  raw <- raw[,sel]
  
  return(raw)
}

###############################################################################################
# Merge function
###############################################################################################

merge_nt <- function(mat1, mat2){
  output <- NULL;
  cur_names <- NULL;
  for(i in 1:ncol(mat1)){
    output <- cbind(output, mat1[,i])
    output <- cbind(output, mat2[,i])
    colnames(output) <- c(cur_names, paste(colnames(mat1)[i], "_Normal", sep=''), paste(colnames(mat1)[i], "_Tumour", sep=''))
    cur_names <- colnames(output)
  }
  return(output)
}

###############################################################################################
# Plot function
###############################################################################################

bp1 <- function(input, title, y_title){
  boxplot.matrix(input, na.rm=TRUE, outline=TRUE, col=bpcol, xaxt="n", xlab="", ylab=y_title, main=title)
  axis(1, labels = FALSE)
  text(x =  seq_along(colnames(input)), y = par("usr")[3]-0.3, srt = 90, adj = 1, labels = colnames(exp.mean.apo), xpd = TRUE)
}

bp <- function(input, title, y_title){
  value <- c()
  cancer <- c()
  group <- c()
  for(i in 1:nrow(input)){
    for(j in colnames(input)){
      value <- c(value, input[i,j])
      cancer <- c(cancer, substr(j, 1, 4))
      group <- c(group, substr(j, 6, 11))
    }
  }
  df <- data.frame(value, cancer, group)
  ggplot(df, aes(x=cancer, y=value, fill=group)) + 
    geom_boxplot()
}
###############################################################################################
# Import expression data
###############################################################################################

for(c in cancers){
  
  # Import data
  f1 <- mcparallel(import_exp(y[1], s[1], g[1], c))
  f2 <- mcparallel(import_exp(y[1], s[2], g[1], c))
  f3 <- mcparallel(import_exp(y[1], s[1], g[2], c))
  f4 <- mcparallel(import_exp(y[1], s[2], g[2], c))
  
  collected <- mccollect(list(f1, f2, f3, f4))
  
  exp.apo.t.raw <- as.matrix(collected[1][[1]])
  exp.apo.n.raw <- as.matrix(collected[2][[1]])
  exp.hyp.t.raw <- as.matrix(collected[3][[1]])
  exp.hyp.n.raw <- as.matrix(collected[4][[1]])
  
  # Append data
  exp.mean.apo.t <- t(rbind.fill.matrix(t(exp.mean.apo.t), t(rowMeans(exp.apo.t.raw, na.rm=TRUE))))
  exp.median.apo.t <- t(rbind.fill.matrix(t(exp.median.apo.t), t(rowMedians(exp.apo.t.raw, na.rm=TRUE))))
  exp.mean.apo.n <- t(rbind.fill.matrix(t(exp.mean.apo.n), t(rowMeans(exp.apo.n.raw, na.rm=TRUE))))
  exp.median.apo.n <- t(rbind.fill.matrix(t(exp.median.apo.n), t(rowMedians(exp.apo.n.raw, na.rm=TRUE))))
  exp.mean.hyp.t <- t(rbind.fill.matrix(t(exp.mean.hyp.t), t(rowMeans(exp.hyp.t.raw, na.rm=TRUE))))
  exp.median.hyp.t <- t(rbind.fill.matrix(t(exp.median.hyp.t), t(rowMedians(exp.hyp.t.raw, na.rm=TRUE))))
  exp.mean.hyp.n <- t(rbind.fill.matrix(t(exp.mean.hyp.n), t(rowMeans(exp.hyp.n.raw, na.rm=TRUE))))
  exp.median.hyp.n <- t(rbind.fill.matrix(t(exp.median.hyp.n), t(rowMedians(exp.hyp.n.raw, na.rm=TRUE))))
  
}

# Row and column names
rownames(exp.mean.apo.t) <- apo.genes
colnames(exp.mean.apo.t) <- cancers
rownames(exp.mean.apo.n) <- apo.genes
colnames(exp.mean.apo.n) <- cancers
rownames(exp.median.apo.t) <- apo.genes
colnames(exp.median.apo.t) <- cancers
rownames(exp.median.apo.n) <- apo.genes
colnames(exp.median.apo.n) <- cancers
rownames(exp.mean.hyp.t) <- hyp.genes
colnames(exp.mean.hyp.t) <- cancers
rownames(exp.mean.hyp.n) <- hyp.genes
colnames(exp.mean.hyp.n) <- cancers
rownames(exp.median.hyp.t) <- hyp.genes
colnames(exp.median.hyp.t) <- cancers
rownames(exp.median.hyp.n) <- hyp.genes
colnames(exp.median.hyp.n) <- cancers

###############################################################################################
# Import correlation data
###############################################################################################

# Mean
f05 <- mcparallel(import_cor(y[1], s[1], m[1], k[1]))
f06 <- mcparallel(import_cor(y[1], s[2], m[1], k[1]))
f07 <- mcparallel(import_cor(y[1], s[1], m[2], k[1]))
f08 <- mcparallel(import_cor(y[1], s[2], m[2], k[1]))

collect_mean <- mccollect(list(f05, f06, f07, f08))

# Median
f09 <- mcparallel(import_cor(y[1], s[1], m[1], k[2]))
f10 <- mcparallel(import_cor(y[1], s[2], m[1], k[2]))
f11 <- mcparallel(import_cor(y[1], s[1], m[2], k[2]))
f12 <- mcparallel(import_cor(y[1], s[2], m[2], k[2]))

collect_median <- mccollect(list(f09, f10, f11, f12))

# PC1
f13 <- mcparallel(import_cor(y[1], s[1], m[1], k[3]))
f14 <- mcparallel(import_cor(y[1], s[2], m[1], k[3]))
f15 <- mcparallel(import_cor(y[1], s[1], m[2], k[3]))
f16 <- mcparallel(import_cor(y[1], s[2], m[2], k[3]))

collect_pc1 <- mccollect(list(f13, f14, f15, f16))

# Assign data
cor.mean.p.n <- collect_mean[2][[1]]
cor.mean.p.t <- collect_mean[1][[1]]

cor.mean.s.n <- collect_mean[4][[1]]
cor.mean.s.t <- collect_mean[3][[1]]

cor.median.p.n <- collect_median[2][[1]]
cor.median.p.t <- collect_median[1][[1]]

cor.median.s.n <- collect_median[4][[1]]
cor.median.s.t <- collect_median[3][[1]]

cor.pc1.p.n <- collect_pc1[2][[1]]
cor.pc1.p.t <- collect_pc1[1][[1]]

cor.pc1.s.n <- collect_pc1[4][[1]]
cor.pc1.s.t <- collect_pc1[3][[1]]

###############################################################################################
# Calculate significance
###############################################################################################

sig <- matrix(nrow=6, ncol=ncol(cor.mean.p.n))
rownames(sig) <- c("Mean P", "Mean S", "Median P", "Median S", "PC1 P", "PC1 S")
colnames(sig) <- colnames(cor.mean.p.n)

for(i in 1:ncol(cor.mean.p.n)){
  sig[1,i] <- wilcox.test(cor.mean.p.n[,i], cor.mean.p.t[,i])$p.value
  sig[2,i] <- wilcox.test(cor.mean.s.n[,i], cor.mean.s.t[,i])$p.value
  sig[3,i] <- wilcox.test(cor.median.p.n[,i], cor.median.p.t[,i])$p.value
  sig[4,i] <- wilcox.test(cor.median.s.n[,i], cor.median.s.t[,i])$p.value
  sig[5,i] <- wilcox.test(cor.pc1.p.n[,i], cor.pc1.p.t[,i])$p.value
  sig[6,i] <- wilcox.test(cor.pc1.s.n[,i], cor.pc1.s.t[,i])$p.value
}

write.table(sig, "apobec/output/boxplots/pv.txt", sep="\t")

###############################################################################################
# Merge normal and tumour data
###############################################################################################

exp.mean.apo <- merge_nt(exp.mean.apo.n, exp.mean.apo.t)
exp.mean.hyp <- merge_nt(exp.mean.hyp.n, exp.mean.hyp.t)

exp.median.apo <- merge_nt(exp.median.apo.n, exp.median.apo.t)
exp.median.hyp <- merge_nt(exp.median.hyp.n, exp.median.hyp.t)

cor.mean.p <- merge_nt(cor.mean.p.n, cor.mean.p.t)
cor.mean.s <- merge_nt(cor.mean.s.n, cor.mean.s.t)

cor.median.p <- merge_nt(cor.median.p.n, cor.median.p.t)
cor.median.s <- merge_nt(cor.median.s.n, cor.median.s.t)

cor.pc1.p <- merge_nt(cor.pc1.p.n, cor.pc1.p.t)
cor.pc1.s <- merge_nt(cor.pc1.s.n, cor.pc1.s.t)

###############################################################################################
# Generate box plots
###############################################################################################

pdf(file="apobec/output/boxplots/exp_mean_apo.pdf", width=w, height=h)
bp(exp.mean.apo, "Expression of APOBEC genes (sample means)", "Gene Expression (normalized)")
dev.off()

pdf(file="apobec/output/boxplots/exp_mean_hyp.pdf", width=w, height=h)
bp(exp.mean.hyp, "Expression of hypoxia signature genes (sample means)", "Gene Expression (normalized)")
dev.off()

pdf(file="apobec/output/boxplots/exp_median_apo.pdf", width=w, height=h)
bp(exp.median.apo, "Expression of APOBEC genes (sample medians)", "Gene Expression (normalized)")
dev.off()

pdf(file="apobec/output/boxplots/exp_median_hyp.pdf", width=w, height=h)
bp(exp.median.hyp, "Expression of hypoxia signature genes (sample medians)", "Gene Expression (normalized)")
dev.off()

pdf(file="apobec/output/boxplots/cor_mean_p.pdf", width=w, height=h)
bp(cor.mean.p, "Correlation between APOBEC genes expression and mean hypoxia signature (Pearson)", "Correlation score")
dev.off()

pdf(file="apobec/output/boxplots/cor_mean_s.pdf", width=w, height=h)
bp(cor.mean.s, "Correlation between APOBEC genes expression and mean hypoxia signature (Spearman)", "Correlation score")
dev.off()

pdf(file="apobec/output/boxplots/cor_median_p.pdf", width=w, height=h)
bp(cor.median.p, "Correlation between APOBEC genes expression and median hypoxia signature (Pearson)", "Correlation score")
dev.off()

pdf(file="apobec/output/boxplots/cor_median_s.pdf", width=w, height=h)
bp(cor.median.s, "Correlation between APOBEC genes expression and median hypoxia signature (Spearman)", "Correlation score")
dev.off()

pdf(file="apobec/output/boxplots/cor_pc1_p.pdf", width=w, height=h)
bp(cor.pc1.p, "Correlation between APOBEC genes expression and 1st p.c. of hypoxia signature (Pearson)", "Correlation score")
dev.off()

pdf(file="apobec/output/boxplots/cor_pc1_s.pdf", width=w, height=h)
bp(cor.pc1.s, "Correlation between APOBEC genes expression and 1st p.c. of hypoxia signature (Spearman)", "Correlation score")
dev.off()

###############################################################################################
# Output
###############################################################################################

writeMat("apobec/output/boxplots/merged_exp_mean_apo.MAT", A=exp.mean.apo)
writeMat("apobec/output/boxplots/merged_exp_mean_hyp.MAT", A=exp.mean.hyp)
writeMat("apobec/output/boxplots/merged_exp_median_apo.MAT", A=exp.median.apo)
writeMat("apobec/output/boxplots/merged_exp_median_hyp.MAT", A=exp.median.hyp)
writeMat("apobec/output/boxplots/merged_cor_mean_p.MAT", A=cor.mean.p)
writeMat("apobec/output/boxplots/merged_cor_mean_s.MAT", A=cor.mean.s)
writeMat("apobec/output/boxplots/merged_cor_median_p.MAT", A=cor.median.p)
writeMat("apobec/output/boxplots/merged_cor_median_s.MAT", A=cor.median.s)
writeMat("apobec/output/boxplots/merged_cor_pc1_p.MAT", A=cor.pc1.p)
writeMat("apobec/output/boxplots/merged_cor_pc1_s.MAT", A=cor.pc1.s)
write.table(as.character(colnames(exp.mean.apo)), "apobec/data/cancers_nt.txt", sep='\t')
