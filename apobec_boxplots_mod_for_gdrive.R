# apobec_boxplots.R
# 4 threads

# Jiachen Liang
# modified by Sakura Maezono

# Generates boxplots with data from apobec_0.3.R and apobec_exp.R
# Compares tumour/normal across cancer types
root <- "G:/My Drive/Jiachen Files/"
setwd(root)

# create new directory
dir.create("G:/My Drive/Jiachen Files/Sakura/output/boxplots/")

# Using Windows?
win <- TRUE

if(win) {
  library(doParallel)
} else {
  library(parallel)
}

library(miscTools)
library(plyr)
library(R.matlab)
library(ggplot2)

###############################################################################################
# Parameters
###############################################################################################

# Cancer types
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", 
             "LUSC", "PRAD", "STAD", "THCA", "UCEC")
# Omit: KIRC, OV, PRAD -- Why?

# Cancer types to compare expression between T/N
#sel <- c("KIRC", "LIHC", "THCA")
sel <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", 
         "LUSC", "PRAD", "STAD", "THCA", "UCEC")

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
apo.genes <- read.table(file  =  "data/tr_a.txt", sep  =  '\t', header  =  TRUE)[,2]
hyp.genes <- read.table(file  =  "data/tr_h.txt", sep  =  '\t', header  =  TRUE)[,2]

# Boxplot colors
bpcol <- c("blue", "red")

# Boxplot dimensions
w <- 15 # width
h <- 10 # height

###############################################################################################
# Incoming data storage
###############################################################################################

exp.mean.apo.t <- matrix(nrow = length(apo.genes), ncol = 0)
rownames(exp.mean.apo.t) <- apo.genes

exp.mean.hyp.t <- matrix(nrow = length(hyp.genes), ncol = 0)
rownames(exp.mean.hyp.t) <- hyp.genes

exp.median.apo.t <- matrix(nrow = length(apo.genes), ncol = 0)
rownames(exp.median.apo.t) <- apo.genes

exp.median.hyp.t <- matrix(nrow = length(hyp.genes), ncol = 0)
rownames(exp.median.hyp.t) <- hyp.genes

exp.mean.apo.n <- matrix(nrow = length(apo.genes), ncol = 0)
rownames(exp.mean.apo.n) <- apo.genes

exp.mean.hyp.n <- matrix(nrow = length(hyp.genes), ncol = 0)
rownames(exp.mean.hyp.n) <- hyp.genes

exp.median.apo.n <- matrix(nrow = length(apo.genes), ncol = 0)
rownames(exp.median.apo.n) <- apo.genes

exp.median.hyp.n <- matrix(nrow = length(hyp.genes), ncol = 0)
rownames(exp.median.hyp.n) <- hyp.genes


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

import_exp <- function(year, state, gene_set, cancer) {
  dir <- paste('data/', year, '/', state, '/', gene_set, '/', cancer, '.txt', sep = '')
  raw <- read.table(file = dir, sep = ' ', header = TRUE, row.names = 1)

  return(raw)
}

import_cor <- function(year, state, method, kind) {
  # going to use what is available for now so instead of "Sakura/output/" use "output" produced by Jiachen
  dir <- paste('output/', year, '/', state, '/', method, '/sig_', kind, '_cor.txt', sep = '')
  raw <- read.table(file = dir, sep = '\t', header = TRUE, row.names = 1)
  
  cancer.names <- colnames(raw)
  raw <- raw[,sel]
  
  return(raw)
}

###############################################################################################
# Merge function
###############################################################################################

merge_nt <- function(mat1, mat2) {
  output <- NULL;
  cur_names <- NULL;
  for(i in 1:ncol(mat1)) {
    output <- cbind(output, mat1[,i])
    output <- cbind(output, mat2[,i])
    colnames(output) <- c(cur_names, paste(colnames(mat1)[i], "_Normal", sep = ''),
                          paste(colnames(mat1)[i], "_Tumour", sep = ''))
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
  require(EnvStats)
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
  plot <- ggplot(df, aes(x = interaction(group, cancer), y = value)) + 
    geom_boxplot(aes(color = group)) +
    scale_color_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
    ggtitle(title) +
    ylab(y_title) + xlab("cancer")
  plot <- plot + facet_grid(.~cancer, scales =  "free_x", switch = "both")
  # remove original x-axis using element_blank()
  plot <- plot + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank()) 
  # add sample size
  plot <- plot + stat_n_text(angle = 90)  
  print(plot)
}


###############################################################################################
# Import expression data
###############################################################################################

for(c in cancers) {
  
  call <- function(input) {
    if (input==1) {
      import_exp(y[1], s[1], g[1], c)
    }else if (input==2) {
      import_exp(y[1], s[2], g[1], c)
    }else if (input==3) {
      import_exp(y[1], s[1], g[2], c)
    }else if (input==4) {
      import_exp(y[1], s[2], g[2], c)
    }
  }
  
  if (win) {
    cl <- makeCluster(4)
    registerDoParallel(cl)
    if (!"OV" %in% cancers) {
      newdflist <- foreach(i = 1:4) %dopar% call(i)
      stopCluster(cl)
      exp.apo.n.raw <- as.matrix(newdflist[2][[1]])
      exp.hyp.n.raw <- as.matrix(newdflist[4][[1]])
      exp.apo.t.raw <- as.matrix(newdflist[1][[1]])
      exp.hyp.t.raw <- as.matrix(newdflist[3][[1]])
    } else {
      newdflist <- foreach(i = c(1,3)) %dopar% call(i)
      stopCluster(cl)
      exp.apo.t.raw <- as.matrix(newdflist[1][[1]])
      exp.hyp.t.raw <- as.matrix(newdflist[2][[1]])
    }

  } else {
    # Import data
    f1 <- mcparallel(import_exp(y[1], s[1], g[1], c))
    f2 <- mcparallel(import_exp(y[1], s[2], g[1], c))
    f3 <- mcparallel(import_exp(y[1], s[1], g[2], c))
    f4 <- mcparallel(import_exp(y[1], s[2], g[2], c))
    
    collected <- mccollect(list(f1, f2, f3, f4))
  }
  

  
  # Append data
  exp.mean.apo.t <- t(rbind.fill.matrix(t(exp.mean.apo.t),
                                        t(rowMeans(exp.apo.t.raw, 
                                                   na.rm = TRUE))))
  exp.median.apo.t <- t(rbind.fill.matrix(t(exp.median.apo.t), 
                                          t(rowMedians(exp.apo.t.raw, 
                                                       na.rm = TRUE))))
  exp.mean.hyp.t <- t(rbind.fill.matrix(t(exp.mean.hyp.t),
                                        t(rowMeans(exp.hyp.t.raw, 
                                                   na.rm = TRUE))))
  exp.median.hyp.t <- t(rbind.fill.matrix(t(exp.median.hyp.t),
                                          t(rowMedians(exp.hyp.t.raw, 
                                                       na.rm = TRUE))))
  if (!"OV" %in% cancers) {
    exp.mean.apo.n <- t(rbind.fill.matrix(t(exp.mean.apo.n), 
                                          t(rowMeans(exp.apo.n.raw,
                                                     na.rm = TRUE))))
    exp.median.apo.n <- t(rbind.fill.matrix(t(exp.median.apo.n),
                                            t(rowMedians(exp.apo.n.raw,
                                                         na.rm = TRUE))))
    exp.mean.hyp.n <- t(rbind.fill.matrix(t(exp.mean.hyp.n),
                                          t(rowMeans(exp.hyp.n.raw,
                                                     na.rm = TRUE))))
    exp.median.hyp.n <- t(rbind.fill.matrix(t(exp.median.hyp.n),
                                            t(rowMedians(exp.hyp.n.raw,
                                                         na.rm = TRUE))))
  } else {
    "skipping... no normal tissues"
  }
  
}

# Row and column names
rownames(exp.mean.apo.t) <- apo.genes
colnames(exp.mean.apo.t) <- cancers
rownames(exp.median.apo.t) <- apo.genes
colnames(exp.median.apo.t) <- cancers
rownames(exp.mean.hyp.t) <- hyp.genes
colnames(exp.mean.hyp.t) <- cancers
rownames(exp.median.hyp.t) <- hyp.genes
colnames(exp.median.hyp.t) <- cancers
if (!"OV" %in% cancers) {
  rownames(exp.mean.apo.n) <- apo.genes
  colnames(exp.mean.apo.n) <- cancers
  rownames(exp.median.apo.n) <- apo.genes
  colnames(exp.median.apo.n) <- cancers
  rownames(exp.mean.hyp.n) <- hyp.genes
  colnames(exp.mean.hyp.n) <- cancers
  rownames(exp.median.hyp.n) <- hyp.genes
  colnames(exp.median.hyp.n) <- cancers
  } else {
  "skipping... no normal tissues"
}

###############################################################################################
# Import correlation data
###############################################################################################

# Mean
call <- function(input) {
  if (input==1) {
    import_cor(y[1], s[1], m[1], k[1])
  }else if (input==2) {
    import_cor(y[1], s[2], m[1], k[1])
  }else if (input==3) {
    import_cor(y[1], s[1], m[2], k[1])
  }else if (input==4) {
    import_cor(y[1], s[2], m[2], k[1])
  }
}

if (win) {
  cl <- makeCluster(4)
  registerDoParallel(cl)
  if (!"OV" %in% cancers) {
    newdflist_mean <- foreach(i = 1:4) %dopar% call(i)
  } else {
    newdflist_mean <- foreach(i = c(1,3)) %dopar% call(i)
  }
  stopCluster(cl)
} else {
  
  f05 <- mcparallel(import_cor(y[1], s[1], m[1], k[1]))
  f06 <- mcparallel(import_cor(y[1], s[2], m[1], k[1]))
  f07 <- mcparallel(import_cor(y[1], s[1], m[2], k[1]))
  f08 <- mcparallel(import_cor(y[1], s[2], m[2], k[1]))
  
  collect_mean <- mccollect(list(f05, f06, f07, f08))
}

# Median
call <- function(input) {
  if (input==1) {
    import_cor(y[1], s[1], m[1], k[2])
  }else if (input==2) {
    import_cor(y[1], s[2], m[1], k[2])
  }else if (input==3) {
    import_cor(y[1], s[1], m[2], k[2])
  }else if (input==4) {
    import_cor(y[1], s[2], m[2], k[2])
  }
}

if (win) {
  cl <- makeCluster(4)
  registerDoParallel(cl)
  if (!"OV" %in% cancers) {
    newdflist_med <- foreach(i = 1:4) %dopar% call(i)
  } else {
    newdflist_med <- foreach(i = c(1,3)) %dopar% call(i)
  }
  stopCluster(cl)
} else {
  
  f09 <- mcparallel(import_cor(y[1], s[1], m[1], k[2]))
  f10 <- mcparallel(import_cor(y[1], s[2], m[1], k[2]))
  f11 <- mcparallel(import_cor(y[1], s[1], m[2], k[2]))
  f12 <- mcparallel(import_cor(y[1], s[2], m[2], k[2]))
  
  collect_median <- mccollect(list(f09, f10, f11, f12))
}

# PC1
call <- function(input) {
  if (input==1) {
    import_cor(y[1], s[1], m[1], k[3])
  }else if (input==2) {
    import_cor(y[1], s[2], m[1], k[3])
  }else if (input==3) {
    import_cor(y[1], s[1], m[2], k[3])
  }else if (input==4) {
    import_cor(y[1], s[2], m[2], k[3])
  }
}

if (win) {
  cl <- makeCluster(4)
  registerDoParallel(cl)
  if (!"OV" %in% cancers) {
    newdflist_PCA <- foreach(i = 1:4) %dopar% call(i)
  } else {
    newdflist_PCA <- foreach(i = c(1,3)) %dopar% call(i)
  }
  stopCluster(cl)
} else {
  
  f13 <- mcparallel(import_cor(y[1], s[1], m[1], k[3]))
  f14 <- mcparallel(import_cor(y[1], s[2], m[1], k[3]))
  f15 <- mcparallel(import_cor(y[1], s[1], m[2], k[3]))
  f16 <- mcparallel(import_cor(y[1], s[2], m[2], k[3]))
  
  collect_pc1 <- mccollect(list(f13, f14, f15, f16))
}


# Assign data
if (win) {
  cor.mean.p.t <- newdflist_mean[1][[1]]
  
  cor.mean.s.t <- newdflist_mean[3][[1]]
  
  cor.median.p.t <- newdflist_med[1][[1]]
  
  cor.median.s.t <- newdflist_med[3][[1]]
  
  cor.pc1.p.t <- newdflist_PCA[1][[1]]
  
  cor.pc1.s.t <- newdflist_PCA[3][[1]]
  
  if (!"OV" %in% cancers) {
    
    cor.mean.p.n <- newdflist_mean[2][[1]]
    
    cor.mean.s.n <- newdflist_mean[4][[1]]
    
    cor.median.p.n <- newdflist_med[2][[1]]
    
    cor.median.s.n <- newdflist_med[4][[1]]
    
    cor.pc1.p.n <- newdflist_PCA[2][[1]]
    
    cor.pc1.s.n <- newdflist_PCA[4][[1]]
    
  } else {
    "skipping... no normal tissues"
  }
  
} else {
  cor.mean.p.t <- collect_mean[1][[1]]
  
  cor.mean.s.t <- collect_mean[3][[1]]
  
  cor.median.p.t <- collect_median[1][[1]]
  
  cor.median.s.t <- collect_median[3][[1]]
  
  cor.pc1.p.t <- collect_pc1[1][[1]]
  
  cor.pc1.s.t <- collect_pc1[3][[1]]
  
  if (!"OV" %in% cancers) {
    
    cor.mean.p.n <- collect_mean[2][[1]]
    
    cor.mean.s.n <- collect_mean[4][[1]]
    
    cor.median.p.n <- collect_median[2][[1]]
    
    cor.median.s.n <- collect_median[4][[1]]
    
    cor.pc1.p.n <- collect_pc1[2][[1]]
    
    cor.pc1.s.n <- collect_pc1[4][[1]]
    
  } else {
    "skipping... no normal tissues"
  }
  
}

###############################################################################################
# Calculate significance
###############################################################################################
if (!"OV" %in% cancers) { 
  sig <- matrix(nrow = 6, ncol = ncol(cor.mean.p.n))
  rownames(sig) <- c("Mean P", "Mean S", "Median P", "Median S", "PC1 P", "PC1 S")
  colnames(sig) <- colnames(cor.mean.p.n)
  # skipping column 3 CESC,  not enough normal tissues n=3, all values are currently NA
  for(i in c(1,2, 4:ncol(cor.mean.p.n))) {
    sig[1,i] <- wilcox.test(as.numeric(cor.mean.p.n[,i]), as.numeric(cor.mean.p.t[,i]))$p.value
    sig[2,i] <- wilcox.test(as.numeric(cor.mean.s.n[,i]), as.numeric(cor.mean.s.t[,i]))$p.value
    sig[5,i] <- wilcox.test(as.numeric(cor.pc1.p.n[,i]), as.numeric(cor.pc1.p.t[,i]))$p.value
    sig[6,i] <- wilcox.test(as.numeric(cor.pc1.s.n[,i]), as.numeric(cor.pc1.s.t[,i]))$p.value
  }
  # for(i in 1:ncol(cor.mean.p.n)) {
  #   sig[3,i] <- wilcox.test(as.numeric(cor.median.p.n[,i]), as.numeric(cor.median.p.t[,i]))$p.value
  #   sig[4,i] <- wilcox.test(as.numeric(cor.median.s.n[,i]), as.numeric(cor.median.s.t[,i]))$p.value
  # }
  
  # skipping column 7 "KIRP" since all values are currently NA
  for(i in c(1,2, 4:6)) {
    sig[3,i] <- wilcox.test(as.numeric(cor.median.p.n[,i]), as.numeric(cor.median.p.t[,i]))$p.value
    sig[4,i] <- wilcox.test(as.numeric(cor.median.s.n[,i]), as.numeric(cor.median.s.t[,i]))$p.value
  }
  # skipping column 14 "UCEC" since all values are currently NA
  for(i in 8:(ncol(cor.mean.p.n)-1)) {
    sig[3,i] <- wilcox.test(as.numeric(cor.median.p.n[,i]), as.numeric(cor.median.p.t[,i]))$p.value
    sig[4,i] <- wilcox.test(as.numeric(cor.median.s.n[,i]), as.numeric(cor.median.s.t[,i]))$p.value
  }
  
  write.table(sig, "Sakura/output/boxplots/pv.txt", sep = "\t")
} else {
  "skipping... no normal tissues"
}
###############################################################################################
# Merge normal and tumour data
###############################################################################################
if (!"OV" %in% cancers) { 
  
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

} else {
  "skipping... no normal tissues"
}
###############################################################################################
# Generate box plots
###############################################################################################
if (!"OV" %in% cancers) { 
  
  pdf(file = "Sakura/output/boxplots/exp_mean_apo.pdf", width = w, height = h)
  bp(exp.mean.apo, "Expression of APOBEC genes (sample means)", "Gene Expression (normalized)")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/exp_mean_hyp.pdf", width = w, height = h)
  bp(exp.mean.hyp, "Expression of hypoxia signature genes (sample means)", "Gene Expression (normalized)")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/exp_median_apo.pdf", width = w, height = h)
  bp(exp.median.apo, "Expression of APOBEC genes (sample medians)", "Gene Expression (normalized)")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/exp_median_hyp.pdf", width = w, height = h)
  bp(exp.median.hyp, "Expression of hypoxia signature genes (sample medians)","Gene Expression (normalized)")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/cor_mean_p.pdf", width = w, height = h)
  bp(cor.mean.p, "Correlation between APOBEC genes expression and mean hypoxia signature (Pearson)", "Correlation score")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/cor_mean_s.pdf", width = w, height = h)
  bp(cor.mean.s, "Correlation between APOBEC genes expression and mean hypoxia signature (Spearman)", "Correlation score")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/cor_median_p.pdf", width = w, height = h)
  bp(cor.median.p, "Correlation between APOBEC genes expression and median hypoxia signature (Pearson)", "Correlation score")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/cor_median_s.pdf", width = w, height = h)
  bp(cor.median.s, "Correlation between APOBEC genes expression and median hypoxia signature (Spearman)", "Correlation score")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/cor_pc1_p.pdf", width = w, height = h)
  bp(cor.pc1.p, "Correlation between APOBEC genes expression and 1st p.c. of hypoxia signature (Pearson)", "Correlation score")
  dev.off()
  
  pdf(file = "Sakura/output/boxplots/cor_pc1_s.pdf", width = w, height = h)
  bp(cor.pc1.s, "Correlation between APOBEC genes expression and 1st p.c. of hypoxia signature (Spearman)", "Correlation score")
  dev.off()
  
  ###############################################################################################
  # Output
  ###############################################################################################
  
  writeMat("Sakura/output/boxplots/merged_exp_mean_apo.MAT", A = exp.mean.apo)
  writeMat("Sakura/output/boxplots/merged_exp_mean_hyp.MAT", A = exp.mean.hyp)
  writeMat("Sakura/output/boxplots/merged_exp_median_apo.MAT", A = exp.median.apo)
  writeMat("Sakura/output/boxplots/merged_exp_median_hyp.MAT", A = exp.median.hyp)
  writeMat("Sakura/output/boxplots/merged_cor_mean_p.MAT", A = cor.mean.p)
  writeMat("Sakura/output/boxplots/merged_cor_mean_s.MAT", A = cor.mean.s)
  writeMat("Sakura/output/boxplots/merged_cor_median_p.MAT", A = cor.median.p)
  writeMat("Sakura/output/boxplots/merged_cor_median_s.MAT", A = cor.median.s)
  writeMat("Sakura/output/boxplots/merged_cor_pc1_p.MAT", A = cor.pc1.p)
  writeMat("Sakura/output/boxplots/merged_cor_pc1_s.MAT", A = cor.pc1.s)
  write.table(as.character(colnames(exp.mean.apo)), "data/cancers_nt.txt", sep = '\t')
} else {
  "skipping... no normal tissues"
}
