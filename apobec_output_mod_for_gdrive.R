# apobec_output.R
# Version 0.1
# 4 threads

# Generates heat maps from output of apobec_0.3.R

# Jiachen Liang
# Buffa Laboratory (Computational Genomics Group)
# Department of Oncology
# University of Oxford

t_start <- Sys.time()

root <- "G:/My Drive/Jiachen Files/"
setwd(root)


library(gplots)


# Using Windows?
win <- TRUE

if(win) {
  library(doParallel)
} else {
  library(parallel)
}

###############################################################################################
# Parameters
###############################################################################################

methods <- c("pearson", "spearman") # correlation methods
margins <- c(8,8) # margins for output image
key <- FALSE # display legend
skey <- TRUE # symmetric key
keysize <- 0.1
fontsize <- 1.2
w <- 5 # width
h <- 5 # height
y <- c("2016", "2015") # data snapshot
s <- c("tumour", "normal") # sample type

thr <- 0.05 # Significance threshold

###############################################################################################
# Main function
###############################################################################################

fun1 <- function(year, state){
  require(gplots)
  
  for(m in methods){
    
    ###########################################################################################
    
    # Import data
    all.cor.mean <- read.table(file = paste("output/", year, "/", state, "/", 
                                            m, "/mean_cor.txt", sep=""),
                               sep = '\t', header = TRUE, row.names = 1)
    all.pv.mean <- read.table(file = paste("output/", year, "/", state, "/", 
                                           m, "/mean_pv.txt", sep=""),
                              sep = '\t', header = TRUE, row.names = 1)
    all.cor.median <- read.table(file = paste("output/", year, "/", state, "/",
                                              m, "/median_cor.txt", sep=""), 
                                 sep = '\t', header = TRUE, row.names = 1)
    all.pv.median <- read.table(file = paste("output/", year, "/", state, "/",
                                             m, "/median_pv.txt", sep=""), 
                                sep = '\t', header = TRUE, row.names = 1)
    all.cor.pc1 <- read.table(file = paste("output/", year, "/", state, "/",
                                           m, "/pc1_cor.txt", sep=""),
                              sep = '\t', header = TRUE, row.names = 1)
    all.pv.pc1 <- read.table(file = paste("output/", year, "/", state, "/",
                                          m, "/pc1_pv.txt", sep=""), 
                             sep = '\t', header = TRUE, row.names = 1)
    
    ###########################################################################################
    
    # Remove NA
    
    for(i in 1:nrow(all.cor.mean)){
      for(j in 1:ncol(all.cor.mean)){
        if(is.na(all.cor.mean[i,j])){
          all.cor.mean[i,j] <- 0
        }
      }
    }
    for(i in 1:nrow(all.cor.median)){
      for(j in 1:ncol(all.cor.median)){
        if(is.na(all.cor.median[i,j])){
          all.cor.median[i,j] <- 0
        }
      }
    }
    for(i in 1:nrow(all.cor.pc1)){
      for(j in 1:ncol(all.cor.pc1)){
        if(is.na(all.cor.pc1[i,j])){
          all.cor.pc1[i,j] <- 0
        }
      }
    }
    
    # Remove insignificant entries
    
    for(i in 1:nrow(all.pv.mean)){
      for(j in 1:ncol(all.pv.mean)){
        if(is.na(all.pv.mean[i,j])){
          all.cor.mean[i,j] <- 0
        }
        else if(all.pv.mean[i,j] > thr){
          all.cor.mean[i,j] <- 0
        }
      }
    }
    for(i in 1:nrow(all.pv.median)){
      for(j in 1:ncol(all.pv.median)){
        if(is.na(all.pv.median[i,j])){
          all.cor.median[i,j] <- 0
        }
        else if(all.pv.median[i,j] > thr){
          all.cor.median[i,j] <- 0
        }
      }
    }
    for(i in 1:nrow(all.pv.pc1)){
      for(j in 1:ncol(all.pv.pc1)){
        if(is.na(all.pv.pc1[i,j])){
          all.cor.pc1[i,j] <- 0
        }
        else if(all.pv.pc1[i,j] > thr){
          all.cor.pc1[i,j] <- 0
        }
      }
    }
    
    ###########################################################################################
    
    # Alphabetical mean
    all.cor.mean <- all.cor.mean[order(rownames(all.cor.mean)),]
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/1_mean_alpha.pdf", sep=""), width=w, height=h)
    hm1 <- heatmap.2(
      data.matrix(all.cor.mean),
      #main=paste("1. Mean -", m),
      Rowv=NA, Colv=NA,
      dendrogram="none",
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    # Alphabetical median
    all.cor.median <- all.cor.median[order(rownames(all.cor.median)),]
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/2_median_alpha.pdf", sep=""), width=w, height=h)
    hm2 <- heatmap.2(
      data.matrix(all.cor.median),
      #main=paste("2. Median -", m),
      Rowv=NA, Colv=NA,
      dendrogram="none",
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    # Alphabetical PC1
    all.cor.pc1 <- all.cor.pc1[order(rownames(all.cor.pc1)),]
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/3_pc1_alpha.pdf", sep=""), width=w, height=h)
    hm3 <- heatmap.2(
      data.matrix(all.cor.pc1),
      #main=paste("3. PC1 -", m),
      Rowv=NA, Colv=NA,
      dendrogram="none",
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    
    ###########################################################################################
    
    # Clustered default mean
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/4_mean_hm2.pdf", sep=""), width=w, height=h)
    hm4 <- heatmap.2(
      data.matrix(all.cor.mean),
      #main=paste("4. Mean (heatmap.2) -", m),
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    # Clustered default median
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/5_median_hm2.pdf", sep=""), width=w, height=h)
    hm5 <- heatmap.2(
      data.matrix(all.cor.median),
      #main=paste("5. Median (heatmap.2) -", m),
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    # Clustered default PC1
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/6_pc1_hm2.pdf", sep=""), width=w, height=h)
    hm6 <- heatmap.2(
      data.matrix(all.cor.pc1),
      #main=paste("6. PC1 (heatmap.2) -", m),
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    
    ###########################################################################################
    
    # Clustered hclust() mean
    mean.d.a <- dist(all.cor.mean, method = "euclidean")
    mean.fit.a <- as.dendrogram(hclust(mean.d.a))
    mean.d.h <- dist(t(all.cor.mean), method = "euclidean")
    mean.fit.h <- as.dendrogram(hclust(mean.d.h))
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/7_mean_hclust.pdf", sep=""), width=w, height=h)
    hm7 <- heatmap.2(
      data.matrix(all.cor.mean),
      #main=paste("7. Mean (hclust) -", m),
      Rowv=mean.fit.a,
      Colv=mean.fit.h,
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    # Clustered hclust() median
    median.d.a <- dist(all.cor.median, method = "euclidean")
    median.fit.a <- as.dendrogram(hclust(median.d.a))
    median.d.h <- dist(t(all.cor.median), method = "euclidean")
    median.fit.h <- as.dendrogram(hclust(median.d.h))
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/8_median_hclust.pdf", sep=""), width=w, height=h)
    hm8 <- heatmap.2(
      data.matrix(all.cor.median),
      #main=paste("8. Median (hclust) -", m),
      Rowv=mean.fit.a,
      Colv=mean.fit.h,
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    # Clustered hclust() PC1
    pc1.d.a <- dist(all.cor.pc1, method = "euclidean")
    pc1.fit.a <- as.dendrogram(hclust(pc1.d.a))
    pc1.d.h <- dist(t(all.cor.pc1), method = "euclidean")
    pc1.fit.h <- as.dendrogram(hclust(pc1.d.h))
    pdf(file=paste("Sakura/output/", year, "/", state, "/", m, "/9_pc1_hclust.pdf", sep=""), width=w, height=h)
    hm9 <- heatmap.2(
      data.matrix(all.cor.pc1),
      #main=paste("9. PC1 (hclust) -", m),
      Rowv=mean.fit.a,
      Colv=mean.fit.h,
      col=redgreen,
      key=key, symkey=skey, keysize=keysize,
      cexRow=fontsize, cexCol=fontsize,
      density.info="none", trace="none",
      margins=margins,
      breaks=seq(-1, 1, length.out=101)
      )
    graphics.off()
  }
  
  ###########################################################################################
  
  # Create new output with only significant values
  
  # Restore NAs
  
  for(i in 1:nrow(all.cor.mean)){
    for(j in 1:ncol(all.cor.mean)){
      if(all.cor.mean[i,j]==0){
        all.cor.mean[i,j] <- NA
      }
    }
  }
  for(i in 1:nrow(all.cor.median)){
    for(j in 1:ncol(all.cor.median)){
      if(all.cor.median[i,j]==0){
        all.cor.median[i,j] <- NA
      }
    }
  }
  for(i in 1:nrow(all.cor.pc1)){
    for(j in 1:ncol(all.cor.pc1)){
      if(all.cor.pc1[i,j]==0){
        all.cor.pc1[i,j] <- NA
      }
    }
  }
  
  for(m in methods){
    # Write to disk
    write.table(all.cor.mean, paste("Sakura/output/", year, "/", state, "/", m, "/sig_mean_cor.txt", sep=""), sep='\t')
    write.table(all.cor.median, paste("Sakura/output/", year, "/", state, "/", m, "/sig_median_cor.txt", sep=""), sep='\t')
    write.table(all.cor.pc1, paste("Sakura/output/", year, "/", state, "/", m, "/sig_pc1_cor.txt", sep=""), sep='\t')
  }
  
  return(state)
}

###############################################################################################
# Parallelization
###############################################################################################

call <- function(input) {
  if (input==1) {
    fun1(y[1], s[1])
  }else if (input==2) {
    fun1(y[1], s[2])
  }
}

if (win) {
  cl <- makeCluster(4)
  registerDoParallel(cl)
  foreach(i = 1:2) %dopar% call(i)
  stopCluster(cl)

} else {
  
  f1 <- mcparallel(fun1(y[1], s[1]))
  f2 <- mcparallel(fun1(y[1], s[2]))
  
  mccollect(list(f1, f2))
}



print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)