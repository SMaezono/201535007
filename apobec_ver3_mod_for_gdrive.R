# apobec.R
# Version 0.3
# 4 threads

# Calculates the correlation between APOBEC gene expression and the hypoxia signature

# Jiachen Liang
# Buffa Laboratory (Computational Genomics Group)
# Department of Oncology
# University of Oxford

t_start <- Sys.time()

root <- "G:/My Drive/Jiachen Files/"
setwd(root)

# Windows or Mac?
win <- TRUE 

library(miscTools)
if(win) {
  library(doParallel)
} else {
  library(parallel)
}
library(plyr)



###############################################################################################
# Parameters
###############################################################################################

# Cancer types
dir_p <- "G:/My Drive/Jiachen Files/data/"
names <- as.character(read.table(paste(dir_p, 'cancer_list.txt',
                                       sep=''), header=T, sep='\t')[,1])
names <- names[-15]

names <- "OV"
# Omit: OV

# Correlation methods used
methods <- c("pearson", "spearman")

# Data snapshots
y <- c("2016", "2015")

# Sample type
s <- c("tumour", "normal")

# Gene names
names.apo <- read.table(file = paste(dir_p, "tr_a.txt", sep=""), sep = '\t',
                        header = TRUE)[,2] # UTP11, AK4, HILPDA, MRGBP, CTSV
names.hyp <- read.table(file = paste(dir_p,"tr_h.txt", sep=""), sep = '\t',
                        header = TRUE)[,2]

# Median expression values (???) --- FILE "2016/median_exp.txt" DOESN'T EXIST
# file changed to "exp_med.txt" for now
medians <- read.table(file=paste(dir_p, "exp_med.txt", sep=""),
                      sep="\t", header=TRUE, row.names=1)

# MEDIANS??? Can't seem to find the right file

###############################################################################################
# Main function
###############################################################################################

fun1 <- function(m, year, state){
  
  count <- 1
  
  # Output matrices
  all.cor.mean <- matrix(nrow=12, ncol=0)
  all.pv.mean <- matrix(nrow=12, ncol=0)
  all.cor.median <- matrix(nrow=12, ncol=0)
  all.pv.median <- matrix(nrow=12, ncol=0)
  all.cor.pc1 <- matrix(nrow=12, ncol=0)
  all.pv.pc1 <- matrix(nrow=12, ncol=0)
  rownames(all.cor.mean) <- names.apo
  rownames(all.pv.mean) <- names.apo
  rownames(all.cor.median) <- names.apo
  rownames(all.pv.median) <- names.apo
  rownames(all.cor.pc1) <- names.apo
  rownames(all.pv.pc1) <- names.apo
  
  for(cancer in names){
    
    print(paste("Starting", m, count, "of 15", "-", cancer))
    
    ###########################################################################################
    
    # Import
    dir.apo <- paste(dir_p, year, "/", state, "/apobec/", cancer, ".txt", sep="")
    dir.hyp <- paste(dir_p, year, "/", state, "/hypoxia/", cancer, ".txt", sep="")
    genes.apo <- read.table(file = dir.apo, sep = '\t', header = TRUE, row.names = 1)
    genes.hyp <- read.table(file = dir.hyp, sep = '\t', header = TRUE, row.names = 1)
    
    #rownames(genes.apo) <- names.apo
    #rownames(genes.hyp) <- names.hyp
    
    genes.apo <- genes.apo[,order(colnames(genes.apo))]
    genes.hyp <- genes.hyp[,order(colnames(genes.hyp))]
    
    thr <- as.numeric(medians[cancer,1])
    
    ###########################################################################################
    
    # # Remove zeroes
    # for(i in 1:nrow(genes.apo)){
    #   for(j in 1:ncol(genes.apo)){
    #     if(genes.apo[i,j]==0){
    #       genes.apo[i,j] <- NA
    #     }
    #   }
    # }
    # for(i in 1:nrow(genes.hyp)){
    #   for(j in 1:ncol(genes.hyp)){
    #     if(genes.hyp[i,j]==0){
    #       genes.hyp[i,j] <- NA
    #     }
    #   }
    # }
    
    # # Remove NaN
    # for(i in 1:nrow(genes.apo)){
    #   for(j in 1:ncol(genes.apo)){
    #     if(is.nan(genes.apo[i,j])){
    #       genes.apo[i,j] <- NA
    #     }
    #   }
    # }
    # for(i in 1:nrow(genes.hyp)){
    #   for(j in 1:ncol(genes.hyp)){
    #     if(is.nan(genes.hyp[i,j])){
    #       genes.hyp[i,j] <- NA
    #     }
    #   }
    # }
    
    to_mask <- c()
    
    for(i in 1:nrow(genes.apo)){
      rm <- mean(as.numeric(genes.apo[i,]), na.rm=TRUE)
      if(is.na(rm)){
        rm <- 0
      }
      if(rm < thr){
        to_mask <- c(to_mask, rownames(genes.apo)[i])
      }
    }
    
    
    ###########################################################################################
    
    # Calculate hypoxia PC1
    genes.hyp.pc <- genes.hyp
    for(i in 1:nrow(genes.hyp.pc)){
      for(j in 1:ncol(genes.hyp.pc)){
        if(is.na(genes.hyp.pc[i,j])){
          genes.hyp.pc[i,j] <- rowMeans(genes.hyp.pc[i,], na.rm=TRUE)
        }
      }
    }
    #pc1 <- princomp(x = t(genes.hyp.pc), scores = TRUE)$scores[,1]
    pc1 <- prcomp(x = t(genes.hyp.pc))$x[,1]
    pc1 <- t(pc1)
    
    # Calculate hypoxia medians
    medians <- data.frame(colMedians(genes.hyp, na.rm=TRUE))
    medians <- t(medians)
    
    # Calculate hypoxia means
    means <- data.frame(colMeans(genes.hyp, na.rm=TRUE))
    means <- t(means)
    
    ###########################################################################################
    
    # Compute correlation
    total.cor.mean <- matrix(nrow=0, ncol=1)
    total.pv.mean <- matrix(nrow=0, ncol=1)
    total.cor.median <- matrix(nrow=0, ncol=1)
    total.pv.median <- matrix(nrow=0, ncol=1)
    total.cor.pc1 <- matrix(nrow=0, ncol=1)
    total.pv.pc1 <- matrix(nrow=0, ncol=1)
    
    for(apo in rownames(genes.apo)){
      
      # Deal with empty rows
      if(rowSums(!is.na(genes.apo[apo,])) > 3){ 
        # select for rows with more than 3 elements (threshold for correlation test)
        # Test: mean
        test <- cor.test(as.numeric(genes.apo[apo,]), as.numeric(means), method=m)
        test.cor.mean <- matrix(as.numeric(test$estimate))
        test.pv.mean <- matrix(test$p.value)
        apo.cor.mean <- matrix(nrow=0, ncol=1)
        apo.pv.mean <- matrix(nrow=0, ncol=1)
        apo.cor.mean <- rbind(apo.cor.mean, test.cor.mean)
        apo.pv.mean <- rbind(apo.pv.mean, test.pv.mean)
        
        # Test: median
        test <- cor.test(as.numeric(genes.apo[apo,]), as.numeric(medians), method=m)
        test.cor.median <- matrix(as.numeric(test$estimate))
        test.pv.median <- matrix(test$p.value)
        apo.cor.median <- matrix(nrow=0, ncol=1)
        apo.pv.median <- matrix(nrow=0, ncol=1)
        apo.cor.median <- rbind(apo.cor.median, test.cor.median)
        apo.pv.median <- rbind(apo.pv.median, test.pv.median)
        
        # Test: PC1
        test <- cor.test(as.numeric(genes.apo[apo,]), as.numeric(pc1), method=m)
        test.cor.pc1 <- matrix(as.numeric(test$estimate))
        test.pv.pc1 <- matrix(test$p.value)
        apo.cor.pc1 <- matrix(nrow=0, ncol=1)
        apo.pv.pc1 <- matrix(nrow=0, ncol=1)
        apo.cor.pc1 <- rbind(apo.cor.pc1, test.cor.pc1)
        apo.pv.pc1 <- rbind(apo.pv.pc1, test.pv.pc1)
      }
      else{
        apo.cor.mean <- NA
        apo.pv.mean <- NA
        apo.cor.median <- NA
        apo.pv.median <- NA
        apo.cor.pc1 <- NA
        apo.pv.pc1 <- NA
      }
      
      # Assemble
      total.cor.mean <- rbind(total.cor.mean, apo.cor.mean)
      total.pv.mean <- rbind(total.pv.mean, apo.pv.mean)
      total.cor.median <- rbind(total.cor.median, apo.cor.median)
      total.pv.median <- rbind(total.pv.median, apo.pv.median)
      total.cor.pc1 <- rbind(total.cor.pc1, apo.cor.pc1)
      total.pv.pc1 <- rbind(total.pv.pc1, apo.pv.pc1)
      
    }
    
    ###########################################################################################
    
    # Assemble
    rownames(total.cor.mean) <- rownames(genes.apo)
    rownames(total.pv.mean) <- rownames(genes.apo)
    rownames(total.cor.median) <- rownames(genes.apo)
    rownames(total.pv.median) <- rownames(genes.apo)
    rownames(total.cor.pc1) <- rownames(genes.apo)
    rownames(total.pv.pc1) <- rownames(genes.apo)
    
    # MASK ROWS
    total.cor.mean[to_mask,] <- NA
    total.pv.mean[to_mask,] <- NA
    total.cor.median[to_mask,] <- NA
    total.pv.median[to_mask,] <- NA
    total.cor.pc1[to_mask,] <- NA
    total.pv.pc1[to_mask,] <- NA
    
    all.cor.mean <- t(rbind.fill.matrix(t(all.cor.mean), t(total.cor.mean)))
    all.pv.mean <- t(rbind.fill.matrix(t(all.pv.mean), t(total.pv.mean)))
    all.cor.median <- t(rbind.fill.matrix(t(all.cor.median), t(total.cor.median)))
    all.pv.median <- t(rbind.fill.matrix(t(all.pv.median), t(total.pv.median)))
    all.cor.pc1 <- t(rbind.fill.matrix(t(all.cor.pc1), t(total.cor.pc1)))
    all.pv.pc1 <- t(rbind.fill.matrix(t(all.pv.pc1), t(total.pv.pc1)))
    
    count <- count + 1
  }
  
  #############################################################################################
  
  colnames(all.cor.mean) <- names
  colnames(all.pv.mean) <- names
  colnames(all.cor.median) <- names
  colnames(all.pv.median) <- names
  colnames(all.cor.pc1) <- names
  colnames(all.pv.pc1) <- names
  
  # Write to disk
  write.table(all.cor.mean, paste("Sakura/output/", year, "/", state, "/", m,
                                  "/mean_cor.txt", sep=""), sep='\t')
  write.table(all.pv.mean, paste("Sakura/output/", year, "/", state, "/", m,
                                 "/mean_pv.txt", sep=""), sep='\t')
  write.table(all.cor.median, paste("Sakura/output/", year, "/", state, "/", m,
                                    "/median_cor.txt", sep=""), sep='\t')
  write.table(all.pv.median, paste("Sakura/output/", year, "/", state, "/", m,
                                   "/median_pv.txt", sep=""), sep='\t')
  write.table(all.cor.pc1, paste("Sakura/output/", year, "/", state, "/", m,
                                 "/pc1_cor.txt", sep=""), sep='\t')
  write.table(all.pv.pc1, paste("Sakura/output/", year, "/", state, "/", m,
                                "/pc1_pv.txt", sep=""), sep='\t')
  
  return(m)
}

###############################################################################################
# Parallelization
###############################################################################################

#f1 <- mcparallel(fun1(methods[1], y[1], s[1]))
#f2 <- mcparallel(fun1(methods[2], y[1], s[1]))
#mccollect(list(f1, f2))

#f <- vector("list", length=length(methods))
#for(i in 1:length(methods)){
#  f[[i]] <- mcparallel(fun1(methods[i], y[1], s[1]))
#}
#mccollect(f)

call <- function(input) {
  if (input==1) {
    fun1(methods[1], y[1], s[1])
  }else if (input==2) {
    fun1(methods[2], y[1], s[1])
  }else if (input==3) {
    fun1(methods[1], y[1], s[2])
  }else if (input==4) {
    fun1(methods[2], y[1], s[2])
  }
}

if (win) {
  cl <- makeCluster(4)
  registerDoParallel(cl)
  if (!"OV" %in% names) {
    foreach(i = 1:4) %dopar% call(i)
  } else {
    foreach(i = 1:2) %dopar% call(i)
  }
  
} else {
  f1 <- mcparallel(fun1(methods[1], y[1], s[1]))
  f2 <- mcparallel(fun1(methods[2], y[1], s[1]))
  f3 <- mcparallel(fun1(methods[1], y[1], s[2]))
  f4 <- mcparallel(fun1(methods[2], y[1], s[2]))

mccollect(list(f1, f2, f3, f4))
}
print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)