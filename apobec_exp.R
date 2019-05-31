# apobec_exp.R
# 4 threads

# Jiachen Liang
# modified by Sakura

# Generates details about expression of APOBEC genes
# Outputs the mean, median, range, and SD for expression of APOBEC genes

t_start <- Sys.time()

# Root directory
#root <- "/Users/jiachen/Downloads/GDAC/"
root <- "G:/My Drive/Jiachen Files/"
setwd(root)

# Using Windows?
win <- TRUE

# Cancer list
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC",
             "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

###############################################################################################
# Parameters
###############################################################################################

names <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD",
           "LUSC", "PRAD", "STAD", "THCA", "UCEC")
# Omit: OV

year <- c("2016", "2015")
state <- c("tumour", "normal")

###############################################################################################
# APOBEC
###############################################################################################

fun_a <- function(year, state){
  for(n in names){
    dir <- paste("Sakura/data/", year, "/", state, "/apobec/", n, ".txt", sep="")
    genes <- read.table(file = dir, sep = '\t', header = TRUE, row.names = 1)
    #names <- read.table(file = "apobec/data/tr_a.txt", sep = '\t', header = TRUE)[,2]
    #rownames(genes) <- names
    
    output <- matrix(nrow=nrow(genes), ncol=6)
    rownames(output) <- rownames(genes)
    colnames(output) <- c('mean', 'median', 'sd', 'min', 'max', 'range')
    
    for(i in rownames(genes)){
      output[i,'mean'] <- mean(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'median'] <- median(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'sd'] <- sd(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'min'] <- min(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'max'] <- max(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'range'] <- output[i,'max'] - output[i,'min']
    }
    write.table(output, paste("Sakura/output/", year, "/", state, "/expression_details_apobec/", n, ".txt", sep=""), sep='\t')
  }
  return(c(year, state, "APOBEC"))
}

###############################################################################################
# Hypoxia
###############################################################################################

fun_h <- function(year, state){
  for(n in names){
    dir <- paste("Sakura/data/", year, "/", state, "/hypoxia/", n, ".txt", sep="")
    genes <- read.table(file = dir, sep = '\t', header = TRUE, row.names = 1)
    #names <- read.table(file = "apobec/data/tr_h.txt", sep = '\t', header = TRUE)[,2]
    #rownames(genes) <- names
    
    output <- matrix(nrow=nrow(genes), ncol=6)
    rownames(output) <- rownames(genes)
    colnames(output) <- c('mean', 'median', 'sd', 'min', 'max', 'range')
    
    for(i in rownames(genes)){
      output[i,'mean'] <- mean(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'median'] <- median(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'sd'] <- sd(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'min'] <- min(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'max'] <- max(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'range'] <- output[i,'max'] - output[i,'min']
    }
    write.table(output, paste("Sakura/output/", year, "/", state, "/expression_details_hypoxia/", n, ".txt", sep=""), sep='\t')
  }
  return(c(year, state, "Hypoxia"))
}

###############################################################################################
# Parallelization
###############################################################################################

call <- function(input){
  if(input==1){
    fun_a(year[1], state[1])
  }else if(input==2){
    fun_h(year[1], state[1])
  }else if(input==3){
    fun_a(year[1], state[2])
  }else if(input==4){
    fun_h(year[1], state[2])
  }
}

if(win){
  cl <- makeCluster(4)
  registerDoParallel(cl)
  foreach(i=1:4) %dopar% call(i)
}else{
  f1 <- mcparallel(fun_a(year[1], state[1]))
  f2 <- mcparallel(fun_h(year[1], state[1]))
  f3 <- mcparallel(fun_a(year[1], state[2]))
  f4 <- mcparallel(fun_h(year[1], state[2]))
  mccollect(list(f1, f2, f3, f4))
}

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)