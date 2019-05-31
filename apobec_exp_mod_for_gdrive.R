# apobec_exp.R
# 4 threads

# Jiachen Liang
# Modified by Sakura Maezono

# Generates details about expression of APOBEC/hypoxia genes
# Outputs the mean, median, range, and SD for expression of APOBEC genes

t_start <- Sys.time()

# Root directory
#root <- "/Users/jiachen/Downloads/GDAC/"
root <- 'G:/My Drive/Jiachen Files/GDAC/'
setwd(root)

# Using Windows?
win <- TRUE

# IMPORTANT: when running cancers without normal samples (e.g. "OV"),
# run them separately

# Cancer list
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", 
             "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

cancers <- "OV"

###############################################################################################
# Parameters
###############################################################################################

names <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", 
           "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

names <- c("OV")


year <- c("2016", "2015")

state <- c("tumour", "normal")

dir_i <- 'G:/My Drive/Jiachen Files/data/'


# create directories
# if year 2016
dir.create('G:/My Drive/Jiachen Files/Sakura/output/2016/')
dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/2016/'

create_dir <- function(state) {
  dir.create(paste(dir_o, state, '/', sep = ""))
  dir.create(paste(dir_o, state, '/', 'expression_details_apobec/', sep = ""))
  dir.create(paste(dir_o, state, '/',  'expression_details_hypoxia/', sep = ""))
  dir.create(paste(dir_o, state, '/',  'pearson/', sep = ""))
  dir.create(paste(dir_o, state, '/',  'spearman/', sep = ""))
}

lapply(state, create_dir)

# retrieve expressions of individual gene signatures
# cancers <- c("CESC", "OV")
# looks like this is using a different file (similar but not the same)
retrieve_indiv_gene_exp <- function(cancer) {
  # retrieve expressions of individual gene signatures (i.e. APOBEC, Hypoxia)
  # retrieve APOBEC gene sig only
  apo.genes <- as.character(read.table(file = paste(dir_i, "tr_a.txt",
                                                    sep = ""), sep = "\t", 
                                       header = TRUE, row.names = 1)$Gene.Symbol)
  apo.genes <- apo.genes[order(apo.genes)]
  # row.names = NULL forces row numbering 
  # to be able to import files with duplicated rownames
  
  if (!cancer == "OV") {
    # APOBEC normal
    exp_apo_n <- read.table(paste(cancer, "_n.txt", sep = ""), sep = '\t', 
                            header = TRUE, row.names = NULL)
    exp_apo_n <- exp_apo_n[(exp_apo_n$row.names %in% apo.genes),]
    rownames(exp_apo_n) <- exp_apo_n$row.names
    exp_apo_n <- exp_apo_n[,-1]
    write.table(exp_apo_n, paste(dir_i,  "2016/normal/apobec/", 
                                 cancer, ".txt", sep = ""))  
  } else {
    print("skipping normal tissues")
  }
  
  # APOBEC tumor
  exp_apo_t <- read.table(paste(cancer, "_t.txt", sep = "") , sep = '\t', 
                          header = TRUE, row.names = NULL)
  exp_apo_t <- exp_apo_t[(exp_apo_t$row.names %in% apo.genes),]
  rownames(exp_apo_t) <- exp_apo_t$row.names
  exp_apo_t <- exp_apo_t[,-1]
  write.table(exp_apo_t, paste(dir_i,  "2016/tumour/apobec/", 
                               cancer, ".txt", sep = ""))   
  
  
  # retrieve APOBEC gene sig only
  hyp.genes <- as.character(read.table(file = paste(dir_i, "tr_h.txt", sep = ""),
                                       sep = "\t", header = TRUE, 
                                       row.names = 1)$Gene.Symbol)
  
  hyp.genes <- hyp.genes[order(hyp.genes)]
  
  if (!cancer == "OV") {
    # hypoxia normal
    exp_hyp_n <- read.table(paste(cancer, "_n.txt", sep = ""), sep = '\t', 
                            header = TRUE, row.names = NULL)
    exp_hyp_n <- exp_hyp_n[(exp_hyp_n$row.names %in% hyp.genes),]
    rownames(exp_hyp_n) <- exp_hyp_n$row.names
    exp_hyp_n <- exp_hyp_n[,-1]
    write.table(exp_hyp_n, paste(dir_i,  "2016/normal/hypoxia/", 
                                 cancer, ".txt", sep = ""))  
  } else {
    print("skipping normal tissues")
  }
  
  # hypoxia tumor
  exp_hyp_t <- read.table(paste(cancer, "_t.txt", sep = "") , sep = '\t', 
                          header = TRUE, row.names = NULL)
  exp_hyp_t <- exp_hyp_t[(exp_hyp_t$row.names %in% hyp.genes),]
  rownames(exp_hyp_t) <- exp_hyp_t$row.names
  exp_hyp_t <- exp_hyp_t[,-1]
  write.table(exp_hyp_t, paste(dir_i,  "2016/tumour/hypoxia/",
                               cancer, ".txt", sep = "")) 
}
  
lapply(cancers, retrieve_indiv_gene_exp)

###############################################################################################
# APOBEC
###############################################################################################

fun_a <- function(year, state) {
  for(n in names) {
    dir <- paste(dir_i, year, "/", state, "/apobec/", n, ".txt", sep = "")
    genes <- read.table(file = dir, sep = ' ', header = TRUE, row.names = 1)
    #names <- read.table(file = "apobec/data/tr_a.txt", sep = '\t', header = TRUE)[,2]
    #rownames(genes) <- names
    
    output <- matrix(nrow=nrow(genes), ncol=6)
    rownames(output) <- rownames(genes)
    colnames(output) <- c('mean', 'median', 'sd', 'min', 'max', 'range')
    
    for(i in rownames(genes)) {
      output[i,'mean'] <- mean(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'median'] <- median(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'sd'] <- sd(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'min'] <- min(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'max'] <- max(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'range'] <- output[i,'max'] - output[i,'min']
    }
    dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/'
    write.table(output, paste(dir_o, year, "/", state, 
                              "/expression_details_apobec/", n, ".txt",
                              sep = ""), sep='\t')
  }
  return(c(year, state, "APOBEC"))
}

###############################################################################################
# Hypoxia
###############################################################################################


fun_h <- function(year, state) {
  for(n in names) {
    dir <- paste(dir_i, year, "/", state, "/hypoxia/", n, ".txt", sep = "")
    genes <- read.table(file = dir, sep = ' ', header = TRUE, row.names = 1)
    #names <- read.table(file = "apobec/data/tr_h.txt", sep = '\t', header = TRUE)[,2]
    #rownames(genes) <- names
    
    output <- matrix(nrow=nrow(genes), ncol=6)
    rownames(output) <- rownames(genes)
    colnames(output) <- c('mean', 'median', 'sd', 'min', 'max', 'range')
    
    for(i in rownames(genes)) {
      output[i,'mean'] <- mean(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'median'] <- median(as.numeric(genes[i,]), na.rm=TRUE)
      output[i,'sd'] <- sd(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'min'] <- min(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'max'] <- max(as.numeric(genes[i,], na.rm=TRUE))
      output[i,'range'] <- output[i,'max'] - output[i,'min']
    }
    dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/'
    write.table(output, paste(dir_o, year, "/", state, 
                              "/expression_details_hypoxia/", n, ".txt",
                              sep = ""), sep='\t')
  }
  return(c(year, state, "Hypoxia"))
}

###############################################################################################
# Parallelization
###############################################################################################

call <- function(input) {
  if (input==1) {
    fun_a(year[1], state[1])
  }else if (input==2) {
    fun_h(year[1], state[1])
  }else if (input==3) {
    fun_a(year[1], state[2])
  }else if (input==4) {
    fun_h(year[1], state[2])
  }
}

if (win) {
  library(doParallel)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  if (!"OV" %in% cancers) {
    foreach(i = 1:4) %dopar% call(i)
  } else {
    foreach(i = 1:2) %dopar% call(i)
  }
  
} else {
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
