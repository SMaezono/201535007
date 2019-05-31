

#setwd("~/OneDrive - OnTheHub - The University of Oxford/Buffa Lab")
setwd("C:\\Users\\Jiachen Liang\\OneDrive - OnTheHub - The University of Oxford")

win <- TRUE


# Import packages
if(win){
  library(doParallel)
}else{
  library(parallel)
}
library(sigQC)

# names = c("dataset1")
# data.matrix = replicate(10, rnorm(20))#random matrix - 10 genes x 20 samples
# mRNA_expr_matrix = list()
# mRNA_expr_matrix[["dataset1"]] = data.matrix
# row.names(mRNA_expr_matrix$dataset1) <- as.character(1:(dim(mRNA_expr_matrix$dataset1)[1]))
# colnames(mRNA_expr_matrix$dataset1) <- as.character(1:(dim(mRNA_expr_matrix$dataset1)[2]))
# 
# #Define the signature
# gene_sigs_list = list()
# signature = "hypoxiaSig"
# gene_sig = c('1', '4', '5')#gene ids
# gene_sigs_list[[signature]] = as.matrix(gene_sig)
# names_sigs = c(signature)


cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", 
             "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
states <- c("t", "n")

main <- function(i){
  library(sigQC)
  s <- states[i]
  mRNA_expr_matrix <- list()
  for(c in cancers){
    input <- read.table(file=paste("C:/Users/Jiachen Liang/Documents/GDAC/", c, "_", s, ".txt", sep=''), sep="\t")
    mRNA_expr_matrix[[c]] <- input
  }
  
  gene_sigs_list <- list()
  
  signature <- "apobec"
  #gene_sig <- read.table(file="apobec/data/apobec_genes.txt", sep='\t')
  gene_sig <- read.table(file="C:\\Users\\Jiachen Liang\\OneDrive - OnTheHub - The University of Oxford\\Buffa Lab\\apobec\\data\\apobec_genes.txt", sep='\t')
  gene_sigs_list[[signature]] = as.matrix(gene_sig)
  
  signature <- "hypoxia"
  #gene_sig <- read.table(file="hypoxia/data/hypoxia_genes.txt", sep='\t')
  gene_sig <- read.table(file="C:\\Users\\Jiachen Liang\\OneDrive - OnTheHub - The University of Oxford\\Buffa Lab\\apobec\\data\\hypoxia_genes.txt", sep='\t')
  gene_sigs_list[[signature]] = as.matrix(gene_sig)
  
  #out_dir <- "apobec/output/sigqc"
  out_dir <- paste("C:\\Users\\Jiachen Liang\\OneDrive - OnTheHub - The University of Oxford\\Buffa Lab\\apobec\\output\\sigqc\\", s, sep='')
  
  make_all_plots(gene_sigs_list = gene_sigs_list, mRNA_expr_matrix = mRNA_expr_matrix, doNegativeControl=FALSE, out_dir = out_dir, showResults=FALSE)
}

if(win){
  cores <- min(detectCores(), length(states))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:length(states)) %dopar% main(i)
  stopCluster(cl)
}else{
  f <- vector("list", length=length(states))
  for(i in 1:length(states)){
    f[[i]] <- mcparallel(main(i))
  }
  mccollect(f)
}

