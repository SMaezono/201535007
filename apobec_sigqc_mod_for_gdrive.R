

#setwd("~/OneDrive - OnTheHub - The University of Oxford/Buffa Lab")
root <- "G:/My Drive/Jiachen Files/"
setwd(root)

win <- TRUE


# Import packages
if(win){
  library(doParallel)
}else{
  library(parallel)
}
library(sigQC)
library(RankProd)

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



states <- c("t", "n")

main <- function(i){
  require(sigQC)
  require(RankProd)
  s <- states[i]
  mRNA_expr_matrix <- list()
  if (s == "t") {
    cancers <- as.character(read.table('data/cancer_list.txt', header=T,
                                       row.names = NULL, sep='\t')[,1])
  } else {
    cancers <- as.character(read.table('data/cancer_list.txt', header=T,
                                       row.names = NULL, sep='\t')[,1])
    cancers <- cancers[-length(cancers)]
  }
 
  for (c in cancers){
    input <- read.table(file=paste("G:/My Drive/Jiachen Files/GDAC/", c, "_",
                                   s, ".txt", sep=''), sep="\t", row.names = NULL)
    # to make all rownames unique
    rownames(input) <- make.names(input$row.names, unique = TRUE)
    input <- input[,-1]
    mRNA_expr_matrix[[c]] <- as.matrix(input)
    gc()

  }
  
  gene_sigs_list <- list()
  
  signature <- "apobec"
  #gene_sig <- read.table(file="apobec/data/apobec_genes.txt", sep='\t')
  gene_sig <- read.table(file="G:/My Drive/Jiachen Files/data/tr_a.txt", sep='\t')
  # only take hugo symbols and no colname
  gene_sig <- as.character(gene_sig[,2])[-1]
  gene_sigs_list[[signature]] = as.character(gene_sig)
  
  signature <- "hypoxia"
  #gene_sig <- read.table(file="hypoxia/data/hypoxia_genes.txt", sep='\t')
  gene_sig <- read.table(file="G:/My Drive/Jiachen Files/data/tr_h 2.txt", sep='\t')
  gene_sig <- as.character(gene_sig[,2])[-1]
  gene_sigs_list[[signature]] = as.character(gene_sig)
  
  #out_dir <- "apobec/output/sigqc"
  out_dir <- paste("G:/My Drive/Jiachen Files/Sakura/output/sigqc/", s, sep='')
  gc()
  make_all_plots(gene_sigs_list = gene_sigs_list,
                 mRNA_expr_matrix = mRNA_expr_matrix,
                 doNegativeControl=FALSE, 
                 out_dir = out_dir, showResults=FALSE)
}

if(win){
  cores <- min(detectCores(), length(states))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:length(states)) %dopar% main(i) 

  stopCluster(cl)
}else{
  f <- vector("list", length=length(states))
  for (i in 1:length(states)){
    f[[i]] <- mcparallel(main(i))
  }
  mccollect(f)
}

