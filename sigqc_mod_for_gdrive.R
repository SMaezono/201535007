library(sigQC)
library(RankProd)
setwd("G:/My Drive/Jiachen Files/")

cancers <- as.character(read.table('data/cancer_list.txt', header = T,
                                   row.names  =  NULL, sep = '\t')[,1])

gene_sig <- read.table(file = "G:/My Drive/Jiachen Files/data/tr_h 2.txt",
                       sep = '\t', header = T)
hypoxia <- as.character(gene_sig[,2])[-1]
hypoxia <- unique(hypoxia)

exp <- vector('list', length = length(cancers))

for (i in 1:length(cancers)) {
  input <- read.table(paste('GDAC/', cancers[i], '_t.txt',
                            sep = ''), header = T, row.names = NULL,
                      sep = '\t')
  rn <- input[,1]
  rn <- make.unique(rn, sep = '.')
  output <- input[,-1]
  rownames(output) <- rn
  output <- output[hypoxia,]
  exp[[i]] <- as.matrix(output)
}

gene_sigs_list <- list()
gene_sigs_list[[1]] <- hypoxia
names(gene_sigs_list) <- 'Buffa52'

mRNA_expr_matrix <- exp
names(mRNA_expr_matrix) <- cancers

names_sigs <- 'Buffa52'
names_datasets <- cancers

dir.create("G:/My Drive/Jiachen Files/Sakura/sigQC/")
out_dir <- paste('G:/My Drive/Jiachen Files/Sakura/sigQC/Run', 
                 "2019", sep = "_") 


make_all_plots(gene_sigs_list, mRNA_expr_matrix, names_sigs = names_sigs,
                 names_datasets = names_datasets, out_dir = out_dir, 
                 doNegativeControl = F)

