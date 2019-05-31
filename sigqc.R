# library(sigQC)

setwd('/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/')

cancers <- as.character(read.table('data/cancer_list.txt', header=T, row.names = NULL, sep='\t')[,1])
hypoxia <- as.character(read.table('data/sign_1.txt', header=T, row.names = NULL, sep='\t')[,1])

exp <- vector('list', length=length(cancers))

for(i in 1:length(cancers)){
  input <- read.table(paste('GDAC/', cancers[i], '_t.txt', sep=''), header=T, row.names=NULL, sep='\t')
  rn <- input[,1]
  rn <- make.unique(rn, sep='.')
  output <- input[,-1]
  rownames(output) <- rn
  exp[[i]] <- output
}

gene_sigs_list <- vector('list', length=1)
gene_sigs_list[[1]] <- as.matrix(hypoxia)
names(gene_sigs_list) <- 'Buffa52'

mRNA_expr_matrix <- exp
names(mRNA_expr_matrix) <- cancers

names_sigs <- 'Buffa52'
names_datasets <- cancers

out_dir <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/sigQC/Run 2018-09-30/'

make_all_plots_2(gene_sigs_list, mRNA_expr_matrix, names_sigs = names_sigs, names_datasets = names_datasets, out_dir = out_dir, doNegativeControl = F)
