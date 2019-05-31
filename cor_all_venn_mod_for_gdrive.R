# cor_all_venn.R

# Generate Venn diagrams for genes and pathways

t_start <- Sys.time()

# Root directory
root <- "G:/My Drive/Jiachen Files/Sakura/output/cor_all/enrichment/"
setwd(root)

# Cancer types
cancers <- c('BRCA')

# Analysis groups
groups <- c('_5', '_A3B', '_hyp')

library(matrixStats)
library(GOplot)
library(VennDiagram)

for(c in cancers){
  
  cancer <- c
  # Import
  mat_5 <- read.table(file=paste('plots/', cancer, '_F_5_exp.txt',
                                 sep=''), sep='\t', row.names=1, header=T)
  mat_a3b <- read.table(file=paste('plots/', cancer, '_F_A3B_exp.txt', 
                                   sep=''), sep='\t', row.names=1, header=T)
  mat_hyp <- read.table(file=paste('plots/', cancer, '_hyp_exp.txt', 
                                   sep=''), sep='\t', row.names=1, header=T)
  # real_hyp <- as.character(read.table(file='/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/sign_1.txt', header=T, row.names=NULL)[,1])
  
  # Remove logFC column
  mat_5 <- mat_5[,1:ncol(mat_5)-1]
  mat_a3b <- mat_a3b[,1:ncol(mat_a3b)-1]
  mat_hyp <- mat_hyp[,1:ncol(mat_hyp)-1]
  
  # Venn diagram for genes
  mat_5_v_mat3b <- intersect(rownames(mat_5), rownames(mat_a3b))
  write.table(mat_5_v_mat3b, "mat_5_v_mat3b_genes.txt")
  mat_5_v_hypoxia <- intersect(rownames(mat_5), rownames(mat_hyp))
  write.table(mat_5_v_hypoxia, "mat_5_v_hypoxia_genes.txt")
  mat3b_v_hypoxia <- intersect(rownames(mat_a3b), rownames(mat_hyp))
  write.table(mat3b_v_hypoxia,"mat3b_v_hypoxia_genes.txt")
  list_exp <- list(
    "Any 5 APOBEC" = rownames(mat_5),
    APOBEC3B = rownames(mat_a3b),
    Hypoxia = rownames(mat_hyp)
  )
  venn_exp <- venn.diagram(list_exp, height = 4000, width = 4000,
                           paste('plots/', cancer, '_venn_exp.tiff', sep=''),
                           main=paste('Overlap of genes -', cancer),
                           fill=c('cadetblue1', 'darksalmon', 'darkolivegreen1'))
  mat_5_v_mat3b <- intersect(colnames(mat_5), colnames(mat_a3b))
  write.table(mat_5_v_mat3b, "mat_5_v_mat3b_pathway.txt")
  mat_5_v_hypoxia <- intersect(colnames(mat_5), colnames(mat_hyp))
  write.table(mat_5_v_hypoxia, "mat_5_v_hypoxia_pathway.txt")
  mat3b_v_hypoxia <- intersect(colnames(mat_a3b), colnames(mat_hyp))
  write.table(mat3b_v_hypoxia, "mat3b_v_hypoxia_pathway.txt")
  
  # Venn diagram for pathways
  list_path <- list(
    "Any 5 APOBEC" = colnames(mat_5),
    APOBEC3B = colnames(mat_a3b),
    Hypoxia = colnames(mat_hyp)
  )
  venn_path <- venn.diagram(list_path, height = 4000, width = 4000,
                            paste('plots/', cancer, '_venn_path.tiff', sep=''),
                            main=paste('Overlap of pathways -', cancer),
                            fill=c('cadetblue1', 'darksalmon', 'darkolivegreen1'))
  
}
