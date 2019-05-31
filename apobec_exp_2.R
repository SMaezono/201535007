# apobec_exp_2.R
# n threads

# Jiachen Liang

# Display the APOBEC gene expressions

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
setwd(root)

# Cancer types
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
#cancers <- 'BRCA'

# APOBEC genes
apo <- read.table(file = "data/tr_a.txt", sep = '\t', header = TRUE)[,2]

# Expression data directory
exp_dir <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/2016/tumour/apobec/'

# Heatmap parameters
w <- 15 # width
h <- 10 # height

library(parallel)
library(matrixStats)
library(ggplot2)

###############################################################################

# Main function

main <- function(c_i){
  
  cancer <- cancers[c_i]
  
  # Import
  dir <- paste(exp_dir, cancer, '.txt', sep='')
  exp <- read.table(dir, header=T, row.names=1, sep='\t')
  
  # Generate expression plots
  plot <- gen_plot(exp, cancer)
  
}

###############################################################################

# Generate plots
gen_plot <- function(exp, cancer){
  
  mrna <- c()
  gene <- c()
  
  for(i in 1:nrow(exp)){
    for(j in 1:ncol(exp)){
      mrna <- c(mrna, exp[i,j])
      gene <- c(gene, rownames(exp)[i])
    }
  }
  
  df <- data.frame(mrna, gene)
  
  plot <- ggplot(df, aes(x=gene, y=mrna)) +
    geom_boxplot() +
    ggtitle(paste("Expression of APOBEC genes in", cancer)) +
    xlab("APOBEC gene") + ylab("mRNA counts (log2)")
  
  filename <- paste('exp/', cancer, '.pdf', sep='')
  pdf(file=filename, width=w, height=h)
  print(plot)
  dev.off()
  
  return(plot)
  
}

###############################################################################

# Parallelization

f <- vector("list", length=length(cancers))
for(i in 1:length(cancers)){
  f[[i]] <- mcparallel(main(i))
}
mccollect(f)

###############################################################################


print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
