# apobec_exp_2.R
# n threads

# Jiachen Liang
# Modified by Sakura Maezono
# Display the APOBEC gene expressions

t_start <- Sys.time()

# Root directory
root <- 'G:/My Drive/Jiachen Files/'
setwd(root)

# Cancer types
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", 
             "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC", "OV")

# cancers <- "OV"

# APOBEC genes
apo <- read.table(file = "data/tr_a.txt", sep = '\t', header = TRUE)[,2]

# Expression data directory
exp_dir <- 'G:/My Drive/Jiachen Files/data/2016/tumour/apobec/'

# if folder doesn't exist
dir.create('G:/My Drive/Jiachen Files/Sakura/exp')

# Heatmap parameters
w <- 15 # width
h <- 10 # height

# Windows or Mac?
win <- TRUE

if(win) {
  library(doParallel)
} else {
  library(parallel)
}

library(matrixStats)
library(ggplot2)

###############################################################################

# Main function

main <- function(c_i){
  
  cancer <- cancers[c_i]
  
  # Import
  dir <- paste(exp_dir, cancer, '.txt', sep='')
  exp <- read.table(dir, header=T, row.names=1, sep=' ')

  
  
  # Generate expression plots
  plot <- gen_plot(exp, cancer)
  
}

###############################################################################

# Generate plots
gen_plot <- function(exp, cancer){
  require(ggplot2)
  require(EnvStats)
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
  plot <- plot + stat_n_text()
  dir_o <- 'G:/My Drive/Jiachen Files/Sakura/'
  filename <- paste(dir_o, 'exp/', cancer, '.pdf', sep='')
  pdf(file=filename, width=w, height=h)
  print(plot)
  dev.off()
  
  return(plot)
  
}

###############################################################################

# Parallelization

if(win) {
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  exp_all <- foreach(i=1:length(cancers)) %dopar% main(i)
  stopCluster(cl)
} else {
  f <- vector("list", length=length(cancers))
  for(i in 1:length(cancers)){
    f[[i]] <- mcparallel(main(i))
  }
  mccollect(f)
}



###############################################################################


print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
