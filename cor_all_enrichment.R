# cor_all_enrichment.R

# Jiachen Liang

# Make circle plots with enrichment data

# Work in progress: integrate expression data
#   Not all data is in gene symbols
#   Must find a way to dynamically convert all to a standard

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/cor_all/enrichment/"
setwd(root)

# Cancer types
cancerc <- 'BRCA'

# Analysis groups
groups <- c('_5', '_A3B', '_F_5', '_F_A3B', '_Adj_5', '_Adj_F_5', '_hyp')

# Limit of pathways to show
path_lim <- 10

# Boxplot dimensions
w <- 15 # width
h <- 17 # height

# Margins
par(mar=c(5,5,5,5))
par(oma=c(5,5,5,5))
par(xpd=TRUE, cex=0.8)

library(parallel)
library(matrixStats)
library(GOplot)

###############################################################################

# Import function

import <- function(cancer, group){
  
  # Read file
  dir <- paste(cancer, group, '.txt', sep='')
  input <- read.table(dir, sep="\t", header=TRUE)
  
  # Filter
  input <- input[(input$Support >= 3),]
  if(nrow(input) > path_lim){
    input <- input[1:path_lim,]
  }
  
  return(input)
}

###############################################################################

# Import expression

import_exp <- function(cancer){
  
  # Read file
  dir <- paste(cancer, '_exp.txt', sep='')
  exp <- read.table(file=dir, sep='\t', header=TRUE, row.names=1)
  
  return(exp)
}

###############################################################################

# Matrix generation functions

# Without expression data
matgen <- function(input){
  
  # List of all genes
  genelist <- c()
  
  # Matrix
  mat <- matrix(ncol=nrow(input), nrow=0)
  colnames(mat) <- input$Items_Details
  
  # Cycle through pathways
  for(i in 1:nrow(input)){
    
    # Get genes in pathway
    genes <- strsplit(as.character(input[i,'Genes']), ',')
    
    for(g in genes[[1]]){
      if(g %in% genelist){
        mat[g,as.character(input[i,"Items_Details"])] <- 1
      }else{
        genelist <- c(genelist, g)
        mat <- rbind(mat, matrix(nrow=1, ncol=nrow(input)))
        rownames(mat) <- genelist
        mat[g,] <- 0
        mat[g,as.character(input[i,"Items_Details"])] <- 1
      }
    }
  }
  
  return(mat)
}

# With expression data
matgen_exp <- function(input, exp){
  
  # Same analysis as above
  no_exp <- as.data.frame(matgen(input))
  
  # Import expression data
  genes <- rownames(no_exp)
  logFC <- as.numeric(rowMeans(exp[genes,], na.rm=T))
  logFC <- logFC-med
  no_exp <- cbind(no_exp, logFC)
  
  return(no_exp)
  
}

###############################################################################

# Plot functions

# Without expression data
plot_out <- function(m, cancer, group){
  if(nrow(m)>70){gs<-2.5;gp<-0.145}else{gs<-3.5;gp<-0.17}
  out <- GOChord(m, nlfc=0, gene.size=gs, gene.space=gp)
  filename <- paste('plots/', cancer, group, '.pdf', sep='')
  pdf(file=filename, width=w, height=h)
  print(out)
  dev.off()
  write.table(m, file=paste('plots/', cancer, group, '.txt', sep=''), sep='\t')
  return(out)
}

# With expressiond data
plot_exp <- function(m, cancer, group){
  if(nrow(m)>70){gs<-2.5;gp<-0.145}else{gs<-3.5;gp<-0.17}
  out <- GOChord(m, gene.order = 'logFC', gene.size=gs, gene.space=gp)
  filename <- paste('plots/', cancer, group, '_exp.pdf', sep='')
  pdf(file=filename, width=w, height=h)
  print(out)
  dev.off()
  write.table(m, file=paste('plots/', cancer, group, '_exp.txt', sep=''), sep='\t')
  return(out)
}

###############################################################################

# Main function

# List of plots to output
plots <- vector("list", length=length(groups))

# Calculate
for(i in 1:length(groups)){
  group <- groups[i]
  raw_data <- import(cancer, group)
  exp <- import_exp(cancer)
  med <- median(as.matrix(exp))
  m <- matgen_exp(raw_data, exp)
  p <- plot_exp(m, cancer, group)
  plots[[i]] <- p
}

###############################################################################


print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)