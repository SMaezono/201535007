# Sakura Maezono # Started: ##------ Mon Apr 01 19:48:29 2019 ------##
# Purpose:
# plot percentage of non-null per signature
# plot 
# IMPORTANT: Only run after mutsig2_mod_for_gdrive2.R


library(tidyverse)
library(doParallel)

cancers <- read.table(file ='G:/My Drive/Jiachen Files/mutsig/cancer_list.txt',
                      sep = '\t', header=F)[,1]
# remove OV for now
cancers <- as.character(cancers)[-9]



nnp_main_fxn <- function(cancer, directory) {
  # directory <- c("output_dcs", "output_dcs_curated")
  percent_non_null <- read.table(file=paste(directory, "/",
                                            cancers[cancer], '_nnp.txt',
                                              sep = ""), sep = '\t')
  colnames(percent_non_null) <- "values"
  # rearrange x by highest to lowest
  require(ggplot2)
  plot_t <- ggplot(data=percent_non_null, aes(x=rownames(percent_non_null),
                                              y=values)) +
    ylim(0, 100) +
    geom_bar(stat = 'identity') +
    ylab("% non-null samples") +
    xlab("") +
    ggtitle(cancers[cancer]) +
    theme(axis.text.x = element_text(angle = 90))
    
  pdf(file=paste(directory, "/", cancers[cancer], '_nnp.pdf', sep = ''),
      width=6, height=4)
  print(plot_t)
  dev.off()
}
  
control_main_fxn <- function(directory) {
  # barplot of control per cancer
  # directory <- c("z_counts", "z_enrichment","z_dcs_all", "z_dcs_cur" )
  require(reshape2)
  z_score <- read.table(file=paste("output_controls/", directory,
                                   '.txt', sep = ""), sep = '\t')
  z_score <- cbind(rownames(z_score), z_score[,-3])
  
  z_score <- melt(z_score)
  colnames(z_score)[c(1:2)] <- c("Cancers", "gene")
  
  # rearrange x by highest to lowest
  require(ggplot2)
  new_filename <- unlist(strsplit(directory, "_"))[2]
  plot_t <- ggplot(data = z_score, aes(x = Cancers, y = value,
                                     fill = gene)) +
    geom_bar(stat = 'identity', position = position_dodge2()) +
    ylab("Z score") +
    ggtitle(paste("Control", new_filename, sep = " "))
  
  pdf(file=paste("output_controls/",new_filename, '.pdf', sep = ""),
      width=6, height=4)
  print(plot_t)
  dev.off()
}

### Main #######################
root <- "G:/My Drive/Jiachen Files/mutsig/"
setwd(root) 

cores <- min(detectCores(), length(cancers))
cl <- makeCluster(cores)
registerDoParallel(cl)
main <- foreach(i = 1:length(cancers)) %dopar% nnp_main_fxn(i, "output_dcs")
stopCluster(cl)
gc()

cores <- min(detectCores(), length(cancers))
cl <- makeCluster(cores)
registerDoParallel(cl)
# remove CESC for now
cancers <- cancers[-3]
main <- foreach(i = 1:length(cancers)) %dopar% nnp_main_fxn(i, "output_dcs_curated")
stopCluster(cl)
gc()


lapply(c("z_counts", "z_enrichment"),control_main_fxn)


