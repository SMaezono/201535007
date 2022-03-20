# Sakura Maezono # Date created: 
# Purpose: 
# stratify samples into TP53 quadrants
# make a violin plot with x-axis as the TP53 quadrants
# and y-axis is APOBEC expression


###############################################################################################
# Associated library functions
###############################################################################################
# Windows or Mac
win <- TRUE

library(matrixStats)
if (win) {
  library(doParallel)
}else{
  library(parallel)
}
library(gplots)
library(ggpubr)
library(EnvStats)

# create folder(s) if it doesn't/ they don't exist

dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/")

# Set working directory
root <- "G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/"
setwd(root)

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'


# APOBEC and TP53 genes
dir_ext <- "G:/My Drive/Jiachen Files/"
apo.genes <- as.character(read.table(file = paste(dir_ext, 
                                                  "data/tr_a.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
# hyp.genes <- as.character(read.table(file = paste(dir_ext,
#                                                   "data/tr_h.txt", sep = ""),
#                                      sep = "\t", header = TRUE, 
#                                      row.names = 1)$Gene.Symbol)

p21.genes <- "TP53"

# Minimum proportion of non-NA samples
thr_na <- 0.8

# "t" for tumour data and "n" for normal data
state <- "t"

dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/Buffa 52/")
dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/APOBEC/")
if (state == "t") {

  dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/Buffa 52/tumour/")
  dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/APOBEC/tumour/")
} else {
  dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/Buffa 52/normal/")
  dir.create("G:/My Drive/Jiachen Files/Sakura/TP53_APOBEC_Stratification/APOBEC/normal/")
}

# Restrict to A3 family?
restrict <- TRUE

# Cancer list
cancers <- as.character(read.table(paste(dir_p, 
                                         'cancer_list.txt',
                                         sep = ''),
                                   header = T, sep = '\t')[,1])

if (state == "n") {
  # "OV" has no normal tissues
  cancers <- cancers[-length(cancers)]
}



###############################################################################
# Import expression matrix from file system
###############################################################################

import_cancer_exp <- function(cancer, state) {
  # cancer = c("BRCA", etc.)
  # state = c("t", "n") 
  
  dir <- paste(dir_e, cancer, '_', state, '.txt', sep = '')
  input <- read.table(dir, sep = '\t', header = T, row.names = NULL)
  
  rn <- input[,1]
  rn <- make.unique(rn, sep = '.')
  
  output <- as.data.frame(input[,-1])
  rownames(output) <- rn
  
  exp <- input[,-1]
  rownames(exp) <- rn
  
  # Fix col names
  names <- colnames(exp)
  for (i in 1:length(names)) {
    split <- strsplit(names[i], '\\.')
    names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
  }
  colnames(exp) <- names
  # output <- data.frame(output)
  # TP53 genes only in order
  p21_t <- as.matrix(output[p21.genes,])
  # # remove rows that only have NAs
  # p21_t <- rm_onlyNArow(p21_t)
  na_filter <-
    (rowCounts(!is.na(p21_t))/ncol(p21_t) >= thr_na)
  p21_t <- p21_t[na_filter,]
  
  # # Replace NA by median of gene
  # med_t <- rowMedians(p21_t, na.rm = T)
  # for (i in 1:nrow(p21_t)) {
  #   for (j in 1:ncol(p21_t)) {
  #     if (is.na(p21_t[i,j])) {
  #       p21_t[i,j] <- med_t[i]
  #     }
  #   }
  # }
  
  input_t <- as.matrix(output[apo.genes,])
  if (restrict) {
    input_t <- input_t[c('APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
                         'APOBEC3F', 'APOBEC3G', 'APOBEC3H'),]
  }
  # input_t <- data.frame(output[order(rownames(input_t)),])
  # remove rows that only have NAs
  # input_t <- rm_onlyNArow(input_t)
  # Filter out genes non-NA in less than x proportion of samples
  na_filter <-  
    (rowCounts(!is.na(input_t))/ncol(input_t)>=thr_na)
  
  input_t <- input_t[na_filter,]
  
  # # Replace NA by median of gene
  # med_t <- rowMedians(input_t, na.rm = T)
  # for(i in 1:nrow(input_t)) {
  #   for(j in 1:ncol(input_t)) {
  #     if(is.na(input_t[i,j])) {
  #       input_t[i,j] <- med_t[i]
  #     }
  #   }
  # }
  
  # matrices of mean & median scores for every patient
  # TP53 gene sig
  hyp_metrics <- as.matrix(p21_t)
  colnames(hyp_metrics) <- p21.genes


  # APOBEC gene sig
  input_t <- t(input_t)
  apo_metrics <- matrix(nrow = length(input_t[,1]), ncol = 2 )
  rownames(apo_metrics) <- rownames(input_t)
  colnames(apo_metrics) <- c("Mean", "Median")
  apo_metrics[,1] <- rowMeans(as.matrix(input_t), na.rm=T)
  apo_metrics[,2] <- rowMedians(as.matrix(input_t), na.rm=T)
  # reorder by mean
  apo_metrics <- apo_metrics[order(apo_metrics[,1]),]
  
  # Output_all, mean and median of gene sig dor each patient
  write.table(hyp_metrics, 
              file = paste(cancer, state,
                           'exp_TP53_genes.txt', 
                           sep = '_'), sep = '\t')
  
  # Output_all, mean and median of gene sig dor each patient
  write.table(apo_metrics, 
              file = paste(cancer, state,
                           '_meanormedian_of_APOBEC_genes.txt', 
                           sep = '_'), sep = '\t')
  # reorder based on median lowest to highest
  hyp_metrics_med <- data.frame(hyp_metrics[order(hyp_metrics[,1]),])
  colnames(hyp_metrics_med) <- "TP53"
  apo_metrics_med <- apo_metrics[order(apo_metrics[,2]),]
  
  cancer_output <- list(hyp_metrics_med, 
                        hyp_metrics_med,
                        data.frame(apo_metrics[,1]),
                        data.frame(apo_metrics_med[,2]),
                        data.frame(p21_t), data.frame(input_t))
  return(cancer_output)
}

Strat_APOBEC_TP53_plot <- function(cancer) {
  # stratify into quartiles (x-axis of the violin plot)
  # graph using a violin plot
  cancer_output <- import_cancer_exp(cancer, state)
  sep_by_quartiles <- function(which_sig) {
    # stratify into quartiles 
    # which_sig becomes x-axis of the violin plot 
    if (which_sig == "TP53") {
      mean <- cbind(rownames(cancer_output[[1]]), cancer_output[[1]])
      median <- cbind(rownames(cancer_output[[2]]), cancer_output[[2]])
      
    } else if (which_sig == "apobec") {
      mean <- cbind(rownames(cancer_output[[3]]), cancer_output[[3]])
      median <- cbind(rownames(cancer_output[[4]]), cancer_output[[4]])
    }
    
    # round() rounds OFF the values by the number of digit desired
    quartile_n <- ceiling(length(mean[,1])*0.25)
    # prints out how many samples were considered
    print(paste(quartile_n, 
                "sample ids were used for the top, upper-middle, lower-middle,",
                "or bottom group based on", which_sig, "gene sig score", "in", cancer,
                sep = " "))
    top_mean <- data.frame(mean[order(mean[,1], decreasing = TRUE),])
    
    top_mean_25 <- data.frame(top_mean[1:quartile_n,]) 
    top_mean_patients <- as.character(top_mean_25[,1])
    
    top_median <- data.frame(median[order(median[,1], 
                                          decreasing = TRUE),])
    top_median_25 <- data.frame(top_median[1:quartile_n,]) 
    top_median_patients <- as.character(top_median_25[,1])
    
    upper_mid_mean <- data.frame(top_mean[c((quartile_n+1):(quartile_n*2)),]) 
    upper_mid_mean_patients <- as.character(upper_mid_mean[,1])
    
    upper_mid_median <- data.frame(top_median[c((quartile_n+1):(quartile_n*2)),]) 
    upper_mid_median_patients <- as.character(upper_mid_median[,1])

    
    bottom_mean <- data.frame(mean[order(mean[,1],
                                      decreasing = FALSE),])
    bottom_mean_25 <- data.frame(bottom_mean[1:quartile_n,]) 
    bottom_mean_patients <- as.character(bottom_mean_25[,1])
    
    bottom_median <- data.frame(median[order(median[,1], 
                                          decreasing = FALSE),])
    bottom_median_25 <- data.frame(bottom_median[1:quartile_n,]) 
    bottom_median_patients <- as.character(bottom_median_25[,1])
    
    lower_mid_mean <- data.frame(bottom_mean[c((quartile_n+1):(quartile_n*2)),]) 
    lower_mid_mean_patients <- as.character(lower_mid_mean[,1])
    
    lower_mid_median <- data.frame(bottom_median[c((quartile_n+1):(quartile_n*2)),]) 
    lower_mid_median_patients <- as.character(lower_mid_median[,1])
    
    # take the values of the other gene signature for their gene expressions
    # this should be the values for y-axis
    
    if (which_sig == "TP53") {
      exp_mean <- cbind(rownames(cancer_output[[3]]), cancer_output[[3]])
      exp_median <- cbind(rownames(cancer_output[[4]]), cancer_output[[4]])
      exp_indiv_gene <- cbind(rownames(cancer_output[[6]]), cancer_output[[6]])
    } else if (which_sig == "apobec") {
      exp_mean <- cbind(rownames(cancer_output[[1]]), cancer_output[[1]])
      exp_median <- cbind(rownames(cancer_output[[2]]), cancer_output[[2]])
      exp_indiv_gene <- cbind(rownames(cancer_output[[5]]), cancer_output[[5]])
    }
    
    inner_fxn <- function(which_metric) {
      if (which_metric == "mean") {
        exp_df <- exp_mean
        top_patients <- top_mean_patients
        upper_mid_patients <- upper_mid_mean_patients
        lower_mid_patients <- lower_mid_mean_patients
        bottom_patients <- bottom_mean_patients
      } else {
        exp_df <- exp_median
        top_patients <- top_median_patients
        upper_mid_patients <- upper_mid_median_patients
        lower_mid_patients <- lower_mid_median_patients
        bottom_patients <- bottom_median_patients
      }
      # metric score of gene sig
      top_exp <- exp_df[top_patients,]
      upper_mid_exp <- exp_df[upper_mid_patients,]
      bottom_exp <- exp_df[bottom_patients,]
      lower_mid_exp <- exp_df[lower_mid_patients,]
      
      # individual score of individual genes of the signature
      top_exp_ind <- exp_indiv_gene[top_patients,]
      upper_mid_exp_ind  <- exp_indiv_gene[upper_mid_patients,]
      bottom_exp_ind  <- exp_indiv_gene[bottom_patients,]
      lower_mid_exp_ind  <- exp_indiv_gene[lower_mid_patients,]
      
      top_exp <- cbind(top_exp, top_exp_ind[,-1])
      upper_mid_exp <- cbind(upper_mid_exp, upper_mid_exp_ind[,-1])
      bottom_exp <- cbind(bottom_exp, bottom_exp_ind[,-1])
      lower_mid_exp <- cbind(lower_mid_exp,lower_mid_exp_ind[,-1])
      
      
      prep_for_graph_fxn <- function(df_exp, which_quartile) {
        # make the colnames the same for combining all of them later
        # add quartile information for violin plot
        require(reshape2)
        colnames(df_exp)[c(1,2)] <- c("Patients", which_metric)
        df_exp <- melt(df_exp, id.vars = c("Patients"))
        df_exp$Quartile <- rep(which_quartile, length.out = length(df_exp[,1]))
        return(df_exp)
      }
      
      top_exp <- prep_for_graph_fxn(top_exp, "top")
      upper_mid_exp <- prep_for_graph_fxn(upper_mid_exp, "upper-mid")
      bottom_exp <- prep_for_graph_fxn(bottom_exp, "bottom")
      lower_mid_exp <- prep_for_graph_fxn(lower_mid_exp, "lower-mid")
      
      final_df_graph <- rbind(bottom_exp, lower_mid_exp, 
                              upper_mid_exp, top_exp)
      if (which_sig == "TP53") {
        folder_name <- "Buffa 52"
      } else if (which_sig == "apobec") {
        folder_name <- "APOBEC"
      }
      if (state == "t") {
        write.table(final_df_graph, 
                    file = paste(folder_name, "/tumour/",
                                 cancer, "_", which_metric, "_", 'stratified_to_quartiles.txt', 
                                 sep = ''), sep = '\t')
      } else {
          write.table(final_df_graph, 
                      file = paste(folder_name, "/normal/",
                                   cancer, "_", which_metric, "_", 'stratified_to_quartiles.txt', 
                                   sep = ''), sep = '\t')
      }

      return(final_df_graph)
    }
    listofdfs <- lapply(c("mean", "median"), inner_fxn)
    return(listofdfs)
  }
  
  TP53_quartiles <- sep_by_quartiles("TP53")
  
  apobec_quartiles <- sep_by_quartiles("apobec")
  
  violin_plot <- function(which_gene, filename, x_label, quadrant_metric) {
    # Violin plots with box plots inside
    # which_gene == c("mean", "median", "APOBEC3A", etc)
    # apobec_quartiles[[1]], apobec_quartiles[[2]])
    # filename <- c("TP53_quartile_mean_apobec_exp",
    # "TP53_quartile_median_apobec_exp", 
    # "apobec_quartile_mean_TP53_exp",
    # "apobec_quartile_median_TP53_exp")
    
    my_comparisons <- list(c("top", "upper-mid"),c("upper-mid", "lower-mid"),
                           c("top", "lower-mid"), c("lower-mid", "bottom"), 
                           c("upper-mid", "bottom"), c("top", "bottom"))
    if (x_label == "TP53" & quadrant_metric == "mean") {
      df <- TP53_quartiles[[1]][TP53_quartiles[[1]]$variable == which_gene,]
    } else if (x_label == "TP53" & quadrant_metric == "median") {
      df <- TP53_quartiles[[2]][TP53_quartiles[[2]]$variable == which_gene,]
    } else if (x_label == "APOBECs" & quadrant_metric == "mean") {
      df <- apobec_quartiles[[1]][apobec_quartiles[[1]]$variable == which_gene,]
    } else if (x_label == "APOBECs" & quadrant_metric == "median") {
      df <- apobec_quartiles[[2]][apobec_quartiles[[2]]$variable == which_gene,]
    }
    
    if (x_label == "TP53") {
      y_label_new <- paste( "APOBEC (", which_gene, ")", sep = "")
    } else {
        y_label_new <- paste( "TP53 (", which_gene, ")", sep = "")
    }
    
    plot <- ggviolin(df, x = "Quartile", y = "value", fill = "Quartile",
                     xlab = paste("Quartile of", x_label, sep = " "),
                     ylab = paste("mRNA expression of", y_label_new, sep = " "),
                     add = "boxplot", add.params = list(fill = "white")) +
      # Add significance levels
      stat_compare_means(comparisons = my_comparisons, label = "p.format") + 
      # Add global the p-value 
      stat_compare_means(label.y.npc = 1) 
    plot <- plot + stat_n_text() 
    pdf(file = paste(cancer, "_", state, "_", filename, "_", which_gene, ".pdf", sep = ""),
        width = 6, height = 6)
    
    print(plot)
    dev.off()
  }
  
  # TP53 quartiles, mean 
  violin_plot("mean", "TP53_quartile_mean_apobec_exp", 
              "TP53", "mean")
  # TP53 quartiles, median 
  violin_plot("median", "TP53_quartile_mean_apobec_exp",
                         "TP53", "median")
  # apobec quartiles, mean 
  violin_plot("mean", "apobec_quartile_mean_TP53_exp",
              "APOBECs", "mean")
  # apobec quartiles, median 
  violin_plot("median", "apobec_quartile_mean_TP53_exp", 
              "APOBECs", "median")
  
  plot_ind_hyp_genes <- function(which_gene) {
    violin_plot(which_gene, "TP53_quartile_mean_apobec_exp", 
                "TP53", "mean")
    violin_plot(which_gene, "TP53_quartile_mean_apobec_exp", 
                "TP53", "median")
    print(paste("plotted violin graphs for", which_gene, sep = " "))
  }
  lapply(colnames(cancer_output[[6]]), plot_ind_hyp_genes)
  
  plot_ind_APOBEC_genes <- function(which_gene) {
    violin_plot(which_gene, "apobec_quartile_mean_TP53_exp",
                "APOBECs", "mean")
    violin_plot(which_gene, "apobec_quartile_median_TP53_exp",
                "APOBECs", "median")
    print(paste("plotted violin graphs for", which_gene, sep = " "))
  }
  lapply(colnames(cancer_output[[5]]), plot_ind_APOBEC_genes)
}

lapply(cancers, Strat_APOBEC_TP53_plot)


root <- "G:/My Drive/Jiachen Files/codes_supp/"
setwd(root)
gc()

# save details of this R script
scriptname <- "strat_into_TP53_quadrants_APOBEC_exp"
savehistory(file = paste(scriptname, ".Rhistory", sep = ""))
sink(file = paste(Sys.Date(), scriptname, ".txt", sep = ""),
     type = c("output", "message"))
# details of the packages and other useful info.
print(sessionInfo()) 
sink(NULL)

# save globalenv/workspace
save.image("strat_into_TP53_quadrants_APOBEC_exp.Rdata")  


  
