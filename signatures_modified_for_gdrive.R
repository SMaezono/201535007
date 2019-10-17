  # signatures.R
  
  # Jiachen Liang
  # 2018
  # modified by Sakura Maezono for gdrive ##------ Tue Feb 26 14:16:59 2019 ------##
  # created subfolder directories 
  # added p-values and sample size for boxplots   
  # switched colors for tumor and normal
  # for big gene signatures (>20), separated the results to several pdf 
  
  # Buffa Laboratory
  # Computational Biology and Integrative Genomics Group
  # Department of Oncology
  # University of Oxford
  
  ###############################################################################################
  
  # Master function to compare a signature with one or many others
  
  # Inputs:
  #   List of cancers
  #   Master signature (sign. 1)
  #   Tested signatures (sign. 2 to n-1)
  #   Housekeeping genes (sign. n)
  #   Expression matrix (normal)
  #   Expression matrix (tumour)
  
  # Outputs:
  #   In each cancer
  #     Boxplot of master signature expression (n/t)
  #       + t-test results
  #     Boxplots of tested signatures expression (n/t)
  #       + t-test results
  #     Heatmap of correlation between all genes in master signature vs all genes
  #     in all tested signatures
  #       + associated rho and p-value tables
  #   Boxplot of mean master signature expression across cancers (n/t)
  #     + t-test results
  #   violin plot of mean master signature expression across cancers (n/t) [NEW 2019/07]
  #   Boxplots of mean tested signatures expression across cancers (n/t)
  #     t-test results
  #   Heatmap of expression of all genes across all signautres, in all samples
  #   Heatmap of expression of all genes across all signatures, in all samples, normalized to housekeeping genes
  #   Heatmap of correlation between the mean expression of master signautre in all samples vs each gene in all tested signatures
  #     + associated rho and p-value tables
  
  #   All heatmaps are set as red for -1 and blue for 1 
  ###############################################################################################
  # Initiation
  ###############################################################################################
  
  t_start <- Sys.time()
  t_last <- Sys.time()
  
  # Using Windows?
  win <- T
  
  
  
  ###############################################################################################
  # Associated library functions
  ###############################################################################################
  
  if (win) {
    library(doParallel)
  } else {
    library(parallel)
  }
  
  # Libraries used
  library(matrixStats)
  library(tidyverse)
  library(gplots)
  # library(ggplot2)
  library(ggpubr) 
  library(ggsignif)
  library(ggpmisc)
  library(EnvStats)
  
  
  ###############################################################################################
  # File locations
  ###############################################################################################
  
  if (win) {
    # Set working directory
    root <- "G:/My Drive/Jiachen Files/data"
    setwd(root)
    
    # Directory: output
    dir.create('G:/My Drive/Jiachen Files/Sakura/sig_output/')
    dir_o <- 'G:/My Drive/Jiachen Files/Sakura/sig_output/'
    # Directory: parameters
    dir_p <- 'G:/My Drive/Jiachen Files/data/'
    
    # Directory: expression matrices
    dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'
  } else {
    # Set working directory
    root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/"
    setwd(root)
    
    # Directory: output
    dir_o <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/signatures_output/'
    
    # Directory: parameters
    dir_p <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/'
    
    # Directory: expression matrices
    dir_e <- '/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/GDAC/'
  }
  
  # create directories [IMPORTANT, files cannot be saved to non-existing
  # directory]
  dir.create(paste(dir_o, 'exp boxplots per signature/', sep = ""))
  dir.create(paste(dir_o, 'exp violin plots per signature/', sep = ""))
  dir.create(paste(dir_o, 'exp boxplots across cancers/', sep = ""))
  dir.create(paste(dir_o, 'exp violin plots across cancers/', sep = ""))
  dir.create(paste(dir_o, 'exp heatmaps/', sep = ""))
  dir.create(paste(dir_o, 'exp heatmaps/norm_hk/', sep = ""))
  dir.create(paste(dir_o, 'cor single gene/', sep = ""))
  dir.create(paste(dir_o, 'cor mean master/', sep = ""))
  dir.create(paste(dir_o, 'cor mean master scatterplots/', sep = ""))
  dir.create(paste(dir_o, 'cor mean master scatterplots/_NT', sep = ""))
  dir.create(paste(dir_o, 'cor mean master/clustered/', sep = ""))
  dir.create(paste(dir_o, 'xCell heatmaps/', sep = ""))
  dir.create(paste(dir_o, 'xCell cor/', sep = ""))
  dir.create(paste(dir_o, 'xCell cor/unclustered/', sep = ""))
  dir.create(paste(dir_o, 'xCell cor all/', sep = ""))
  dir.create(paste(dir_o, 'xCell cor all/unclustered/', sep = ""))
  dir.create(paste(dir_o, 'heatmaps master individual genes/', sep = ""))
  dir.create(paste(dir_o, 
                   'heatmaps master individual genes/signature clustered/',
                   sep = ""))
  dir.create(paste(dir_o, 
                   'heatmaps master individual genes/full clustered/',
                   sep = ""))
  
  ###############################################################################################
  # Parameters
  ###############################################################################################
  
  # Number of tested signatures (default: 4)
  nts <- 21
  
  
  # Box/violin plot dimensions
  w <- 20 # width default: 8
  h <- 10 # height default: 6
  
  # Heatmap parameters
  # heatmap_colors <- colorRampPalette(c("red", "white", "blue"))(n = 100)
  margins <- c(5,5) # margins for output image
  key <- FALSE # display legend # figure margins too large if TRUE
  skey <- TRUE # symmetric key
  keysize <- 1
  fontsize <- 0.8
  hm_w <- 16 # width
  hm_h <- 14 # height
  
  # Significance threshold
  thr <- 0.05
  
  # Using means (T) vs medians(F)
  um <- T
  
  # Minimum proportion of non-NA samples
  thr_na <- 0.80
  
  # remove rows that are only NA
  rm_onlyNArow <- function(df) df[rowSums(is.na(df)) != ncol(df),]
  
  rm_onlyNAcol <- function(df) df[,colSums(is.na(df)) != nrow(df)]
  
  
  
  ###############################################################################################
  # Import expression matrix from file system
  ###############################################################################################
  
  exp_import <- function(cancer, state) {
    
    dir <- paste(dir_e, cancer, '_', state, '.txt', sep = '')
    input <- read.table(dir, sep = '\t', header = T, row.names = NULL)
    
    rn <- input[,1]
    rn <- make.unique(rn, sep = '.')
    
    output <- data.frame(input[,-1])
    rownames(output) <- rn
    
    # Fix col names
    names <- colnames(output)
    for (i in 1:length(names)) {
      split <- strsplit(names[i], '\\.')
      names[i] <- paste(split[[1]][1], split[[1]][2], 
                        split[[1]][3], sep = '.')
    }
    colnames(output) <- names
    
    return(output)
    
  }
  
  ###############################################################################################
  # Expression box/violin plots for every gene in a given signature
  ###############################################################################################
  
  ebp_sign <- function(cancer_index) {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    
    # Box/violin plot dimensions
    w <- 20 # width default: 8
    h <- 10 # height default: 6
    # get data
    cancer <- cancers[cancer_index]
    if (("BRCA1_mut" %in% cancers|"OV" %in% cancers|
         "BRCA CDKN1_mut" %in% cancers|
         "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
         "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
         "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
         "UCEC CDKN1_mut" %in% cancers) & cancer_index == length(cancers)) {
      print("skipping... no normal tissues")
    } else if ("OV p53_mut" %in% cancers) {
      print("skipping... no normal tissues")
    } else if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
      print("skipping... no normal tissues")
    } else if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & cancer_index == 3) {
      exp_n <- exp_n_all[[2]]
    } else if (("KIRP p53_mut" %in% cancers) & cancer_index %in%  c(2,12)) {
      print("skipping... no normal tissues")
    } else {
      exp_n <- exp_n_all[[cancer_index]]
    }
    
    exp_t <- exp_t_all[[cancer_index]]
    
    # mean expression per sample (for next step in analysis)
    if (("BRCA1_mut" %in% cancers|"BRCA CDKN1_mut" %in% cancers|
         "OV" %in% cancers|"BLCA CDKN1_mut" %in% cancers|
         "LUAD CDKN1_mut" %in% cancers|
         "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
         "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
         "UCEC CDKN1_mut" %in% cancers) & cancer_index == length(cancers)) {
      print("skipping... no normal tissues")
    } else if ("OV p53_mut" %in% cancers) {
      print("skipping... no normal tissues")
    } else if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
      print("skipping... no normal tissues")
    } else if (("KIRP p53_mut" %in% cancers) & cancer_index %in% c(2, 12)) {
      print("skipping... no normal tissues")
    } else {
      meps_n <- matrix(nrow = length(signatures), ncol = ncol(exp_n))
    }
    meps_t <- matrix(nrow = length(signatures), ncol = ncol(exp_t))
    rn <- c()
    gc()
    # for each signature
    for (i in 1:length(signatures)) {
      
      # signature name
      name <- colnames(signatures[[i]])
      
      # subsize the expression matrices
      if (("BRCA1_mut" %in% cancers|"OV" %in% cancers|"BRCA CDKN1_mut" %in% cancers| 
           "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
           "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
           "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & cancer_index == length(cancers)) {
        print("skipping... no normal tissues")
      } else if ("OV p53_mut" %in% cancers) {
        print("skipping... no normal tissues")
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
        print("skipping... no normal tissues")
      } else if ("KIRP p53_mut" %in% cancers & cancer_index %in% c(2,12)) {
        print("skipping... no normal tissues")
      } else {
        exp_n_sub <- exp_n[as.character(signatures[[i]][,1]),]
      }
      
      exp_t_sub <- exp_t[as.character(signatures[[i]][,1]),]
      
      # build data frame
      values <- c()
      genes <- c()
      status <- c()
      
      if (("OV" %in% cancers & cancer_index == length(cancers))|
          "OV p53_mut" %in% cancers) {
        print("skipping... no normal tissues")
        for (j in 1:nrow(exp_t_sub)) {
          for (k in 1:ncol(exp_t_sub)) {
            values <- c(values, exp_t_sub[j,k])
            genes <- c(genes, rownames(exp_t_sub)[j])
            status <- c(status, 'Tumour')
          }
        }
      } else if ("KIRP p53_mut" %in% cancers & cancer_index %in%  c(2,12)) {
        print("skipping... no normal tissues")
        for (j in 1:nrow(exp_t_sub)) {
          for (k in 1:ncol(exp_t_sub)) {
            values <- c(values, exp_t_sub[j,k])
            genes <- c(genes, rownames(exp_t_sub)[j])
            status <- c(status, 'Tumour')
          }
        }
      } else if (("BRCA1_mut" %in% cancers|
                  "BRCA CDKN1_mut" %in% cancers|
                  "BLCA CDKN1_mut" %in% cancers|
                  "LUAD CDKN1_mut" %in% cancers|
                  "COAD CDKN1_mut" %in% cancers|
                  "LUSC CDKN1_mut" %in% cancers|
                  "PRAD CDKN1_mut" %in% cancers|
                  "STAD CDKN1_mut" %in% cancers |
                  "UCEC CDKN1_mut" %in% cancers) & cancer_index %in% c(4:length(cancers))) {
        print("skipping... no normal tissues")
        for (j in 1:nrow(exp_t_sub)) {
          for (k in 1:ncol(exp_t_sub)) {
            values <- c(values, exp_t_sub[j,k])
            genes <- c(genes, rownames(exp_t_sub)[j])
            status <- c(status, 'Tumour')
          }
        }
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
        print("skipping... no normal tissues")
        for (j in 1:nrow(exp_t_sub)) {
          for (k in 1:ncol(exp_t_sub)) {
            values <- c(values, exp_t_sub[j,k])
            genes <- c(genes, rownames(exp_t_sub)[j])
            status <- c(status, 'Tumour')
          }
        }
      } else {
        for (j in 1:nrow(exp_n_sub)) {
          for (k in 1:ncol(exp_n_sub)) {
            values <- c(values, exp_n_sub[j,k])
            genes <- c(genes, rownames(exp_n_sub)[j])
            status <- c(status, 'Normal')
          }
        }
        for (j in 1:nrow(exp_t_sub)) {
          for (k in 1:ncol(exp_t_sub)) {
            values <- c(values, exp_t_sub[j,k])
            genes <- c(genes, rownames(exp_t_sub)[j])
            status <- c(status, 'Tumour')
          }
        }
      }
      df <- data.frame(values, genes, status)
      # remove rows with NA values 
      df <- df[!is.na(df$values),]
      ###########################################################################
      if (length(df[,1]) > 1) {
        if (i %in% c(4:15, 19)) {
          generalfontsize <- 30
          } else {
            generalfontsize <- 20
          }

        # output boxplots 
        plot_out <- ggplot(df, aes(x = interaction(status, genes), y = values)) +
          geom_boxplot(aes(color = status)) +
          ggtitle(paste('Expression of', name, 'genes in', cancer)) +
          xlab("Gene") + ylab("mRNA count (log2, normalized)") +
          scale_color_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
          theme_classic() +  theme(text = element_text(size = generalfontsize),
                                   axis.text.x = element_text(angle = 45, hjust = 1)) 
        
        
        # free means scales vary across all facets, free_x is across rows and
        # y is across col (so the bars do not have to be in different position in x)
        
        plot_out <- plot_out + facet_grid(.~genes, scales = "free_x", switch = "both")
        # remove original x-axis using element_blank()
        plot_out <- plot_out + theme_classic() +  theme(axis.text.x = element_blank(),
                                                        axis.ticks.x = element_blank(),
                                                        axis.text.y = element_text(size = 35),
                                                        axis.title = element_text(size = 35)) 
        
        # add significance layer
        if (("BRCA1_mut" %in% cancers & cancer_index|
            "BRCA CDKN1_mut" %in% cancers|"OV" %in% cancers|
             "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
             "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
             "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
             "UCEC CDKN1_mut" %in% cancers) & cancer_index %in% c(4:length(cancers))) {
          plot_out <- plot_out
        } else if ("OV p53_mut" %in% cancers) {
          plot_out <- plot_out
        } else if (("CESC CDKN1_mut" %in% cancers|
                      "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
          plot_out <- plot_out
        } else if ("KIRP p53_mut" %in% cancers & cancer_index %in% c(2,12)) {
          plot_out <- plot_out
        } else {
          plot_out <- plot_out + geom_signif(comparisons = combn(levels(df$status), 
                                                                 2, simplify = F))
        }
        
        
        # add sample size
        plot_out <- plot_out + stat_n_text(angle = 90) 
        
        # # adjust the width of the graph (if need be)
        # if (length(signatures[[i]][,1]) > 40) {
        #   w <- 20 
        # }
        
        pdf(file = paste(dir_o, 'exp boxplots per signature/', cancer, '_',
                         name, '.pdf', sep = ''), width = w, height = h)
        print(plot_out)
        dev.off()
        
        # Boxplots by genes in multiple facets 
        # t.test is parametric and wilcox.test is non-parametric
        
        # short.panel.labs if TRUE simplifies labels
        plot_out2_fxn <- function(df, numberofpages) {
          # numberofpages = " " for all gene sig in one page,
          # "1" for first 20 genes, "2" for the next 20 and so on 
          plot_out2 <- ggpar(ggboxplot(df, x = "status", y = "values", 
                                       color = "status", 
                                       palette = c(Normal = "darkturquoise", 
                                                   Tumour = "deeppink"), 
                                       line.color = "gray", line.size = 0.4, 
                                       facet.by = "genes", 
                                       short.panel.labs = TRUE) 
                             + labs(x = "genes", y =  "mRNA count (log2, normalized)", 
                                    title = paste('Expression of', name, 'genes in', cancer)), 
                             font.title = c(generalfontsize, "bold", "black"),
                             font.x = c(generalfontsize,"bold", "black"), 
                             font.y = c(generalfontsize,"bold", "black"),
                             xtickslab.rt = 0, 
                             legend = "right") 
          
          # parametric 
          # use p.format but if the p.signif is used, ns is non-sig and
          # * means sig p value 
          plot_out2_1 <- plot_out2 + 
            stat_compare_means(method = "t.test", label = "p.format",
                               paired = FALSE, label.x.npc = "centre",
                               label.y.npc = 0.95) 
          # add sample size
          plot_out2_1 <- plot_out2_1 + stat_n_text() 
          pdf(file = paste(dir_o, 'exp boxplots per signature/', cancer, '_',
                           name, 'ttest', numberofpages, '.pdf', sep = ''), 
              width = w, height = h)
          print(plot_out2_1)
          dev.off()
          
          # non-parametric
          plot_out2_2 <- plot_out2 + 
            stat_compare_means(method = "wilcox.test", label = "p.format",
                               paired = FALSE, label.x.npc = "centre", 
                               label.y.npc = 0.95)
          plot_out2_2 <- plot_out2_2 + stat_n_text() 
          print(plot_out2_2)
          pdf(file = paste(dir_o, 'exp boxplots per signature/', cancer, '_',
                           name, 'wtest', numberofpages, '.pdf', sep = ''), width = w, height = h)
          print(plot_out2_2)
          dev.off()
        }
        df <- df[order(df$genes),]
        gene_sig <- unique(df$genes)
        if (length(gene_sig) < 20) {
          plot_out2_fxn(df, "1")
        } else if (length(gene_sig) > 20 & length(gene_sig) < 41) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          rem_df <- df[df$genes %in% gene_sig[c(21:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(rem_df, "2")
        } else if (length(gene_sig) > 40 & length(gene_sig) < 61) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          rem_df <- df[df$genes %in% gene_sig[c(41:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(rem_df, "3")
        } else if (length(gene_sig) > 60 & length(gene_sig) < 81) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          rem_df <- df[df$genes %in% gene_sig[c(61:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(rem_df, "4")
        } else if (length(gene_sig) > 80 & length(gene_sig) < 101) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          new_df4 <- df[df$genes %in% gene_sig[c(61:80)],]
          rem_df <- df[df$genes %in% gene_sig[c(81:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(new_df4, "4")
          plot_out2_fxn(rem_df, "5")
        } else if (length(gene_sig) > 100 & length(gene_sig) < 121) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          new_df4 <- df[df$genes %in% gene_sig[c(61:80)],]
          new_df5 <- df[df$genes %in% gene_sig[c(81:100)],]
          rem_df <- df[df$genes %in% gene_sig[c(101:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(rem_df4, "4")
          plot_out2_fxn(rem_df5, "5")
          plot_out2_fxn(rem_df, "6")
        } else {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          new_df4 <- df[df$genes %in% gene_sig[c(61:80)],]
          new_df5 <- df[df$genes %in% gene_sig[c(81:100)],]
          new_df6 <- df[df$genes %in% gene_sig[c(101:120)],]
          new_df7 <- df[df$genes %in% gene_sig[c(121:140)],]
          new_df8 <- df[df$genes %in% gene_sig[c(141:160)],]
          new_df9 <- df[df$genes %in% gene_sig[c(161:180)],]
          rem_df <- df[df$genes %in% gene_sig[c(181:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(new_df4, "4")
          plot_out2_fxn(new_df5, "5")
          plot_out2_fxn(new_df6, "6")
          plot_out2_fxn(new_df7, "7")
          plot_out2_fxn(new_df8, "8")
          plot_out2_fxn(new_df9, "9")
          plot_out2_fxn(rem_df, "10")
        }
        ###########################################################
        ### CONINUE 
        # + geom_jitter(height = 0, width = 0.1)
        # output violin plots 
        plot_out <- ggplot(df, aes(x = interaction(status, genes), y = values)) +
          geom_violin(aes(color = status)) + geom_jitter(aes(color = status),
                                                         height = 0, width = 0.1) +
          ggtitle(paste('Expression of', name, 'genes in', cancer)) +
          xlab("Gene") + ylab("mRNA count (log2, normalized)") +
          scale_color_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
          theme_classic() + theme(text = element_text(size = generalfontsize),
                                  axis.text.x = element_text(angle = 45, hjust = 1)) 
        # free means scales vary across all facets, free_x is across rows and
        # y is across col (so the bars do not have to be in different position in x)
        
        
        plot_out <- plot_out + facet_grid(.~genes, scales =  "free_x", switch = "both")
        # remove original x-axis using element_blank()
        plot_out <- plot_out + theme_classic() + theme(axis.text.x = element_blank(),
                                                       axis.ticks.x = element_blank(),
                                                       axis.text.y = element_text(size = 35),
                                                       axis.title = element_text(size = 35)) 
        
        # add significance layer
        if (("BRCA1_mut" %in% cancers|
             "BRCA CDKN1_mut" %in% cancers|"OV" %in% cancers|
             "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
             "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
             "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
             "UCEC CDKN1_mut" %in% cancers) & cancer_index %in% c(4:length(cancers))) {
          plot_out <- plot_out
        } else if ("OV p53_mut" %in% cancers) {
          plot_out <- plot_out
        } else if (("CESC CDKN1_mut" %in% cancers|
                    "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
          plot_out <- plot_out
        } else if ("KIRP p53_mut" %in% cancers & cancer_index %in% c(2,12)) {
          plot_out <- plot_out
        } else {
          plot_out <- plot_out + geom_signif(comparisons = combn(levels(df$status), 
                                                                 2, simplify = F))
        }
        
        # add sample size
        plot_out <- plot_out + stat_n_text(angle = 90) 
        
        # # adjust the width of the graph (if need be)
        # if (length(signatures[[i]][,1]) > 40) {
        #   w <- 20 
        # }
        
        pdf(file = paste(dir_o, 'exp violin plots per signature/', cancer, '_',
                         name, '.pdf', sep = ''), width = w, height = h)
        print(plot_out)
        dev.off()
        
        # violin plots by genes in multiple facets 
        # t.test is parametric and wilcox.test is non-parametric
        
        # short.panel.labs if TRUE simplifies labels
        plot_out2_fxn <- function(df, numberofpages) {
          # numberofpages = " " for all gene sig in one page,
          # "1" for first 20 genes, "2" for the next 20 and so on 
          plot_out2 <- ggpar(ggviolin(df, x = "status", y = "values", 
                                      color = "status", 
                                      palette = c(Normal = "darkturquoise",
                                                  Tumour = "deeppink"), 
                                      line.color = "gray", line.size = 0.4, 
                                      facet.by = "genes", add = "jitter", size = 0.1,
                                      short.panel.labs = TRUE) 
                             + labs(x = "genes", y =  "mRNA count (log2, normalized)", 
                                    title = paste('Expression of', name, 
                                                  'genes in', cancer)), 
                             font.title = c(generalfontsize, "bold", "black"),
                             font.x = c(generalfontsize,"bold", "black"), 
                             font.y = c(generalfontsize,"bold", "black"),
                             xtickslab.rt = 0, 
                             legend = "right") 
          
          # parametric 
          # use p.format but if the p.signif is used, ns is non-sig and
          # * means sig p value 
          plot_out2_1 <- plot_out2 + 
            stat_compare_means(method = "t.test", label = "p.format",
                               paired = FALSE, label.x.npc = "centre",
                               label.y.npc = 0.95) 
          # add sample size
          plot_out2_1 <- plot_out2_1 + stat_n_text() 
          pdf(file = paste(dir_o, 'exp violin plots per signature/', cancer, '_',
                           name, 'ttest', numberofpages, '.pdf', sep = ''), 
              width = w, height = h)
          print(plot_out2_1)
          dev.off()
          
          # non-parametric
          plot_out2_2 <- plot_out2 + 
            stat_compare_means(method = "wilcox.test", label = "p.format",
                               paired = FALSE, label.x.npc = "centre", 
                               label.y.npc = 0.95)
          plot_out2_2 <- plot_out2_2 + stat_n_text() 
          print(plot_out2_2)
          pdf(file = paste(dir_o, 'exp violin plots per signature/', cancer, '_',
                           name, 'wtest', numberofpages, '.pdf', sep = ''),
              width = w, height = h)
          print(plot_out2_2)
          dev.off()
        }
        df <- df[order(df$genes),]
        gene_sig <- unique(df$genes)
        if (length(gene_sig) < 20) {
          plot_out2_fxn(df, "1")
        } else if (length(gene_sig) > 20 & length(gene_sig) < 41) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          rem_df <- df[df$genes %in% gene_sig[c(21:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(rem_df, "2")
        } else if (length(gene_sig) > 40 & length(gene_sig) < 61) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          rem_df <- df[df$genes %in% gene_sig[c(41:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(rem_df, "3")
        } else if (length(gene_sig) > 60 & length(gene_sig) < 81) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          rem_df <- df[df$genes %in% gene_sig[c(61:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(rem_df, "4")
        } else if (length(gene_sig) > 80 & length(gene_sig) < 101) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          new_df4 <- df[df$genes %in% gene_sig[c(61:80)],]
          rem_df <- df[df$genes %in% gene_sig[c(81:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(new_df4, "4")
          plot_out2_fxn(rem_df, "5")
        } else if (length(gene_sig) > 100 & length(gene_sig) < 121) {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          new_df4 <- df[df$genes %in% gene_sig[c(61:80)],]
          new_df5 <- df[df$genes %in% gene_sig[c(81:100)],]
          rem_df <- df[df$genes %in% gene_sig[c(101:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(rem_df4, "4")
          plot_out2_fxn(rem_df5, "5")
          plot_out2_fxn(rem_df, "6")
        } else {
          new_df <- df[df$genes %in% gene_sig[c(1:20)],]
          new_df2 <- df[df$genes %in% gene_sig[c(21:40)],]
          new_df3 <- df[df$genes %in% gene_sig[c(41:60)],]
          new_df4 <- df[df$genes %in% gene_sig[c(61:80)],]
          new_df5 <- df[df$genes %in% gene_sig[c(81:100)],]
          new_df6 <- df[df$genes %in% gene_sig[c(101:120)],]
          new_df7 <- df[df$genes %in% gene_sig[c(121:140)],]
          new_df8 <- df[df$genes %in% gene_sig[c(141:160)],]
          new_df9 <- df[df$genes %in% gene_sig[c(161:180)],]
          rem_df <- df[df$genes %in% gene_sig[c(181:length(gene_sig))],]
          plot_out2_fxn(new_df, "1")
          plot_out2_fxn(new_df2, "2")
          plot_out2_fxn(new_df3, "3")
          plot_out2_fxn(new_df4, "4")
          plot_out2_fxn(new_df5, "5")
          plot_out2_fxn(new_df6, "6")
          plot_out2_fxn(new_df7, "7")
          plot_out2_fxn(new_df8, "8")
          plot_out2_fxn(new_df9, "9")
          plot_out2_fxn(rem_df, "10")
        }
        ###########################################################
      }
      
      # perform t-test for each gene
      if (("OV" %in% cancers & cancer_index == length(cancers))|
           "OV p53_mut" %in% cancers) {
         print("skipping... no normal tissues") 
        
        # get mean expression per sample (for next step in analysis)
        if (um) {
          meps_t[i,] <- colMeans(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } else {
          meps_t[i,] <- colMedians(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } 
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
        print("skipping... no normal tissues") 
        
        # get mean expression per sample (for next step in analysis)
        if (um) {
          meps_t[i,] <- colMeans(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } else {
          meps_t[i,] <- colMedians(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } 
      } else if ("KIRP p53_mut" %in% cancers & cancer_index %in% c(2,12)) {
        print("skipping... no normal tissues") 
        
        # get mean expression per sample (for next step in analysis)
        if (um) {
          meps_t[i,] <- colMeans(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } else {
          meps_t[i,] <- colMedians(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } 
      } else if (("BRCA1_mut" %in% cancers|
                   "BRCA CDKN1_mut" %in% cancers|
                   "BLCA CDKN1_mut" %in% cancers|
                   "LUAD CDKN1_mut" %in% cancers|
                   "COAD CDKN1_mut" %in% cancers|
                   "LUSC CDKN1_mut" %in% cancers|
                   "PRAD CDKN1_mut" %in% cancers|
                   "STAD CDKN1_mut" %in% cancers|
                   "UCEC CDKN1_mut" %in% cancers) & cancer_index %in% c(4:length(cancers))) {
        print("skipping... no normal tissues") 
        # get mean expression per sample (for next step in analysis)
        if (um) {
          meps_t[i,] <- colMeans(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } else {
          meps_t[i,] <- colMedians(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        }
        
      } else {
        pv <- matrix(nrow = nrow(exp_n_sub), ncol = 1)
        for (j in 1:nrow(exp_n_sub)) {
          if ((length(na.omit(as.numeric(exp_n_sub[j,]))) < 3) | 
              (length(na.omit(as.numeric(exp_t_sub[j,]))) < 3)) {
            # filter out genes that have less than 3 values
            pv[j,1] <- NA
          } else {
            pv[j,1] <- t.test(exp_n_sub[j,], exp_t_sub[j,])$p.value
          }
        }
        rownames(pv) <- rownames(exp_n_sub)
        colnames(pv) <- 'p-value (t test)'
        write.table(pv, file = paste(dir_o, 'exp boxplots per signature/', cancer, 
                                     '_', name, '_', 'pv.txt', sep = ''), sep = '\t')
        # get mean expression per sample (for next step in analysis)
        if (um) {
          meps_n[i,] <- colMeans(exp_n_sub, na.rm = T)
          meps_t[i,] <- colMeans(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        } else {
          meps_n[i,] <- colMedians(exp_n_sub, na.rm = T)
          meps_t[i,] <- colMedians(exp_t_sub, na.rm = T)
          rn <- c(rn, name)
        }
      }
    }
    
    # output mean expression per sample for all signatures (for next step in analysis)
    if (("OV" %in% cancers & cancer_index == length(cancers))|
        "OV p53_mut" %in% cancers) {
      rownames(meps_t) <- rn
      colnames(meps_t) <- colnames(exp_t)
      meps <- list(meps_t)
    } else if (("CESC CDKN1_mut" %in% cancers|
           "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
      rownames(meps_t) <- rn
      colnames(meps_t) <- colnames(exp_t)
      meps <- list(meps_t)
    } else if ("KIRP p53_mut" %in% cancers & cancer_index == c(2,12)) {
      rownames(meps_t) <- rn
      colnames(meps_t) <- colnames(exp_t)
      meps <- list(meps_t)
    } else if (("BRCA1_mut" %in% cancers|
                 "BRCA CDKN1_mut" %in% cancers|
                 "BLCA CDKN1_mut" %in% cancers| 
                 "LUAD CDKN1_mut" %in% cancers|
                 "COAD CDKN1_mut" %in% cancers|
                 "LUSC CDKN1_mut" %in% cancers|
                 "PRAD CDKN1_mut" %in% cancers|
                 "STAD CDKN1_mut" %in% cancers|
                 "UCEC CDKN1_mut" %in% cancers) & cancer_index %in% c(4:length(cancers))) {
      rownames(meps_t) <- rn
      colnames(meps_t) <- colnames(exp_t)
      meps <- list(meps_t)
    } else {
      rownames(meps_n) <- rn
      colnames(meps_n) <- colnames(exp_n)
      rownames(meps_t) <- rn
      colnames(meps_t) <- colnames(exp_t)
      meps <- list(meps_n, meps_t)
    }
    
    return(meps)
  }
  
  ###############################################################################################
  # Expression box/violin plots across cancers
  ###############################################################################################
  
  ebp_cancers <- function(s_i) {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################

    # signature name
    name <- colnames(signatures[[s_i]])
    
    # build data frame
    # value
    v <- c()
    # cancer
    c <- c()
    state <- c()
    
    # for each cancer
    for (i in 1:length(meps_all)) {
      # for each sample in signature
      if (("OV" %in% cancers|"BRCA CDKN1_mut" %in% cancers|
           "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
           "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
           "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & i == length(cancers)) {
        print("skipping because no normal sample available")
        for (j in 1:ncol(meps_all[[i]][[1]])) {
          v <- c(v, meps_all[[i]][[1]][s_i, j])
          c <- c(c, cancers[i])
          state <- c(state, 'Tumour')
        }
      } else if ("OV p53_mut" %in% cancers) {
        print("skipping because no normal sample available")
        for (j in 1:ncol(meps_all[[i]][[1]])) {
          v <- c(v, meps_all[[i]][[1]][s_i, j])
          c <- c(c, cancers[i])
          state <- c(state, 'Tumour')
        }
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & i %in% c(2,4)) {
        print("skipping because no normal sample available")
        for (j in 1:ncol(meps_all[[i]][[1]])) {
          v <- c(v, meps_all[[i]][[1]][s_i, j])
          c <- c(c, cancers[i])
          state <- c(state, 'Tumour')
        }
      } else if ("KIRP p53_mut" %in% cancers & i %in% c(2,12)) {
        print("skipping because no normal sample available")
        for (j in 1:ncol(meps_all[[i]][[1]])) {
          v <- c(v, meps_all[[i]][[1]][s_i, j])
          c <- c(c, cancers[i])
          state <- c(state, 'Tumour')
        }
      } else {
        for (j in 1:ncol(meps_all[[i]][[1]])) {
          v <- c(v, meps_all[[i]][[1]][s_i, j])
          c <- c(c, cancers[i])
          state <- c(state, 'Normal')
        }
        
        for (j in 1:ncol(meps_all[[i]][[2]])) {
          v <- c(v, meps_all[[i]][[2]][s_i, j])
          c <- c(c, cancers[i])
          state <- c(state, 'Tumour')
        }
        
      }
    }
    
    df <- data.frame(v, c, state)
    
    # change font size according to how many cancers
    if (length(meps_all) < 5) {
        generalfontsize <- 30
    } else {
        generalfontsize <- 20
    }
  
    
    # output box plot
    if (um) {
      plot_out <- ggplot(df, aes(x = interaction(state, c), y = v)) +
        geom_boxplot(aes(color = state)) +
        ggtitle(paste('Expression of', name, 'genes across cancers')) +
        xlab("Cancer") + ylab("mean mRNA count (log2, normalized)") +
        scale_color_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
        theme_classic() + theme(text = element_text(size = generalfontsize),
              axis.text.x = element_text(angle = 45, hjust = 1))
      
    
      plot_out <- plot_out + facet_grid(.~c, scales =  "free_x", switch = "both")
      # remove original x-axis using element_blank()
      plot_out <- plot_out + theme_classic() + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(), 
                                   axis.text.y = element_text(size = 35),
                                   axis.title = element_text(size = 35)) 
      
      # add significance layer
      if (("BRCA1_mut" %in% cancers|"BRCA CDKN1_mut" %in% cancers|
           "OV" %in% cancers|"BLCA CDKN1_mut" %in% cancers|
           "LUAD CDKN1_mut" %in% cancers|"COAD CDKN1_mut" %in% cancers|
           "LUSC CDKN1_mut" %in% cancers|"PRAD CDKN1_mut" %in% cancers|
           "STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & i %in% c(4:length(cancers))) {
        plot_out <- plot_out
      } else if ("OV p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & i %in% c(2,4)) {
        plot_out <- plot_out
      } else if ("KIRP p53_mut" %in% cancers) {
        plot_out <- plot_out
        } else {
        plot_out <- plot_out + geom_signif(comparisons = combn(levels(df$state), 
                                                               2, simplify = F))
      }
      
      # add sample size
      plot_out <- plot_out + stat_n_text(angle = 90)  
      pdf(file = paste(dir_o, 'exp boxplots across cancers/', name, 
                       '.pdf', sep = ''),
          width = w, height = h)
      print(plot_out)
      dev.off()
      # Boxplots by genes in multiple facets 
      # t.test is parametric and wilcox.test is non-parametric
      
      # short.panel.labs if TRUE simplifies labels
      plot_out2 <- ggpar(ggboxplot(df, x = "state", y = "v", 
                                   color = "state", 
                                   palette = c(Normal = "darkturquoise", Tumour = "deeppink"), 
                                   line.color = "gray", line.size = 0.4, 
                                   facet.by = "c", 
                                   short.panel.labs = TRUE) 
                         + labs(x = "Cancer patient groups",
                                y = "mean mRNA count (log2, normalized)", 
                                title = paste('Expression of', name, 
                                              'genes across Cancer patient groups')), 
                         font.title = c(generalfontsize, "bold", "black"),
                         font.x = c(generalfontsize,"bold", "black"), 
                         font.y = c(generalfontsize,"bold", "black"),
                         xtickslab.rt = 0, 
                         legend = "right") 
      
      # parametric 
      # use p.format but if the p.signif is used, ns is non-sig and
      # * means sig p v 
      plot_out2_1 <- plot_out2 + stat_compare_means(method = "t.test", 
                                                    label = "p.format",
                                                    paired = FALSE, 
                                                    label.x.npc = "centre",
                                                    label.y.npc = 0.95) 
      # add sample size
      plot_out2_1 <- plot_out2_1 + stat_n_text() 
      pdf(file = paste(dir_o, 'exp boxplots across cancers/', 
                       name, 'ttest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_1)
      dev.off()
      # non-parametric
      plot_out2_2 <- plot_out2 + stat_compare_means(method = "wilcox.test", 
                                                    label = "p.format",
                                                    paired = FALSE,
                                                    label.x.npc = "centre", 
                                                    label.y.npc = 0.95)
      plot_out2_2 <- plot_out2_2 + stat_n_text() 
      print(plot_out2_2)
      pdf(file = paste(dir_o, 'exp boxplots across cancers/', name,
                       'wtest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_2)
      dev.off()
    } else {
      plot_out <- ggplot(df, aes(x = interaction(state, c), y = v)) +
        geom_boxplot(aes(color = state)) +
        ggtitle(paste('Expression of', name, 'genes across cancers')) +
        xlab("Cancer") + ylab("median mRNA count (log2, normalized)") +
        scale_fill_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
        theme_classic() + theme(text = element_text(size = generalfontsize),
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      plot_out <- plot_out + facet_grid(.~c, scales = "free_x", switch = "both")
      # remove original x-axis using element_blank()
      plot_out <- plot_out + theme_classic() + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   axis.text.y = element_text(size = 35),
                                   axis.title = element_text(size = 35)) 
      
      # add significance layer
      if (("BRCA1_mut" %in% cancers|"BRCA CDKN1_mut" %in% cancers|
           "OV" %in% cancers|
           "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
           "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
           "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & i %in% c(4:length(cancers))) {
        plot_out <- plot_out
      } else if ("OV p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & i %in% c(2,4)) {
        plot_out <- plot_out
      } else if ("KIRP p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else {
        plot_out <- plot_out + geom_signif(comparisons = combn(levels(df$state), 
                                                               2, simplify = F))
      }
      # add sample size
      plot_out <- plot_out + stat_n_text(angle = 90)  
      pdf(file = paste(dir_o, 'exp boxplots across cancers/', name, 
                       '.pdf', sep = ''),
          width = w, height = h)
      print(plot_out)
      dev.off()
      # Boxplots by genes in multiple facets 
      # t.test is parametric and wilcox.test is non-parametric
      
      # short.panel.labs if TRUE simlifies labels
      plot_out2 <- ggpar(ggboxplot(df, x = "state", y = "v", 
                                   color = "state", 
                                   palette = c(Normal = "darkturquoise",
                                               Tumour = "deeppink"), 
                                   line.color = "gray", line.size = 0.4, 
                                   facet.by = "c", 
                                   short.panel.labs = TRUE) 
                         + labs(x = "Cancer patient groups",
                                y = "median mRNA count (log2, normalized)", 
                                title = paste('Expression of', name, 
                                              'genes across Cancer patient groups')), 
                         font.title = c(generalfontsize, "bold", "black"),
                         font.x = c(generalfontsize,"bold", "black"), 
                         font.y = c(generalfontsize,"bold", "black"),
                         xtickslab.rt = 0, 
                         legend = "right") 
      
      # parametric 
      # use p.format but if the p.signif is used, ns is non-sig and
      # * means sig pv 
      plot_out2_1 <- plot_out2 + stat_compare_means(method = "t.test", 
                                                    label = "p.format",
                                                    paired = FALSE,
                                                    label.x.npc = "centre",
                                                    label.y.npc = 0.95) 
      # add sample size
      plot_out2_1 <- plot_out2_1 + stat_n_text() 
      pdf(file = paste(dir_o, 'exp boxplots across cancers/',
                       name, 'ttest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_1)
      dev.off()
      # non-parametric
      plot_out2_2 <- plot_out2 + stat_compare_means(method = "wilcox.test",
                                                    label = "p.format",
                                                    paired = FALSE, 
                                                    label.x.npc = "centre", 
                                                    label.y.npc = 0.95)
      plot_out2_2 <- plot_out2_2 + stat_n_text() 
      print(plot_out2_2)
      pdf(file = paste(dir_o, 'exp boxplots across cancers/',
                       name, 'wtest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_2)
      dev.off()
    }
    
    # output violin plot
    if (um) {
      plot_out <- ggplot(df, aes(x = interaction(state, c), y = v)) +
        geom_violin(aes(color = c)) + geom_jitter(aes(color = state),
                                                      height = 0, width = 0.1) +
        ggtitle(paste('Expression of', name, 'genes across cancers')) +
        xlab("Cancer") + ylab("mean mRNA count (log2, normalized)") +
        scale_color_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
        theme_classic() + theme(text = element_text(size = generalfontsize),
                      axis.text.x = element_text(angle = 45, hjust = 1))
      
      
      
      plot_out <- plot_out + facet_grid(.~c, scales =  "free_x", switch = "both")
      # remove original x-axis using element_blank()
      plot_out <- plot_out + theme_classic() + theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.text.y = element_text(size = 35),
                                           axis.title = element_text(size = 35)) 
      
      # add significance layer
      if (("BRCA1_mut" %in% cancers|
          "BRCA CDKN1_mut" %in% cancers|
           "BLCA CDKN1_mut" %in% cancers|
          "LUAD CDKN1_mut" %in% cancers|
           "COAD CDKN1_mut" %in% cancers|
          "LUSC CDKN1_mut" %in% cancers|
           "PRAD CDKN1_mut" %in% cancers|
          "STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & i %in% c(4:length(cancers))) {
        plot_out <- plot_out
      } else if ("OV p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & i %in% c(2,4)) {
        plot_out <- plot_out
      } else if ("KIRP p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else {
        plot_out <- plot_out + geom_signif(comparisons = combn(levels(df$state), 
                                                               2, simplify = F))
      }
      
      # add sample size
      plot_out <- plot_out + stat_n_text(angle = 90)  
      pdf(file = paste(dir_o, 'exp violin plots across cancers/', name, 
                       '.pdf', sep = ''),
          width = w, height = h)
      print(plot_out)
      dev.off()
      # violin plots by genes in multiple facets 
      # t.test is parametric and wilcox.test is non-parametric
      
      # short.panel.labs if TRUE simplifies labels
      plot_out2 <- ggpar(ggviolin(df, x = "state", y = "v", 
                                  color = "state", 
                                  palette = c(Normal = "darkturquoise",
                                              Tumour = "deeppink"), 
                                  line.color = "gray", line.size = 0.4, 
                                  add = "jitter", size = 0.1,
                                  facet.by = "c", 
                                  short.panel.labs = TRUE) 
                         + labs(x = "Cancer patient groups",
                                y =  "mean mRNA count (log2, normalized)", 
                                title = paste('Expression of', name, 
                                              'genes across Cancer patient groups')), 
                         font.title = c(generalfontsize, "bold", "black"),
                         font.x = c(generalfontsize,"bold", "black"), 
                         font.y = c(generalfontsize,"bold", "black"),
                         xtickslab.rt = 0, 
                         legend = "right") 
      
      # parametric 
      # use p.format but if the p.signif is used, ns is non-sig and
      # * means sig p v 
      plot_out2_1 <- plot_out2 + stat_compare_means(method = "t.test", 
                                                    label = "p.format",
                                                    paired = FALSE, 
                                                    label.x.npc = "centre",
                                                    label.y.npc = 0.95) 
      # add sample size
      plot_out2_1 <- plot_out2_1 + stat_n_text() 
      pdf(file = paste(dir_o, 'exp violin plots across cancers/', 
                       name, 'ttest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_1)
      dev.off()
      # non-parametric
      plot_out2_2 <- plot_out2 + stat_compare_means(method = "wilcox.test", 
                                                    label = "p.format",
                                                    paired = FALSE,
                                                    label.x.npc = "centre", 
                                                    label.y.npc = 0.95)
      plot_out2_2 <- plot_out2_2 + stat_n_text() 
      print(plot_out2_2)
      pdf(file = paste(dir_o, 'exp violin plots across cancers/', name,
                       'wtest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_2)
      dev.off()
    } else {
      plot_out <- ggplot(df, aes(x = interaction(state, c), y = v)) +
        geom_violin(aes(color = state)) + geom_jitter(aes(color = state),
                                                      height = 0, width = 0.1) +
        ggtitle(paste('Expression of', name, 'genes across cancers')) +
        xlab("Cancer") + ylab("median mRNA count (log2, normalized)") +
        scale_fill_manual(values = c(Normal = "darkturquoise", Tumour = "deeppink")) +
        theme_classic() + theme(text = element_text(size = generalfontsize),
                      axis.text.x = element_text(angle = 45, hjust = 1))
     
      plot_out <- plot_out + facet_grid(.~c, scales = "free_x", switch = "both")
      # remove original x-axis using element_blank()
      plot_out <- plot_out + theme_classic() + theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank()) 
      
      # add significance layer
      if (("BRCA1_mut" %in% cancers|
           "BRCA CDKN1_mut" %in% cancers |
           "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
           "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
           "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & i %in% c(4:length(cancers))) {
        plot_out <- plot_out
      } else if ("OV p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else if (("CESC CDKN1_mut" %in% cancers|
                  "THCA CDKN1_mut" %in% cancers) & i %in% c(2,4)) {
        plot_out <- plot_out
      } else if ("KIRP p53_mut" %in% cancers) {
        plot_out <- plot_out
      } else {
        plot_out <- plot_out + geom_signif(comparisons = combn(levels(df$state), 
                                                               2, simplify = F))
      }
      # add sample size
      plot_out <- plot_out + stat_n_text(angle = 90)  
      pdf(file = paste(dir_o, 'exp violin plots across cancers/', name, 
                       '.pdf', sep = ''),
          width = w, height = h)
      print(plot_out)
      dev.off()
      # violin plots by genes in multiple facets 
      # t.test is parametric and wilcox.test is non-parametric
      
      # short.panel.labs if TRUE simlifies labels
      plot_out2 <- ggpar(ggviolin(df, x = "state", y = "v", 
                                  color = "state", 
                                  palette = c(Normal = "darkturquoise", 
                                              Tumour = "deeppink"), 
                                  line.color = "gray", line.size = 0.4, 
                                  add = "jitter", size = 0.1,
                                  facet.by = "c", 
                                  short.panel.labs = TRUE) 
                         + labs(x = "Cancer patient groups",
                                y =  "median mRNA count (log2, normalized)", 
                                title = paste('Expression of', name, 
                                              'genes across Cancer patient groups')), 
                         font.title = c(generalfontsize, "bold", "black"),
                         font.x = c(generalfontsize,"bold", "black"), 
                         font.y = c(generalfontsize,"bold", "black"),
                         xtickslab.rt = 0, 
                         legend = "right") 
      
      # parametric 
      # use p.format but if the p.signif is used, ns is non-sig and
      # * means sig p v 
      plot_out2_1 <- plot_out2 + stat_compare_means(method = "t.test", 
                                                    label = "p.format",
                                                    paired = FALSE,
                                                    label.x.npc = "centre",
                                                    label.y.npc = 0.95) 
      # add sample size
      plot_out2_1 <- plot_out2_1 + stat_n_text() 
      pdf(file = paste(dir_o, 'exp violin plots across cancers/',
                       name, 'ttest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_1)
      dev.off()
      # non-parametric
      plot_out2_2 <- plot_out2 + stat_compare_means(method = "wilcox.test",
                                                    label = "p.format",
                                                    paired = FALSE, 
                                                    label.x.npc = "centre", 
                                                    label.y.npc = 0.95)
      plot_out2_2 <- plot_out2_2 + stat_n_text() 
      print(plot_out2_2)
      pdf(file = paste(dir_o, 'exp violin plots across cancers/',
                       name, 'wtest.pdf', sep = ''),
          width = w, height = h)
      print(plot_out2_2)
      dev.off()
    }
    
    if (("OV" %in% cancers|"BRCA CDKN1_mut" %in% cancers|
         "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
         "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
         "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
         "UCEC CDKN1_mut" %in% cancers) & i == length(cancers)) {
      # perform t-test for each gene
      print("skipping... no normal samples")
    } else if ("OV p53_mut" %in% cancers) {
      # perform t-test for each gene
      print ("skipping... no normal samples")
    } else if ("BRCA1_mut" %in% cancers & i %in% c(4:6)) {
      print("skipping... <3 normal samples")
    } else if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & i %in% c(2,4)) {
      print ("skipping... no normal samples")
    } else if ("KIRP p53_mut" %in% cancers) {
      print ("skipping... no normal samples")
    } else {
      # perform t-test for each gene
      pv <- matrix(nrow = length(meps_all), ncol = 1)
      for (i in 1:length(meps_all)) {
        if ((length(na.omit(as.numeric(meps_all[[i]][[1]][s_i,]))) < 3)|
            (length(na.omit(as.numeric(meps_all[[i]][[2]][s_i,]))) < 3)) {
          # filter out < 3 values
          pv[i,1] <- NA
        } else {
          pv[i,1] <- t.test(meps_all[[i]][[1]][s_i,], meps_all[[i]][[2]][s_i,])$p.value
        }
      }
      rownames(pv) <- cancers
      
      colnames(pv) <- 'p-value (t test)'
      
      write.table(pv, file = paste(dir_o, 'exp boxplots across cancers/', name, 
                                   'pv.txt', sep = ''), sep = '\t')
      
      model <- aov(formula = v ~ c + state, data = df)
      tukey <- anova(model)
      
      write.table(tukey, file = paste(dir_o, 'exp boxplots across cancers/', name, 
                                      'anova.txt', sep = ''), sep = '\t') 
    }
    
  }
  
  ###############################################################################################
  # Generate expression heatmap
  ###############################################################################################
  
  exp_hm <- function(cancer_index){
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    # combine all signatures into one list
    combined_sign <- c()
    for (i in 1:length(signatures)) {
      for (j in 1:nrow(signatures[[i]])) {
        combined_sign <- c(combined_sign, as.character(signatures[[i]][j,1]))
      }
    }
    
    # subset data
    cancer <- cancers[cancer_index]
    exp_t <- exp_t_all[[cancer_index]]
    exp_t_sub <- exp_t[combined_sign,]
    
    # remove NAs
    exp_t_sub[is.na(exp_t_sub)] <- 0
    
    # heatmap limits
    med <- exp_med[cancer_index, 1]
    ul <- 2*med
    ll <- 0
    
    # plot heatmap
    pdf(file = paste(dir_o, 'exp heatmaps/', cancer, '.pdf', sep = ''),
        width = hm_w, height = hm_h)
    hm <- heatmap.2(
      data.matrix(exp_t_sub),
      # main=paste("1. Mean -", m),
      Rowv = NA,
      dendrogram = "none",
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      labCol = NA,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = margins,
      breaks = seq(ll, ul, length.out = 101)
    )
    dev.off()
    graphics.off()
    # normalize for housekeeping genes
    hk_mean <- t(colMeans(exp_t_sub[as.character(signatures[[length(signatures)]][,1]),],
                          na.rm = T))
    exp_t_sub_norm <- exp_t_sub
    for (i in 1:ncol(exp_t_sub_norm)) {
      exp_t_sub_norm[,i] <- exp_t_sub_norm[,i] - hk_mean[i]
    }
    
    # heatmap limits
    ul <- 10
    ll <- -10
    
    # plot heatmap
    
    pdf(file = paste(dir_o, 'exp heatmaps/norm_hk/', cancer, '.pdf', sep = ''),
        width = hm_w, height = hm_h)
    hm <- heatmap.2(
      data.matrix(exp_t_sub_norm),
      #main=paste("1. Mean -", m),
      Rowv = NA,
      dendrogram = "none",
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      labCol = NA,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = margins,
      breaks = seq(ll, ul, length.out = 101)
    )
    dev.off()
    graphics.off()
  }
  
  ###############################################################################################
  # Correlation heatmaps (single genes)
  ###############################################################################################
  
  cor_hm <- function(cancer_index) {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    # combine all tested signatures into one list
    combined_sign <- c()
    # for (i in 2:length(signatures)){
    #   for (j in 1:nrow(signatures[[i]])){
    #     combined_sign <- c(combined_sign, as.character(signatures[[i]][j,1]))
    #   }
    # }
    # OVERRIDE: DO NOT INCLUDE OTHER SIGNATURES
    combined_sign <- as.character(signatures[[2]][,1]) 
    # master signature
    master_sign <- as.character(signatures[[1]][,1])
    
    # subset data
    cancer <- cancers[cancer_index]
    exp_t <- exp_t_all[[cancer_index]]
    exp_t_sub <- exp_t[combined_sign,] # tested signatures
    exp_t_master <- exp_t[master_sign,] # master signature
    
    if (("OV" %in% cancers & cancer_index == length(cancers))|
         "OV p53_mut" %in% cancers) {
      print("no normal tissues with mutation")
    } else if (("BRCA1_mut" %in% cancers|"BRCA CDKN1_mut" %in% cancers| 
                "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
                "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
                "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
                "UCEC CDKN1_mut" %in% cancers) & cancer_index %in% c(4:6)) {
      print("<3 normal tissues with mutation")
    } else if ("KIRP p53_mut" %in% cancers) {
      print ("skipping... no normal samples")
    } else if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
      print("<3 normal tissues with mutation")
    } else {
      if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & cancer_index == 3) {
        exp_n <- data.frame(exp_n_all[[2]])
      } else {
        exp_n <- data.frame(exp_n_all[[cancer_index]])
      }
      
      exp_n_sub <- data.frame(exp_n[combined_sign,]) # tested signatures
      exp_n_master <- data.frame(exp_n[master_sign,]) # master signature
      # correlation analysis (spearman)
      table_cor <- matrix(nrow = nrow(exp_t_master), ncol = nrow(exp_t_sub))
      table_pv <- matrix(nrow = nrow(exp_t_master), ncol = nrow(exp_t_sub))
      
      for (i in 1:nrow(exp_t_master)) {
        for (j in 1:nrow(exp_t_sub)) {
          result = tryCatch({
            test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(exp_t_sub[j,]),
                             method = 'spearman')
            table_cor[i,j] <- test$estimate
            table_pv[i,j] <- test$p.value
          }, error = function(e) {
            table_cor[i,j] <- NA
            table_pv[i,j] <- NA
          }, finally = {
            
          })
        }
      }
      
      # remove insignificant values
      table_cor[table_pv >= 0.05] <- 0
      
      # naming
      rownames(table_cor) <- master_sign
      colnames(table_cor) <- combined_sign
      rownames(table_pv) <- master_sign
      colnames(table_pv) <- combined_sign
      
      # plot heatmap
      pdf(file = paste(dir_o, 'cor single gene/', cancer, '.pdf', sep = ''), 
          width = 12, height = 4)
      hm <- heatmap.2(
        data.matrix(table_cor),
        #main=paste("1. Mean -", m),
        Rowv = NA,
        col = bluered,
        key = key, symkey = skey, keysize = keysize,
        cexRow = fontsize, cexCol = fontsize,
        density.info = "none", trace = "none",
        margins = margins,
        breaks = seq(-1, 1, length.out = 101)
      )
      dev.off()
      graphics.off()
      
      # write tables
      write.table(table_cor, file = paste(dir_o, 'cor single gene/', cancer,
                                          '_cor.txt', sep = ''), sep = '\t')
      write.table(table_pv, file = paste(dir_o, 'cor single gene/', cancer, 
                                         '_pv.txt', sep = ''), sep = '\t')
    }
    
    # get mean of signatures (for downstream analysis)
    table_means <- matrix(nrow = length(signatures) - 1, ncol = ncol(exp_t_sub))
    for (i in 2:length(signatures)) {
      if (um) {
        table_means[i - 1,] <- colMeans(exp_t[as.character(signatures[[i]][,1]),])
      } else {
        table_means[i - 1,] <- colMedians(exp_t[as.character(signatures[[i]][,1]),])
      }
    }
    
    rownames(table_means) <- sig_names[-1]
    na_filter <-
      (rowCounts(!is.na(table_means))/ncol(table_means) >= thr_na)
    table_means <- table_means[na_filter,]
    new_sig_names <- rownames(table_means)
    
    if (("OV" %in% cancers|
         "BRCA CDKN1_mut" %in% cancers|
         "BLCA CDKN1_mut" %in% cancers|
         "LUAD CDKN1_mut" %in% cancers|
         "COAD CDKN1_mut" %in% cancers|
         "LUSC CDKN1_mut" %in% cancers|
         "PRAD CDKN1_mut" %in% cancers|
         "STAD CDKN1_mut" %in% cancers|
         "UCEC CDKN1_mut" %in% cancers) & cancer_index == length(cancers)) {
      print("no normal tissue samples")
    } else if ("OV p53_mut" %in% cancers) {
      print("skipping... no normal tissues")
    } else if ("KIRP p53_mut" %in% cancers) {
      print ("skipping... no normal samples")
    } else if (("CESC CDKN1_mut" %in% cancers|
                "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
      print("skipping... no normal tissues")
    } else if ("BRCA1_mut" %in% cancers & cancer_index %in% c(4:length(cancers))) {
      print("<2 normal tissue samples")
    } else {
      table_means_n <- matrix(nrow = length(signatures) - 1, ncol = ncol(exp_n_sub))
      for (i in 2:length(signatures)) {
        if (um) {
          table_means_n[i - 1,] <- colMeans(exp_n[as.character(signatures[[i]][,1]),])
        } else {
          table_means_n[i - 1,] <- colMedians(exp_n[as.character(signatures[[i]][,1]),])
        }
      }
    }
    
    # correlate means with master signature
    table_cor_means <- matrix(nrow = nrow(exp_t_master), ncol = nrow(table_means))
    table_pv_means <- matrix(nrow = nrow(exp_t_master), ncol = nrow(table_means))
    
    for (i in 1:nrow(exp_t_master)) {
      for (j in 1:nrow(table_means)) {
        
        result = tryCatch({
          test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(table_means[j,]),
                                    method = 'spearman')
          table_cor_means[i,j] <- test$estimate
          table_pv_means[i,j] <- test$p.value
        }, error = function(e) {
          table_cor_means[i,j] <- NA
          table_pv_means[i,j] <- NA
        }, finally = {
          
        })
      }
    }
  
    # Generate scatter plots (tumour vs normal cancer cells)
    for (i in 1:nrow(exp_t_master)) {
      
      out_dir <- paste(dir_o, 'cor mean master scatterplots/', cancer, '-', 
                       as.character(signatures[[1]][,1])[i], '_cor.pdf', sep = '')
      a <- as.numeric(exp_t_master[i,])
      h <- as.numeric(table_means[1,])
      df <- data.frame(a, h)
      plot <- ggplot(df, aes(x = a, y = h)) +
        geom_point(alpha = 0.6, colour = "deeppink") +
        # draw a regression line
        geom_smooth(colour = "black", fill = "white", method = "lm",
                      se = FALSE, formula =  y~x) +
        ggtitle(paste(as.character(signatures[[1]][,1])[i], 'vs Hypoxia -',
                      cancer, '(Tumour)')) +
        xlab(as.character(signatures[[1]][,1])[i]) +
        ylab('mean Hypoxia signature') +
        theme_classic() + theme(text = element_text(size = 12))
      plot <- plot  + stat_cor(method = "spearman", label.x.npc = "middle") 
  
      pdf(file = out_dir, width = 6, height = 4)
      
      print(plot)
      dev.off()
    }
    
    for (i in 1:nrow(exp_t_master)) {
      out_dir <- paste(dir_o, 'cor mean master scatterplots/_NT/', cancer, '-',
                       as.character(signatures[[1]][,1])[i], '_cor_NT.pdf', sep = '')
      if (("OV" %in% cancers|
           "BRCA CDKN1_mut" %in% cancers|
           "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
           "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
           "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
           "UCEC CDKN1_mut" %in% cancers) & cancer_index == length(cancers)) {
          print("no normal tissue samples")
          a <- c(as.numeric(exp_t_master[i,]))
          h <- c(as.numeric(table_means[1,]))
          tissue <- rep('tumour', length(as.numeric(exp_t_master[i,])))
          df <- data.frame(a, h, tissue)
  
        } else if ("OV p53_mut" %in% cancers) {
          print("no normal tissue samples")
          a <- c(as.numeric(exp_t_master[i,]))
          h <- c(as.numeric(table_means[1,]))
          tissue <- rep('tumour', length(as.numeric(exp_t_master[i,])))
          df <- data.frame(a, h, tissue)
        
        } else if ("BRCA1_mut" %in% cancers & cancer_index %in% c(4:length(cancers))) {
          print("no normal tissue samples")
          a <- c(as.numeric(exp_t_master[i,]))
          h <- c(as.numeric(table_means[1,]))
          tissue <- rep('tumour', length(as.numeric(exp_t_master[i,])))
          df <- data.frame(a, h, tissue)
        
        } else if (("CESC CDKN1_mut" %in% cancers|
                     "THCA CDKN1_mut" %in% cancers) & cancer_index %in% c(2,4)) {
          print("no normal tissue samples")
          a <- c(as.numeric(exp_t_master[i,]))
          h <- c(as.numeric(table_means[1,]))
          tissue <- rep('tumour', length(as.numeric(exp_t_master[i,])))
          df <- data.frame(a, h, tissue)
          
        } else if ("KIRP p53_mut" %in% cancers) {
          print("no normal tissue samples")
          a <- c(as.numeric(exp_t_master[i,]))
          h <- c(as.numeric(table_means[1,]))
          tissue <- rep('tumour', length(as.numeric(exp_t_master[i,])))
          df <- data.frame(a, h, tissue)
          
        } else {
          a <- c(as.numeric(exp_t_master[i,]), as.numeric(exp_n_master[i,]))
          h <- c(as.numeric(table_means[1,]), as.numeric(table_means_n[1,]))
          tissue  <- c(rep('tumour', length(as.numeric(exp_t_master[i,]))), 
                     rep('normal', length(as.numeric(exp_n_master[i,]))))
          df <- data.frame(a, h, tissue)
        
      }
      
      plot <- ggplot(df, aes(x = a, y = h, color = tissue)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = c(normal = "darkturquoise", tumour = "deeppink")) +
        # draw a regression line
        geom_smooth(colour = "black", fill = "white", method = "lm",
                    se = FALSE, formula =  y~x) +
        ggtitle(paste(as.character(signatures[[1]][,1])[i], 'vs Hypoxia -', 
                      cancer, '(Tumour)')) +
        xlab(as.character(signatures[[1]][,1])[i]) +
        ylab('mean Hypoxia signature') + theme_classic() +
        theme(text = element_text(size = 12))
      plot <- plot  + stat_cor(method = "spearman", label.x.npc = "middle")
      pdf(file = out_dir, width = 6, height = 4)
      print(plot)
      dev.off()
    }
    
    # remove insignificant values
    table_cor_means[table_pv_means >= 0.05] <- 0
    
    # get signature names
    sign_names <- c()
    
    for (i in 2:length(signatures)) {
      sign_names <- c(new_sig_names, colnames(signatures[[i]]))
    }
    
    # naming
    rownames(table_cor_means) <- master_sign
    colnames(table_cor_means) <- new_sig_names
    rownames(table_pv_means) <- master_sign
    colnames(table_pv_means) <- new_sig_names
    
    if (length(colnames(table_cor_means)) < nts) {
      sig_names2 <- sig_names[-1]
      too_many_NAs_genes <- sig_names2[!sig_names2 %in% new_sig_names] 
      NA_df <- matrix(nrow = length(master_sign), ncol = length(too_many_NAs_genes))
      colnames(NA_df) <- too_many_NAs_genes
      table_cor_means <- cbind(table_cor_means, NA_df)
      # reorder columns based on original order of signatures
      table_cor_means <- table_cor_means[,sig_names2]
      table_pv_means <- cbind(table_pv_means, NA_df)
      # reorder columns based on original order of signatures
      table_pv_means <- table_pv_means[,sig_names2]
    }
    
    # output
    out <- list(table_cor_means, table_pv_means)
    return(out)
    
  }
  
  ###############################################################################################
  # Correlation heatmaps (signature means)
  ###############################################################################################
  
  cor_mean_hm <- function(mean_cor) {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    
    for (i in 2:length(signatures)) {
      
      sign_name <- colnames(signatures[[i]])
      
      # assemble matrix
      mat_cor <- matrix(nrow = nrow(signatures[[1]]), ncol = length(cancers))
      mat_pv <- matrix(nrow = nrow(signatures[[1]]), ncol = length(cancers))
      for (j in 1:length(cancers)) {
        mc <- mean_cor[[j]]
        mc_cor <- mc[[1]]
        mc_pv <- mc[[2]]
        mat_cor[,j] <- mc_cor[, i - 1]
        mat_pv[,j] <- mc_pv[, i - 1]
      }
      
      mat_cor <- data.frame(mat_cor)
      mat_pv <- data.frame(mat_pv)
      # names
      rownames(mat_cor) <- as.character(signatures[[1]][,1])
      colnames(mat_cor) <- cancers
      rownames(mat_pv) <- as.character(signatures[[1]][,1])
      colnames(mat_pv) <- cancers
      
      
      if (length(mat_cor[is.na(matrix(mat_cor))]) < length(mat_cor) & length(data.frame(mat_cor)) >= 2) {
        # plot heatmap
        pdf(file = paste(dir_o, 'cor mean master/', sign_name, '.pdf', sep = ''), 
            width = 12, height = 10)
        hm <- heatmap.2(
          data.matrix(mat_cor),
          #main=paste("1. Mean -", m),
          Rowv = NA, Colv = NA,
          dendrogram = "none",
          col = bluered,
          key = key, symkey = skey, keysize = keysize,
          cexRow = fontsize, cexCol = fontsize,
          density.info = "none", trace = "none",
          margins = margins,
          breaks = seq(-1, 1, length.out = 101)
        )
        dev.off()
        graphics.off()
        # write file
        write.table(mat_cor, file = paste(dir_o, 'cor mean master/',
                                          sign_name, '_cor.txt', sep = ''), sep = '\t')
        mat_cor <- rm_onlyNAcol(mat_cor)
        mat_pv <- rm_onlyNAcol(mat_pv)
        
        if (length(mat_cor[is.na(matrix(mat_cor))]) < length(mat_cor) &  length(data.frame(mat_cor)) >= 2) {
          # plot heatmap (clustered)
          pdf(file = paste(dir_o, 'cor mean master/clustered/', sign_name, '.pdf', 
                           sep = ''), width = 12, height = 10)
          hm <- heatmap.2(
            data.matrix(mat_cor),
            Colv = NA,
            #main=paste("1. Mean -", m),
            col = bluered,
            key = key, symkey = skey, keysize = keysize,
            cexRow = fontsize, cexCol = fontsize,
            density.info = "none", trace = "none",
            margins = margins,
            breaks = seq(-1, 1, length.out = 101)
          )
          dev.off()
          graphics.off()
          # write file
          write.table(mat_pv, file = paste(dir_o, 'cor mean master/', 
                                           sign_name, '_pv.txt', sep = ''), sep = '\t')
        }
        
      } else {
        "Skip plotting since there is no values to plot"
      }
    }
    
    
  }
  
  ###############################################################################################
  # xCell
  ###############################################################################################
  
  import_xcell <- function() {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    xc_raw <- read.table(file = paste(dir_p, 'xCell_TCGA_RSEM.txt', sep = ''), 
                         sep = '\t', header = T, row.names = 1)
    names <- colnames(xc_raw)
    for (i in 1:length(names)) {
      split <- strsplit(names[i], '\\.')
      names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep='.')
    }
    colnames(xc_raw) <- names
    xc <- list('vector', length = length(cancers))
    for (i in 1:length(cancers)) {
      cn <- colnames(exp_t_all[[i]])
      common <- intersect(cn, colnames(xc_raw))
      entry <- xc_raw[,common]
      xc[[i]] <- entry
    }
    return(xc)
  }
  
  xcell <- function(cancer_index) {
    
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    
    cancer <- cancers[cancer_index]
    
    # plot heatmap
    
    pdf(file = paste(dir_o, 'xCell heatmaps/', cancer, '.pdf', sep = ''), 
        width = 12, height = 12)
    hm <- heatmap.2(
      data.matrix(xc[[cancer_index]]),
      #main=paste("1. Mean -", m),
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = margins,
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    # row means for further analysis
    if (um) {
      xc_rm <- rowMeans(xc[[cancer_index]])
    } else {
      xc_rm <- rowMedians(xc[[cancer_index]])
    }
    
    # Prepare for correlation analysis
    cancer <- cancers[cancer_index]
    exp_t <- exp_t_all[[cancer_index]]
    common <- intersect(colnames(xc[[cancer_index]]), colnames(exp_t))
    exp_t <- exp_t[,common]
    
    # APOBEC vs xCell
    exp_t_sub <- exp_t[as.character(signatures[[1]][,1]),]
    exp_t_master <- xc[[cancer_index]]
    
    # correlation analysis (spearman)
    table_cor <- matrix(nrow = nrow(exp_t_master), ncol = nrow(exp_t_sub))
    table_pv <- matrix(nrow = nrow(exp_t_master), ncol = nrow(exp_t_sub))
    for (i in 1:nrow(exp_t_master)){
      for (j in 1:nrow(exp_t_sub)){
        result = tryCatch({
          test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(exp_t_sub[j,]),
                           method = 'spearman')
          table_cor[i,j] <- test$estimate
          table_pv[i,j] <- test$p.value
        }, error = function(e) {
          table_cor[i,j] <- NA
          table_pv[i,j] <- NA
        }, finally = {
          
        })
      }
    }
    
    # remove insignificant values
    table_cor[table_pv >= 0.05] <- 0
    
    # naming
    rownames(table_cor) <- rownames(xc[[cancer_index]])
    colnames(table_cor) <- rownames(exp_t_sub)
    rownames(table_pv) <- rownames(xc[[cancer_index]])
    colnames(table_pv) <- rownames(exp_t_sub)
    
    # plot heatmap
    
    pdf(file = paste(dir_o, 'xCell cor/', cancer, '_', sig_names[1], '.pdf',
                     sep = ''), width = 6, height = 16)
    hm <- heatmap.2(
      data.matrix(table_cor),
      #main=paste("1. Mean -", m),
      Rowv = F, Colv = F,
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    pdf(file = paste(dir_o, 'xCell cor/unclustered/', cancer, '_', sig_names[1], 
                     '.pdf', sep = ''), width = 6, height = 16)
    hm <- heatmap.2(
      data.matrix(table_cor),
      #main=paste("1. Mean -", m),
      Rowv = F, Colv = F,
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    
    # write tables
    write.table(table_cor, file = paste(dir_o, 'xCell cor/', cancer, '_', 
                                        sig_names[1], '_cor.txt', sep = ''), sep = '\t')
    write.table(table_pv, file = paste(dir_o, 'xCell cor/', cancer, '_',
                                       sig_names[1], '_pv.txt', sep = ''), sep = '\t')
    
    # get mean of signatures (for downstream analysis)
    table_means <- matrix(nrow = 2, ncol = ncol(exp_t_sub))
    if (um){
      table_means[1,] <- colMeans(exp_t_sub, na.rm=T)
    } else {
      table_means[1,] <- colMedians(exp_t_sub, na.rm=T)
    }
    
    # Hypoxia vs xCell
    exp_t_sub <- exp_t[as.character(signatures[[2]][,1]),]
    exp_t_master <- xc[[cancer_index]]
    
    # correlation analysis (spearman)
    table_cor <- matrix(nrow = nrow(exp_t_master), ncol = nrow(exp_t_sub))
    table_pv <- matrix(nrow = nrow(exp_t_master), ncol = nrow(exp_t_sub))
    for (i in 1:nrow(exp_t_master)){
      for (j in 1:nrow(exp_t_sub)){
        result = tryCatch({
          test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(exp_t_sub[j,]),
                           method = 'spearman')
          table_cor[i,j] <- test$estimate
          table_pv[i,j] <- test$p.value
        }, error = function(e) {
          table_cor[i,j] <- NA
          table_pv[i,j] <- NA
        }, finally = {
          
        })
      }
    }
    
    # remove insignificant values
    table_cor[table_pv >= 0.05] <- 0
    
    # naming
    rownames(table_cor) <- rownames(xc[[cancer_index]])
    colnames(table_cor) <- rownames(exp_t_sub)
    rownames(table_pv) <- rownames(xc[[cancer_index]])
    colnames(table_pv) <- rownames(exp_t_sub)
    
    # plot heatmap
    pdf(file = paste(dir_o, 'xCell cor/', cancer, '_', sig_names[2], '.pdf',
                     sep = ''), width = 12, height = 16)
    hm <- heatmap.2(
      data.matrix(table_cor),
      #main=paste("1. Mean -", m),
      Rowv = F, Colv = F,
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    pdf(file = paste(dir_o, 'xCell cor/unclustered/', cancer, '_', sig_names[2],
                     '.pdf', sep = ''), width = 12, height = 16)
    hm <- heatmap.2(
      data.matrix(table_cor),
      #main=paste("1. Mean -", m),
      Rowv = F, Colv = F,
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    # write tables
    write.table(table_cor, file = paste(dir_o, 'xCell cor/', cancer, '_', 
                                        sig_names[2], '_cor.txt', sep = ''), sep = '\t')
    write.table(table_pv, file = paste(dir_o, 'xCell cor/', cancer, '_',
                                       sig_names[2], '_pv.txt', sep = ''), sep = '\t')
    
    # get mean of signatures (for downstream analysis)
    if (um) {
      table_means[2,] <- colMeans(exp_t_sub, na.rm = T)
    } else {
      table_means[2,] <- colMedians(exp_t_sub, na.rm = T)
    }
    
    # correlate means with master signature
    table_cor_means <- matrix(nrow = nrow(exp_t_master), ncol = nrow(table_means))
    table_pv_means <- matrix(nrow = nrow(exp_t_master), ncol = nrow(table_means))
    for (i in 1:nrow(exp_t_master)) {
      for (j in 1:nrow(table_means)) {
        result = tryCatch({
          test <- cor.test(as.numeric(exp_t_master[i,]), as.numeric(table_means[j,]),
                           method = 'spearman')
          table_cor_means[i,j] <- test$estimate
          table_pv_means[i,j] <- test$p.value
        }, error = function(e) {
          table_cor_means[i,j] <- NA
          table_pv_means[i,j] <- NA
        }, finally = {
          
        })
      }
    }
    
    
    
    # remove insignificant values from means
    table_cor_means[table_pv_means >= 0.05] <- NA
    
    # naming means
    rownames(table_cor_means) <- rownames(xc[[cancer_index]])
    rownames(table_pv_means) <- rownames(xc[[cancer_index]])
    colnames(table_cor_means) <- c(sig_names[1], sig_names[2])
    colnames(table_pv_means) <- c(sig_names[1], sig_names[2])
    
    out <- list(table_cor_means, table_pv_means, xc_rm)
    return(out)
    
  }
  
  xcell_all <- function(xc_out) {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    
    # get means data from previous analysis
    score_means <- matrix(nrow = length(xc_out[[1]][[3]]), ncol = length(xc_out))
    for (i in 1:length(xc_out)) {
      score_means[,i] <- xc_out[[i]][[3]]
    }
    rownames(score_means) <- rownames(xc[[1]])
    colnames(score_means) <- cancers
    
    # remove insignificant values
    score_means[is.na(score_means)] <- 0
    
    # file naming
    
    file.name <- c()
    if (um) {
      file.name <- 'xCell cor all/mean_xcell_scores.pdf'
    } else {
      file.name <- 'xCell cor all/median_xcell_scores.pdf'
    }
    
    # graph
    pdf(file = paste(dir_o, file.name, sep = ''), width =6, height = 18)
    hm <- heatmap.2(
      score_means,
      #main=paste("1. Mean -", m),
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    score_means <- matrix(nrow = length(as.numeric(xc_out[[1]][[1]][,1])),
                          ncol = length(xc_out))
    for (i in 1:length(xc_out)) {
      score_means[,i] <- as.numeric(xc_out[[i]][[1]][,1])
    }
    rownames(score_means) <- rownames(xc[[i]])
    colnames(score_means) <- cancers
    score_means[is.na(score_means)] <- 0
    pdf(file = paste(dir_o, 'xCell cor all/cor_', 
                     sig_names[1], '.pdf', sep = ''),
        width =6, height = 18)
    hm <- heatmap.2(
      score_means,
      #main=paste("1. Mean -", m),
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    pdf(file = paste(dir_o, 'xCell cor all/unclustered/cor_', sig_names[1], 
                     '.pdf', sep = ''), width =6, height = 18)
    hm <- heatmap.2(
      score_means,
      #main=paste("1. Mean -", m),
      Rowv = F, Colv = F,
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace="none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out=101)
    )
    dev.off()
    graphics.off()
    score_means <- matrix(nrow = length(as.numeric(xc_out[[1]][[1]][,2])), 
                          ncol = length(xc_out))
    for (i in 1:length(xc_out)){
      score_means[,i] <- as.numeric(xc_out[[i]][[1]][,2])
    }
    rownames(score_means) <- rownames(xc[[i]])
    colnames(score_means) <- cancers
    score_means[is.na(score_means)] <- 0
    pdf(file = paste(dir_o, 'xCell cor all/cor_', sig_names[2], 
                     '.pdf', sep = ''), width = 6, height = 18)
    hm <- heatmap.2(
      score_means,
      #main=paste("1. Mean -", m),
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
    pdf(file = paste(dir_o, 'xCell cor all/unclustered/cor_', sig_names[2],
                     '.pdf', sep = ''), width = 6, height = 18)
    hm <- heatmap.2(
      score_means,
      #main=paste("1. Mean -", m),
      Rowv = F, Colv = F,
      col = bluered,
      key = key, symkey = skey, keysize = keysize,
      cexRow = fontsize, cexCol = fontsize,
      density.info = "none", trace = "none",
      margins = c(10,20),
      breaks = seq(-1, 1, length.out = 101)
    )
    dev.off()
    graphics.off()
  }
  
  ###############################################################################################
  # Heatmaps of correlation of every gene in master signature with others + xCell
  ###############################################################################################
  
  
  hm_master <- function(master_index) {
    require(doParallel)
    require(matrixStats)
    require(gplots)
    require(ggplot2)
    require(ggpubr) 
    require(ggpmisc)
    require(EnvStats)
    ###################
    
    this.gene <- as.character(signatures[[1]][master_index,1])
    
    # for each signature
    for (i in 2:length(signatures)) {
      
      # matrix to store correlation values
      mat.cor <- matrix(nrow = nrow(signatures[[i]]) + 1, ncol = length(cancers) + 1)
      mat.pv <- matrix(nrow = nrow(signatures[[i]]) + 1, ncol = length(cancers) + 1)
      if (um) {
        rownames(mat.cor) <- c(as.character(signatures[[i]][,1]), 'Signature mean')
        rownames(mat.pv) <- c(as.character(signatures[[i]][,1]), 'Signature mean')
        colnames(mat.cor) <- c(cancers, 'Pan-cancer mean')
        colnames(mat.pv) <- c(cancers, 'Pan-cancer mean')
      } else {
        rownames(mat.cor) <- c(as.character(signatures[[i]][,1]), 'Signature median')
        rownames(mat.pv) <- c(as.character(signatures[[i]][,1]), 'Signature median')
        colnames(mat.cor) <- c(cancers, 'Pan-cancer median')
        colnames(mat.pv) <- c(cancers, 'Pan-cancer median')
      }
      
      # for each cancer
      for (cancer_index in 1:length(cancers) + 1) {
        
        # get the expression for this master gene in one cancer
        # exception: the last index refers to across cancers
        exp_this.gene <- NULL
        if (cancer_index <= length(cancers)) {
          exp_this.gene <- as.numeric(exp_t_all[[cancer_index]][this.gene,])
        } else {
          exp_this.gene <- as.numeric(exp_t_full[this.gene,])
        }
        
        # for each gene in the signature
        for (j in 1:nrow(signatures[[i]])) {
          
          # get the expression for the gene
          sig.gene <- as.character(signatures[[i]][j,1])
          if (cancer_index <= length(cancers)) {
            exp_sig.gene <- as.numeric(exp_t_all[[cancer_index]][sig.gene,])
          } else {
            exp_sig.gene <- as.numeric(exp_t_full[sig.gene,])
          }
          
          # correlate the expression of master gene with signature gene
          result = tryCatch({
            test <- cor.test(exp_this.gene, exp_sig.gene, method = 'spearman')
            mat.cor[j,cancer_index] <- test$estimate
            mat.pv[j,cancer_index] <- test$p.value
          }, error = function(e) {
            mat.cor[j,cancer_index] <- NA
            mat.pv[j, cancer_index] <- NA
          }, finally = {
            
          })
        }
        
        # mean of the signature
        sigmean <- c()
        if (um) {
          if (cancer_index <= length(cancers)) { # individual cancers
            sigmean <- colMeans(exp_t_all[[cancer_index]][as.character(signatures[[i]][,1]),])
          } else { 
            # across cancers
            sigmean <- colMeans(exp_t_full[as.character(signatures[[i]][,1]),])
          }
        } else {
          if (cancer_index <= length(cancers)) {
            sigmean <- colMedians(exp_t_all[[cancer_index]][as.character(signatures[[i]][,1]),])
          } else {
            sigmean <- colMedians(exp_t_full[as.character(signatures[[i]][,1]),])
          }
        }
        sigmean <- as.numeric(sigmean)
        
        df <- data.frame(exp_this.gene, sigmean)
        
        df <- rm_onlyNArow(df)
        df <- df[rowSums(is.na(df)) == 2,]
        if (length(df[,1]) == 0) {
          print("skip correlation plot since there is not enough samples")
        } else {
          # correlate the expressions using the mean of the signature
          
          test <- cor.test(exp_this.gene, sigmean, method = 'spearman')
          mat.cor[j + 1,cancer_index] <- test$estimate
          mat.pv[j + 1,cancer_index] <- test$p.value
        }
        
        
        
      }
      
      # filter out p>0.05
      mat.cor.filtered <- mat.cor
      mat.cor.filtered[mat.pv > thr] <- 0
      mat.cor.filtered[is.na(mat.cor.filtered)] <- 0
      
      # plot heatmaps (unclustered)
      pdf(file = paste(dir_o, 'heatmaps master individual genes/',
                       this.gene, '_', sig_names[i], '.pdf', sep = ''),
          width = 6, height = 18)
      hm <- heatmap.2(
        mat.cor.filtered,
        Rowv = F, Colv = F,
        col = redgreen,
        key = key, symkey = skey, keysize = keysize,
        cexRow = fontsize, cexCol = fontsize,
        density.info = "none", trace = "none",
        margins = c(10,10),
        breaks = seq(-1, 1, length.out = 101)
      )
      dev.off()
      graphics.off()
      # plot heatmaps (signature clustered)
      pdf(file = paste(dir_o, 'heatmaps master individual genes/signature clustered/',
                       this.gene, '_', sig_names[i], '.pdf', sep = ''), width = 6, height = 18)
      hm <- heatmap.2(
        mat.cor.filtered,
        Colv = NA,
        col = redgreen,
        key = key, symkey = skey, keysize = keysize,
        cexRow = fontsize, cexCol = fontsize,
        density.info = "none", trace = "none",
        margins = c(10,10),
        breaks = seq(-1, 1, length.out = 101)
      )
      dev.off()
      graphics.off()
      # plot heatmaps (fully clustered)
      pdf(file = paste(dir_o, 'heatmaps master individual genes/full clustered/',
                       this.gene, '_', sig_names[i], '.pdf', sep = ''), 
          width = 6, height = 18)
      hm <- heatmap.2(
        mat.cor.filtered,
        col = redgreen,
        key = key, symkey = skey, keysize = keysize,
        cexRow = fontsize, cexCol = fontsize,
        density.info = "none", trace = "none",
        margins = c(10,10),
        breaks = seq(-1, 1, length.out = 101)
      )
      dev.off()
      graphics.off()
      # output tables
      write.table(mat.cor, file = paste(dir_o, 'heatmaps master individual genes/',
                                        this.gene, '_', sig_names[i], '_cor.txt', sep = ''))
      write.table(mat.pv, file = paste(dir_o, 'heatmaps master individual genes/',
                                       this.gene, '_', sig_names[i], '_pv.txt', sep = ''))
      
    }
    
  }

##############################################################################################
# Main function
###############################################################################################

if (win) { 
  ### WINDOWS ##########################################################################
  
  # Import parameters
  print('Importing parameters')
  cancers <- as.character(read.table(paste(dir_p,
                                           'cancer_list (ALL cancers).txt',
                                           sep = ''),
                                     header = T, sep = '\t')[,1])
  
  signatures <- vector('list', nts + 1)
  signatures[[1]] <- read.table(paste(dir_p, 'sign_master.txt', sep = ''), 
                                header = T, sep = '\t')
  sig_names <- as.character(colnames(signatures[[1]]))
  
  for (i in 1:nts) {
    signatures[[i + 1]] <- read.table(paste(dir_p, 'sign_', i, '.txt', sep = ''),
                                      header = T, sep = '\t')
    sig_names <- c(sig_names, as.character(colnames(signatures[[i + 1]])))
  }
  
  
  # gene signature to associate the master sig with (corr test)
  # by default it's hypoxia
  # which_desired_sig_vs_master <- as.character(signatures[[2]][,1]) 
  
  which_desired_sig_vs_master <- as.character(signatures[[2]][,1]) 
  
  # Import expression matrices (parallel)
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Importing expression matrices')
  cores <- min(detectCores(), length(cancers))
  if ("BRCA1_mut" %in% cancers|"OV" %in% cancers|
      "BRCA CDKN1_mut" %in% cancers|"BLCA CDKN1_mut" %in% cancers| 
      "LUAD CDKN1_mut" %in% cancers|"COAD CDKN1_mut" %in% cancers|
      "LUSC CDKN1_mut" %in% cancers|"PRAD CDKN1_mut" %in% cancers|
      "STAD CDKN1_mut" %in% cancers|"UCEC CDKN1_mut" %in% cancers) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    exp_n_all <- foreach(i = 1:(length(cancers) - 1)) %dopar% exp_import(cancers[i], 'n')
    stopCluster(cl)
  } else if ("CESC CDKN1_mut" %in% cancers| "THCA CDKN1_mut" %in% cancers) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    exp_n_all <- foreach(i = c(1,3)) %dopar% exp_import(cancers[i], 'n')
    stopCluster(cl)
  } else if ("KIRP p53_mut" %in% cancers) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    if ("KIRP p53_mut" == cancers[2]) {
      exp_n_all <- foreach(i = 1:(length(cancers) - 1)) %dopar% exp_import(cancers[i], 'n')
    } else {
      exp_n_all[c(1:5)] <- foreach(i = c(1:5)) %dopar% exp_import(cancers[i], 'n')
      exp_n_all[6] <- data.frame(matrix(data = NA, nrow = 20531, ncol = 1))
      exp_n_all[c(7:11)] <- foreach(i = c(7:11)) %dopar% exp_import(cancers[i], 'n')
      exp_n_all[12] <- data.frame(matrix(data = NA, nrow = 20531, ncol = 1))
      exp_n_all[c(13:21)] <- foreach(i = c(13:21)) %dopar% exp_import(cancers[i], 'n')
      exp_n_all[22] <- data.frame(matrix(data = NA, nrow = 20531, ncol = 1))
      exp_n_all[c(23:24)] <- foreach(i = c(23:(length(cancers) - 1))) %dopar% exp_import(cancers[i], 'n')
    }
    
    stopCluster(cl)
  } else if ("OV p53_mut" %in% cancers) {
    print ("skipping...no normal tissues")
  } else {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    exp_n_all <- foreach(i = 1:length(cancers)) %dopar% exp_import(cancers[i], 'n')
    stopCluster(cl)
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  exp_t_all <- foreach(i = 1:length(cancers)) %dopar% exp_import(cancers[i], 't')
  stopCluster(cl)
  gc()
  
  # Import xCell data
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Importing xCell data')
  xc <- import_xcell()
  gc()
  
  # Rearrange all samples across cancers in one matrix
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Pre-processing data')
  exp_t_full <- exp_t_all[[1]]
  for (i in 2:length(exp_t_all)) {
    exp_t_full <- cbind(exp_t_full, exp_t_all[[i]])
  }
  gc()
  
  # Median expression of all cancers
  exp_med <- matrix(nrow = length(cancers), ncol = 2)
  if ("BRCA1_mut" %in% cancers|"OV" %in% cancers|
      "BRCA CDKN1_mut" %in% cancers| 
      "BLCA CDKN1_mut" %in% cancers|
      "LUAD CDKN1_mut" %in% cancers|
      "COAD CDKN1_mut" %in% cancers|
      "LUSC CDKN1_mut" %in% cancers|
      "PRAD CDKN1_mut" %in% cancers|
      "STAD CDKN1_mut" %in% cancers|
      "UCEC CDKN1_mut" %in% cancers) {
    for (i in 1:(length(cancers) - 1)) {
      median_n <- exp_n_all[[i]]
      median_n[is.na(median_n)] <- 0
      exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)
    }
    
  } else if ("CESC CDKN1_mut" %in% cancers|"THCA CDKN1_mut" %in% cancers) {
    for (i in c(1,2)) {
      median_n <- exp_n_all[[i]]
      median_n[is.na(median_n)] <- 0
      if (i == 1) {
        exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)
      } else if (i == 2) {
        exp_med[3,2] <- median(as.matrix(median_n), na.rm = T)
      }
      
    }
  } else if ("KIRP p53_mut" %in% cancers) {
    if ("KIRP p53_mut" == cancers[2]) {
      for (i in 1:(length(cancers) - 1)) {
        median_n <- exp_n_all[[i]]
        median_n[is.na(median_n)] <- 0
        exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)
      }
    } else {
      for (i in c(1:5, 7:11, 13:21, 23:(length(cancers) - 1))) {
        median_n <- exp_n_all[[i]]
        median_n[is.na(median_n)] <- 0
        exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)
      }
    }
    
  } else if ("OV p53_mut" %in% cancers) {
    print("skipping... no normal tissues")
    } else {
      for (i in 1:length(cancers)) {
      median_n <- exp_n_all[[i]]
      median_n[is.na(median_n)] <- 0
      exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)}
    }
  for (i in 1:length(cancers)) {
    median_t <- exp_t_all[[i]]
    median_t[is.na(median_t)] <- 0
    exp_med[i,1] <- median(as.matrix(median_t), na.rm = T)
  }
  
  rownames(exp_med) <- cancers
  # colnames switched 
  colnames(exp_med) <- c('Tumour', 'Normal')
  write.table(exp_med, file = 'exp_med.txt', sep = '\t')
  gc()
 
  # Generate expression box plots/violin for every gene in all signatures, in all cancers
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating expression box/violin plots for genes in every signature')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  meps_all <- foreach(i = c(1:length(cancers))) %dopar% ebp_sign(i)
  
  stopCluster(cl)
  gc()
 
  # Generate expression boxplots/violin across cancers
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating expression box/violin plots across cancers')
  cores <- min(detectCores(), length(signatures))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i = 1:length(signatures)) %dopar% ebp_cancers(i)
  stopCluster(cl)
  gc()
# Generate expression heatmaps
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating expression heat maps')
  if ("BRCA1_mut" %in% cancers|"OV" %in% cancers|
      "BRCA CDKN1_mut" %in% cancers|
      "BLCA CDKN1_mut" %in% cancers|"LUAD CDKN1_mut" %in% cancers|
      "COAD CDKN1_mut" %in% cancers|"LUSC CDKN1_mut" %in% cancers|
      "PRAD CDKN1_mut" %in% cancers|"STAD CDKN1_mut" %in% cancers|
      "UCEC CDKN1_mut" %in% cancers) {
    cores <- min(detectCores(), length(cancers))
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    foreach(i = 1:(length(cancers) - 1)) %dopar% exp_hm(i)
    stopCluster(cl)
  } else if (("CESC CDKN1_mut" %in% cancers|
              "THCA CDKN1_mut" %in% cancers)) {
    cores <- min(detectCores(), length(cancers))
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    foreach(i = c(1,3)) %dopar% exp_hm(i)
    stopCluster(cl)
    
  } else if ("KIRP p53_mut" %in% cancers) {
    if ("KIRP p53_mut" == cancers[2]) {
      cores <- min(detectCores(), length(cancers))
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      foreach(i = 1:(length(cancers) - 1)) %dopar% exp_hm(i)
      stopCluster(cl)
    } else {
      cores <- min(detectCores(), length(cancers))
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      foreach(i = c(1:5, 7:11, 13:21, 23:(length(cancers) - 1))) %dopar% exp_hm(i)
      stopCluster(cl)
      }
    } else {
    cores <- min(detectCores(), length(cancers))
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    foreach(i = 1:length(cancers)) %dopar% exp_hm(i)
    stopCluster(cl)
    gc()
  }
  # Generate single gene correlation heatmaps
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating single gene correlation heatmaps')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  mean_cor <- foreach(i = 1:length(cancers)) %dopar% cor_hm(i)
  stopCluster(cl)
  gc()
  
  # Generate signature mean correlation heatmaps
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating signature mean correlation heatmaps')
  cor_mean_hm(mean_cor)
  gc()
  # xCell
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('xCell')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  xc_out <- foreach(i = 1:length(cancers)) %dopar% xcell(i)
  stopCluster(cl)
  xcell_all(xc_out)
  
  # Generate heatmaps of correlation of every gene in master signature with others
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating heatmaps of correlation of every gene in master signature with others')
  cores <- min(detectCores(), length(cancers))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i = 1:nrow(signatures[[1]])) %dopar% hm_master(i)
  stopCluster(cl)
  
  # Terminate
  print(Sys.time() - t_last)
  gc()
} else { ##### Mac OS ###########################################################################
  
  # Import parameters
  print('Importing parameters')
  cancers <- as.character(read.table(paste(dir_p, 'cancer_list.txt', sep = ''), 
                                     header = T, sep = '\t')[,1])
  signatures <- vector('list', nts+1)
  signatures[[1]] <- read.table(paste(dir_p, 'sign_master.txt', sep = ''),
                                header = T, sep = '\t')
  sig_names <- as.character(colnames(signatures[[1]]))
  for (i in 1:nts) {
    signatures[[i + 1]] <- read.table(paste(dir_p, 'sign_', i, '.txt', sep = ''), 
                                      header = T, sep = '\t')
    sig_names <- c(sig_names, as.character(colnames(signatures[[i + 1]])))
  }
  
  # Import expression matrices (parallel)
  print(Sys.time()- t_last)
  t_last <- Sys.time()
  print('Importing expression matrices')
  f <- vector("list", length = length(cancers))
  for (i in 1:length(cancers)) {
    f[[i]] <- mcparallel(exp_import(cancers[i], 'n'))
  }
  exp_n_all <- mccollect(f) # all normal exp matrices
  f <- vector("list", length = length(cancers))
  for (i in 1:length(cancers)) {
    f[[i]] <- mcparallel(exp_import(cancers[i], 't'))
  }
  exp_t_all <- mccollect(f) # tumour exp matrices
  
  # Import xCell data
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Importing xCell data')
  xc <- import_xcell()
  
  # Rearrange all samples across cancers in one matrix
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Pre-processing data')
  exp_t_full <- exp_t_all[[1]]
  for (i in 2:length(exp_t_all)) {
    exp_t_full <- cbind(exp_t_full, exp_t_all[[i]])
  }
  
  # Median expression of all cancers
  exp_med <- matrix(nrow = length(cancers), ncol = 2)
  for (i in 1:length(cancers)) {
    median_n <- exp_n_all[[i]]
    median_n[is.na(median_n)] <- 0
    median_t <- exp_t_all[[i]]
    median_t[is.na(median_t)] <- 0
    exp_med[i,1] <- median(as.matrix(median_t), na.rm = T)
    exp_med[i,2] <- median(as.matrix(median_n), na.rm = T)
  }
  rownames(exp_med) <- cancers
  colnames(exp_med) <- c('Normal', 'Tumour')
  write.table(exp_med, file = 'data/exp_med.txt', sep = '\t')
  
  # Generate expression boxplots for every gene in all signatures, in all cancers
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating expression box/violin plots for genes in every signature')
  f <- vector("list", length = length(cancers))
  for (i in 1:length(cancers)){
    f[[i]] <- mcparallel(ebp_sign(i))
  }
  meps_all <- mccollect(f)
  
  # Generate expression boxplots across cancers
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating expression box/violin plots across cancers')
  f <- vector("list", length = length(signatures))
  for (i in 1:length(signatures)){
    f[[i]] <- mcparallel(ebp_cancers(i))
  }
  mccollect(f)
  
  # Generate expression heatmaps
  print(Sys.time() - t_last)
  t_last <- Sys.time()
  print('Generating expression heat maps')
  f <- vector("list", length = length(cancers))
  for (i in 1:length(cancers)) {
    f[[i]] <- mcparallel(exp_hm(i))
  }
  mccollect(f)
  
  # Generate single gene correlation heatmaps
  print(Sys.time()- t_last)
  t_last <- Sys.time()
  print('Generating single gene correlation heatmaps')
  f <- vector("list", length = length(cancers))
  for (i in 1:length(cancers)) {
    f[[i]] <- mcparallel(cor_hm(i))
  }
  mean_cor <- mccollect(f)
  
  # Generate signature mean correlation heatmaps
  print(Sys.time()- t_last)
  t_last <- Sys.time()
  print('Generating signature mean correlation heatmaps')
  cor_mean_hm(mean_cor)
  
  # xCell
  print(Sys.time()- t_last)
  t_last <- Sys.time()
  print('xCell')
  f <- vector("list", length = length(cancers))
  for (i in 1:length(cancers)){
    f[[i]] <- mcparallel(xcell(i))
  }
  xc_out <- mccollect(f)
  xcell_all(xc_out)
  
  # Generate heatmaps of correlation of every gene in master signature with others
  print(Sys.time()- t_last)
  t_last <- Sys.time()
  print('Generating heatmaps of correlation of every gene in master signature with others')
  f <- vector("list", length = nrow(signatures[[1]]))
  for (i in 1:nrow(signatures[[1]])){
    f[[i]] <- mcparallel(hm_master(i))
  }
  mccollect(f)
  
  # Terminate
  print(Sys.time() - t_last)
}

###############################################################################################
# Termination
###############################################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)

setwd("G:/My Drive/Jiachen Files/codes_supp")
# save details of this R script
scriptname <- "signatures_modified_for_gdrive"
savehistory(file = paste(scriptname, ".Rhistory", sep = ""))
sink(file = paste(Sys.Date(), scriptname, ".txt", sep = ""),
     type = c("output", "message"))
# details of the packages and other useful info.
print(sessionInfo()) 
sink(NULL)

# save globalenv/workspace
save.image("signatures_modified_for_gdrive.Rdata")  

