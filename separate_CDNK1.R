# separate_subtypes.R
# Modified by Sakura Maezono ##------ Wed Mar 13 13:17:54 2019 ------##
# Separates cancer expression matrix by CDKN1A mutation status

setwd('G:/My Drive/Jiachen Files/GDAC/')

# HNSC, KIRP has no CDKN1 mutation
cancer <- c('BRCA', 'BLCA', 'CESC', 'COAD', 'LUAD', 'LUSC', 'LIHC', 'STAD',
            'KIRC','PRAD', 'OV', 'UCEC', 'THCA')
sep_subtypes <- function(cancercode) {
  if (cancercode == 'OV') {
    print ('No normal samples')
  } else {
    # Import normal
    input <- read.table(paste(cancercode,'n.txt', sep = "_"), sep = '\t',
                        header = T, row.names = NULL)
    rn <- input[,1]
    rn <- make.unique(rn, sep = '.')
    output <- input[,-1]
    rownames(output) <- rn
    names <- colnames(output)
    for (i in 1:length(names)) {
      split <- strsplit(names[i], '\\.')
      names[i] <- paste(split[[1]][1], split[[1]][2], split[[1]][3],
                        sep = '.')
    }
    colnames(output) <- names
    exp_n <- output
  }
  # Import tumour
  input <- read.table(paste(cancercode,'t.txt', sep = "_"), sep = '\t',
                      header = T, row.names = NULL)
  rn <- input[,1]
  rn <- make.unique(rn, sep = '.')
  output <- input[,-1]
  rownames(output) <- rn
  names <- colnames(output)
  for (i in 1:length(names)) {
    split <- strsplit(names[i], '\\.')
    names[i] <- paste(split[[1]][1], split[[1]][2],
                      split[[1]][3], sep = '.')
  }
  colnames(output) <- names
  exp_t <- output
  
  # Import subtype data
  
  dir_sub <- paste('G:/My Drive/Jiachen Files/data/', cancercode,
                   '_CDKN1.txt', sep='')
  subtypes <- read.table(dir_sub, sep = '\t', header = T, row.names = NULL)
  
  # if (cancercode == 'HNSC') {
  #   subtypes <- subtypes[,c(1,2,5)]
  #   colnames(subtypes) <- c('STUDY_ID', 'SAMPLE_ID', 'CDKN1A')
  # }
  
  mut <- as.character(subtypes$SAMPLE_ID[!is.na(subtypes$CDKN1A)])
  wt <- as.character(subtypes$SAMPLE_ID[is.na(subtypes$CDKN1A)])
  
  # Fix names
  mut <- substr(mut, 1, 12)
  wt <- substr(wt, 1, 12)
  mut <- gsub('-', '.', mut)
  wt <- gsub('-', '.', wt)
  
  # Separate matrices
  if (cancercode == 'OV') {
    print ('No normal samples')
  } else {
    mut_n <- exp_n[,intersect(colnames(exp_n), mut)]
    wt_n <- exp_n[,intersect(colnames(exp_n), wt)]
  }
    mut_t <- exp_t[,intersect(colnames(exp_t), mut)]
    wt_t <- exp_t[,intersect(colnames(exp_t), wt)]
  
  # Export
    if (cancercode == 'OV') {
      print ('No normal samples')
    } else {
  write.table(mut_n, file = paste(cancercode, 'CDKN1_mut_n.txt', sep = ' '),
              sep = '\t')
  write.table(wt_n, file = paste(cancercode, 'CDKN1_wt_n.txt', sep = ' '),
              sep = '\t')
    }
  write.table(mut_t, file = paste(cancercode, 'CDKN1_mut_t.txt', sep = ' '),
              sep = '\t')
  write.table(wt_t, file = paste(cancercode, 'CDKN1_wt_t.txt', sep = ' '), 
              sep = '\t')
}


lapply(cancer[1:length(cancer)],sep_subtypes)

sep_subtypes('CESC')
