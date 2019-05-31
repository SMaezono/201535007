# assembly.R

# Jiachen Liang
# 2018
# modified by Sakura Maezono

# IMPORTANT: Missing stromal, Immune, Epithelial information and
# Values of the hypoxia signature from the same patient is slighty different



# Buffa Laboratory
# Computational Biology and Integrative Genomics Group
# Department of Oncology
# University of Oxford

###############################################################################################

# Function to combine multiple forms of data from TCGA
# Currently implemented for BRCA only

# Inputs:
#   List of signatures
#   List of subtypes
#   TCGA expression data (tumour)
#   TCGA subtype data
#   Clinical data

# Outputs:
#   Data frame with all data

###############################################################################################
# Initiation
###############################################################################################

t_start <- Sys.time()

###############################################################################################
# Associated library functions
###############################################################################################

library(matrixStats)
library(survival)

###############################################################################################
# File locations
###############################################################################################

# Set working directory
root <- "G:/My Drive/Jiachen Files/data/"
setwd(root)

###############################################################################################
# Parameters
###############################################################################################

# Number of tested signatures
nts <- 4

###############################################################################################
# Main function
###############################################################################################

# Import parameters
print('Importing parameters')
subtypes <- read.table('subtypes.txt', header = F, sep = '\t')
signatures <- vector('list', nts+1)
signatures[[1]] <- read.table('sign_master.txt', header = T, sep = '\t')
for(i in 1:nts){
  signatures[[i+1]] <- read.table(paste('sign_', i, '.txt', sep = ''), header = T, sep = '\t')
}

# Import expression matrix
filename = paste("G:/My Drive/Jiachen Files/GDAC/", "BRCA", "_", "t", ".txt", sep = "")
exp <- read.table(filename, sep = '\t', header = TRUE, row.names = NULL)
exp_rownames <- exp[,1]
exp <- exp[,-1]
exp_rn_unique <- make.unique(exp_rownames, sep = '.')
rownames(exp) <- exp_rn_unique

# Subset expression matrix
combined_sign <- c()
for(i in 1:length(signatures)){
  for(j in 1:nrow(signatures[[i]])){
    combined_sign <- c(combined_sign, as.character(signatures[[i]][j,1]))
  }
}
exp_sub <- exp[combined_sign,]

# Fix expression matrix names
exp_sample_names <- colnames(exp_sub)
exp_sample_names_fixed <- c()
for(i in 1:length(exp_sample_names)){
  split <- strsplit(exp_sample_names[i], '\\.')
  left <- paste(split[[1]][1], split[[1]][2], split[[1]][3], sep = '.')
  exp_sample_names_fixed <- c(exp_sample_names_fixed, left)
}
colnames(exp_sub) <- exp_sample_names_fixed

# Import clinical data
clin <- read.table('BRCA_clinical.txt', sep = '\t', header = T, row.names = 1)

# Import selected clinical variables
clin_var <- as.character(read.table('selected_clinical.txt', sep = '\t', header = F)[,1])

# Subset clinical data
clin_sub <- clin[clin_var,]

# Import subtype data
subtypes <- read.table('patient_annotation.txt', sep = '\t', header = T, row.names = 1)

# Import selected subtypes
sel_sub <- as.character(read.table('subtypes.txt', sep = '\t', header = F)[,1])

# Select for relevant subtypes
subtypes_sub <- subtypes[subtypes$PAM50.mRNA %in% sel_sub, ]

# Import list of mutation files
mut_manifest <- as.character(read.table('BRCA_mutations/MANIFEST.txt', 
                                        sep = ' ', header = F, row.names = NULL)[,2])
for(i in 1:length(mut_manifest)){
  n <- strsplit(mut_manifest[i], '\\.')[[1]][1]
  nx <- paste(strsplit(n, '-')[[1]][1], strsplit(n, '-')[[1]][2], strsplit(n, '-')[[1]][3], sep = '.')
  mut_manifest[i] <- nx
}

# Get common samples
samples_exp <- colnames(exp_sub)
samples_clin <- colnames(clin_sub)
samples_sub <- rownames(subtypes_sub)
samples_common <- intersect(intersect(samples_exp, samples_clin), intersect(samples_sub, mut_manifest))

# Get signature names
sign_names <- c()
for(i in 1:length(signatures)){
  sign_names <- c(sign_names, colnames(signatures[[i]]))
}

# Define parameters to present
params <- c(
  sign_names,
  'Subtype',
  as.character(clin_var),
  'Mutated APOBEC3A',
  'Mutated APOBEC3B',
  't',
  'e'
  )

# Construct initial matrix
mat <- matrix(nrow = length(samples_common), ncol = length(params))
rownames(mat) <- samples_common
colnames(mat) <- params
mat <- as.data.frame(mat)

# Fill expression
for(s in rownames(mat)){
  for(i in 1:length(signatures)){
    sign <- as.character(signatures[[i]][,1])
    mean_sign <- mean(exp_sub[sign,s], na.rm = T)
    mat[s,i] <- mean_sign
  }
}

# Fill subtype
for(s in rownames(mat)){
  mat[s,'Subtype'] <- as.character(subtypes_sub[s,'PAM50.mRNA'])
}

# Fill clinical
for(s in rownames(mat)){
  mat[s,nts+3] <- as.character(clin_sub[1,s])
  mat[s,nts+4] <- as.numeric(as.character(clin_sub[2,s]))
  mat[s,nts+5] <- as.numeric(as.character(clin_sub[3,s]))
  mat[s,nts+6] <- as.numeric(as.character(clin_sub[4,s]))
  mat[s,nts+7] <- as.character(clin_sub[5,s])
  mat[s,nts+8] <- as.character(clin_sub[6,s])
  mat[s,nts+9] <- as.character(clin_sub[7,s])
  mat[s,nts+10] <- as.character(clin_sub[8,s])
  if(!is.na(mat[s,nts+7])){
    if(mat[s,nts+7] == 'positive'){
      mat[s,nts+7] <- T
    }else{
      mat[s,nts+7] <- F
    }
  }
  if(!is.na(mat[s,nts+8])){
    if(mat[s,nts+8] == 'positive'){
      mat[s,nts+8] <- T
    }else{
      mat[s,nts+8] <- F
    }
  }
}

# Fill mutation data
for(s in rownames(mat)){
  sx <- paste(strsplit(s, '\\.')[[1]][1], strsplit(s, '\\.')[[1]][2], strsplit(s, '\\.')[[1]][3], sep = '-')
  mut_raw <- read.table(file = paste('BRCA_mutations/', sx, '-01.maf.txt', sep = ''), 
                        sep = '\t', fill = T, header = T, row.names  =  NULL, quote  =  "")
  mut_genes <- as.character(mut_raw$Hugo_Symbol)
  if('APOBEC3A' %in% mut_genes){
    mat[s,'Mutated APOBEC3A'] <- T
  }else{
    mat[s,'Mutated APOBEC3A'] <- F
  }
  if('APOBEC3B' %in% mut_genes){
    mat[s,'Mutated APOBEC3B'] <- T
  }else{
    mat[s,'Mutated APOBEC3B'] <- F
  }
}

# Fill time and event variables
for(s in rownames(mat)){
  if(mat[s,nts+3] == 'alive'){
    mat[s,'e'] <- F
    mat[s,'t'] <- mat[s,nts+5]
  }else if(mat[s,nts+3] == 'dead'){
    mat[s,'e'] <- T
    mat[s,'t'] <- mat[s,nts+4]
  }
}

# xCell

# Roberts mut sign

# Alexandrov mut sign

# Write table
write.table(mat, file = 'assembled_exp_clin_Sakura.txt', sep = '\t')

# coxph(Surv(t, e) ~ breast_carcinoma_estrogen_receptor_status, mat)
# survfit(Surv(t,e) ~ breast_carcinoma_estrogen_receptor_status, mat)

###############################################################################################
# Termination
###############################################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)
