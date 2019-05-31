# Sakura Maezono
##------ Mon Mar 25 16:17:15 2019 ------##
# Retrieve germline mutations from TCGA from GDC legacy

# Libraries used

# for modifying tables etc
library(tidyverse) 
library(TCGAbiolinks)


# All clinical annotations
clinical_data <- read.csv(paste("BRCA", "anno.csv", sep = "_"), header = F)
clinical_data <- clinical_data[,-1]
colnames(clincal_data) <- clinical_data[1,]
