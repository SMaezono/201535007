# What is this script for? 

setwd("G:/My Drive/Jiachen Files/")

root <- "G:/My Drive/Jiachen Files/GDAC/TCGA/"

dir <- paste(root, "BRCA", ".txt", sep="")

y.genes <- c('BPY2', 'DAZ1', 'DAZ2', 'DDX3Y', 'PRKY', 'RBMY1A1', 'SRY', 
             'TSPY1', 'USP9Y', 'UTY', 'ZFY')

input <- read.table(dir, sep="\t", header = TRUE, row.names=1)
input <- input[-1,]

good <- rownames(input)
for(i in 1:length(good)){
  good[i] <- strsplit(good[i], "\\|")[[1]][1]
}

index.y <- match(y.genes, good)
y <- input[index.y,]
rownames(y) <- y.genes

y <- matrix(as.numeric(as.matrix(y)), nrow=nrow(y), ncol=ncol(y))
rownames(y) <- y.genes

y <- log2(y)

y[is.infinite(y)] <- NA

med_na <- median(y, na.rm=TRUE)
mean_na <- mean(y, na.rm=TRUE)
sd_na <- sd(y, na.rm=TRUE)

med_na
mean_na
sd_na
mean_na+3*sd_na

y[is.na(y)] <- 0

med_zero <- median(y)
mean_zero <- mean(y)
sd_zero <- sd(y)
mean_zero + 3*sd_zero

med_zero
mean_zero

pos <- mean(y[y>0])
pos
