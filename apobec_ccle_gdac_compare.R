# apobec_ccle_gdac_compare.R

# Jiachen Liang

# Correlate CCLE and TCGA data for APOBEC expression

t_start <- Sys.time()

# Root directory
root <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/output/ccle/"
setwd(root)

# TCGA directory
tcga_dir <- "/Users/jiachen/OneDrive - OnTheHub - The University of Oxford/Buffa Lab/apobec/data/stats/"

library(matrixStats)

###############################################################################

# Import

mean_bl <- read.table(file='mean_bl.txt', sep='\t', header=TRUE, row.names=1)
mean_fc <- read.table(file='mean_fc.txt', sep='\t', header=TRUE, row.names=1)
median_bl <- read.table(file='median_bl.txt', sep='\t', header=TRUE, row.names=1)
median_fc <- read.table(file='median_fc.txt', sep='\t', header=TRUE, row.names=1)

prad_t_a <- read.table(file=paste(tcga_dir, 'PRAD_t_a.txt', sep=''), sep='\t', header=TRUE, row.names=1)
brca_t_a <- read.table(file=paste(tcga_dir, 'BRCA_t_a.txt', sep=''), sep='\t', header=TRUE, row.names=1)
coad_t_a <- read.table(file=paste(tcga_dir, 'COAD_t_a.txt', sep=''), sep='\t', header=TRUE, row.names=1)

# Arrange
# format: prad/brca/coad _ mean/median _ bl/fc

prad_mean_bl <- matrix(nrow=nrow(mean_bl), ncol=2)
brca_mean_bl <- matrix(nrow=nrow(mean_bl), ncol=2)
coad_mean_bl <- matrix(nrow=nrow(mean_bl), ncol=2)

prad_mean_bl[,1] <- mean_bl[,1]
brca_mean_bl[,1] <- mean_bl[,2]
coad_mean_bl[,1] <- mean_bl[,4]
prad_mean_bl[,2] <- prad_t_a[,1]
brca_mean_bl[,2] <- brca_t_a[,1]
coad_mean_bl[,2] <- coad_t_a[,1]

prad_median_bl <- matrix(nrow=nrow(median_bl), ncol=2)
brca_median_bl <- matrix(nrow=nrow(median_bl), ncol=2)
coad_median_bl <- matrix(nrow=nrow(median_bl), ncol=2)

prad_median_bl[,1] <- median_bl[,1]
brca_median_bl[,1] <- median_bl[,2]
coad_median_bl[,1] <- median_bl[,4]
prad_median_bl[,2] <- prad_t_a[,2]
brca_median_bl[,2] <- brca_t_a[,2]
coad_median_bl[,2] <- coad_t_a[,2]

prad_mean_fc <- matrix(nrow=nrow(mean_fc), ncol=2)
brca_mean_fc <- matrix(nrow=nrow(mean_fc), ncol=2)
coad_mean_fc <- matrix(nrow=nrow(mean_fc), ncol=2)

prad_mean_fc[,1] <- mean_fc[,1]
brca_mean_fc[,1] <- mean_fc[,2]
coad_mean_fc[,1] <- mean_fc[,4]
prad_mean_fc[,2] <- prad_t_a[,1]
brca_mean_fc[,2] <- brca_t_a[,1]
coad_mean_fc[,2] <- coad_t_a[,1]

prad_median_fc <- matrix(nrow=nrow(median_fc), ncol=2)
brca_median_fc <- matrix(nrow=nrow(median_fc), ncol=2)
coad_median_fc <- matrix(nrow=nrow(median_fc), ncol=2)

prad_median_fc[,1] <- median_fc[,1]
brca_median_fc[,1] <- median_fc[,2]
coad_median_fc[,1] <- median_fc[,4]
prad_median_fc[,2] <- prad_t_a[,2]
brca_median_fc[,2] <- brca_t_a[,2]
coad_median_fc[,2] <- coad_t_a[,2]

# Output matrices

cor.matrix.p <- matrix(nrow=3, ncol=4)
rownames(cor.matrix.p) <- c('Prostate', 'Breast', 'Colon')
colnames(cor.matrix.p) <- c('FC (mean)', 'FC (median)', 'BL (mean)', 'BL (median)')
cor.matrix.s <- matrix(nrow=3, ncol=4)
rownames(cor.matrix.s) <- c('Prostate', 'Breast', 'Colon')
colnames(cor.matrix.s) <- c('FC (mean)', 'FC (median)', 'BL (mean)', 'BL (median)')

pv.matrix.p <- matrix(nrow=3, ncol=4)
rownames(pv.matrix.p) <- c('Prostate', 'Breast', 'Colon')
colnames(pv.matrix.p) <- c('FC (mean)', 'FC (median)', 'BL (mean)', 'BL (median)')
pv.matrix.s <- matrix(nrow=3, ncol=4)
rownames(pv.matrix.s) <- c('Prostate', 'Breast', 'Colon')
colnames(pv.matrix.s) <- c('FC (mean)', 'FC (median)', 'BL (mean)', 'BL (median)')

# PRAD

cor.matrix.p[1,1] <- cor.test(as.numeric(prad_mean_fc[,1]), as.numeric(prad_mean_fc[,2]), method='pearson')$estimate
pv.matrix.p[1,1] <- cor.test(as.numeric(prad_mean_fc[,1]), as.numeric(prad_mean_fc[,2]), method='pearson')$p.value
cor.matrix.s[1,1] <- cor.test(as.numeric(prad_mean_fc[,1]), as.numeric(prad_mean_fc[,2]), method='spearman')$estimate
pv.matrix.s[1,1] <- cor.test(as.numeric(prad_mean_fc[,1]), as.numeric(prad_mean_fc[,2]), method='spearman')$p.value

cor.matrix.p[1,2] <- cor.test(as.numeric(prad_median_fc[,1]), as.numeric(prad_median_fc[,2]), method='pearson')$estimate
pv.matrix.p[1,2] <- cor.test(as.numeric(prad_median_fc[,1]), as.numeric(prad_median_fc[,2]), method='pearson')$p.value
cor.matrix.s[1,2] <- cor.test(as.numeric(prad_median_fc[,1]), as.numeric(prad_median_fc[,2]), method='spearman')$estimate
pv.matrix.s[1,2] <- cor.test(as.numeric(prad_median_fc[,1]), as.numeric(prad_median_fc[,2]), method='spearman')$p.value

cor.matrix.p[1,3] <- cor.test(as.numeric(prad_mean_bl[,1]), as.numeric(prad_mean_bl[,2]), method='pearson')$estimate
pv.matrix.p[1,3] <- cor.test(as.numeric(prad_mean_bl[,1]), as.numeric(prad_mean_bl[,2]), method='pearson')$p.value
cor.matrix.s[1,3] <- cor.test(as.numeric(prad_mean_bl[,1]), as.numeric(prad_mean_bl[,2]), method='spearman')$estimate
pv.matrix.s[1,3] <- cor.test(as.numeric(prad_mean_bl[,1]), as.numeric(prad_mean_bl[,2]), method='spearman')$p.value

cor.matrix.p[1,4] <- cor.test(as.numeric(prad_median_bl[,1]), as.numeric(prad_median_bl[,2]), method='pearson')$estimate
pv.matrix.p[1,4] <- cor.test(as.numeric(prad_median_bl[,1]), as.numeric(prad_median_bl[,2]), method='pearson')$p.value
cor.matrix.s[1,4] <- cor.test(as.numeric(prad_median_bl[,1]), as.numeric(prad_median_bl[,2]), method='spearman')$estimate
pv.matrix.s[1,4] <- cor.test(as.numeric(prad_median_bl[,1]), as.numeric(prad_median_bl[,2]), method='spearman')$p.value

# BRCA

cor.matrix.p[2,1] <- cor.test(as.numeric(brca_mean_fc[,1]), as.numeric(brca_mean_fc[,2]), method='pearson')$estimate
pv.matrix.p[2,1] <- cor.test(as.numeric(brca_mean_fc[,1]), as.numeric(brca_mean_fc[,2]), method='pearson')$p.value
cor.matrix.s[2,1] <- cor.test(as.numeric(brca_mean_fc[,1]), as.numeric(brca_mean_fc[,2]), method='spearman')$estimate
pv.matrix.s[2,1] <- cor.test(as.numeric(brca_mean_fc[,1]), as.numeric(brca_mean_fc[,2]), method='spearman')$p.value

cor.matrix.p[2,2] <- cor.test(as.numeric(brca_median_fc[,1]), as.numeric(brca_median_fc[,2]), method='pearson')$estimate
pv.matrix.p[2,2] <- cor.test(as.numeric(brca_median_fc[,1]), as.numeric(brca_median_fc[,2]), method='pearson')$p.value
cor.matrix.s[2,2] <- cor.test(as.numeric(brca_median_fc[,1]), as.numeric(brca_median_fc[,2]), method='spearman')$estimate
pv.matrix.s[2,2] <- cor.test(as.numeric(brca_median_fc[,1]), as.numeric(brca_median_fc[,2]), method='spearman')$p.value

cor.matrix.p[2,3] <- cor.test(as.numeric(brca_mean_bl[,1]), as.numeric(brca_mean_bl[,2]), method='pearson')$estimate
pv.matrix.p[2,3] <- cor.test(as.numeric(brca_mean_bl[,1]), as.numeric(brca_mean_bl[,2]), method='pearson')$p.value
cor.matrix.s[2,3] <- cor.test(as.numeric(brca_mean_bl[,1]), as.numeric(brca_mean_bl[,2]), method='spearman')$estimate
pv.matrix.s[2,3] <- cor.test(as.numeric(brca_mean_bl[,1]), as.numeric(brca_mean_bl[,2]), method='spearman')$p.value

cor.matrix.p[2,4] <- cor.test(as.numeric(brca_median_bl[,1]), as.numeric(brca_median_bl[,2]), method='pearson')$estimate
pv.matrix.p[2,4] <- cor.test(as.numeric(brca_median_bl[,1]), as.numeric(brca_median_bl[,2]), method='pearson')$p.value
cor.matrix.s[2,4] <- cor.test(as.numeric(brca_median_bl[,1]), as.numeric(brca_median_bl[,2]), method='spearman')$estimate
pv.matrix.s[2,4] <- cor.test(as.numeric(brca_median_bl[,1]), as.numeric(brca_median_bl[,2]), method='spearman')$p.value

# COAD

cor.matrix.p[3,1] <- cor.test(as.numeric(coad_mean_fc[,1]), as.numeric(coad_mean_fc[,2]), method='pearson')$estimate
pv.matrix.p[3,1] <- cor.test(as.numeric(coad_mean_fc[,1]), as.numeric(coad_mean_fc[,2]), method='pearson')$p.value
cor.matrix.s[3,1] <- cor.test(as.numeric(coad_mean_fc[,1]), as.numeric(coad_mean_fc[,2]), method='spearman')$estimate
pv.matrix.s[3,1] <- cor.test(as.numeric(coad_mean_fc[,1]), as.numeric(coad_mean_fc[,2]), method='spearman')$p.value

cor.matrix.p[3,2] <- cor.test(as.numeric(coad_median_fc[,1]), as.numeric(coad_median_fc[,2]), method='pearson')$estimate
pv.matrix.p[3,2] <- cor.test(as.numeric(coad_median_fc[,1]), as.numeric(coad_median_fc[,2]), method='pearson')$p.value
cor.matrix.s[3,2] <- cor.test(as.numeric(coad_median_fc[,1]), as.numeric(coad_median_fc[,2]), method='spearman')$estimate
pv.matrix.s[3,2] <- cor.test(as.numeric(coad_median_fc[,1]), as.numeric(coad_median_fc[,2]), method='spearman')$p.value

cor.matrix.p[3,3] <- cor.test(as.numeric(coad_mean_bl[,1]), as.numeric(coad_mean_bl[,2]), method='pearson')$estimate
pv.matrix.p[3,3] <- cor.test(as.numeric(coad_mean_bl[,1]), as.numeric(coad_mean_bl[,2]), method='pearson')$p.value
cor.matrix.s[3,3] <- cor.test(as.numeric(coad_mean_bl[,1]), as.numeric(coad_mean_bl[,2]), method='spearman')$estimate
pv.matrix.s[3,3] <- cor.test(as.numeric(coad_mean_bl[,1]), as.numeric(coad_mean_bl[,2]), method='spearman')$p.value

cor.matrix.p[3,4] <- cor.test(as.numeric(coad_median_bl[,1]), as.numeric(coad_median_bl[,2]), method='pearson')$estimate
pv.matrix.p[3,4] <- cor.test(as.numeric(coad_median_bl[,1]), as.numeric(coad_median_bl[,2]), method='pearson')$p.value
cor.matrix.s[3,4] <- cor.test(as.numeric(coad_median_bl[,1]), as.numeric(coad_median_bl[,2]), method='spearman')$estimate
pv.matrix.s[3,4] <- cor.test(as.numeric(coad_median_bl[,1]), as.numeric(coad_median_bl[,2]), method='spearman')$p.value

# Write

write.table(cor.matrix.p, file='cor_pearson.txt', sep='\t')
write.table(cor.matrix.s, file='cor_spearman.txt', sep='\t')
write.table(pv.matrix.p, file='pv_pearson.txt', sep='\t')
write.table(pv.matrix.s, file='pv_spearman.txt', sep='\t')

# Scatter plots

rownames(prad_mean_bl) <- rownames(mean_bl)
colnames(prad_mean_bl) <- c('CCLE', 'TCGA')
prad_mean_bl <- as.data.frame(prad_mean_bl)
plot <- ggplot(prad_mean_bl, aes(x=CCLE, y=TCGA)) +
  geom_point()
pdf(file="prad_mean_bl.pdf", width=12, height=9)
print(plot)
dev.off()

rownames(brca_mean_bl) <- rownames(mean_bl)
colnames(brca_mean_bl) <- c('CCLE', 'TCGA')
brca_mean_bl <- as.data.frame(brca_mean_bl)
plot <- ggplot(brca_mean_bl, aes(x=CCLE, y=TCGA)) +
  geom_point() +
  geom_text(aes(label=rownames(brca_mean_bl)),hjust=0, vjust=0)
pdf(file="brca_mean_bl.pdf", width=12, height=9)
print(plot)
dev.off()

rownames(coad_mean_bl) <- rownames(mean_bl)
colnames(coad_mean_bl) <- c('CCLE', 'TCGA')
coad_mean_bl <- as.data.frame(coad_mean_bl)
plot <- ggplot(coad_mean_bl, aes(x=CCLE, y=TCGA)) +
  geom_point() +
  geom_text(aes(label=rownames(brca_mean_bl)),hjust=0, vjust=0)
pdf(file="coad_mean_bl.pdf", width=12, height=9)
print(plot)
dev.off()

rownames(prad_mean_fc) <- rownames(mean_fc)
colnames(prad_mean_fc) <- c('CCLE', 'TCGA')
prad_mean_fc <- as.data.frame(prad_mean_fc)
plot <- ggplot(prad_mean_fc, aes(x=CCLE, y=TCGA)) +
  geom_point() +
  geom_text(aes(label=rownames(brca_mean_bl)),hjust=0, vjust=0)
pdf(file="prad_mean_fc.pdf", width=12, height=9)
print(plot)
dev.off()

rownames(brca_mean_fc) <- rownames(mean_fc)
colnames(brca_mean_fc) <- c('CCLE', 'TCGA')
brca_mean_fc <- as.data.frame(brca_mean_fc)
plot <- ggplot(brca_mean_fc, aes(x=CCLE, y=TCGA)) +
  geom_point() +
  geom_text(aes(label=rownames(brca_mean_bl)),hjust=0, vjust=0)
pdf(file="brca_mean_fc.pdf", width=12, height=9)
print(plot)
dev.off()

rownames(coad_mean_fc) <- rownames(mean_fc)
colnames(coad_mean_fc) <- c('CCLE', 'TCGA')
coad_mean_fc <- as.data.frame(coad_mean_fc)
plot <- ggplot(coad_mean_fc, aes(x=CCLE, y=TCGA)) +
  geom_point() +
  geom_text(aes(label=rownames(brca_mean_bl)),hjust=0, vjust=0)
pdf(file="coad_mean_fc.pdf", width=12, height=9)
print(plot)
dev.off()

# Filter out insignificant values
cor.matrix.p[pv.matrix.p>=0.05] <- NA
cor.matrix.s[pv.matrix.s>=0.05] <- NA

# Write again
write.table(cor.matrix.p, file='cor_pearson_filtered.txt', sep='\t')
write.table(cor.matrix.s, file='cor_spearman_filtered.txt', sep='\t')

###############################################################################

print("----- DONE -----")

t_end <- Sys.time()
t_diff <- t_end - t_start

print(t_diff)