# Sakura Maezono # Created on: ##------ Mon Jul 01 00:33:03 2019 ------##

# Purpose: To perform gene enrichment analyses using DAVID
###############################################################################################
# Associated library functions
###############################################################################################
# Windows or Mac
win <- TRUE

if (win) {
  library(doParallel)
}else{
  library(parallel)
}
library(gplots)
library(RDAVIDWebService)

###############################################################################################
# File locations
###############################################################################################

# create folder(s) if it doesn't/ they don't exist

# Set working directory
root <- "G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/"
setwd(root)

# Directory: parameters
dir_p <- 'G:/My Drive/Jiachen Files/data/'

# Directory: expression matrices
dir_e <- 'G:/My Drive/Jiachen Files/GDAC/'

# Directory: for output
dir_o <- 'G:/My Drive/Jiachen Files/Sakura/output/cor_all_vs_pca/'


###############################################################################################
# In-house functions
###############################################################################################
# latest biomaRt annotation  
getBM_hsapiens_ensembl <- function(orig_annotation, listofgenes) {
  #reannotate specified gene signature to the latest biomaRt
  
  require(biomaRt)
  affy_en <- getBM(attributes = c(orig_annotation,"entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                   filters = orig_annotation,
                   values = c(as.character(listofgenes)), 
                   mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))
  
  # order based on origininal annotation
  affy_en <- affy_en[order(affy_en[,1]),]
  
  # "select()" function from biomaRt messes with the 
  # dplyr::select so it is best to detach this package
  detach("package:biomaRt", unload = TRUE)
  
  return(affy_en)
}

reannotate_hngc_to_entrez_fxn <- function(genelist, up_or_down) {
  # reannotate the ensembl ids from WeiChen
  # (for BRCA, PRAD, PAAD, COAD) using biomaRt

  ensembl_entrez <- getBM_hsapiens_ensembl("entrezgene", genelist)
  biomart_SF_eg <- na.omit(ensembl_entrez)

  write.csv(biomart_SF_eg, paste("Annotations",cancer,
                                 "SF_WeiChen.csv", sep = "_"))

  # match ensembl ids
  matching_affyid <- intersect(genelist,
                               unique(biomart_SF_eg$ensembl_gene_id))
  print(paste(length(matching_affyid),
              " of the hugo symbols have other latest biomaRt annotations"))

  # detemine the ensembl ids missing
  missing_affy_id <- genelist[!genelist %in% biomart_SF_eg$ensembl_gene_id]

  if (length(missing_affy_id) > 0) {
    # reannotate to hgnc symbol
    print(paste(missing_affy_id,
                "is/are not present in the latest biomaRt"))
  } else {
    print(paste("All ensembl ids from WeiChen have the latest biomaRt"))
  }

  # only keep the entrez and ensembl id columns
  biomart_SF_eg <- biomart_SF_eg[c(2:4)]

  # only keep the unique rows
  biomart_SF_eg <- unique(biomart_SF_eg)

  # order based on entrez symbol
  biomart_SF_eg <- biomart_SF_eg[order(biomart_SF_eg$entrezgene),]
  rownames(biomart_SF_eg) <- biomart_SF_eg$entrezgene

  # SF entrez gene version
  SF_entrez_mRNA <- as.character(unique(biomart_SF_eg$entrezgene))
  print(paste(length(SF_entrez_mRNA), "entrez ids"))

  # unique ensembl ids version
  SF_en_mRNA <- as.character(unique(biomart_SF_eg$ensembl_gene_id.1))
  print(paste(length(SF_en_mRNA), "ensembl ids"))

  # unique gene symbol version
  SF_symbol_mRNA <- as.character(unique(biomart_SF_eg$hgnc_symbol))
  print(paste(length(SF_symbol_mRNA), "hgnc symbols"))

  return(biomart_SF_eg)
}


###############################################################################
# Import expression matrix from file system
###############################################################################
# APOBEC and hypoxia genes
dir_ext <- "G:/My Drive/Jiachen Files/data/"
apo.genes <- as.character(read.table(file = paste(dir_ext, 
                                                  "tr_a.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)
apo.genes <- apo.genes[order(apo.genes)]
hyp.genes <- as.character(read.table(file = paste(dir_ext,
                                                  "tr_h.txt", sep = ""),
                                     sep = "\t", header = TRUE, 
                                     row.names = 1)$Gene.Symbol)

# "t" for tumour data and "n" for normal data
state <- "t"

# Cancer list
cancers <- as.character(read.table(paste(dir_p, 
                                         'cancer_list.txt',
                                         sep = ''),
                                   header = T, sep = '\t')[,1])

if (state == "n") {
  # "OV" has no normal tissues
  cancers <- cancers[-length(cancers)]
}


### In house functions ####

# remove rows that are only NA
rm_onlyNArow <- function(df) df[rowSums(is.na(df)) != ncol(df),]
rm_onlyNAcol <- function(df) df[colSums(is.na(df)) != nrow(df),]
rm_onlyNANrow <- function(df) df[rowSums(is.nan(df)) != ncol(df),]
rm_row_without_sig <- function(df) df[rowSums(df > 0.05) != ncol(df),]


# Create a DAVIDWebService object connected to David, using your registration email.
# To register, go to: http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.
david <- DAVIDWebService$new(email = 's1535007@u.tsukuba.ac.jp',
                             url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

# Define foreground and background gene lists.
# The foreground list should be contained within the background list.
cancer <- "BRCA"
# gene_sig <- c("vs_hypoxia_sig", "vs_all_genes")
gene_sig <- "vs_all_genes"


myForegroundGenes <- read.table(file = paste(dir_o, cancer, "_", 
                                             "tumour", "_", gene_sig,
                                             "_",  "cor_to_PC1.txt",
                                             sep = ""), sep = '\t')

if (length(myForegroundGenes[,1])>3000) {
  myForegroundGenes <- as.character(myForegroundGenes[c(1:3000),1])
} else {
  myForegroundGenes <- as.character(myForegroundGenes[,1])
}


                            
# egids <- select(org.Hs.eg.db, myForegroundGenes, "ENTREZID", "SYMBOL")[,2]                           
myBackgroundGenes <- myForegroundGenes

FG <- addList(david, myForegroundGenes, idType="REFSEQ_MRNA",
              listName="isClass", listType="Gene")
BG <- addList(david, myBackgroundGenes, idType="REFSEQ_MRNA",
              listName="all", listType="Background")



# Inspect FG and BG to see the proportion of genes recognized by David, and those that are unmapped.
FG
BG

# Inspect "david" object to see the gene lists selected as foreground and background.
david

# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)

# Print functional annotation clustering to file (limited to 3000 genes).
getClusterReportFile(david, "FuncAnnotClust.tsv")