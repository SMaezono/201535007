# download.R

# Jiachen Liang

# Downloads data sets from GDAC

# WARNING: Unzip function currently does not work on Windows
# Use 7zip to extract txt files, and put them in the root directory before running gdac_preprocess.R

# Using Windows?
win <- TRUE

cancers <- c("BLCA", "BRCA", "CESC", "COAD", "HNSC", "KIRC", "KIRP", "LIHC",
             "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC", "OV")

# Root directory
#root <- "/Users/jiachen/Downloads/GDAC/"
root <- 'G:/My Drive/Jiachen Files/data/'
setwd(root)

for(c in cancers){
  path <- paste(
    "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/",
    c,
    "/20160128/gdac.broadinstitute.org_",
    c,
    ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz",
    sep = ""
    )
  gz_dir <- paste(root, c, ".tar.gz", sep="")
  tar_dir <- paste(root, c, ".tar", sep="")
  if(win){
    download.file(url = path, destfile = gz_dir)
  }else{
    download.file(url = path, destfile = tar_dir)
  }
  if(!win){
    file <- paste(
      c,
      ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
      sep="")
    untar(tar_dir, files=file, exdir=root)
    file.rename(
      file,
      paste('TCGA/', c, ".txt", sep="")
    )
  }
}