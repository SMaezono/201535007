whichSignatures <- function(tumor.ref = NA, 
                           sample.id, 
                           signatures.ref = signatures.nature2013, 
                           associated = c(), 
                           signatures.limit = NA, 
                           # do sensitivity test (0.05,0.07, etc); original 0.06
                           signature.cutoff = 0.07,
                           contexts.needed = FALSE, 
                           tri.counts.method = "default") {
  
  if(class(tumor.ref) == 'matrix') {
    stop(paste('Input tumor.ref needs to be a data frame or location of input text file', sep = ''))
  }
  
  if(exists("tumor.ref", mode = "list")){
    tumor     <- tumor.ref
    if(contexts.needed == TRUE){
      tumor   <- getTriContextFraction(mut.counts.ref = tumor, trimer.counts.method = tri.counts.method) 
    }
  } else {
    if(file.exists(tumor.ref)){
      tumor   <- utils::read.table(tumor.ref, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
      if(contexts.needed == TRUE){
        tumor <- getTriContextFraction(tumor, trimer.counts.method = tri.counts.method) 
      }
    } else {
      print("tumor.ref is neither a file nor a loaded data frame")
    }
  }  
  
  if (missing(sample.id) && nrow(tumor) == 1) {
    sample.id = rownames(tumor)[1]
  }
  # Take patient id given
  tumor <- as.matrix(tumor)
  if(!sample.id %in% rownames(tumor)){
    stop(paste(sample.id, " not found in rownames of tumor.ref", sep = ''))
  }
  tumor <- subset(tumor, rownames(tumor) == sample.id)
  if(round(rowSums(tumor), digits = 1) != 1){
    stop(paste('Sample: ', sample.id, ' is not normalized\n', 'Consider using "contexts.needed = TRUE"', sep = ' '))
  }
  
  #Read in Stratton signatures file  
  if(exists("signatures.ref", mode = "list")){
    signatures   <- signatures.ref
  } else {
    if(file.exists(signatures.ref)){
      signatures <- read.table(signatures.ref, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE, row.names=1)
    } else {
      print("signatures.ref is neither a file nor a loaded data frame") 
    }  
  }
  
  signatures    <- as.matrix(signatures)
  original.sigs <- signatures
  
  # Check column names are formatted correctly
  if(length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))){
    colnames(tumor) <- changeColumnNames(colnames(tumor))
    if(length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))){
      stop("Check column names on input file")
    }
  }
  
  # Ensure that columns in tumor match the order of those in signatures
  tumor <- tumor[,colnames(signatures), drop = FALSE]
  
  #Take a subset of the signatures
  if(!is.null(associated)){
    signatures <- signatures[rownames(signatures) %in% associated, ]
  }
  
  if(is.na(signatures.limit)){
    signatures.limit <- nrow(signatures)
  }
  
  # Remove signatures from possibilities if they have a "strong" peak not seen in the tumor sample
  zero.contexts   <- colnames(tumor)[tumor < 0.01]
  corr.sigs       <- which(signatures[,zero.contexts] >= 0.2, arr.ind = T)
  signatures      <- signatures[which(!rownames(signatures) %in% rownames(corr.sigs)),,drop = FALSE]
  #print(paste(rownames(corr.sigs), " not considered in the analysis.", sep = ""))
  
  #Set the weights matrix to 0
  weights         <- matrix(0, nrow = nrow(tumor), ncol = nrow(signatures), dimnames = list(rownames(tumor), rownames(signatures)))
  
  seed            <- findSeed(tumor, signatures)
  weights[seed]   <- 1
  w               <- weights*10
  
  error_diff      <- Inf
  error_threshold <- 1e-3
  
  num <- 0
  while(error_diff > error_threshold){
    num        <- num + 1
    #print(num)
    error_pre  <- getError(tumor, signatures, w)
    #print(w)
    w          <- updateW_GR(tumor, signatures, w, signatures.limit = signatures.limit)
    error_post <- getError(tumor, signatures, w)
    #print(paste("old_error = ", error_pre, sep = ""))
    #print(paste("new_error = ", error_post, sep = ""))
    #print(w)
    error_diff <- (error_pre-error_post)/error_pre
  }
  
  weights <- w/sum(w)
  unknown <- 0
  
  ## filtering on a given threshold value (0.06 default)
  weights[weights < signature.cutoff ] <- 0
  unknown <- 1 - sum(weights)
  
  product <- weights %*% signatures
  diff    <- tumor - product
  
  x       <- matrix(data = 0, nrow = 1, ncol = nrow(original.sigs), dimnames = list(rownames(weights), rownames(original.sigs)))
  x       <- data.frame(x)
  x[colnames(weights)] <- weights
  weights <- x
  
  out        <- list(weights, tumor, product, diff, unknown)
  names(out) <- c("weights", "tumor", "product", "diff", "unknown")
  return(out)
  
}
