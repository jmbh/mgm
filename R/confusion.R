

confusion <- function(tg, # True Graph
                      eg, # Estimated Graph
                      undirected = TRUE, 
                      diagonal = FALSE
                      ) 
  {
  
  
  # ----- Input Checks -----
  
  # Square Matrix
  if(nrow(tg)!=ncol(tg)) stop('True graph has to be represented by a square matrix.')
  if(nrow(eg)!=ncol(eg)) stop('Estimated graph has to be represented by a square matrix.')
  
  # Zero-One
  if(sum(tg %in% c(0,1)) != nrow(tg)^2) stop('The true graph has to be represented by a matrix containing only 0s and 1s.')
  if(sum(eg %in% c(0,1)) != nrow(eg)^2) stop('The estimated graph has to be represented by a matrix containing only 0s and 1s.')
  
  
  if(undirected) {
    
    # Symmetric
    if(!isSymmetric(tg)) stop('True graph has to be symmetric. For directed graphs choose undirected=FALSE.')
    if(!isSymmetric(eg)) stop('Estimated graph has to be symmetric. For directed graphs choose undirected=FALSE.')
    
  }
  
  # ----- Compute Measures -----
  
  if(undirected) {
    
    # Take upper triangular entries of both matrices
    tg_utr <- tg[upper.tri(tg)]
    eg_utr <- eg[upper.tri(eg)]
    
    sen <- mean(eg_utr[which(tg_utr == 1)] == 1) # Sensitivity
    pre <- mean(tg_utr[which(eg_utr == 1)] == 1) # Precision
    acc <- mean(eg_utr == tg_utr) # Accuracy
    spe <- mean(eg_utr[which(tg_utr == 0)] == 0) # Specificity
    
    
  } else {
    
    # Turn matrix into vector
    if(diagonal) {
      diag_tg <- diag(tg)
      diag_eg <- diag(eg)
    } else {
      diag_tg <- NULL
      diag_eg <- NULL
    }
    
    tg_vec <- c(tg[upper.tri(tg)], tg[lower.tri(tg)], diag_tg)
    eg_vec <- c(eg[upper.tri(eg)], eg[lower.tri(eg)], diag_eg)
    
    
    sen <- mean(eg_vec[which(tg_vec == 1)] == 1) # Sensitivity
    pre <- mean(tg_vec[which(eg_vec == 1)] == 1) # Precision
    acc <- mean(eg_vec == tg_vec) # Accuracy
    spe <- mean(eg_vec[which(tg_vec == 0)] == 0) # Specificity
    
    
  }
  
  
  outlist <- list('sensitivity'= sen, 
                  'precision' = pre, 
                  'accuracy' = acc, 
                  'specificity' = spe, 
                  'call'=list('tg' = tg,
                              'eg' = eg,
                              'undirected' = undirected,
                              'diagonal' = diagonal))
  
  return(outlist)
  
  
  
} # eoF

