

# Gives identical interactions (e.g. 1-2-3 and 3-1-2) the same flag

# Input: Matrix with row = number of interactions, col = involved variables
# Output: numeric indicator vector



FlagSymmetric <- function(x) {
  
  vec_sim <- rep(NA, nrow(x))
  ind_ord <- ncol(x)
  
  counter <- 1
  
  for(i in 1:nrow(x)) {
    
    if(is.na(vec_sim[i])) {
      
      vec_sim[i] <- counter
      
      for(j in (i+1):nrow(x)) {
        
        if( (i+1) > nrow(x) ) next # in case of very few interactions
        
        ind <- x[j, ] %in% x[i, ]
        if(sum(ind)==ind_ord) vec_sim[j] <- counter
        
      }
      
      counter <- counter + 1
      
    }
  }
  return(vec_sim)
}



FlagSymmetricFast <- function(x) {
  
  x <- data.frame(x)
  flag <- as.numeric(factor(apply(x, 1, function(x) paste0(sort(x), collapse = "-"))))
  flag_consec <- as.numeric(factor(flag, levels = unique(flag))) # I think unnecessary
  
  return(flag_consec)
}



