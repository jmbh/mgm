

plot.mgm <- function(x, y = NULL, ...) {
  
  # x = bootstrap object
  
  if('rs' %in% class(x)) {
    
    
    # ---------------- For MGM ----------------
    if('mgm' %in% class(x)) {
      
      arguments <- list(...)
      
      # A) For MGM
      
      # Select those edges which are nonzero in at least one of the bootstrap samples
      ind <- unlist(x$edgeNonZero)
      m_values <- do.call(rbind, x$edgeWeights)[ind,]
      v_order <- abs(rowMeans(m_values))
      
      # Get Edge Names
      v_names <- unlist(x$edgeNames)
      
      # Order decreasingly by mean edge weight
      m_values <- m_values[order(v_order, decreasing=FALSE), ]
      v_names <- v_names[order(v_order, decreasing=FALSE)]
      t_m_values <- t(m_values)
      colnames(t_m_values) <- v_names
      
      # check whether any edges are stimated nonzero at all
      if( sum(unlist(x$edgeNonZero)) == 0) {
        
        warning('No Edges were estimated nonzero in any of the resampled datasets.')
        
      } else {
        
        # Plotting
        par(mar=c(3,7,2,1))
        boxplot(t_m_values, horizontal=TRUE, xaxt='n',  yaxt='n')
        axis(2, 1:length(v_names), las=2, cex.axis=.8, labels = v_names)
        axis(1)
        
      }
      
    } 
    
    # ---------------- For VAR ----------------
    
    if('var' %in% class(x)) {
      
    }
    
    
    
  }
  
  
  
} # EoF