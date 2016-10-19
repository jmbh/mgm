

summary.rs <- function(object, ZeroEdges = FALSE, ...) 
  
{
  
  # ++++++++++ Output for input: Bootstrap Object ++++++++++
  
  if(!('rs' %in% class(object))) stop('Please provide the output of rs.mgm() as an input to this function.')
  
  # --- For MGM ---
  
  if('mgm' %in% class(object)) {
    
    # compute summary measures of bootstrapped edge weights
    l_measures <- lapply(object$edgeWeights, function(ew) {
      round(c(mean(ew), 
              median(ew), 
              quantile(ew, probs = c(.05, .25, .75, .95))),2)
    })
    
    # Put edge IDs and summary measures together
    out <- cbind(do.call(rbind, object$edgeIDs), 
                 do.call(rbind, l_measures))
    colnames(out) <- c('Edge A', 'Edge B', 'Mean', 'Median', '5% Q', '25% Q', '75% Q', '95% Q')
    v_order <- out[,'Mean']
    
    # Order from large to small mean Edge weights
    out <- out[order(v_order, decreasing=TRUE), ]
    
    # check whether any edges are stimated nonzero at all
    if( sum(unlist(object$edgeNonZero)) == 0) {
      warning('No Edges were estimated nonzero in any of the resampled datasets.')
    } else {
      
      # subset nonzero Edges
      if(ZeroEdges==FALSE) {
        ind_nonzero <- unlist(object$edgeNonZero)[order(v_order, decreasing=TRUE)]
        return(out[ind_nonzero>0,])
      } else {
        return(out)
      }
    }
      
      # --- For VAR ---
    } else {
   
      # compute summary measures of bootstrapped edge weights
      l_measures <- lapply(object$edgeWeights, function(ew) {
        round(c(mean(ew), 
                median(ew), 
                quantile(ew, probs = c(.05, .25, .75, .95))),2)
      })
      
      # Put edge IDs and summary measures together
      out <- cbind(do.call(rbind, object$edgeIDs), 
                   do.call(rbind, l_measures))
      colnames(out) <- c('Edge A', 'Edge B', 'Mean', 'Median', '5% Q', '25% Q', '75% Q', '95% Q')
      v_order <- out[,'Mean']
      
      # Order from large to small mean Edge weights
      out <- out[order(v_order, decreasing=TRUE), ]
      
      # check whether any edges are stimated nonzero at all
      if( sum(unlist(object$edgeNonZero)) == 0) {
        warning('No Edges were estimated nonzero in any of the resampled datasets.')
      } else {
        
        # subset nonzero Edges
        if(ZeroEdges==FALSE) {
          ind_nonzero <- unlist(object$edgeNonZero)[order(v_order, decreasing=TRUE)]
          return(out[ind_nonzero>0,])
        } else {
          return(out)
        }
      }
      
    } # end of MGM/VAR if statement

  
} # EoF