
print.mgm <- function(x, ...) {
  
  if('tv.mgm' %in% class(x) | 'tv.mgm' %in% class(x)) {
    
    if('tv.mgm' %in% class(x)) { 
      mc <- 'Time-varying Mixed Graphical Model' 
    } else {
      mc <- 'Time-varying Vector Autoregressive Model' 
    }
    
    cat('mgm fit-object',
        '\n\nModel class: ', mc, 
        '\nNodes: ' , length(x$call$type),  
        '\nTimesteps: ' , x$call$tsteps)
    
  } else {
    
    if('var' %in% class(x)) { 
      mc <- 'Vector Autoregressive Model' 
    } else {
      mc <- 'Mixed Graphical Model' 
    }
    
    cat('mgm fit-object', 
    '\n\nModel class: ', mc, 
        '\nNodes: ' , length(x$call$type))
    
  }
  
}
