
print.mgm <- function(x, ...) {
  
  # bootstrap object YES
  if('boot' %in% class(x)) {
    
    if('var' %in% class(x)) { 
      mc <- 'Bootstrapped Mixed Vector Autoregressive Models'
    } else {
      mc <- 'Bootstrapped Mixed Graphical Models'
    }
    
    cat('mgm fit-object', 
        '\n\nModel class: ', mc, 
        '\nNodes: ' , length(x$call$type),
        '\nBootstrap Samples: ' , length(x$B_samples))
    
    # bootstrap object NO
  } else {
    
    if('tv.mgm' %in% class(x) | 'tv.mgm' %in% class(x)) {
      
      if('tv.mgm' %in% class(x)) { 
        mc <- 'Time-varying Mixed Graphical Model' 
      } else {
        mc <- 'Time-varying Mixed Vector Autoregressive Model' 
      }
      
      cat('mgm fit-object',
          '\n\nModel class: ', mc, 
          '\nNodes: ' , length(x$call$type),  
          '\nTimesteps: ' , x$call$tsteps)
      
    } else {
      
      if('var' %in% class(x)) { 
        mc <- 'Mixed Vector Autoregressive Model' 
      } else {
        mc <- 'Mixed Graphical Model' 
      }
      
      cat('mgm fit-object', 
          '\n\nModel class: ', mc, 
          '\nNodes: ' , length(x$call$type))
      
    }
  }
  
} # eoF
