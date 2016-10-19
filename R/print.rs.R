

print.rs <- function(x, ...) {
  
  # bootstrap object YES
  if('rs' %in% class(x)) {
    
    if('var' %in% class(x)) { 
      mc <- 'Resampled Mixed Vector Autoregressive Models'
    } else {
      mc <- 'Resampled Mixed Graphical Models'
    }
    
    cat('mgm fit-object', 
        '\n\nModel class: ', mc, 
        '\nNodes: ' , length(x$call$type),
        '\nIterations: ' , x$call_rs$B,
        '\nReplace: ' , x$call_rs$replace,
        '\nSample Size: ' , x$call_rs$N)
    
    # bootstrap object NO
  } 
  
}