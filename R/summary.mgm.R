

summary.mgm <- function(object, data = NULL, ...) 
  
{
  
  # ---------- Loop over Time Steps ----------  
  
  out_list <- list()
  
  # stationary or time varying?
  if('tv.mgm' %in% class(object) | 'tv.var' %in% class(object)) {
    tsteps <- object$call$tsteps
  } else {
    tsteps <- 1
  }
  
  
  # compute nodewise errors
  if(!is.null(data)) {
  l_errors <- predict.mgm(object, data)  
  l_errors <- lapply(l_errors, function(x) x$error)
  }
  
  for(ts in 1:tsteps) {
    
    if(tsteps>1) { # for time varying
      call <- object$call
      node.models <- object$t.models[[ts]]$node.models
    } else { # for stationary
      call <- object$call
      node.models <- object$node.models
    }
    
    type <- call$type
    nNodes <- length(type)
    
    # ---------- compute measures to report ----------
    
    l_tau <- list()
    for(v in 1:nNodes) l_tau[[v]] <- node.models[[v]]$threshold
    l_lambda <- list()
    for(v in 1:nNodes) l_lambda[[v]] <- node.models[[v]]$lambda
    l_EBIC <- list()
    for(v in 1:nNodes) l_EBIC[[v]] <- node.models[[v]]$EBIC
    
    
    # ---------- Make nice dataframe for save/print ----------
    
    df_out <- data.frame(matrix(NA, nNodes, 1))
    colnames(df_out) <- 'Variable'
    
    # variable lable
    df_out$Variable <- 1:nNodes
    df_out$Type <- type
    
    # degree
    if(tsteps>1) {wadj <- object$wadj[,,ts] } else {wadj <- object$wadj}
    adj <- wadj; adj[adj!=0] <- 1
    
    if('var' %in% class(object)) {
    diag(adj) <- 0
    df_out$degree.in <- colSums(adj)
    df_out$degree.out <- rowSums(adj)
    } else {
    df_out$degree <- colSums(adj)
    }
    
    # fit parameters
    df_out$Lambda <- round(unlist(l_lambda),3)
    df_out$Threshold <- round(unlist(l_tau),3)
    df_out$EBIC <- round(unlist(l_EBIC),3) 
    
    # add errors to data frame
    df_out$Error <- l_errors[[ts]]$Error
    df_out$ErrorType <- l_errors[[ts]]$ErrorType 

    out_list[[ts]] <- df_out
    
  } # end for: timesteps
  
  if(tsteps==1) {
    return(out_list[[1]])
  } else {
    return(out_list)
  }
  
  
} # EoF


