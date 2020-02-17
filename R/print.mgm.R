

print.mgm <- function(x,
                      ...)
  
  
  
{
  
  model_classes <- c('Mixed Graphical Model (MGM)',
                     'mixed Vector Autoregressive (mVAR) model',
                     'Time-varying Mixed Graphical Model (tv-MGM)',
                     'Time-varying mixed Vector Autoregressive (tv-mVAR) model')
  
  # ---------- print for fit objects ----------
  
  if(!('predicted' %in% class(x)) & !('bwSelect' %in% class(x))) {
    
    if('core' %in% class(x)) {
      
      msg_basic <- paste0('mgm fit-object',
                          '\n\nModel class: ', model_classes[1],
                          '\nOrder: ' , x$call$k,
                          '\nNodes: ' , length(x$call$type))
      
      
      if(!is.null(x$call$moderators))  {
        
        x$call$k <- 3
        
        msg_basic <- paste0('mgm fit-object',
                            '\n\nModel class: ', model_classes[1],
                            '\nOrder: ' , x$call$k,
                            '\nNodes: ' , length(x$call$type)) # call again with updated x$call$k
        
        
        if(is.matrix(x$call$moderators)) mod_text <- "Custom specification" else mod_text <- paste(x$call$moderators, collapse = ", ")

        msg_basic <- paste0(msg_basic,
                            paste0('\nModerators: ' , mod_text)) 
        
    
      }
        
      
      if(!is.null(x$call$condition)) {
        
        nCond <- length(x$call$condition)
        names <- names(x$call$condition)
        msgCond <- paste0(names, "=", unlist(x$call$condition))
        msg_basic <- paste0(msg_basic, "\nFixed: ", paste(msgCond, collapse = ", "))
        
      }
      
      cat(msg_basic)

      
    } # end if: basic mgm object?
    
    
    
    
    
    
    if('mvar' %in% class(x)) {
      
      n_incl <- sum(x$call$data_lagged$included == TRUE)
      n_exp <- sum(x$call$data_lagged$included == FALSE)
      n <- n_incl + n_exp
      perc <- n_incl / n
      perc <- round(perc, 4) * 100
      
      cat('mgm fit-object',
          '\n\nModel class: ', model_classes[2],
          '\nLags: ' , x$call$lags,
          '\nRows included in VAR design matrix: ' , n_incl ,'/', n, '(', perc, '%)',
          '\nNodes: ' , length(x$call$type))
    }
    
    
    
    if('tvmgm' %in% class(x)) {
      
      if(!is.null(x$call$moderators)) {
        
        x$call$k <- 3
        
        cat('mgm fit-object',
            '\n\nModel class: ', model_classes[3],
            '\nOrder: ' , x$call$k,
            '\nNodes: ' , length(x$call$type),
            '\nModerators: Variable ' , x$call$moderators,
            '\nEstimation points: ' , length(x$call$estpoints))
        
      } else {
        
        cat('mgm fit-object',
            '\n\nModel class: ', model_classes[3],
            '\nOrder: ' , x$call$k,
            '\nNodes: ' , length(x$call$type),
            '\nEstimation points: ' , length(x$call$estpoints))
        
      }
      
      
      
    }
    
    
    if('tvmvar' %in% class(x)) {
      
      
      n_incl <- sum(x$call$data_lagged$included == TRUE)
      n_exp <- sum(x$call$data_lagged$included == FALSE)
      n <- n_incl + n_exp
      perc <- n_incl / n
      perc <- round(perc, 4) * 100
      
      
      cat('mgm fit-object',
          '\n\nModel class: ', model_classes[4],
          '\nLags: ' , x$call$lags,
          '\nRows included in VAR design matrix: ' , n_incl ,'/', n, '(', perc, '%)',
          '\nNodes: ' , length(x$call$type),
          '\nEstimation points: ' , length(x$call$estpoints))
    }
    
  }
  
  
  # ---------- print for prediction object ----------
  
  if('predicted' %in% class(x)) {
    
    if('mgm' %in% class(x)) mc <- model_classes[1]
    if('mvar' %in% class(x)) mc <- model_classes[2]
    if('tvmgm' %in% class(x)) mc <- model_classes[3]
    if('tvmvar' %in% class(x)) mc <- model_classes[4]
    
    
    cat('mgm prediction-object',
        '\n\nModel class: ', mc,
        '\nError Types:', paste(names(x$call$errorCon), names(x$call$errorCat)))
    
  }
  
  
  # ---------- print for bwSelect object ----------
  
  
  if('bwSelect' %in% class(x)) {
    
    if('mgm' %in% class(x)) mc <- model_classes[3]
    if('mvar' %in% class(x)) mc <- model_classes[4]
    
    cat('mgm bandwidth selection-object',
        '\n\nModel class: ', mc,
        '\nBandwith path: ', paste(x$call$bwSeq),
        '\nNumber of Folds: ', paste(x$call$bwFolds),
        '\nFoldsize: ', paste(x$call$bwFoldsize),
        '\nOptimal Bandwidth: ', x$call$bwSeq[which.min(x$meanError)])
    
  }
  
  
  
  
} # eoF
