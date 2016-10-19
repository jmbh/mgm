


predict.mgm <- function(object, 
                        data, 
                        variables='all', 
                        error.continuous = 'RMSE', # or 'VarExpl'
                        error.categorical = 'CorrectClass', # or 'CorrectClassNorm'
                        ...) {
  
  # ---------- input checks ----------
  
  if(object$call$method!='glm') stop("Predictions can only be computed for estimation method='glm'.")
  if(ncol(data)!=length(object$call$type)) stop('Please provide a dataset matching the one that was used for estimation.')
  
  # define relevant variables
  if(variables[1]=='all') {
    Nseq <- 1:ncol(data)
  } else if(class(variables)=='character') {
    Nseq <- which(colnames(data) %in% variables);
  } else if (sum(!(round(variables)==variables))==0) # check for integers
  { 
    Nseq <- variables
  } else { 
    stop("Please provide the a vector of column indices or column names as input for the argument 'variables' " )
  }
  
  # Assign column names
  if(is.null(colnames(data))) colnames(data)  <- 1:ncol(data)
  
  # ---------- Loop over Time Points ----------
  
  out_list <- list()
  
  if('tv.mgm' %in% class(object) | 'tv.var' %in% class(object)) {
    tsteps <- object$call$tsteps
  } else {
    tsteps <- 1
  }
  
  for(ts in 1:tsteps) {
    
    if(tsteps>1) { # for time varying
      call <- object$call
      node.models <- object$t.models[[ts]]$node.models
      weights_raw <- object$t.models[[ts]]$call$weights
      weights <- weights_raw / sum(weights_raw)
    } else { # for stationary
      call <- object$call
      node.models <- object$node.models
      weights <- rep(1/nrow(data), nrow(data))
      if('var' %in% class(object)) weights <- weights[-1]
    }
    
    # ---------- Loop over Variables & collect Errors ----------
    
    error_list <- list()
    pred_list <- list()
    pred_list_prob <- list()
    
    if(('var' %in% class(object))==FALSE) {
      
      # ----- for mgm objects ----
      type <- call$type
      nNodes <- length(type)
      nCases <- nrow(data)
      for(sc in which(type=='g')) data[,sc] <- scale(data[,sc]) # scale gaussians
      
      for(v in Nseq) {
        
        # Make model matrix (compute dummies for categoricals)
        data_m <- data.frame(data)
        for(cat in which(type=='c')) data_m[,cat] <- as.factor(data_m[,cat])
        form <- as.formula(paste(colnames(data_m)[v],"~ (.)"))
        X <- model.matrix(form, data= data_m)[,-1]   
        
        if(type[v]=='c') {
          
          ## Prediction Categorical
          coefs <- node.models[[v]]$coefs
          K <- length(coefs)
          Potentials <- matrix(NA, nCases, K)
          
          # loop over categories & compute potentials
          
          for(k in 1:K) {
            # Compute exp(potentials)
            Potentials[,k] <- exp(coefs[[k]][1] + X %*% matrix(coefs[[k]][-1][1:ncol(X)], nrow=length(coefs[[k]][-1][1:ncol(X)])))
          }
          
          # compute category-probabilities
          Probabilities <- Potentials[,1:K] / rowSums(Potentials[,1:K])
          pred_class_id <-  apply(Probabilities, 1, which.max) # classify
          pred_list[[v]] <- PredictedCat <- sort(unique(data[,v]))[pred_class_id]
          pred_list_prob[[v]] <- Probabilities
          error_list[[v]] <- sum(weights*(PredictedCat==data[,v])) # proportion correctly classified
          if(error.categorical == 'CorrectClassNorm') {
            tb <- table(data[,v])
            norm_constant <- max(tb/sum(tb)) # normalize by maximum additional accuracy that can be predicted beyond the larger marginal frequency
            error_list[[v]] <- (error_list[[v]]-norm_constant) / (1-norm_constant)
          }
          
        } else {
          ## Prediction Continuous
          # predicitions
          coefs <- as.numeric(node.models[[v]]$coefs) # get coefficients
          pred_list[[v]] <- preds <-  coefs[1] + X %*% matrix(coefs[-1][1:ncol(X)], nrow=length(coefs[-1][1:ncol(X)])) # predict
          pred_list_prob[[v]] <- preds
          if(error.continuous=='RMSE') {
            error_list[[v]] <- sqrt(sum(weights*(preds-as.numeric(data[,v]))^2) ) # compute RMSE
          } else {
            error_list[[v]] <- 1 - var(preds-data[,v]) / var(data[,v])
          }
        }
        
      } # end of variable loop 'mgm'
      
    } else {
      
      # ----- for VAR objects ----
      type <- call$type
      type2 <- c(type, type)
      nNodes <- length(type)
      nCases <- nrow(data)
      for(sc in which(type[1:nNodes]=='g')) data[,sc] <- scale(data[,sc]) # scale gaussians
      
      for(v in 1:nNodes) {
        
        # Make model matrix (compute dummies for categoricals)
        data_lag1 <- data.frame(VARreshape(as.matrix(data)))
        for(cat in which(type2=='c')) data_lag1[,cat] <- as.factor(data_lag1[,cat])
        data_lag1_cut <-   data_lag1[,-c((1:nNodes)[-v])]
        form <- as.formula(paste(colnames(data_lag1_cut)[1],"~ (.)"))
        X <- model.matrix(form, data= data_lag1_cut)[,-1]   
        
        if(type[v]=='c') {
          
          ## Prediction Categorical
          coefs <- node.models[[v]]$coefs
          K <- length(coefs)
          Potentials <- matrix(NA, nCases-1, K+1)
          
          # loop over categories & compute potentials
          for(k in 1:K) {
            # Compute exp(potentials)
            Potentials[,k] <- exp(coefs[[k]][1] + X %*% matrix(coefs[[k]][-1][1:ncol(X)], nrow=length(coefs[[k]][-1][1:ncol(X)])))
          }
          
          # compute category-probabilities
          Potentials[,K+1] <- sum(Potentials[,1:K])
          Probabilities <- Potentials[,1:K] / rowSums(Potentials[,1:K])
          pred_class_id <-  apply(Probabilities, 1, which.max) # classify
          pred_list[[v]] <- PredictedCat <- sort(unique(data[-1,v]))[pred_class_id]
          pred_list_prob[[v]] <- Probabilities
          error_list[[v]] <- sum(weights*(PredictedCat==data[-1,v])) # proportion correctly classified
          if(error.categorical == 'CorrectClassNorm') {
            tb <- table(data[-1,v])
            norm_constant <- max(tb/sum(tb)) # normalize by maximum additional accuracy that can be predicted beyond the larger marginal frequency
            error_list[[v]] <- (error_list[[v]]-norm_constant) / (1-norm_constant)
          }
        } else {
          ## Prediction Continuous
          # predicitions
          coefs <- as.numeric(node.models[[v]]$coefs) # get coefficients
          pred_list[[v]] <- preds <-  coefs[1] + X %*% matrix(coefs[-1][1:ncol(X)], nrow=length(coefs[-1][1:ncol(X)])) # predict
          pred_list_prob[[v]] <- preds
          if(error.continuous=='RMSE') {
            error_list[[v]] <- sqrt(sum(weights*(preds-as.numeric(data[-1,v]))^2) ) # compute RMSE
          } else {
            error_list[[v]] <- 1 - var(preds-data[-1,v]) / var(data[-1,v])
          }
        }
        
      } # end of variable loop 'var'
      
    } # end if: mgm vs var
    
    
    # ---------- output ----------
    
    # Define error type
    eType <- rep(NA, nNodes)
    
    eType[type[1:nNodes]=='c'] <- error.categorical
    eType[type[1:nNodes]!='c'] <- error.continuous
    
    error_out <- data.frame(matrix(NA, length(Nseq), 3))
    colnames(error_out) <- c("Variable", 'Error', "ErrorType")
    error_out$Variable <- colnames(data)[Nseq]
    error_out$Error <- round(unlist(error_list),3)
    error_out$ErrorType <- eType[Nseq]
    
    # collapse predictions into matrix
    predmat <- do.call(cbind, pred_list)
    
    out_list[[ts]] <- list('error'=error_out, 'pred'=predmat, 'pred_prob'=pred_list_prob)
    
  } # end for: timesteps
  
  if(tsteps==1) {
    return(out_list[[1]])
  } else {
    return(out_list)
  }
  
  
} # EoF
