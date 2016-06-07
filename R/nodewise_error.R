

nodewise_error <- function(fitobj, data) {
  
  # ---------- input checks ----------
  
  if(fitobj$call$method!='glm') stop("Nodewise error can only be computed for method='glm'.")
  
  # ---------- Loop over Time Points ----------
  
  out_list <- list()
  
  if('tv.mgm' %in% class(fitobj) | 'tv.var' %in% class(fitobj)) {
    tsteps <- fitobj$call$tsteps
  } else {
    tsteps <- 1
  }
  
  for(ts in 1:tsteps) {
    
  if(tsteps>1) { # for time varying
    call <- fitobj$call
    node.models <- fitobj$t.models[[ts]]$node.models
    weights_raw <- fitobj$t.models[[ts]]$call$weights
    weights <- weights_raw / sum(weights_raw)
  } else { # for stationary
    call <- fitobj$call
    node.models <- fitobj$node.models
    weights <- rep(1/nrow(data), nrow(data))
    if('var' %in% class(fitobj)) weights <- weights[-1]
  }
  
  # ---------- Loop over Variables & collect Errors ----------
  
  error_list <- list()
  
  if(('var' %in% class(fitobj))==FALSE) {
    
    # ----- for mgm objects ----
    type <- call$type
    nNodes <- length(type)
    nCases <- nrow(data)
    for(sc in which(type=='g')) data[,sc] <- scale(data[,sc]) # scale gaussians
    
    for(v in 1:nNodes) {
      
      # Make model matrix (compute dummies for categoricals)
      data_m <- data.frame(data)
      for(cat in which(type=='c')) data_m[,cat] <- as.factor(data_m[,cat])
      form <- as.formula(paste(colnames(data_m)[v],"~ (.)"))
      X <- model.matrix(form, data= data_m)[,-1]   
      
      if(type[v]=='c') {
        
        ## Prediction Categorical
        coefs <- node.models[[v]]$coefs
        K <- length(coefs)
        Potentials <- matrix(NA, nCases, K+1)
        
        # loop over categories & compute potentials
        
        for(k in 1:K) {
          # Compute exp(potentials)
          Potentials[,k] <- exp(coefs[[k]][1] + X %*% matrix(coefs[[k]][-1][1:ncol(X)], nrow=length(coefs[[k]][-1][1:ncol(X)])))
        }
        
        # compute category-probabilities
        Potentials[,K+1] <- sum(Potentials[,1:K])
        Probabilities <- Potentials[,1:K] / rowSums(Potentials[,1:K])
        PredictedCat <- apply(Probabilities, 1, which.max) # classify
        error_list[[v]] <- sum(weights*(PredictedCat==data[,v])) # proportion correctly classified
        
      } else {
        ## Prediction Continuous
        # predicitions
        coefs <- as.numeric(node.models[[v]]$coefs) # get coefficients
        preds <-  coefs[1] + X %*% matrix(coefs[-1][1:ncol(X)], nrow=length(coefs[-1][1:ncol(X)])) # predict
        error_list[[v]] <- sqrt(sum(weights*(preds-as.numeric(data[,1]))^2) ) # compute RMSE
        
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
      data_lag1 <- data.frame(f_VARreshape(as.matrix(data)))
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
        PredictedCat <- apply(Probabilities, 1, which.max) # classify
        error_list[[v]] <- sum(weights*(PredictedCat==data[-1,v])) # proportion correctly classified
        
      } else {
        ## Prediction Continuous
        # predicitions
        coefs <- as.numeric(node.models[[v]]$coefs) # get coefficients
        preds <-  coefs[1] + X %*% matrix(coefs[-1][1:ncol(X)], nrow=length(coefs[-1][1:ncol(X)])) # predict
        error_list[[v]] <- sqrt(sum(weights*(preds-as.numeric(data[-1,v]))^2)) # compute RMSE
        
      }
      
    } # end of variable loop 'var'
    
  } # end if: mgm vs var

  # ---------- output ----------
  
  # Define variable names
  if(is.null(colnames(data))) colnames(data)  <- 1:nNodes
  
  # Define error type
  eType <- rep(NA, nNodes)
  eType[type[1:nNodes]=='c'] <- '%Correct'
  eType[type[1:nNodes]!='c'] <- 'RMSE'
  
  error_out <- data.frame(matrix(NA, nNodes, 3))
  colnames(error_out) <- c("Variable", 'Error', "ErrorType")
  error_out$Variable <- colnames(data)
  error_out$Error <- round(unlist(error_list),3)
  error_out$ErrorType <- eType
  
  out_list[[ts]] <- error_out
  
  } # end for: timesteps
  
  if(tsteps==1) {
    return(out_list[[1]])
  } else {
    return(out_list)
  }

} # EoF

