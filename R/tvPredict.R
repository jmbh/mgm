

tvPredict <- function(data, # time series data 
                      type, # necessary for prediction: continuous or categorical
                      test, # indicator vector: training time points
                      VAR, 
                      model # time-varying model
) 
  
  
{
  
  # ----- Input Checks -----
  
  
  
  
  # ----- Aux Variables -----
  
  nTest <- length(test)
  p <- ncol(data)
  n <- nrow(data)
  
  # For now Fix, later one can choose.
  error.categorical = 'CorrectClass'
  error.continuous = 'VarExpl'
  
  # model <- l_bwModels[[bw]]
  
  
  
  
  # ----- Compute Predictions -----
  
  data_pred <- matrix(NA, nTest, p) # Storage for predicted values
  
  # -- For Contemporaneous Models -- (this code is very similar to predict.mgm() )
  if(!VAR) {
    
    l_predm <- list() # storage for all predictions at each tv model location
    
    for(i in 1:nTest) {
      
      l_preds <- list()
      
      # Select model for time point i
      node.models <- model$t.models[[i]]$node.models
      
      for(sc in which(type=='g')) data[,sc] <- scale(data[,sc]) # scale gaussians
      
      for(v in 1:p) {
        
        # Make model matrix (compute dummies for categoricals)
        data_m <- data.frame(data)
        for(cat in which(type=='c')) data_m[,cat] <- as.factor(data_m[,cat])
        form <- as.formula(paste(colnames(data_m)[v],"~ (.)"))
        X <- model.matrix(form, data= data_m)[,-1]   
        
        if(type[v]=='c') {
          
          ## Prediction Categorical
          coefs <- node.models[[v]]$coefs
          K <- length(coefs)
          Potentials <- matrix(NA, n, K)
          
          # loop over categories & compute potentials
          
          for(k in 1:K) {
            # Compute exp(potentials)
            Potentials[,k] <- exp(coefs[[k]][1] + X %*% matrix(coefs[[k]][-1][1:ncol(X)], nrow=length(coefs[[k]][-1][1:ncol(X)])))
          }
          
          # compute category-probabilities
          Probabilities <- Potentials[,1:K] / rowSums(Potentials[,1:K])
          pred_class_id <- apply(Probabilities, 1, which.max) # Classify
          l_preds[[v]] <- sort(unique(data[,v]))[pred_class_id]
          
        } else {
          ## Prediction Continuous
          coefs <- as.numeric(node.models[[v]]$coefs) # get coefficients
          l_preds[[v]] <-  coefs[1] + X %*% matrix(coefs[-1][1:ncol(X)], nrow=length(coefs[-1][1:ncol(X)])) # predict
        }
        
      } # end of variable loop 'mgm'
      
      # Collect everything
      l_predm[[i]] <- do.call(cbind, l_preds) # All predictions
      data_pred[i,] <- l_predm[[i]][test[i], ] # Only predictions on test set
      
    }

    
    # -- For VAR Models --
  } else {
    
    
  }
  
  
  # ----- Compute Errors -----
  
  PF <- matrix(NA, nTest, p)
  
  for(v in 1:p) {
    
    if(type[v]=='c') {
      
      error <- (data_pred[,v] == data[test, v]) * 1 # 0/1 loss
      if(error.categorical == 'CorrectClassNorm') {
        tb <- table(data[,v])
        norm_constant <- max(tb/sum(tb)) # normalize by maximum additional accuracy that can be predicted beyond the larger marginal frequency
        error <- (error-norm_constant) / (1-norm_constant)
      } 
      PF[, v] <- error

    } else {
      
      preds <- data_pred[, v]
      if(error.continuous=='RMSE') {
        error <- sqrt(sum((preds - as.numeric(data[test,v]))^2)/ nTest) # compute RMSE
      } else {
        error <- 1 - var(preds-data[test,v]) / var(data[test,v])
      }
      PF[, v] <- error
      
    }
    
    
  }
  
  
  # ----- Output -----
  
  outlist <- list('per' = PF,
                  'perVar' = colMeans(PF),
                  'perMean' = mean(PF))
  
  return(outlist)
  
}



