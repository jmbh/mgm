

predictCore_stat <- function(object,
                             data, 
                             consec = NULL)
  
  
{
  
  
  # ----- Compute Aux Variables -----
  
  
  cobj <- class(object)[2]
  nNodes <- ncol(data)
  nCases <- nrow(data)
  call <- object$call
  type <- call$type
  level <- call$level
  k <- call$k
  n_lags <- length(call$lags)
  if(cobj == "mvar") max_lags <- max(call$lags)
  
  
  # Create outlist and storage
  predCoreObj <- list('pred' = vector('list', length = nNodes),
                      'prob' = vector('list', length = nNodes),
                      'true' = NULL,
                      'included' = NULL)
  
  
  # Some Generic Data preparation
  data <- data.frame(data)
  colnames(data)[1:nNodes] <- paste("V", 1:nNodes, '.', sep = "")
  for(sc in which(type=='g')) data[,sc] <- scale(data[,sc]) # scale gaussians
  data_df <- data
  for(sc2 in which(type=='c')) data_df[,sc2] <- as.factor(data_df[,sc2]) # categoricals as factors
  
  nodeModels <- object$nodemodels
  d <- object$call$k - 1
  
  # ----- A.1) mgm -----
  
  if(cobj == 'core') {
    
    for(v in 1:nNodes) {
      
      # -- Define model matrix --
      
      # ----- Construct Design Matrix -----
      
      X_standard <- X <- ModelMatrix_standard(data = data,
                                              type = type,
                                              d = d, 
                                              v = v, 
                                              moderators = object$call$moderators)
      
      if(object$call$overparameterize) {
        
        X_over <- ModelMatrix(data = data, # fix that input, that's stupid
                              type = type,
                              level = level,
                              labels = colnames(data),
                              d = d, 
                              moderators = object$call$moderators,
                              v = v)
        
        X <- X_over
        
      } # end if: overparameterize?
      
      data_df <- as.data.frame(data_df)
      y <- as.numeric(data[, v])
      
      
      if(type[v]=='c') {
        
        ## Prediction Categorical
        coefs <- nodeModels[[v]]$model
        n_cat <- length(coefs)
        Potentials <- matrix(NA, nCases, n_cat)
        
        # loop over categories & compute potentials
        for(cat in 1:n_cat) Potentials[,cat] <- exp(coefs[[cat]][1] + X %*% matrix(coefs[[cat]][-1], nrow=length(coefs[[cat]][-1])))
        
        # compute category-probabilities
        Probabilities <- Potentials[,1:n_cat] / rowSums(Potentials[,1:n_cat])
        pred_class_id <-  apply(Probabilities, 1, which.max) # classify
        predCoreObj$pred[[v]] <- sort(unique(data[,v]))[pred_class_id] # matches the ordering of predicted categories in glmnet(), which also orders with increasing integers
        predCoreObj$prob[[v]] <- Probabilities
        
        
      } else {
        
        ## Prediction Continuous (same for Gauss and Pois, because in both cases the natural parameter is equal to the expectation)
        # predicitions
        coefs <- as.numeric(nodeModels[[v]]$model) # get coefficients
        predCoreObj$pred[[v]] <- coefs[1] + X %*% matrix(coefs[-1], nrow=length(coefs[-1])) # predict
        
      }
      
    } # end of variable loop 'mgm'
    
    # Save ground truth
    predCoreObj$true <- data
    
  }
  
  
  
  # ----- A.2) mvar -----
  
  if(cobj == 'mvar') {
    
    # Prepare Data (already cuts max(lags) first observations to compute design matrix)
    data_lagged <- lagData(data = data, 
                           lags = object$call$lags, 
                           consec = consec)
    
    predCoreObj$included <- data_lagged$included # this specifies additionally, whether measurements are successive after the max(lags) measurement
    
    data_response <- data_lagged$data_response
    l_data_lags <- data_lagged$l_data_lags
    data_response <- apply(data_response, 2, as.numeric) # to avoid confusion with labels of categories if there are factors
    
    data_response <- data_response[predCoreObj$included, ]
    l_data_lags <- lapply(l_data_lags, function(z) z <- z[predCoreObj$included, ])
    
    nCases <- nrow(data_response)
    
    for(v in 1:nNodes) {
      
      
      # ----- Create VAR Design Matrix -----
      
      # append response with predictors
      y <- data_response[, v] # response variable v
      data_input_MM <- do.call(cbind, l_data_lags)
      data_input_MM <- as.data.frame(data_input_MM)
      
      # turn categoricals into factors for model.matrix()
      for(i in which(type=='c')) data_input_MM[, i] <- as.factor(data_input_MM[, i])  
      
      # need to have response and predictors in one dataframe for model.matrix()
      data_v <- cbind(y, data_input_MM) 
      
      # Dummy coding
      form <- as.formula('y ~ (.)')
      
      # Construct standard design matrix (to get number of parameters for tau threshold)
      X_standard <- model.matrix(form, data = data_v)[, -1] # delete intercept (added by glmnet later)
      
      if(object$call$overparameterize) {
        
        # Compute augmented type and level vectors
        type_aug <- rep(type, n_lags)
        level_aug <- rep(level, n_lags)
        
        # Construct over-parameterized design matrix
        X_over <- ModelMatrix(data = data_input_MM,
                              type = type_aug,
                              level = level_aug,
                              labels = colnames(data_input_MM),
                              d = 1, 
                              v = NULL)
        X <- X_over
        
      } else {
        
        X <- X_standard
        
      }
      
      # ---- Prediction -----
      
      if(type[v]=='c') {
        
        ## Prediction Categorical
        coefs <- nodeModels[[v]]$model
        n_cat <- length(coefs)
        Potentials <- matrix(NA, nCases, n_cat)
        
        # loop over categories & compute potentials
        for(cat in 1:n_cat) Potentials[,cat] <- exp(coefs[[cat]][1] + X %*% matrix(coefs[[cat]][-1], nrow=length(coefs[[cat]][-1])))
        
        # compute category-probabilities
        Probabilities <- Potentials[, 1:n_cat] / rowSums(Potentials[, 1:n_cat])
        pred_class_id <-  apply(Probabilities, 1, which.max) # classify
        predCoreObj$pred[[v]] <- sort(unique(data[,v]))[pred_class_id] # matches the ordering of predicted categories in glmnet(), which also orders with increasing integers
        predCoreObj$prob[[v]] <- Probabilities
        
      } else {
        
        ## Prediction Continuous (same for Gauss and Pois, because in both cases the natural parameter is equal to the expectation)
        # predicitions
        coefs <- as.numeric(nodeModels[[v]]$model) # get coefficients
        predCoreObj$pred[[v]] <- coefs[1] + X %*% matrix(coefs[-1], nrow=length(coefs[-1])) # predict
        
      }
      
    } # end loop: v
    
    
    # Save true (reduced) data
    predCoreObj$true <- data_response
    
  } # end if: mvar?
  
  
  # ----- Return -----
  
  return(predCoreObj)
  
  
} # eoF





