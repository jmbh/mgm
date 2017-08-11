

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
  
  
  # ----- A.1) mgm -----
  
  if(cobj == 'core') {
    
    for(v in 1:nNodes) {
      
      # -- Define model matrix --
      
      if(call$overparameterize) {
        
        # Construct over-parameterized design matrix
        X <- ModelMatrix(data = data[, -v],
                         type = type[-v],
                         level = level[-v],
                         labels = colnames(data)[-v],
                         d = k - 1)
        
        y <- as.numeric(data[, v])
        
        
      } else {
        
        if (call$k == 2){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
        } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^", call$k - 1)) }
        
        X <- model.matrix(form, data=data_df)[, -1] # delete intercept (added by glmnet later)
        y <- as.numeric(data[, v])
        
      }
      
      
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
    
    # if consec is used, then reduce dataset to usable observations
    ind_included_wo_begin <- data_lagged$included #[-c(1:n_lags)] # the weights-vector has alread length nrow-max(lags)
    
    data_response <- data_response[ind_included_wo_begin, ]
    l_data_lags <- lapply(l_data_lags, function(z) z <- z[ind_included_wo_begin, ])
    
    # browser()
    
    nCases <- nrow(data_response)
    
    for(v in 1:nNodes) {
      
      # ---- Create design matrix -----
      
      # append response with predictors
      y <- data_response[,v] # response variable v
      data_v <- cbind(y, do.call(cbind, l_data_lags)) # combine
      data <- data_v[, -1]
      
      if(call$overparameterize) {
        
        # Compute augmented type and level vectors
        type_aug <- rep(type, n_lags)
        level_aug <- rep(level, n_lags)
        
        # Construct over-parameterized design matrix
        X <- ModelMatrix(data = data,
                         type = type_aug,
                         level = level_aug,
                         labels = colnames(data),
                         d = 1)
        
      } else {
        
        # Construct standard design matrix
        form <- as.formula('y ~ (.)')
        data_df <- data
        type_aug <- rep(type, n_lags)
        for(j in which(type_aug=='c')) data_df[, j] <- as.factor(data_df[, j])
        X <- model.matrix(form, data = data_df)[, -1] # delete intercept (added by glmnet later)
        
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
  
  # browser()
  
  
  # ----- Return -----
  
  return(predCoreObj)
  
  
} # eoF





