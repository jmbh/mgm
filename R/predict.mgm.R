

predict.mgm <- function(object, # One of the four mgm objects
                        data, # data in same format as used for fitting
                        errorCon, # specifying error or providing function for continuous variables
                        errorCat, # specifying error or providing function for categorical variables
                        tvMethod, # 'weighted' vs. 'closestModel'
                        consec,
                        ...)
  
  
{
  
  
  # ----- Fill in defaults ----- 
  
  if(missing(consec)) consec <- NULL
  
  
  # ----- Compute Aux Variables -----
  
  cobj <- class(object)[2]
  type <- object$call$type
  p <- ncol(data)
  n <- nrow(data)
  
  if(cobj %in% c('code', 'tvmgm'))  n_pred <- n  
  if(cobj %in% c('mvar', 'tvmvar'))  n_pred <- n - max(object$call$lags)
  
  
  
  level <- object$call$level
  
  # Scale Gaussians if scale==TRUE
  if(object$call$scale) {
    for(j in which(type=='g')) data[, j] <- scale(data[, j])
  }
  
  # ---------- Input Checks ----------
  
  # Checks for all situations
  
  if(!(cobj %in% c('core', 'mvar', 'tvmgm', 'tvmvar'))) stop('Please provide a mgm fit object.')
  
  
  # Checkes for stationary
  
  if(cobj %in% c('core', 'mvar')) {
    
    # is the model object saved?
    if(is.null(object$nodemodels)) stop(paste0('Prediction is only possible if the model-object is saved. See ?', cobj))
    
    # does data have same number of columns as the fitting data?
    if(!(ncol(data) == length(object$call$type))) stop('The data used for prediction has to have the same number of variables as the data used for estimation.')
    
  }
  
  
  # Checks for stationary/time-varying mVAR models
  if(cobj %in% c('tvmvar', 'mvar')) {
    
    # Does consec match nrow(data) ?
    if(!is.null(consec)) if(length(consec) != nrow(data)) stop("The length of consec has to be equal to the number of rows of data")
    
  }
  
  # Checks only for time-varying objects
  if(cobj %in% c('tvmgm', 'tvmvar')) {
    
    if(missing(tvMethod)) stop('Specify the type of error for the time-varying model, see ?tvmgm or ?tvmvar!')
    
    # Has the provided time-series data the same lenght as the time-series used for fitting?
    if(cobj == 'tvmvar') {
      if(nrow(data)-max(object$tvmodels[[1]]$call$lags) != length(object$tvmodels[[1]]$call$weights)) {
        stop('For time-varying models the time-series used for prediction has to have the same lenght as the time-series used for estimation.')
      }
    }
    
    if(cobj == 'tvmgm') {
      if(nrow(data) != length(object$tvmodels[[1]]$call$weights)) {
        stop('For time-varying models the time-series used for prediction has to have the same lenght as the time-series used for estimation.')
      }
    }
    
    # does data have same number of columns as the fitting data?
    if(!(ncol(data) == length(object$tvmodels[[1]]$call$type))) stop('The data used for prediction has to have the same number of variables as the data used for estimation.')
    
    # Is the model object saved?
    if(is.null(object$tvmodels[[1]]$nodemodels)) stop(paste0('Prediction is only possible if the model-object is saved. See ?', cobj))
    
    
  } else {
    
    tvMethod <- 'NA'
    
  }
  
  
  
  # ----- Fill in defaults -----
  
  # todo: if error types missing, don't compute errors!
  if(missing(errorCon)) {
    if(any(c('g', 'p') %in% type)) {
      errorCon <- c('R2', 'RMSE')
    } else {
      errorCon <- NULL
    }
  }
  
  
  if(missing(errorCat)) {
    if('c' %in% type) {
      errorCat <- c('CC', 'nCC')
    } else {
      errorCat <- NULL
    }
  }
  
  
  # ----- Define Error Functions -----
  
  # Proportion Variance Explained
  f_error_R2 <- function(true, pred, weights = NULL) {
    
    if(is.null(weights)) weights <- rep(1, length(true))
    
    res <- true - pred
    R2 <- 1 - wtd.var(res, weights)/wtd.var(true, weights)
    return(R2)
  }
  
  # Root Mean Squared Error
  f_error_RMSE <- function(true, pred, weights = NULL) {
    
    if(is.null(weights)) weights <- rep(1, length(true))
    
    weights_norm <- weights / sum(weights)
    RMSE <- sqrt(sum(weights_norm * (true - pred)^2))
    
    return(RMSE)
  }
  
  # Proportion Correct Classification / Accuracy
  f_error_CC <- function(true, pred, weights = NULL) {
    
    if(is.null(weights)) weights <- rep(1, length(true))
    weights_norm <- weights / sum(weights)
    
    CC <- sum((true == pred)*weights_norm)
    return(CC)
  }
  
  
  # Normalized Accuracy (Accuracy full model minus Accuracy intercept model)
  f_error_nCC <- function(true, pred, weights = NULL) {
    
    if(is.null(weights)) weights <- rep(1, length(true))
    weights_norm <- weights / sum(weights)
    
    CC <- sum((true == pred)*weights_norm)
    tb <- table(true)
    norm_constant <- max(tb/sum(tb)) # normalize by maximum additional accuracy that can be predicted beyond the larger marginal frequency
    nCC <- (CC - norm_constant) / (1 - norm_constant)
    return(nCC)
  }
  
  
  
  
  # Create list for loss functions
  
  # Continuous
  l_errorCon <- list()
  if(class(errorCon)=='character') {
    
    if('R2' %in% errorCon & !('RMSE' %in% errorCon)) l_errorCon <- list('R2' = f_error_R2)
    if('RMSE' %in% errorCon & !('R2' %in% errorCon)) l_errorCon <- list('RMSE' = f_error_RMSE)
    if(all(c('R2', 'RMSE') %in% errorCon)) l_errorCon <- list('RMSE' = f_error_RMSE, 'R2' = f_error_R2)
    
  } else if(class(errorCon)=='list') {
    
    l_errorCon <- errorCon
    
    for(e in 1:length(errorCon)) {
      if(l_errorCon[[e]] == 'RMSE') l_errorCon[[e]] <- f_error_RMSE
      if(l_errorCon[[e]] == 'R2') l_errorCon[[e]] <- f_error_R2
    }
    
  } else if(is.null(errorCon)) {
    
    l_errorCon <- NULL
    
  } else {
    stop('errorCon must either be a list of supported error functions or a list of supported and/or custom error functions.')
  }
  
  
  
  # Categorical
  l_errorCat <- list()
  if(class(errorCat)=='character') {
    
    if('CC' %in% errorCat & !('nCC' %in% errorCat)) l_errorCat <- list('CC' = f_error_CC)
    if('nCC' %in% errorCat & !('CC' %in% errorCat)) l_errorCat <- list('nCC' = f_error_nCC)
    if(all(c('CC', 'nCC') %in% errorCat)) l_errorCat <- list('CC' = f_error_CC, 'nCC' = f_error_nCC)
    
  } else if(class(errorCat)=='list') {
    
    for(e in 1:length(errorCat)) {
      if(l_errorCat[[e]] == 'CC') l_errorCat[[e]] <- f_error_CC
      if(l_errorCat[[e]] == 'nCC') l_errorCat[[e]] <- f_error_nCC
    }
    
  } else if(is.null(errorCat)) {
    
    l_errorCat <- NULL
    
  } else {
    stop('errorCat must either be a list of supported error functions or a list of supported and/or custom error functions.')
  }
  
  
  
  # ----- Create Prediction Object -----
  
  predObj <- list('call' = NULL,
                  'predicted' = NULL,
                  'probabilties' = NULL,
                  'true' = NULL,
                  'errors' = NULL)
  
  
  # Copy the Call
  
  predObj$call <- list('object' = object,
                       'data' = data,
                       'errorCon' = l_errorCon,
                       'errorCat' = l_errorCat,
                       'tvMethod' = tvMethod)
  
  
  
  
  # ---------- A) Stationary Models ----------
  
  if(cobj %in% c('core', 'mvar')) {
    
    # ----- Make Predictions -----
    
    corePred <- predictCore_stat(object = object,
                                 data = data,
                                 consec = consec)
    
    m_pred <- do.call(cbind, corePred$pred) # Collapse predictions in matrix
    
    preds <- m_pred
    probs <- corePred$prob
    true <- corePred$true
    
  } # end if: stationary?
  
  
  # ---------- B) Time-varying Models ----------
  
  if(cobj %in% c('tvmgm', 'tvmvar')) {
    
    if(cobj == 'tvmgm') cobj_ep <- c('mgm', 'core')
    if(cobj == 'tvmvar') cobj_ep <- c('mgm', 'mvar')
    
    n_estpoints <- length(object$call$estpoints)
    
    
    # ++++++++++ Prediction Option B.1: 'weighted' ++++++++++
    
    l_errors_ep_con <- l_errors_ep_cat <- l_preds <- l_probs <- l_true <- l_weights <- vector('list', length = n_estpoints) # Storage
    
    if(tvMethod == 'weighted') {
      
      for(ep in 1:n_estpoints) {
        
        # ----- Make Predictions -----
        
        
        object_ep <- object$tvmodels[[ep]]
        class(object_ep) <- cobj_ep
        
        corePred <- l_preds[[ep]] <- predictCore_stat(object = object_ep,
                                                      data = data, 
                                                      consec = consec)
        
        if(!is.null(consec)) n_pred <- sum(corePred$included)
        
        # Save separate for output:
        l_preds[[ep]] <- do.call(cbind, corePred$pred)
        l_probs[[ep]] <- corePred$prob
        l_true[[ep]] <- corePred$true
        
        # Save weights (subset by consec, if provided)
        wo <- object_ep$call$weights
        if(!is.null(consec)) wo <- wo[corePred$included]
        l_weights[[ep]] <- wo
        
      } # end for: estpoints
      
      
      # ----- Collapse Predictions via Weighting -----
      
      # add up weights
      m_weights <- do.call(cbind, l_weights)
      
      v_weights <- rowSums(m_weights)
      
      p_ind_con <- which(type != 'c')
      p_ind_cat <- which(type == 'c')
      
      
      # Storage
      l_w_predict_cat <- list()
      a_w_predict_con <- array(NA, dim = c(n_pred, length(p_ind_con), n_estpoints))
      
      
      # --- Continuous Variables ---
      
      if(length(p_ind_con) > 0) {
        
        for(ep in 1:n_estpoints) {
          
          object_ep <- object$tvmodels[[ep]]
          
          # For continuous variables (just convex combination of values)
          a_w_predict_con[, , ep] <- apply(l_preds[[ep]][,p_ind_con], 2, function(x) {
            w_pred_col <- x*m_weights[, ep]
            return(w_pred_col)
          })
          
        }
        
        # Continuous: sum up and scale back
        m_w_predict_con1 <- apply(a_w_predict_con, 1:2, sum)
        m_w_predict_con2 <- apply(m_w_predict_con1, 2, function(x)  {
          
          cb <- cbind(x, v_weights)
          apply(cb, 1, function(x)  x[1] / x[2])
          
        })
        
      } # if: any continuous variables?
      
      
      # --- Categorical Variables ---
      
      if(length(p_ind_cat) > 0) {
        
        # For categorical variables (convex combination of probabilities)
        l_w_predict_cat <- list()
        cat_counter <- 1
        
        # for all categorical variables
        for(v in p_ind_cat) {
          
          a_w_predict_cat <- array(NA, dim = c(n_pred, level[v], n_estpoints))
          
          # for all estimation points
          for(ep in 1:n_estpoints) {
            
            a_w_predict_cat[, , ep] <- l_probs[[ep]][[v]]
            a_w_predict_cat[, , ep] <- apply(a_w_predict_cat[, , ep], 2, function(x) x * m_weights[, ep])
            
          }
          
          # one array for each categorical variable
          l_w_predict_cat[[cat_counter]] <- a_w_predict_cat        
          cat_counter <- cat_counter + 1
          
        }
        
        # Get category labels
        l_cat_labels <- list()
        for(i in 1:length(p_ind_cat)) l_cat_labels[[i]] <- sort(unique(data[, p_ind_cat[i]]))
        
        # Sum up and scale back probabilities
        l_probs_agg <- lapply(l_w_predict_cat, function(x) {
          p_agg1 <- apply(x, 1:2, sum)
          p_agg2 <- apply(p_agg1, 2, function(x) x / v_weights)
          
        })
        
        # Predict category
        l_cat_predicted <-  list()
        for(i in 1:length(p_ind_cat)) {
          l_cat_predicted[[i]] <- l_cat_labels[[i]][apply(l_probs_agg[[i]], 1, which.max)]
        }
        m_cat_predicted <- do.call(cbind, l_cat_predicted)
        
      } # if: any categorical variables?
      
      
      # --- Putting all back together ---
      
      m_preds_final <- matrix(NA, n_pred, p)
      if(length(p_ind_con) > 0)  m_preds_final[, p_ind_con] <- m_w_predict_con2
      if(length(p_ind_cat) > 0) m_preds_final[, p_ind_cat] <- m_cat_predicted
      
      # Bring in uniform output format
      preds <- m_preds_final
      if(length(p_ind_cat) > 0) probs <- l_probs_agg else probs <- NULL
      true <- corePred$true
      
      
      
      # ----- Calculate Errors -----
      
      
      # I think can now be done as in stationary case; maybe simplify code, s.t. error computation is outside the tv/stationary method
      
      
      # m_pred <- do.call(cbind, corePred$pred) # Collapse predictions in matrix
      # v_error <- v_errortype <- rep(NA, p) # Create Storage
      # 
      # 
      # # Errors Continuous
      # l_errors_con <- list()
      # if(!is.null(l_errorCon)) {
      #   for(e in 1:length(l_errorCon)) {
      #     v_errors <- rep(NA, p)
      #     for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](true = corePred$true[, j],
      #                                                                       pred = m_pred[, j],
      #                                                                       weights = object$tvmodels[[ep]]$call$weights)
      #     l_errors_con[[e]] <- v_errors
      #   }
      #   names(l_errors_con) <- names(l_errorCon)
      # }
      # l_errors_ep_con[[ep]] <- l_errors_con
      # 
      # # Errors Categorical
      # l_errors_cat <- list()
      # if(!is.null(l_errorCat)) {
      #   for(e in 1:length(l_errorCat)) {
      #     v_errors <- rep(NA, p)
      #     for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](true = corePred$true[, j],
      #                                                                       pred = m_pred[, j],
      #                                                                       weights = object$tvmodels[[ep]]$call$weights)
      #     l_errors_cat[[e]] <- v_errors
      #   }
      #   names(l_errors_cat) <- names(l_errorCat)
      # }
      # l_errors_ep_cat[[ep]] <- l_errors_cat
      
      
      
      # Aggregate: across estimation points - yes!
      
      
    } # end if: tvMethod weighted?
    
    
    
    # ++++++++++ Prediction Option B.1: 'closestModel' ++++++++++
    
    if(tvMethod == 'closestModel') {
      
      l_preds <- l_true <- l_probs <- vector('list', length = n_estpoints) # Storage
      
      # ----- Make Predictions -----
      
      for(ep in 1:n_estpoints) {
        
        object_ep <- object$tvmodels[[ep]]
        class(object_ep) <- cobj_ep
        
        # we fit each model on all data to avoid problems with scaling and overlap with predictors because of the VAR model
        corePred <- predictCore_stat(object = object_ep,
                                     data = data, 
                                     consec = consec) # only data closest to model at current time point
        
        
        # --- Determine closest Model for all time points ---
        
        # only for first ep, because invariant from 1:n_estpoints
        if(ep == 1) {
          
          v_timepoints <- object$call$timepoints # for both tvmgm and tvmvar
          n_timepoints <- length(v_timepoints)
          
          m_assign <- matrix(NA,
                             nrow = n_timepoints,
                             ncol = n_estpoints)
          for(tp in 1:n_timepoints) m_assign[tp, ] <- abs(v_timepoints[tp] - object$call$estpointsNorm)
          
          v_assign <- apply(m_assign, 1, which.min)
          
          # if consec is defined, subset v_assign accordingly
          if(!is.null(consec)) v_assign <- v_assign[corePred$included]
        }
        
        
        # Copy relevant predictions, defined by v_assign
        if('core' %in% cobj_ep) ind_batch <- v_assign == ep
        if('mvar' %in% cobj_ep) ind_batch <- v_assign[-c(1:max(object$call$lags))] == ep
        
        # if consec is defined, only output predictions for valid rows
        l_preds[[ep]] <- do.call(cbind, corePred$pred)[ind_batch,]
        l_true[[ep]] <- corePred$true[ind_batch,]
        l_probs[[ep]] <- lapply(corePred$prob, function(x) x[ind_batch,])
        
      } # end for: estpoints
      
      # Combine predictions
      preds <- do.call(rbind, l_preds)
      probs <- l_probs
      true <- do.call(rbind, l_true)
      
      
      # ----- Calculate Errors (as in stationary case, actually) -----
      
      # # Errors Continuous
      # l_errors_con <- list()
      # if(!is.null(l_errorCon)) {
      #   for(e in 1:length(l_errorCon)) {
      #     v_errors <- rep(NA, p)
      #     for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](TrueAll[, j], PredAll[, j])
      #     l_errors_con[[e]] <- v_errors
      #   }
      #   names(l_errors_con) <- names(l_errorCon)
      # }
      # 
      # # Errors Categorical
      # l_errors_cat <- list()
      # if(!is.null(l_errorCat)) {
      #   for(e in 1:length(l_errorCat)) {
      #     v_errors <- rep(NA, p)
      #     for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](TrueAll[, j], PredAll[, j])
      #     l_errors_cat[[e]] <- v_errors
      #   }
      #   names(l_errors_cat) <- names(l_errorCat)
      # }
      
      
    } # end if: tvMethod closestModel?
    
    
  } # end if: time-varying?
  
  
  # ---------- Compute Nodewise Errors ----------
  
  # input: true and preds!
  
  # Errors Continuous
  l_errors_con <- list()
  if(!is.null(l_errorCon)) {
    for(e in 1:length(l_errorCon)) {
      v_errors <- rep(NA, p)
      for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](true[, j], preds[, j])
      l_errors_con[[e]] <- v_errors
    }
    names(l_errors_con) <- names(l_errorCon)
  }
  
  
  # Errors Categorical
  l_errors_cat <- list()
  if(!is.null(l_errorCat)) {
    for(e in 1:length(l_errorCat)) {
      v_errors <- rep(NA, p)
      for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](true[, j], preds[, j])
      l_errors_cat[[e]] <- v_errors
    }
    names(l_errors_cat) <- names(l_errorCat)
  }
  
  
  
  
  
  # ---------- Compute Error Table for Output ----------
  
  # This commented out code creates an array instead a matrix, with error measures for different time points
  
  # if(cobj %in% c('tvmgm', 'tvmvar') & tvMethod=='weighted') {
  #   
  #   # ... if time-varying & tvMethod==weighted: array
  #   
  #   if(all(c('CC', 'nCC', 'CCmarg') %in% errorCat)) exCCmarg <- 1 else exCCmarg <- 0
  #   
  #   ea <- array(0, dim = c(p, 1+length(l_errorCon)+length(l_errorCat)+exCCmarg, n_estpoints))
  #   ea[, 1,] <- 1:p
  #   if(length(l_errors_cat) == 0) names_cat <- NULL else names_cat <- paste0('Error.', names(l_errors_cat))
  #   if(length(l_errors_con) == 0) names_con <- NULL else names_con <- paste0('Error.', names(l_errors_con))
  #   
  #   if(length(l_errors_con) > 0) for(ep in 1:n_estpoints) ea[,2:(1+length(l_errorCon)), ep] <- round(do.call(cbind, l_errors_ep_con[[ep]]),3)
  #   if(length(l_errors_cat) > 0) for(ep in 1:n_estpoints) ea[,(2+length(l_errorCon)):(length(l_errorCat)+length(l_errorCon)+1), ep] <- round(do.call(cbind, l_errors_ep_cat[[ep]]),3)
  #   
  #   if(all(c('CC', 'nCC', 'CCmarg') %in% errorCat)) dimnames(ea)[[2]] <- c('Variable', names_con, names_cat, 'CCmarg') else dimnames(ea)[[2]] <- c('Variable', names_con, names_cat)
  #   
  #   
  #   # Intercept performance for categoricals
  #   if(all(c('CC', 'nCC', 'CCmarg') %in% errorCat)) {
  #     
  #     for(ep in 1:n_estpoints) {
  #       CCmarg <- rep(NA, p)
  #       CC <- ea[ , which(colnames(ea) =='Error.CC'), ep]
  #       nCC <- ea[ , which(colnames(ea) =='Error.nCC'), ep]
  #       for(j in which(type=='c')) CCmarg[j] <- (nCC[j] - CC[j]) / (nCC[j] - 1)
  #       
  #       ea[, dim(ea)[2], ep] <- round(CCmarg, 3)
  #     }
  #   }
  #   
  #   predObj$errors <- ea
  #   
  #   
  # } else {
  
  
  # get colnames if available, otherwise fill in 1:p
  if(is.null(colnames(data))) {
    cnames <- 1:p
  } else {
    cnames <- colnames(data)
  }
  
  # create matrix with errords depending on specified errors
  ea <- as.data.frame(matrix(0, nrow = p, ncol = 1+length(l_errorCon)+length(l_errorCat)))
  ea[, 1] <- cnames
  
  # Get column names of Errors
  if(length(l_errors_cat) == 0) names_cat <- NULL else names_cat <- paste0('Error.', names(l_errors_cat))
  if(length(l_errors_con) == 0) names_con <- NULL else names_con <- paste0('Error.', names(l_errors_con))
  
  colnames(ea) <- c('Variable', names_con, names_cat)
  
  if(length(l_errors_con) > 0) ea[,2:(1+length(l_errorCon))] <- as.numeric(round(do.call(cbind, l_errors_con),3))
  if(length(l_errors_cat) > 0)  ea[,(2+length(l_errorCon)):(length(l_errorCat)+length(l_errorCon)+1)] <- as.numeric(round(do.call(cbind, l_errors_cat),3))
  
  # Intercept performance for categoricals
  if(all(c('CC', 'nCC', 'CCmarg') %in% errorCat)) {
    
    CCmarg <- rep(NA, p)
    CC <- ea[ , which(colnames(ea) =='Error.CC')]
    nCC <- ea[ , which(colnames(ea) =='Error.nCC')]
    for(j in which(type=='c')) CCmarg[j] <- (nCC[j] - CC[j]) / (nCC[j] - 1)
    
    ea <- cbind(ea, round(CCmarg, 3))
    dimnames(ea)[[2]][dim(ea)[2]] <- 'CCmarg'
  }
  
  predObj$errors <- ea
  

  
  # }
  
  
  # ---------- Prepare Remaining Output ----------
  
  predObj$predicted <- preds
  predObj$probabilties <- probs
  predObj$true <- true
  
  # add included vector for VAR models
  if(cobj %in% c('mvar', 'tvmvar')) predObj$included <- corePred$included else predObj$included <- NULL
  
  
  # ---------- Output ----------
  
  class(predObj) <- c(cobj, 'predicted')
  
  return(predObj)
  
  
} # eoF





