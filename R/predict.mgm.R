
# Used within predict.mgm to compute nice output error table

f_makeErrorTable <- function(data, 
                             l_errorCon,
                             l_errorCat,
                             l_errors_con,
                             l_errors_cat,
                             p,
                             errorCat,
                             errorCon,
                             type)
  
{
  
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
  if(length(l_errors_cat) == 0) names_cat <- NULL else names_cat <- paste0('', names(l_errors_cat))
  if(length(l_errors_con) == 0) names_con <- NULL else names_con <- paste0('', names(l_errors_con))
  
  colnames(ea) <- c('Variable', names_con, names_cat)
  
  if(length(l_errors_con) > 0) ea[,2:(1+length(l_errorCon))] <- as.numeric(round(do.call(cbind, l_errors_con),3))
  if(length(l_errors_cat) > 0)  ea[,(2+length(l_errorCon)):(length(l_errorCat)+length(l_errorCon)+1)] <- as.numeric(round(do.call(cbind, l_errors_cat),3))
  
  # Intercept performance for categoricals
  if(all(c('CC', 'nCC', 'CCmarg') %in% errorCat)) {
    
    CCmarg <- rep(NA, p)
    CC <- ea[ , which(colnames(ea) =='CC')]
    nCC <- ea[ , which(colnames(ea) =='nCC')]
    for(j in which(type=='c')) CCmarg[j] <- (nCC[j] - CC[j]) / (nCC[j] - 1)
    
    ea <- cbind(ea, round(CCmarg, 3))
    dimnames(ea)[[2]][dim(ea)[2]] <- 'CCmarg'
  }
  
  return(ea)
  
}





predict.mgm <- function(object, # One of the four mgm objects
                        data, # data in same format as used for fitting
                        errorCon, # specifying error or providing function for continuous variables
                        errorCat, # specifying error or providing function for categorical variables
                        tvMethod, # 'weighted' vs. 'closestModel'
                        consec,
                        beepvar, 
                        dayvar,
                        ...)
  
  
{
  
  
  
  # ----- Check for NAs in the data ----- 
  
  if(any(is.na(data))) stop("The data contains missing values. Please provide a complete data set without missing values.")
  
  # ----- Fill in defaults ----- 
  
  if(missing(consec)) consec <- NULL
  if(missing(beepvar)) beepvar <- NULL
  if(missing(dayvar)) dayvar <- NULL
  
  # ----- Compute consec argument (copied from mvar() ) -----
  
  # Input checks (can only specify consec OR beepvar and dayvar)
  
  if(!is.null(consec) & !is.null(beepvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  if(!is.null(consec) & !is.null(dayvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  
  if(!is.null(dayvar)) if(is.null(beepvar)) stop("Argument beepvar not specified.")
  if(!is.null(beepvar)) if(is.null(dayvar)) stop("Argument dayvar not specified.")
  
  if(!is.null(beepvar) & !is.null(dayvar)) {
    
    consec <- beepday2consec(beepvar = beepvar,
                             dayvar = dayvar)
    
  } # if: specification of consecutiveness via beepvar and dayvar
  
  
  # ----- Compute Aux Variables -----
  
  cobj <- class(object)[2]
  type <- object$call$type
  p <- ncol(data)
  n <- nrow(data)
  
  if(cobj %in% c('code', 'tvmgm'))  n_pred <- n  
  if(cobj %in% c('mvar', 'tvmvar'))  {

    max_lags <- max(object$call$lags)
    n_pred <- n - max_lags
  } 
  
  
  
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
      if(nrow(data) != length(object$tvmodels[[1]]$call$weights)) {
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
      errorCat <- c('CC', 'nCC', 'CCmarg')
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
  if(inherits(errorCon, "character")) { # this is the same as class(errorCon) == "character", but this seems not to be allowed anymore on CRAN [Jult, 2022]
    
    if('R2' %in% errorCon & !('RMSE' %in% errorCon)) l_errorCon <- list('R2' = f_error_R2)
    if('RMSE' %in% errorCon & !('R2' %in% errorCon)) l_errorCon <- list('RMSE' = f_error_RMSE)
    if(all(c('R2', 'RMSE') %in% errorCon)) l_errorCon <- list('RMSE' = f_error_RMSE, 'R2' = f_error_R2)
    
  } else if(inherits(errorCon, "list")) {
    
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
  if(inherits(errorCat, "character")) {
    
    if('CC' %in% errorCat & !('nCC' %in% errorCat)) l_errorCat <- list('CC' = f_error_CC)
    if('nCC' %in% errorCat & !('CC' %in% errorCat)) l_errorCat <- list('nCC' = f_error_nCC)
    if(all(c('CC', 'nCC') %in% errorCat)) l_errorCat <- list('CC' = f_error_CC, 'nCC' = f_error_nCC)
    
  } else if(inherits(errorCat, "list")) {
    
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
                  'errors' = NULL, 
                  'tverrors' = NULL)
  
  
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
        
        # if(!is.null(consec)) 
        if(cobj == 'tvmvar') n_pred <- sum(corePred$included) else { n_pred <- n}
        
        # Save separate for output:
        l_preds[[ep]] <- do.call(cbind, corePred$pred)
        l_probs[[ep]] <- corePred$prob
        l_true[[ep]] <- corePred$true
        
        # Save weights (subset by consec, if provided, if mVAR model)
        wo <- object_ep$call$weights
        if(cobj == 'tvmvar') wo <- wo[corePred$included] 
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
          a_w_predict_con[, , ep] <- apply(l_preds[[ep]][, p_ind_con], 2, function(x) {
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
  
  
  # ---------- Compute time-varying Nodewise Errors ----------
  
  if(tvMethod == "weighted") {
    
    l_tverrors <- list()
    
    for(ep in 1:n_estpoints) {

      object_ep <- object$tvmodels[[ep]]
      
      # --- Compute weights (needed if new data is used, otherwise I could use weights_design from the model object) ---
      
      # Get time vector as basis for weighting
      n_pred <- nrow(true)
      if(is.null(object$call$timepoints)) {
        timepoints <- seq(0, 1, length = n_pred)
      } else {
        
        a_w_predict_con[, , ep]
        
        # Adapt timepoints to time-varying mVAR/MGM      
        if(cobj == 'tvmvar') {
          timepoints <- object$call$timepoints[corePred$included]
        } else {
          timepoints <- object$call$timepoints
        }

      }
      
      # Compute weighting as function of estimation point and abndwidth
      weights <- dnorm(x = timepoints, 
                       mean = object$call$estpoints[ep], 
                       sd = object$call$bandwidth)
              

      # --- Errors Continuous ---
      
      l_errors_con_tv <- list()
      # Are any error functions specified?
      if(!is.null(l_errorCon)) {
        # Loop over error functions
        for(e in 1:length(l_errorCon)) {
          v_errors <- rep(NA, p)
          # Loop over p variables
          for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](true = true[, j], 
                                                                            pred = l_preds[[ep]][, j], 
                                                                            weights = weights)
          l_errors_con_tv[[e]] <- v_errors
        }
        names(l_errors_con_tv) <- names(l_errorCon)
      }
      
      
      # --- Errors Categorical ---
      
      l_errors_cat_tv <- list()
      # Are any error functions specified?
      if(!is.null(l_errorCat)) {
        # Loop over error functions
        for(e in 1:length(l_errorCat)) {
          v_errors <- rep(NA, p)
          # Loop over p variables
          for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](true = true[, j], 
                                                                            pred = preds[, j],
                                                                            weights = weights)
          l_errors_cat_tv[[e]] <- v_errors
        }
        names(l_errors_cat_tv) <- names(l_errorCat)
      }
      
      
      l_tverrors[[ep]] <- f_makeErrorTable(data, 
                                           l_errorCon,
                                           l_errorCat,
                                           l_errors_con_tv,
                                           l_errors_cat_tv,
                                           p,
                                           errorCat,
                                           errorCon, 
                                           type)
      
      
      
    } # end for: est points
    
    
    
    predObj$tverrors <- l_tverrors
    
    
  } # end if: tvMethod weighted?
  
  
  # Make error table for overall errors
  predObj$errors <- f_makeErrorTable(data, 
                                     l_errorCon,
                                     l_errorCat,
                                     l_errors_con,
                                     l_errors_cat,
                                     p,
                                     errorCat,
                                     errorCon,
                                     type)
  
  
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





