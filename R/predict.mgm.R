

predict.mgm <- function(object, # One of the four mgm objects
                        data, # data in same format as used for fitting
                        errorCon, # specifying error or providing function for continuous variables
                        errorCat, # specifying error or providing function for categorical variables
                        tvErrorType, # 'weighted' vs. 'closestModel'
                        ...)


{


  # ----- Compute Aux Variables -----

  cobj <- class(object)
  type <- object$call$type
  p <- ncol(data)


  # ---------- Input Checks ----------

  if(!(cobj %in% c('mgm', 'mvar', 'tvmgm', 'tvmvar'))) stop('Please provide a mgm fit object.')


  if(cobj %in% c('mgm', 'mvar')) {

    # is the model object saved?
    if(is.null(object$nodemodels)) stop(paste0('Prediction is only possible if the model-object is saved. See ?', cobj))

    # does data have same number of columns as the fitting data?
    if(!(ncol(data) == length(object$call$type))) stop('The data used for prediction has to have the same number of variables as the data used for estimation.')

  }


  if(cobj %in% c('tvmgm', 'tvmvar')) {

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


  }


  # Checks only for time-varying objects
  if(cobj %in% c('tvmgm', 'tvmvar')) {

    if(missing(tvErrorType)) stop('Specify the type of error for the time-varying model, see ?tvmgm or ?tvmvar!')

  } else {

    tvErrorType <- 'NA'

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

    if('CC' %in% errorCon & !('nCC' %in% errorCon)) l_errorCat <- list('CC' = f_error_CC)
    if('nCC' %in% errorCon & !('CC' %in% errorCon)) l_errorCat <- list('nCC' = f_error_nCC)
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
                       'tvErrorType' = tvErrorType)



  # ---------- A) Stationary Models ----------

  if(cobj %in% c('mgm', 'mvar')) {

    # ----- Make Predictions -----

    corePred <- predictCore_stat(object = object,
                                 data = data)


    # ----- Calculate Errors -----

    m_pred <- do.call(cbind, corePred$pred) # Collapse predictions in matrix

    # Errors Continuous
    l_errors_con <- list()
    if(!is.null(l_errorCon)) {
      for(e in 1:length(l_errorCon)) {
        v_errors <- rep(NA, p)
        for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](corePred$true[, j], m_pred[, j])
        l_errors_con[[e]] <- v_errors
      }
      names(l_errors_con) <- names(l_errorCon)
    }


    # Errors Categorical
    l_errors_cat <- list()
    if(!is.null(l_errorCat)) {
      for(e in 1:length(l_errorCat)) {
        v_errors <- rep(NA, p)
        for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](corePred$true[, j], m_pred[, j])
        l_errors_cat[[e]] <- v_errors
      }
      names(l_errors_cat) <- names(l_errorCat)
    }


    # So we have consistent variable names
    l_preds <- m_pred
    l_probs <- corePred$prob
    l_true <- corePred$true

  } # end if: stationary?







  # ---------- B) Time-varying Models ----------

  if(cobj %in% c('tvmgm', 'tvmvar')) {

    if(cobj == 'tvmgm') cobj_ep <- 'mgm'
    if(cobj == 'tvmvar') cobj_ep <- 'mvar'

    n_estpoints <- length(object$call$estpoints)


    # ++++++++++ Prediction Option B.1: 'weighted' ++++++++++

    l_errors_ep_con <- l_errors_ep_cat <- l_preds <- l_probs <- l_true <- vector('list', length = n_estpoints) # Storage

    if(tvErrorType == 'weighted') {

      for(ep in 1:n_estpoints) {

        # ----- Make Predictions -----

        object_ep <- object$tvmodels[[ep]]
        class(object_ep) <- cobj_ep

        corePred <- l_preds[[ep]] <- predictCore_stat(object = object_ep,
                                                      data = data)

        # Save separate for output:
        l_preds[[ep]] <- do.call(cbind, corePred$pred)
        l_probs[[ep]] <- corePred$prob
        l_true[[ep]] <- corePred$true


        # ----- Calculate Errors -----

        m_pred <- do.call(cbind, corePred$pred) # Collapse predictions in matrix
        v_error <- v_errortype <- rep(NA, p) # Create Storage


        # Errors Continuous
        l_errors_con <- list()
        if(!is.null(l_errorCon)) {
          for(e in 1:length(l_errorCon)) {
            v_errors <- rep(NA, p)
            for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](true = corePred$true[, j],
                                                                              pred = m_pred[, j],
                                                                              weights = object$tvmodels[[ep]]$call$weights)
            l_errors_con[[e]] <- v_errors
          }
          names(l_errors_con) <- names(l_errorCon)
        }
        l_errors_ep_con[[ep]] <- l_errors_con

        # Errors Categorical
        l_errors_cat <- list()
        if(!is.null(l_errorCat)) {
          for(e in 1:length(l_errorCat)) {
            v_errors <- rep(NA, p)
            for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](true = corePred$true[, j],
                                                                              pred = m_pred[, j],
                                                                              weights = object$tvmodels[[ep]]$call$weights)
            l_errors_cat[[e]] <- v_errors
          }
          names(l_errors_cat) <- names(l_errorCat)
        }
        l_errors_ep_cat[[ep]] <- l_errors_cat


      } # end for: estpoints

      # Aggregate: across estimation points - no!

    } # end if: tvErrorType weighted?



    # ++++++++++ Prediction Option B.1: 'closestModel' ++++++++++

    if(tvErrorType == 'closestModel') {

      l_preds <- l_true <- l_probs <- vector('list', length = n_estpoints) # Storage

      # ----- Determine closest Model for all time points -----

      v_timepoints <- object$call$timepoints # for both tvmgm and tvmvar
      n_timepoints <- length(v_timepoints)

      m_assign <- matrix(NA,
                         nrow = n_timepoints,
                         ncol = n_estpoints)
      for(tp in 1:n_timepoints) m_assign[tp, ] <- abs(v_timepoints[tp] - object$call$estpointsNorm)

      v_assign <- apply(m_assign, 1, which.min)


      # ----- Make Predictions -----

      for(ep in 1:n_estpoints) {

        # determine closest time points

        object_ep <- object$tvmodels[[ep]]
        class(object_ep) <- cobj_ep

        # we fit each model on all data to avoid problems with scaling and overlap with predictors because of the VAR model
        codePred <- predictCore_stat(object = object_ep,
                                     data = data) # only data closest to model at current time point

        # Copy relevant predictions, defined by v_assign
        if(cobj_ep == 'mgm') ind_batch <- v_assign==ep
        if(cobj_ep == 'mvar') ind_batch <- v_assign[-c(1:max(object$call$lags))]==ep

        l_preds[[ep]] <- do.call(cbind, codePred$pred)[ind_batch,]
        l_true[[ep]] <- codePred$true[ind_batch,]
        l_probs[[ep]] <- lapply(codePred$prob, function(x) x[ind_batch,])

      } # end for: estpoints


      # Combine predictions
      PredAll <- do.call(rbind, l_preds)
      TrueAll <- do.call(rbind, l_true)


      # ----- Calculate Errors (as in stationary case, actually) -----

      # Errors Continuous
      l_errors_con <- list()
      if(!is.null(l_errorCon)) {
        for(e in 1:length(l_errorCon)) {
          v_errors <- rep(NA, p)
          for(j in 1:p)  if(type[j] != 'c')  v_errors[j] <- l_errorCon[[e]](TrueAll[, j], PredAll[, j])
          l_errors_con[[e]] <- v_errors
        }
        names(l_errors_con) <- names(l_errorCon)
      }

      # Errors Categorical
      l_errors_cat <- list()
      if(!is.null(l_errorCat)) {
        for(e in 1:length(l_errorCat)) {
          v_errors <- rep(NA, p)
          for(j in 1:p)  if(type[j] == 'c')  v_errors[j] <- l_errorCat[[e]](TrueAll[, j], PredAll[, j])
          l_errors_cat[[e]] <- v_errors
        }
        names(l_errors_cat) <- names(l_errorCat)
      }

      # So we have consistent variable names
      l_true <- TrueAll
      l_preds <- PredAll


    } # end if: tvErrorType closestModel?


  } # end if: time-varying?


  # ---------- Compute Error Matrix/Array ----------


  if(cobj %in% c('tvmgm', 'tvmvar') & tvErrorType=='weighted') {

    # ... if time-varying & tvErrorType==weighted: array

    ea <- array(0, dim = c(p, 1+length(l_errorCon)+length(l_errorCat), n_estpoints))
    ea[, 1,] <- 1:p
    if(length(l_errors_cat) == 0) names_cat <- NULL else names_cat <- paste0('Error.', names(l_errors_cat))
    if(length(l_errors_con) == 0) names_con <- NULL else names_con <- paste0('Error.', names(l_errors_con))
    dimnames(ea)[[2]] <- c('Variable', names_con, names_cat)

    if(length(l_errors_con) > 0) for(ep in 1:n_estpoints) ea[,2:(1+length(l_errorCon)), ep] <- round(do.call(cbind, l_errors_ep_con[[ep]]),3)
    if(length(l_errors_cat) > 0) for(ep in 1:n_estpoints) ea[,(2+length(l_errorCon)):(length(l_errorCat)+length(l_errorCon)+1), ep] <- round(do.call(cbind, l_errors_ep_cat[[ep]]),3)
    predObj$errors <- ea

  } else {

    # ... if stationary or if time-varying & tvErrorType==closestModel: matrix
    ea <- array(0, dim = c(p, 1+length(l_errorCon)+length(l_errorCat)))
    ea[, 1] <- 1:p
    if(length(l_errors_cat) == 0) names_cat <- NULL else names_cat <- paste0('Error.', names(l_errors_cat))
    if(length(l_errors_con) == 0) names_con <- NULL else names_con <- paste0('Error.', names(l_errors_con))
    dimnames(ea)[[2]] <- c('Variable', names_con, names_cat)

    # browser()

    if(length(l_errors_con) > 0) ea[,2:(1+length(l_errorCon))] <- round(do.call(cbind, l_errors_con),3)
    if(length(l_errors_cat) > 0)  ea[,(2+length(l_errorCon)):(length(l_errorCat)+length(l_errorCon)+1)] <- round(do.call(cbind, l_errors_cat),3)
    predObj$errors <- ea

  }


  # ---------- Prepare Remaining Output ----------

  predObj$predicted <- l_preds
  predObj$probabilties <- l_probs
  predObj$true <- l_true


  # ---------- Output ----------

  class(predObj) <- c(cobj, 'predicted')

  return(predObj)


} # eoF





