

tvmvar <- function(data,         # n x p data matrix
                   type,         # p vector indicating the type of each variable
                   level,        # p vector indivating the levels of each variable
                   timepoints,
                   estpoints,    # vector of points, where model should be estimated, on the scale [0:(n-max(lag))], do not need to be integers
                   bandwidth,    # choice of bandwidth
                   ...           # arguments passed to mvar
)
  
  
{
  
  
  # -------------------- Input Checks -------------------
  
  args <- list(...)
  lags <- args$lags
  
  # ----- Fill in Defaults -----
  
  if(is.null(args$lambdaSeq)) args$lambdaSeq <- NULL
  if(is.null(args$lambdaSel)) args$lambdaSel <- 'EBIC'
  if(is.null(args$lambdaFolds)) args$lambdaFolds <- 10
  if(is.null(args$lambdaGam)) args$lambdaGam <- .25
  if(is.null(args$alphaSeq)) args$alphaSeq <- 1
  if(is.null(args$alphaSel)) args$alphaSel <- 'CV'
  if(is.null(args$alphaFolds)) args$alphaFolds <- 10
  if(is.null(args$alphaGam)) args$alphaGam <- .25
  if(is.null(args$threshold)) args$threshold <- 'HW'
  if(is.null(args$method)) args$method <- 'glm'
  if(is.null(args$binarySign)) args$binarySign <- FALSE
  if(is.null(args$scale)) args$scale <- TRUE
  
  if(is.null(args$consec)) args$consec <- NULL
  
  if(is.null(args$verbatim)) args$verbatim <- FALSE
  if(is.null(args$warnings)) args$warnings <- FALSE
  if(is.null(args$saveModels)) args$saveModels <- TRUE
  if(is.null(args$saveData)) args$saveData <- FALSE
  if(is.null(args$pbar)) args$pbar <- pbar <-  TRUE
  if(is.null(args$overparameterize)) args$overparameterize <- FALSE
  if(is.null(args$signInfo)) args$signInfo <- TRUE
  
  if(missing(level)) level <- NULL
  
  if(args$verbatim) args$pbar <- FALSE
  if(args$verbatim) args$warnings <- FALSE
  
  # Switch all warnings off
  if(!args$warnings) {
    oldw <- getOption("warn")
    options(warn = -1)
  }
  
  
  # ----- Create Output Object -----
  
  tvmvar_object <- list('call' = NULL,
                        'wadj' = NULL,
                        'signs' = NULL,
                        'edgecolor' = NULL,
                        'intercepts' = NULL,
                        'tvmodels' = list())
  
  
  # ----- Copy the Call -----
  
  tvmvar_object$call <- list('data' = NULL,
                             'type' = type,
                             'level' = level,
                             'timepoints' = NULL,
                             'estpoints' = estpoints,
                             'estpointsNorm' = NULL,
                             'bandwidth' = bandwidth,
                             'lags' = args$lags,
                             'consec' = args$consec,
                             'lambdaSeq' = args$lambdaSeq,
                             'lambdaSel' = args$lambdaSel,
                             'lambdaFolds' = args$lambdaFolds,
                             'lambdaGam' = args$lambdaGam,
                             'alphaSeq' = args$alphaSeq,
                             'alphaSel' = args$alphaSel,
                             'alphaFolds' = args$alphaFolds,
                             'alphaGam' = args$alphaGam,
                             'threshold' = args$threshold,
                             'method' = args$method,
                             'binarySign' = args$binarySign,
                             'scale' = args$scale,
                             'verbatim' = args$verbatim,
                             'warnings' = args$warnings,
                             'saveModels' = args$saveModels,
                             'saveData' = args$saveData,
                             'pbar' = args$pbar,
                             'overparameterize' = args$overparameterize,
                             'signInfo' = args$signInfo)
  
  if(args$saveData) tvmvar_object$call$data <- data
  
  
  # ----- Compute Aux Variables -----
  
  n <- nrow(data)
  n_var <- n - max(lags) # this is how many rows there are after transforming the data
  p <- ncol(data)
  
  
  # ----- Checks on estpoints and bandwidth -----
  
  if(any(estpoints < 0)) stop('Estimation have to bepositive')
  if(any(estpoints>n_var)) stop('Estimation points have to be on scale [0, n - max(lags)]')
  if(bandwidth <= 0) stop('The bandwidth parameter has to be strictly positive')
  
  
  # -------------------- Compute Weightings -------------------
  
  # Define time vector: if not provided, assume equally spaced time points
  if(missing(timepoints)) {
    timevec <- seq(0, 1, length = n_var)
  } else {
    # normalize to [0,1]
    timepoints <- timepoints[-(1:max(lags))] # delete first x rows that have to be exluded by definition of VAR model
    timepoints <- timepoints - min(timepoints)
    timevec <- timepoints / max(timepoints)
  }
  tvmvar_object$call$timepoints <- timevec
  
  
  # Normalize time estimation points to interval [0,1]
  estpoints_norm <- estpoints / n_var
  tvmvar_object$call$estpointsNorm <- estpoints_norm
  no_estpoints <- length(estpoints_norm)
  
  # Compute weights foe each
  l_weights <- list()
  for(i in 1:no_estpoints) {
    l_weights[[i]] <- dnorm(timevec, mean = estpoints_norm[i], sd = bandwidth)
    # Normalize to [x,1]
    l_weights[[i]] <- l_weights[[i]] / max(l_weights[[i]])
    
    # If tvmvar is used within bwSelect: set weights to zero for test samples
    if(!is.null(args$zero_weights)) {
      # Sanity
      if(length(l_weights[[i]]) != length(args$zero_weights)) stop('Zero weights vector does not have same length as weights vector!')
      # Set to zero
      l_weights[[i]] <- l_weights[[i]] * args$zero_weights
    }
  } # end for:i (estpoints)
  
  
  # -------------------- Loop over Estpoints -------------------
  
  # Progress bar
  if(args$pbar==TRUE) pb <- txtProgressBar(min = 0, max = no_estpoints, initial = 0, char="-", style = 3)
  
  # Storage
  l_mvar_models <- list()
  
  for(i in 1:no_estpoints) {
    
    l_mvar_models[[i]] <- mvar(data = data,
                               type = type,
                               level = level,
                               lambdaSeq = args$lambdaSeq,
                               lambdaSel = args$lambdaSel,
                               lambdaFolds = args$lambdaFolds,
                               lambdaGam = args$lambdaGam,
                               alphaSeq = args$alphaSeq,
                               alphaSel = args$alphaSel,
                               alphaFolds = args$alphaFolds,
                               alphaGam = args$alphaGam,
                               lags = args$lags,
                               consec = args$consec,
                               weights = l_weights[[i]],
                               threshold = args$threshold,
                               method = args$method,
                               binarySign = args$binarySign,
                               scale = args$scale,
                               verbatim = args$verbatim,
                               pbar = FALSE,
                               warnings = args$warnings,
                               saveModels = args$saveModels,
                               saveData = args$saveData,
                               overparameterize = args$overparameterize,
                               signInfo = FALSE) # to avoid msg for each model
    
    # Update Progress Bar
    if(args$pbar==TRUE) setTxtProgressBar(pb, i)
    
  } # End for: timepoints
  
  # Save into output list
  tvmvar_object$tvmodels <- l_mvar_models
  
  
  # -------------------- Process Output -------------------
  # Restructure output to make time-varying model easier accessible
  
  n_lags <- length(lags)
  a_wadj <- a_signs <- a_edgecolor <- array(dim=c(p, p, n_lags, no_estpoints))

  for(lag in 1:n_lags) {
    
    for(i in 1:no_estpoints) {
      
      a_wadj[ , , lag, i] <- l_mvar_models[[i]]$wadj[,,lag]
      a_signs[ , , lag, i] <- l_mvar_models[[i]]$signs[,,lag]
      a_edgecolor[ , , lag, i] <- l_mvar_models[[i]]$edgecolor[,,lag]
      
    }
    
  }
  
  # Save intercepts
  l_intercepts <- list()
  for(i in 1:no_estpoints) l_intercepts[[i]] <- l_mvar_models[[i]]$intercepts
  tvmvar_object$intercepts <- l_intercepts
  
  # Save in output list
  tvmvar_object$wadj <- a_wadj
  tvmvar_object$signs <- a_signs
  tvmvar_object$edgecolor <- a_edgecolor
  
  
  # Compute effectively used Sample size (relative to n)
  Ne <- lapply(l_weights, sum)
  tvmvar_object$Ne <- unlist(Ne)
  
  
  # Copy lagged Data if saveData == TRUE (take from first time point, because always the same)
  if(args$saveData) tvmvar_object$call$data_lagged <-  l_mvar_models[[1]]$call$data_lagged
  
  
  # -------------------- Output -------------------
  
  class(tvmvar_object) <- c('mgm', 'tvmvar')
  
  if(pbar) {
    if(signInfo) cat('\nNote that the sign of parameter estimates is stored separately; see ?tvmvar')    
  } else {
    if(signInfo) cat('Note that the sign of parameter estimates is stored separately; see ?tvmvar')    
  }
  
  
  return(tvmvar_object)
  
  
}


