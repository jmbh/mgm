

tvmgm <- function(data,         # n x p data matrix
                  type,         # p vector indicating the type of each variable
                  level,        # p vector indivating the levels of each variable
                  timepoints,
                  estpoints,
                  bandwidth,
                  ...           # mgm() arguments
)

{
  if (nrow(data) == length(estpoints)) message("mgm-parallel")
  
  # -------------------- Input Checks -------------------

  args <- list(...)


  # ----- Compute Aux Variables -----

  n <- nrow(data)
  p <- ncol(data)
  no_estpoints <- length(estpoints)


  # ----- Fill in Defaults -----

  if(is.null(args$lambdaSeq)) args$lambdaSeq <- NULL
  if(is.null(args$lambdaSel)) args$lambdaSel <- 'EBIC'
  if(is.null(args$lambdaFolds)) args$lambdaFolds <- 10
  if(is.null(args$lambdaGam)) args$lambdaGam <- .25
  if(is.null(args$alphaSeq)) args$alphaSeq <- 1
  if(is.null(args$alphaSel)) args$alphaSel <- 'CV'
  if(is.null(args$alphaFolds)) args$alphaFolds <- 10
  if(is.null(args$alphaGam)) args$alphaGam <- .25
  if(is.null(args$k)) args$k <- 2
  if(is.null(args$moderators)) args$moderators <- NULL
  if(is.null(args$ruleReg)) args$ruleReg <- 'AND'
  if(is.null(args$threshold)) args$threshold <- 'HW'
  if(is.null(args$method)) args$method <- 'glm'
  if(is.null(args$binarySign)) args$binarySign <- FALSE
  if(is.null(args$scale)) args$scale <- TRUE

  if(is.null(args$verbatim)) args$verbatim <- FALSE
  if(is.null(args$warnings)) args$warnings <- FALSE
  if(is.null(args$saveModels)) args$saveModels <- TRUE
  if(is.null(args$saveData)) args$saveData <- FALSE
  if(is.null(args$pbar)) args$pbar <- pbar <-  TRUE
  if(is.null(args$overparameterize)) args$overparameterize <- FALSE
  if(is.null(args$thresholdCat)) if(args$overparameterize) args$thresholdCat <- TRUE else args$thresholdCat <- TRUE # always better
  if(is.null(args$signInfo)) args$signInfo <-  TRUE

  if(missing(level)) level <- NULL

  if(args$verbatim) args$pbar <- FALSE
  if(args$verbatim) args$warnings <- FALSE

  # Switch all warnings off
  if(!args$warnings) {
    oldw <- getOption("warn")
    options(warn = -1)
  }

  # Define time vector: if not provided, assume equally spaced time points
  if(missing(timepoints)) {
    timevec <- seq(0,1, length = n)
  } else {
    # normalize to [0,1]
    timepoints <- timepoints - min(timepoints)
    timevec <- timepoints / max(timepoints)
  }


  # ----- Storage: Create empty mgm object -----

  tvmgmobj <- list('call' = NULL,
                   'pairwise' = list('wadj' = NULL,
                                     'signs' = NULL,
                                     'edgecolor'= NULL),
                   'interactions' = list('indicator' = NULL,
                                      'weights' = NULL,
                                      'signs' = NULL),
                   'intercepts' = NULL,
                   'tvmodels' = list())


  # ----- Copy the Call -----

  tvmgmobj$call <- list('data' = NULL,
                        'type' = type,
                        'level' = level,
                        'timepoints' = timevec,
                        'estpoints' = estpoints,
                        'estpointsNorm' = NULL,
                        'bandwidth' = bandwidth,
                        'lambdaSeq' = args$lambdaSeq,
                        'lambdaSel' = args$lambdaSel,
                        'lambdaFolds' = args$lambdaFolds,
                        'lambdaGam' = args$lambdaGam,
                        'alphaSeq' = args$alphaSeq,
                        'alphaSel' = args$alphaSel,
                        'alphaFolds' = args$alphaFolds,
                        'alphaGam' = args$alphaGam,
                        'k' = args$k,
                        "moderators" = args$moderators, 
                        'ruleReg' = args$ruleReg,
                        'threshold' = args$threshold,
                        'method' = args$method,
                        'binarySign' = args$binarySign,
                        'scale' = args$scale,
                        'verbatim' = args$verbatim,
                        'pbar' = args$pbar,
                        'warnings' = args$warnings,
                        'saveModels' = args$saveModels,
                        'saveData' = args$saveData,
                        'overparameterize' = args$overparameterize,
                        "thresholdCat" = args$thresholdCat,
                        'signInfo' = args$signInfo)


  if(args$saveData) tvmgmobj$call$data <- data


  # ----- Checks on timepoints, estpoints and bandwidth -----

  if(any(estpoints < 0 | estpoints > 1)) stop('Estimation points have to be specified on the unit interval [0,1].')
  if(bandwidth <= 0) stop('The bandwidth parameter has to be strictly positive')


  # -------------------- Compute Weightings -------------------
  # Note: slightly different from tvmvar: time points don't have to be equally spaced
  # time measurements can be entered via 'timepoints' argument

  estpoints_norm <- estpoints 
  tvmgmobj$call$estpointsNorm <- estpoints_norm



  # Compute weights foe each
  l_weights <- list()
  for(i in 1:no_estpoints) {
    l_weights[[i]] <- dnorm(timevec, mean = estpoints_norm[i], sd = bandwidth)
    l_weights[[i]] <- l_weights[[i]] / max(l_weights[[i]]) # normalize to [x,1]

    # If tvmgm is used within bwSelect: set weights to zero for test samples
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
  l_tvmgm_models <- list()

  # for(i in 1:no_estpoints) {
  l_tvmgm_models <- foreach::`%dopar%`(
    foreach::foreach(i = 1:no_estpoints, .packages = "mgm"),
    {
      # l_tvmgm_models[[i]] <-
      mgm(data = data,
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
          k = args$k,
          ruleReg = args$ruleReg,
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
    # if(args$pbar==TRUE) setTxtProgressBar(pb, i)

  }) # End for: timepoints

  # Save into output list
  tvmgmobj$tvmodels <- l_tvmgm_models

  
  # -------------------- Process Output -------------------

  # Output:
  # pairwise and factorgraph: 3d arrays
  # interactions and intercepts: lists

  # col_factorgraph <- ncol(l_FG[[1]]$weightedgraph) # all the same, so i = 1 does the job

  # Create Storage
  a_pairwise_wadj <- a_pairwise_signs <- a_pairwise_edgecolor <- array(dim = c(p, p, no_estpoints))
  # a_factorgraph_graph <- a_factorgraph_signs <- a_factorgraph_edgecolor <-array(dim = c(col_factorgraph, col_factorgraph, no_estpoints))
  l_factorgraph_order <- l_factorgraph_nodetype <- l_interactions_indicator <- l_interactions_weights <- l_interactions_weightsAgg <- l_interactions_signs <- l_intercepts <- list()

  # Loop over timepoints and restructure
    for(i in 1:no_estpoints) {

      # Pairwise
      a_pairwise_wadj[,,i] <- tvmgmobj$tvmodels[[i]]$pairwise$wadj
      a_pairwise_signs[,,i] <- tvmgmobj$tvmodels[[i]]$pairwise$signs
      a_pairwise_edgecolor[,,i] <- tvmgmobj$tvmodels[[i]]$pairwise$edgecolor

      # Factorgraph Raw
      l_interactions_indicator[[i]] <- tvmgmobj$tvmodels[[i]]$interactions$indicator
      l_interactions_weightsAgg[[i]] <- tvmgmobj$tvmodels[[i]]$interactions$weightsAgg
      l_interactions_weights[[i]] <- tvmgmobj$tvmodels[[i]]$interactions$weights
      l_interactions_signs[[i]] <- tvmgmobj$tvmodels[[i]]$interactions$signs

      # Intercepts
      l_intercepts[[i]] <- tvmgmobj$tvmodels[[i]]$intercepts

    }


  # Fill into output list
  tvmgmobj$pairwise$wadj <- a_pairwise_wadj
  tvmgmobj$pairwise$signs <- a_pairwise_signs
  tvmgmobj$pairwise$edgecolor <- a_pairwise_edgecolor
  tvmgmobj$interactions$indicator <- l_interactions_indicator
  tvmgmobj$interactions$weights <- l_interactions_weights
  tvmgmobj$interactions$weightsAgg <- l_interactions_weightsAgg
  tvmgmobj$interactions$signs <- l_interactions_signs
  tvmgmobj$intercepts <- l_intercepts
  
  
  # Compute effectively used Sample size (relative to n)
  Ne <- lapply(l_weights, sum)
  tvmgmobj$Ne <- unlist(Ne)
  


  # -------------------- Output -------------------

  # Save Node Models 
  if(!args$saveModels) {
    tvmgmobj$nodemodels <- NULL
  }

  # intercepts
  # tvmgmobj$intercepts <- list_Intercepts

  # Switch warings back on
  if(!args$warnings) {
    options(warn = oldw)
  }

  if(args$pbar) {
    if(args$signInfo) cat('\nNote that the sign of parameter estimates is stored separately; see ?tvmgm')    
  } else {
    if(args$signInfo) cat('Note that the sign of parameter estimates is stored separately; see ?tvmgm')    
  }
  

  # Assign class
  class(tvmgmobj) <- c('mgm', 'tvmgm')

  return(tvmgmobj)


} # eoF

