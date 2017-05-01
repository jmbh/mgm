

tvmgm <- function(data,         # n x p data matrix
                  type,         # p vector indicating the type of each variable
                  level,        # p vector indivating the levels of each variable
                  timepoints,
                  estpoints,
                  bandwidth,
                  ...           # mgm() arguments
)

{

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
  if(is.null(args$overparameterize)) args$overparameterize <-  TRUE
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
                   'factorgraph' = list('graph' = NULL,
                                        'signs' = NULL,
                                        'edgecolor' = NULL,
                                        'order' = NULL,
                                        'nodetype' = NULL),
                   'rawfactor' = list('indicator' = NULL,
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
                        'ruleReg' = args$ruleReg,
                        'threshold' = args$threshold,
                        'method' = args$method,
                        'binarySign' = args$binarySign,
                        'verbatim' = args$verbatim,
                        'pbar' = args$pbar,
                        'warnings' = args$warnings,
                        'saveModels' = args$saveModels,
                        'saveData' = args$saveData,
                        'overparameterize' = args$overparameterize,
                        'signInfo' = args$signInfo)


  if(args$saveData) tvmgmobj$call$data <- data


  # ----- Checks on timepoints, estpoints and bandwidth -----

  if(any(estpoints < 0)) stop('Estimation have to bepositive')
  if(any(estpoints>n)) stop('Estimation points have to be on scale [0, n]')
  if(bandwidth <= 0) stop('The bandwidth parameter has to be strictly positive')


  # -------------------- Compute Weightings -------------------
  # Note: slightly different from tvmvar: time points don't have to be equally spaced
  # time measurements can be entered via 'timepoints' argument


  # Normalize time estimation points to interval [0,1]
  estpoints_norm <- estpoints / n # cannot be normalized, because we don't know the end point
  tvmgmobj$call$estpointsNorm <- estpoints_norm



  # Compute weights foe each
  l_weights <- list()
  for(i in 1:no_estpoints) {
    l_weights[[i]] <- dnorm(timevec, mean = estpoints_norm[i], sd = bandwidth)

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

  for(i in 1:no_estpoints) {

    l_tvmgm_models[[i]] <- mgm(data = data,
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
                               signInfo = args$signInfo)


    # Update Progress Bar
    if(args$pbar==TRUE) setTxtProgressBar(pb, i)

  } # End for: timepoints

  # Save into output list
  tvmgmobj$tvmodels <- l_tvmgm_models


  # -------------------- Create time-varying Factorgraph -------------------

  d <- args$k - 1

  # Obtain set of all interactions that occur at least once in the time series
  l_all_indicators <- list()
  l_unique_indicators <- vector('list', length = d)
  l_unique_indicators <- lapply(l_unique_indicators, function(x) vector('list', length=no_estpoints))
  for(i in 1:no_estpoints)  {
    l_all_indicators[[i]] <- tvmgmobj$tvmodels[[i]]$rawfactor$indicator
    for(ii in 1:d) l_unique_indicators[[ii]][[i]] <- l_all_indicators[[i]][[ii]]
  }

  # get unique set of interactions for each order
  l_unique_indicators <- lapply(l_unique_indicators, function(x) {
    x_comb <- do.call(rbind, x)
    x_comb[!duplicated(x_comb),]
  })


  # Call time-varying DrawFG function:
  l_FG <- list()
  for(i in 1:no_estpoints)  {
    l_FG[[i]] <- DrawFGtv(list_ind = tvmgmobj$tvmodels[[i]]$rawfactor$indicator,
                        list_weights = tvmgmobj$tvmodels[[i]]$rawfactor$weightsAgg,
                        list_signs = tvmgmobj$tvmodels[[i]]$rawfactor$signs,
                        list_ind_all = l_unique_indicators,
                        p = p)
  }


  # -------------------- Process Output -------------------

  # Output:
  # pairwise and factorgraph: 3d arrays
  # rawfactor and intercepts: lists

  col_factorgraph <- ncol(l_FG[[1]]$weightedgraph) # all the same, so i = 1 does the job

  # Create Storage
  a_pairwise_wadj <- a_pairwise_signs <- a_pairwise_edgecolor <- array(dim = c(p, p, no_estpoints))
  a_factorgraph_graph <- a_factorgraph_signs <- a_factorgraph_edgecolor <-array(dim = c(col_factorgraph, col_factorgraph, no_estpoints))
  l_factorgraph_order <- l_factorgraph_nodetype <- l_rawfactor_indicator <- l_rawfactor_weights <- l_rawfactor_signs <- l_intercepts <- list()

  # Loop over timepoints and restructure
    for(i in 1:no_estpoints) {

      # Pairwise
      a_pairwise_wadj[,,i] <- tvmgmobj$tvmodels[[i]]$pairwise$wadj
      a_pairwise_signs[,,i] <- tvmgmobj$tvmodels[[i]]$pairwise$signs
      a_pairwise_edgecolor[,,i] <- tvmgmobj$tvmodels[[i]]$pairwise$edgecolor

      # Factorgraph
      a_factorgraph_graph[,,i] <- l_FG[[i]]$weightedgraph # Fill in from extended Factor graph (see above, DrawFGtv())
      a_factorgraph_signs[,,i] <- l_FG[[i]]$signs
      a_factorgraph_edgecolor[,,i] <- l_FG[[i]]$signcolor
      l_factorgraph_order[[i]] <- l_FG[[i]]$order
      l_factorgraph_nodetype[[i]] <- l_FG[[i]]$nodetype

      # Factorgraph Raw
      l_rawfactor_indicator[[i]] <- tvmgmobj$tvmodels[[i]]$rawfactor$indicator
      l_rawfactor_weights[[i]] <- tvmgmobj$tvmodels[[i]]$rawfactor$weightsAgg
      l_rawfactor_signs[[i]] <- tvmgmobj$tvmodels[[i]]$rawfactor$signs

      # Intercepts
      l_intercepts[[i]] <- tvmgmobj$tvmodels[[i]]$intercepts

    }


  # Fill into output list

  tvmgmobj$pairwise$wadj <- a_pairwise_wadj
  tvmgmobj$pairwise$signs <- a_pairwise_signs
  tvmgmobj$pairwise$edgecolor <- a_pairwise_edgecolor

  tvmgmobj$factorgraph$graph <- a_factorgraph_graph
  tvmgmobj$factorgraph$signs <- a_factorgraph_signs
  tvmgmobj$factorgraph$edgecolor <- a_factorgraph_edgecolor
  tvmgmobj$factorgraph$order <- l_factorgraph_order
  tvmgmobj$factorgraph$nodetype <- l_factorgraph_nodetype

  tvmgmobj$rawfactor$indicator <- l_rawfactor_indicator
  tvmgmobj$rawfactor$weightsAgg <- l_rawfactor_weights
  tvmgmobj$rawfactor$signs <- l_rawfactor_signs

  tvmgmobj$intercepts <- l_intercepts


  # -------------------- Output -------------------

  # Save Node Models and extracted raw factors?
  if(!args$saveModels) {
    tvmgmobj$nodemodels <- NULL
    tvmgmobj$rawfactor <- NULL
  }

  # intercepts
  # tvmgmobj$intercepts <- list_Intercepts

  # Switch warings back on
  if(!args$warnings) {
    options(warn = oldw)
  }

  if(args$signInfo) cat('Note that the sign of parameter estimates is stored separately; see ?tvmgm / ?mgm')

  # Assign class
  class(tvmgmobj) <- c('tvmgm')

  return(tvmgmobj)


} # eoF

