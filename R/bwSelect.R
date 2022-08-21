

bwSelect <- function(data,
                     type,
                     level,
                     bwSeq,
                     bwFolds, # number of training rounds (or folds in CV)
                     bwFoldsize, # should default to n_var / bwFolds (then proper stratified CV, but very slow for big datasets)
                     modeltype,
                     pbar,
                     ... # arguments for mvar/mgm, passed through by tvmar
)


{
  
  # -------------------- Input Checks -------------------

  args <- list(...)
  
  if(is.null(args$mgm_par)) args$mgm_par <- TRUE
  if(is.null(args$tvmgm_par)) args$tvmgm_par <- FALSE
  
  # Input checks specific for mvar()
  if(modeltype == 'mvar') {
    if(is.null(args$lags)) stop('No lags specified')
    if(any(args$lags<1)) stop('Specified lags have to be in {1, 2, ..., n -1}')
    if(any(round(args$lags)!=args$lags)) stop('Specified lags have to be positive integers')
  }

  if(any(bwSeq < 0)) stop('All bandwidth values have to be positive.')

  if(missing(bwFolds)) stop('No number of folds specified.')
  if(bwFolds<0) stop('Number of folds has to be a positive integer')
  if(bwFolds!=round(bwFolds)) stop('Number of folds has to be a positive integer')

  if(missing(bwFoldsize)) stop('No fold size specified.')
  if(bwFoldsize<0) stop('Fold size has to be a positive integer')
  if(bwFolds!=round(bwFolds)) stop('Fold size has to be a positive integer')

  # Input checks specific for mvar() / mgm() are passed down to these functions


  # ----- Fill in Defaults -----
  
  # browser()
  
  if(missing(bwFolds)) bwFolds <- NULL
  if(missing(bwFoldsize)) bwFoldsize <- NULL
  if(missing(pbar)) pbar <- TRUE
  if(is.null(args$saveData)) args$saveData <- FALSE
  if(is.null(args$k)) args$k <- k <- 2
  if(is.null(args$overparameterize)) args$overparameterize <- FALSE


  # ----- Compute Aux Variables -----

  n <- nrow(data)
  if(modeltype == 'mvar') n_var <- n - max(args$lags) # this is how many rows there are after transforming the data
  p <- ncol(data)

  n_bw <- length(bwSeq) # length of bandwidth path

  # ----- Checks -----

  if(missing(bwSeq)) stop('No bandwidth sequence specified')
  if(any(bwSeq <= 0)) stop('All entries in the bandwidth sequence have to be strictly positive')


  # ----- Create Output Structure -----

  bwSelectObj <- list('call' = NULL,
                      'bwModels' = NULL,
                      'fullErrorFolds' = NULL,
                      'fullError' = NULL,
                      'meanError' = NULL,
                      'testsets' = NULL,
                      'zeroweights' = NULL)


  # ----- Copy The Call -----

  bwSelectObj$call <- list('data' = NULL,
                           'type' = type,
                           'level' = level,
                           'bwSeq' = bwSeq,
                           'bwFolds' = bwFolds,
                           'bwFoldsize' = bwFoldsize,
                           'modeltype' = modeltype,
                           'pbar' = pbar)


  if(args$saveData) bwSelectObj$call$data <- data



  # ----- Progress Bar Business -----

  # if(is.null(args$pbar)) args$pbar <- TRUE
  # Set up Progress bar
  if(pbar==TRUE) pb_bw <- txtProgressBar(min = 0, max = n_bw, initial = 0, char="-", style = 3)

  # -------------------- Estimate Path -------------------

  l_bw_models <- vector('list', length = n_bw) # Storage
  l_performance <- list()


  # !!!! IF EBIC REACTIVATED: needs argument 'nEstpoints' !!!!

  # if(bwSel == 'EBIC') {
  #
  #   warning('The EBIC does NOT select the optimal bandwidth. It is only implemented for research purposes. Use cross validation instead.')
  #
  #   estpoints <- seq(0, n_var, length = nEstpoints)
  #
  #   # Set up Storage
  #   l_performance <- vector('list', length = n_bw)
  #   dummy_estpoints <- vector('list', length = nEstpoints)
  #   dummy_p <- vector('list', length = p)
  #   dummy_estpoints <- lapply(dummy_estpoints, function(x) dummy_p)
  #   l_performance <- lapply(l_performance, function(x) dummy_estpoints)
  #
  #
  #   for(i in 1:n_bw) {
  #
  #     l_bw_models[[i]] <- tvmvar(data = data,
  #                                type = type,
  #                                level = level,
  #                                estpoints = estpoints,
  #                                bandwidth = bwSeq[i],
  #                                pbar,
  #                                ...)
  #
  #
  #     # Get the p x estpoints x n_bw EBICs out of the tvmvar objects
  #     for(v in 1:p) {
  #       for(j in 1:nEstpoints) {
  #         l_performance[[i]][[j]][[v]] <- l_bw_models[[i]]$tvmodels[[j]]$nodeModels[[v]]$EBIC
  #       }
  #     }
  #
  #
  #     # Update Progress Bar
  #     if(args$pbar==TRUE) setTxtProgressBar(pb_bw, i)
  #
  #   }
  #
  #   # ----- Evaluate -----
  #
  #   # Collapse across nodes
  #   l_perf_time <- lapply(l_performance, function(x) {
  #     lapply(x, function(y) mean(unlist(y)))
  #   })
  #
  #   # Collapse across nodes & time points
  #   l_perf_mean <- lapply(l_perf_time, function(x) mean(unlist(x)))
  #
  #
  # }

  # -------------------- Estimate Path: VAR model -------------------


  if(modeltype == 'mvar') {

    # Define Testset in each fold:  we take every kth object in some distance and move the sequence by 1 in each fold
    l_testsets <- list()
    l_zero_weights <- list()
    for(fold in 1:bwFolds) {
      l_testsets[[fold]] <- round(seq(fold, n_var - bwFolds + fold, length = bwFoldsize)) # nice!
      l_zero_weights[[fold]] <- rep(1, n_var)
      l_zero_weights[[fold]] [l_testsets[[fold]]] <- 0
    }


    # Storage
    l_performance <- list()
    l_bw_models <- list()
    
    for(i in 1:n_bw) {

      # Cross validation scheme
      l_foldModels <- list()
      l_foldPerform <- list()

      for(fold in 1:bwFolds) {

        # Fit model on training dataset fold
        l_foldModels[[fold]] <- tvmvar(data = data, # full data, test cases are deleted via setting weights to zero
                                       type = type,
                                       level = level,
                                       estpoints = l_testsets[[fold]] / n_var, # since we use [0,1] interval since mgm 1.2-3
                                       bandwidth = bwSeq[i],
                                       pbar = FALSE, # switch off progress bar in tvmvar
                                       zero_weights = l_zero_weights[[fold]], # vector indicating zero weights for test cases
                                       saveModels = TRUE, # otherwise we can't make any predictions
                                       signInfo = FALSE,
                                       ...)
    
      

        # Make Predictions at test-locations
        l_foldPerform[[fold]] <- bwSelPredict(data = data,
                                              type = type,
                                              level = level,
                                              obj = l_foldModels[[fold]],
                                              test = l_testsets[[fold]],
                                              lags = args$lags,
                                              modeltype = modeltype,
                                              overparameterize = args$overparameterize,
                                              consec = l_foldModels[[fold]]$call$consec,
                                              k = k,
                                              ...)

      }

      l_bw_models[[i]] <- l_foldModels

      # Compute mean performance per Fold on maximal level: time x nodes
      # Fill in Performance
      l_performance[[i]] <- l_foldPerform


      # Update Progress Bar
      if(pbar==TRUE) setTxtProgressBar(pb_bw, i)

    }

  }


  # -------------------- Estimate Path: MGM -------------------


  if(modeltype == 'mgm') {

    # Define Testset in each fold:  we take every kth object in some distance and move the sequence by 1 in each fold
    l_testsets <- list()
    l_zero_weights <- list()
    for(fold in 1:bwFolds) {
      l_testsets[[fold]] <- round(seq(fold, n - bwFolds + fold, length = bwFoldsize)) # now normalized
      l_zero_weights[[fold]] <- rep(1, n)
      l_zero_weights[[fold]] [l_testsets[[fold]]] <- 0
    }

    # Storage
    l_performance <- list()
    l_bw_models <- list()

    for(i in 1:n_bw) {

      # Cross validation scheme
      l_foldModels <- list()
      l_foldPerform <- list()

      for(fold in 1:bwFolds) {

        # Fit model on training dataset fold
        l_foldModels[[fold]] <- tvmgm(data = data, # full data, test cases are deleted via setting weights to zero
                                    type = type,
                                    level = level,
                                    estpoints = l_testsets[[fold]] / n, # since we use [0,1] interval since mgm 1.2-3
                                    bandwidth = bwSeq[i],
                                    pbar = FALSE, # switch off progress bar in tvmvar
                                    zero_weights = l_zero_weights[[fold]], # vector indicating zero weights for test cases
                                    saveModels = TRUE, # otherwise we can't make predictions
                                    signInfo = FALSE,
                                    mgm_par = args$mgm_par,
                                    tvmgm_par = args$tvmgm_par,
                                    ...)

        

        # Make Predictions at test-locations
        l_foldPerform[[fold]] <- bwSelPredict(data = data,
                                              type = type,
                                              level = level,
                                              obj = l_foldModels[[fold]],
                                              test = l_testsets[[fold]],
                                              lags = args$lags,
                                              modeltype = modeltype,
                                              overparameterize = args$overparameterize,
                                              k = args$k,
                                              ...)

      }

      l_bw_models[[i]] <- l_foldModels

      # Compute mean performance per Fold on maximal level: time x nodes
      # Fill in Performance
      l_performance[[i]] <- l_foldPerform

      # Update Progress Bar
      if(pbar==TRUE) setTxtProgressBar(pb_bw, i)

    }

  }



  # ----- Evaluate: Take mean across folds -----

  # Collapse across folds
  l_perf_mean <- lapply(l_performance, function(x) {
    mean(unlist(lapply(x, function(y) y$error_mean)))
  })

  l_perf_all <- lapply(l_performance, function(x) {
    dum <- lapply(x, function(y) y$errors)
    n_dum <- length(dum)
    Reduce('+', dum) / n_dum
  })





  # -------------------- Output -------------------

  # Estimated Models
  if(is.null(args$saveModels)) {
    bwSelectObj$bwModels <- NULL
  } else if(args$saveModels) {
    bwSelectObj$bwModels <- l_bw_models
  } else {
    bwSelectObj$bwModels <- NULL
  }

  # Testsets
  bwSelectObj$testsets <- l_testsets
  bwSelectObj$zeroweights <- l_zero_weights

  # Error Output
  bwSelectObj$fullErrorFolds <- l_performance
  bwSelectObj$fullError <- l_perf_all
  bwSelectObj$meanError <- unlist(l_perf_mean)
  names(bwSelectObj$meanError) <- bwSeq


  class(bwSelectObj) <- c(modeltype, 'bwSelect')

  return(bwSelectObj)


} # eoF








