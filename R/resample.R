
resample <- function(object, # one of the four mgm model objects (mgm, mvar, tvmgm, tvmvar) 
                     data, # the data
                     nB,
                     seeds,
                     blocks, # type of resampling for time-series data
                     quantiles,
                     pbar, # progress bar
                     verbatim,
                     ...)
  
{
  
  
  # ----- Fill in defaults -----
  
  if(missing(pbar)) pbar <- TRUE
  if(missing(verbatim)) verbatim <- FALSE
  if(missing(seeds)) seeds <- 1:nB
  if(missing(blocks)) blocks <- 10
  if(missing(quantiles)) quantiles <- c(0.05, .95)
  
  # --- Define aux variables ---
  
  n <- nrow(data)
  p <- length(object$call$type)
  o_call <- object$call
  
  
  # --- Create output object ---
  
  call <- list('object' = object,
               'data' = data,
               'nB' = nB,
               'blocks' = blocks,
               'pbar' = pbar)
  
  outlist <- list('call' = call,
                  'bootParameters' = NULL,
                  'bootQuantiles' = NULL,
                  'models' = NULL,
                  "Times" = rep(NA, nB),
                  "totalTime" = NULL)
  
  
  # ----- General input checks -----
  
  if(missing(nB)) stop('No number of bootstrap samples provided.')
  if(!is.logical(pbar)) stop('The argument pbar can only be TRUE/FALSE.')    
  if(any(quantiles < 0) | any(quantiles > 1)) stop("Specified quantiles have to be in [0,1].")
  
  if(!is.null(object$call$k)) if(object$call$k > 2) warning("Summary statistics in the output are only implemented for pairwise mgms/tvmgms.")
  
  if(!any(c("core", "mvar", "tvmvar", "tvmgm") %in% class(object))) stop("Please provide an mgm object as input.")
  
  
  # -------------------- Resampling for mgm-objects -------------------------------------------------------------------------------- 
  
  if("core" %in% class(object)) {
    
    # --- Define bootstrap samples ---
    
    l_ind <- list()
    for(b in 1:nB) {
      
      set.seed(seeds[b])
      
      l_ind[[b]] <- sample(1:n, 
                           size = n, 
                           replace = TRUE)
      
    }
    
    
    # --- Estimate mgms on bootstrap samples ---
    
    l_b_models <- list() # Storage
    
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nB, initial=0, char="-", style = 3) # initialize progress bar
    
    for(b in 1:nB) {
      
      set.seed(seeds[b]) # for cross validation
      if(verbatim) print(paste0("Seed = ", seeds[b]))
      tt <- proc.time()[3]
      
      l_b_models[[b]] <- mgm(data = data[l_ind[[b]], ],
                             type = o_call$type,
                             level = o_call$level,
                             lambdaSeq = o_call$lambdaSeq,
                             lambdaSel = o_call$lambdaSel,
                             lambdaFolds = o_call$lambdaFolds,
                             lambdaGam = o_call$lambdaGam,
                             alphaSeq = o_call$alphaSeq,
                             alphaSel = o_call$alphaSel,
                             alphaFolds = o_call$alphaFolds,
                             alphaGam = o_call$alphaGam,
                             k = o_call$k,
                             moderators = o_call$moderators,
                             ruleReg = o_call$ruleReg,
                             weights = o_call$weights[l_ind[[b]]], # just copying, no other vector describing rows
                             threshold = o_call$threshold,
                             method = o_call$method,
                             binarySign = o_call$binarySign,
                             scale = o_call$scale,
                             verbatim = o_call$verbatim,
                             pbar = FALSE,
                             warnings = o_call$warnings,
                             saveModels = o_call$saveModels,
                             saveData = o_call$saveData,
                             overparameterize = o_call$overparameterize,
                             signInfo = FALSE) # to avoid msg for each model
      
      outlist$Times[b] <- proc.time()[3] - tt
      if(pbar==TRUE) setTxtProgressBar(pb, b)
      
    }
    
    outlist$models <- l_b_models
    
    
  } #end: if mgm
  
  
  # -------------------- Resampling for tvmgm-objects --------------------------------------------------------------------------------
  
  if("tvmgm" %in% class(object)) {
    
    # --- Define bootstrap samples ---
    
    # Break data into blocks of equal time-duration (unsing timepoints vector)
    Qt <- quantile(object$call$timepoints, probs = seq(0, 1, length = blocks + 1))
    ind_blocks <- cut(x = object$call$timepoints, 
                      breaks = Qt,
                      labels = FALSE)
    ind_blocks[1] <- 1 # for some reason the first entry is NA
    
    
    # Storage
    l_ind <- list() 
    
    for(b in 1:nB) {
      
      set.seed(seeds[b])
      
      l_ind_blocks <- list()
      
      
      for(b2 in 1:blocks) {
        
        # Take bootstrap samples within blocks  
        block_b2 <- which(ind_blocks == b2)
        l_ind_blocks[[b2]] <- sample(x = block_b2, 
                                     size = length(block_b2), 
                                     replace = TRUE)
        
      }
      
      l_ind[[b]] <- unlist(l_ind_blocks)
      
    } # end: loop bootstraps
    
    
    # --- Estimate mgms on bootstrap samples ---
    
    l_b_models <- list() # Storage
    
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nB, initial=0, char="-", style = 3) # initialize progress bar
    
    for(b in 1:nB) {
      
      set.seed(seeds[b])
      if(verbatim) print(paste0("Seed = ", seeds[b]))
      tt <- proc.time()[3]
      
      l_b_models[[b]] <- tvmgm(data = data[l_ind[[b]],],
                               type = o_call$type,
                               level = o_call$level,
                               timepoints = o_call$timepoints[l_ind[[b]]], # just copying
                               estpoints = o_call$estpoints,
                               bandwidth = o_call$bandwidth,
                               # mgm arguments passed via (...):
                               lambdaSeq = o_call$lambdaSeq,
                               lambdaSel = o_call$lambdaSel,
                               lambdaFolds = o_call$lambdaFolds,
                               lambdaGam = o_call$lambdaGam,
                               alphaSeq = o_call$alphaSeq,
                               alphaSel = o_call$alphaSel,
                               alphaFolds = o_call$alphaFolds,
                               alphaGam = o_call$alphaGam,
                               k = o_call$k,
                               ruleReg = o_call$ruleReg,
                               weights = o_call$weights[l_ind[[b]]], # just copying
                               threshold = o_call$threshold,
                               method = o_call$method,
                               binarySign = o_call$binarySign,
                               scale = o_call$scale,
                               verbatim = o_call$verbatim,
                               pbar = FALSE, # we use pbar over bootstrap samples instead
                               warnings = o_call$warnings,
                               saveModels = o_call$saveModels,
                               saveData = o_call$saveData,
                               overparameterize = o_call$overparameterize,
                               signInfo = FALSE)
      
      outlist$Times[b] <- proc.time()[3] - tt
      if(pbar==TRUE) setTxtProgressBar(pb, b)
      
      
    } # end: loop bootstraps
    
    outlist$models <- l_b_models
    
    
  } #end if: tvmgm?
  
  
  
  
  # -------------------- Resampling for mvar-objects --------------------------------------------------------------------------------
  
  if("mvar" %in% class(object)) {
    
    # --- Define bootstrap samples ---
    
    # Compute design matrix to find out how many rows it has
    if(!is.null(o_call$beepvar) & !is.null(o_call$dayvar)) {
      o_call$consec <- beepday2consec(beepvar = o_call$beepvar,
                                      dayvar = o_call$dayvar)
    } # if: specification of consecutiveness via beepvar and dayvar
    
    data_lagged <- lagData(data = data,
                           lags = o_call$lags,
                           consec = o_call$consec)
    
    n_design <- sum(data_lagged$included) 
    ind_valid_rows <-  (1:n)[data_lagged$included]
    
    # Take bootstrap sample from rows in design matrix
    l_ind <- list()
    
    for(b in 1:nB) {
      
      set.seed(seeds[b])
      
      l_ind[[b]] <- sample(x = ind_valid_rows, 
                           size = n_design, 
                           replace = TRUE)
      
    }
    
    
    # --- Estimate mvar on bootstrap samples ---
    
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nB, initial=0, char="-", style = 3) # initialize progress bar
    
    l_b_models <- list()
    
    for(b in 1:nB) {
      
      set.seed(seeds[b])
      if(verbatim) print(paste0("Seed = ", seeds[b]))
      tt <- proc.time()[3]
      
      l_b_models[[b]] <- mvar(data = data, # not changed here because changed via boot_ind inside of mvar()
                              type = o_call$type,
                              level = o_call$level,
                              lambdaSeq = o_call$lambdaSeq,
                              lambdaSel = o_call$lambdaSel,
                              lambdaFolds = o_call$lambdaFolds,
                              lambdaGam = o_call$lambdaGam,
                              alphaSeq = o_call$alphaSeq,
                              alphaSel = o_call$alphaSel,
                              alphaFolds = o_call$alphaFolds,
                              alphaGam = o_call$alphaGam,
                              lags = o_call$lags,
                              consec = o_call$consec,
                              # dayvar = o_call$dayvar, 
                              # beepvar = o_call$beepvar, # already comes through consec
                              weights = o_call$weights,
                              threshold = o_call$threshold,
                              method = o_call$method,
                              binarySign = o_call$binarySign,
                              scale = o_call$scale,
                              verbatim = o_call$verbatim,
                              pbar = FALSE,
                              warnings = o_call$warnings,
                              saveModels = o_call$saveModels,
                              saveData = o_call$saveData,
                              overparameterize = o_call$overparameterize,
                              signInfo = FALSE, # to avoid msg for each model
                              # extra arguments to do bootstrap passed via (...)
                              bootstrap = TRUE,
                              boot_ind = l_ind[[b]]) 
      
      outlist$Times[b] <- proc.time()[3] - tt
      if(pbar==TRUE) setTxtProgressBar(pb, b)
      
    } # end loop: bootstraps
    
    outlist$models <- l_b_models
    
  } # end: if: mvar?
  
  
  
  
  # -------------------- Resampling for tvmvar-objects --------------------------------------------------------------------------------
  
  if("tvmvar" %in% class(object)) {
    
    # --- Define bootstrap samples ---
    
    # Compute design matrix to find out how many rows it has
    if(!is.null(o_call$beepvar) & !is.null(o_call$dayvar)) {
      o_call$consec <- beepday2consec(beepvar = o_call$beepvar,
                                      dayvar = o_call$dayvar)
    } # if: specification of consecutiveness via beepvar and dayvar; otherwise: already specified
    
    data_lagged <- lagData(data = data, 
                           lags = o_call$lags, 
                           consec = o_call$consec)
    
    # Make use of time-vector to be able to have block-bootstrap on actual time scale  
    timepoints_design <- o_call$timepoints[data_lagged$included]

    # Break data into blocks of equal time-duration (unsing timepoints vector)
    Qt <- quantile(timepoints_design, probs = seq(0, 1, length = blocks + 1))
    ind_blocks <- cut(x = timepoints_design,  # important: indices in the design matrix
                      breaks = Qt,
                      labels = FALSE)
    ind_blocks[1] <- 1 # for some reason the first entry is NA
    
    # Rows included in the design matrix
    ind_valid_rows <- (1:nrow(data))[data_lagged$included]
    
    # Storage
    l_ind <- list() 
    
    for(b in 1:nB) {
      
      set.seed(seeds[b])
      
      l_ind_blocks <- list()
      
      for(b2 in 1:blocks) {
        
        # Take bootstrap samples within blocks defined above from rows that are included (ind_valid_rows)
        block_b2 <- ind_valid_rows[which(ind_blocks == b2)] 
        l_ind_blocks[[b2]] <- sample(x = block_b2, 
                                     size = length(block_b2), 
                                     replace = TRUE)
        
      }
      
      l_ind[[b]] <- unlist(l_ind_blocks)
      
    } # end: loop bootstraps
    

    # --- Estimate tvmvar on bootstrap samples ---
    
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nB, initial=0, char="-", style = 3) # initialize progress bar
    
    l_b_models <- list()
    
    for(b in 1:nB) {
      
      # browser()
      
      set.seed(seeds[b])
      if(verbatim) print(paste0("Seed = ", seeds[b]))
      tt <- proc.time()[3]
      
      l_b_models[[b]] <- tvmvar(data = data, # not changed here because changed via boot_ind inside of mvar() / tvmvar()
                                type = o_call$type, 
                                level = o_call$level, 
                                timepoints = o_call$timepoints, 
                                estpoints = o_call$estpoints,
                                bandwidth = o_call$bandwidth,
                                lambdaSeq = o_call$lambdaSeq,
                                lambdaSel = o_call$lambdaSel,
                                lambdaFolds = o_call$lambdaFolds,
                                lambdaGam = o_call$lambdaGam,
                                alphaSeq = o_call$alphaSeq,
                                alphaSel = o_call$alphaSel,
                                alphaFolds = o_call$alphaFolds,
                                alphaGam = o_call$alphaGam,
                                lags = o_call$lags,
                                consec = o_call$consec,
                                # dayvar = o_call$dayvar, 
                                # beepvar = o_call$beepvar,  # already comes through consec
                                weights = o_call$weights,
                                threshold = o_call$threshold,
                                method = o_call$method,
                                binarySign = o_call$binarySign,
                                scale = o_call$scale,
                                verbatim = o_call$verbatim,
                                pbar = FALSE,
                                warnings = o_call$warnings,
                                saveModels = o_call$saveModels,
                                saveData = o_call$saveData,
                                overparameterize = o_call$overparameterize,
                                signInfo = FALSE,
                                # extra arguments for bootstrap
                                bootstrap = TRUE,
                                boot_ind = l_ind[[b]])
      
      outlist$Times[b] <- proc.time()[3] - tt
      if(pbar==TRUE) setTxtProgressBar(pb, b)
      
    } # end loop: bootstraps
    
    outlist$models <- l_b_models
    
    
    
  } # end if: tvmvar?
  
  
  # ----------------------------------------------------------------------------------------
  # ----- Summarize parameters in arrays ---------------------------------------------------
  # ----------------------------------------------------------------------------------------
  
  # Compute Aux variables
  
  model_obj <- object
  obj_class <- class(model_obj)[2]
  
  # -------------------- 1) mvar --------------------
  
  if(obj_class == 'mvar') {
    
    p <- length(model_obj$call$type)
    nlags <- length(model_obj$call$lags)
    nquantiles <- length(quantiles)
    
    # Collect all estimates in one array
    collect_array <- collect_array_sign <- array(NA, dim = c(p, p, nlags, nB))
    for(b in 1:nB) collect_array[, , , b] <- outlist$models[[b]]$wadj
    for(b in 1:nB) collect_array_sign[, , , b] <- outlist$models[[b]]$signs
    
    # add sign
    collect_array_wS <- collect_array
    ind_negative <- which(collect_array_sign == -1, arr.ind = TRUE)
    collect_array_wS[ind_negative] <- collect_array_wS[ind_negative] * -1
    
    # Compute quantiles
    quantile_array <- apply(collect_array_wS, 1:3, function(x) quantile(x, probs = quantiles))
    quantile_array_res <- array(dim = c(p, p, nlags, nquantiles))
    for(qu in 1:nquantiles) quantile_array_res[, , , qu] <- quantile_array[qu, , , ]
    
    # # Compute median
    # median_array <- apply(collect_array, 1:3, function(x) median(x))
    
    outlist$bootParameters <- collect_array_wS
    outlist$bootQuantiles <- quantile_array_res
    
  }
  
  
  # -------------------- 2) tvmvar --------------------
  
  if(obj_class == 'tvmvar') {
    
    p <- length(model_obj$call$type)
    nlags <- length(model_obj$call$lags)
    n_estpoints <- length(model_obj$call$estpoints)
    nquantiles <- length(quantiles)
    
    # Collect all estimates in one array
    collect_array <- collect_array_sign <- array(NA, dim = c(p, p, nlags, n_estpoints, nB))
    for(b in 1:nB) collect_array[, , , , b] <- outlist$models[[b]]$wadj
    for(b in 1:nB) collect_array_sign[, , , , b] <- outlist$models[[b]]$signs 
    
    # add sign
    collect_array_wS <- collect_array
    ind_negative <- which(collect_array_sign == -1, arr.ind = TRUE)
    collect_array_wS[ind_negative] <- collect_array_wS[ind_negative] * -1
    
    # Compute quantiles
    quantile_array <- apply(collect_array, 1:4, function(x) quantile(x, probs = quantiles))
    quantile_array_res <- array(dim = c(p, p, nlags, n_estpoints, nquantiles))
    for(qu in 1:nquantiles) quantile_array_res[, , , , qu] <- quantile_array[qu, , , , ]
    
    # # Compute median
    # median_array <- apply(collect_array, 1:3, function(x) median(x))
    
    outlist$bootParameters <- collect_array_wS
    outlist$bootQuantiles <- quantile_array_res
    
  }
  
  
  # -------------------- 3) pairwise mgms --------------------
  
  if(obj_class == 'core') {
    
    p <- length(model_obj$call$type)
    nquantiles <- length(quantiles)
    
    
    ## Collect all estimates
    collect_array <- collect_array_sign <- array(NA, dim = c(p, p, nB))
    for(b in 1:nB) collect_array[, , b] <- outlist$models[[b]]$pairwise$wadj
    for(b in 1:nB) collect_array_sign[, , b] <- outlist$models[[b]]$pairwise$signs
    
    # add sign
    collect_array_wS <- collect_array
    ind_negative <- which(collect_array_sign == -1, arr.ind = TRUE)
    collect_array_wS[ind_negative] <- collect_array_wS[ind_negative] * -1
    
    # Compute quantiles
    quantile_array <- apply(collect_array_wS, 1:2, function(x) quantile(x, probs = quantiles))
    quantile_array_res <- array(dim = c(p, p, nquantiles))
    for(qu in 1:nquantiles) quantile_array_res[, , qu] <- quantile_array[qu, , ]
    
    outlist$bootParameters <- collect_array_wS
    outlist$bootQuantiles <- quantile_array_res
    
  }
  
  # -------------------- 4) pairwise tvmgms --------------------
  
  if(obj_class == 'tvmgm') {
    
    p <- length(model_obj$call$type)
    nquantiles <- length(quantiles)
    n_estpoints <- length(model_obj$call$estpoints)
    
    ## Collect all estimates
    collect_array <- collect_array_sign <- array(NA, dim = c(p, p, n_estpoints, nB))
    for(b in 1:nB) collect_array[, , , b] <- outlist$models[[b]]$pairwise$wadj
    for(b in 1:nB) collect_array_sign[, , , b] <- outlist$models[[b]]$pairwise$signs
    
    # add sign
    collect_array_wS <- collect_array
    ind_negative <- which(collect_array_sign == -1, arr.ind = TRUE)
    collect_array_wS[ind_negative] <- collect_array_wS[ind_negative] * -1
    
    # Compute quantiles
    quantile_array <- apply(collect_array_wS, 1:3, function(x) quantile(x, probs = quantiles))
    quantile_array_res <- array(dim = c(p, p, n_estpoints, nquantiles))
    for(qu in 1:nquantiles) quantile_array_res[, , , qu] <- quantile_array[qu, , , ]
    
    outlist$bootParameters <- collect_array_wS
    outlist$bootQuantiles <- quantile_array_res
    
  }
  
  
  
  # ----- Return outlist -----
  
  outlist$totalTime <- sum(outlist$Times)
  class(outlist) <- c('resample', class(object)[2])
  
  return(outlist)
  
  
} # eoF






