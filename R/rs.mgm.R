

rs.mgm <- function(..., B = 100, N = NULL, replace = TRUE, pbar=TRUE) {
  
  
  # --------- Checking & Preparing Input ----------
  
  arguments <- list(...) # get mgmfit arguments from (...) construct
  
  # # Input Option A: provide fit object
  # c_a <- class(arguments[[1]])
  # if('mgm' %in% c_a) {
  #   arguments <- arguments[[1]]$call
  #   arguments$VAR <- FALSE
  # }
  # if('var' %in% c_a) {
  #   arguments <- arguments[[1]]$call
  #   arguments$VAR <- TRUE
  # }
  
  # Input Option B: provide mgmfit/var.mgm arguments seperately
  
  # Check whether minimum required input is provided
  if(is.null(arguments$data)) stop('Specify the data, see ?mgmfit / ?var.mgm')
  if(is.null(arguments$type)) stop('Specify the type of the variables , see ?mgmfit / ?var.mgm')
  if(is.null(arguments$lev)) stop('Specify the levels of the variables, see ?mgmfit / ?var.mgm')
  if(is.null(arguments$VAR)) stop('Specify whether a MGM or a mixed VAR models hould be estimated, see ?mgmfit / ?var.mgm')
  
  # For other arguments check whether specified, if not assign default values
  if(is.null(arguments$lambda.sel)) arguments$lambda.sel <- 'EBIC'
  if(is.null(arguments$folds)) arguments$folds <- 10
  if(is.null(arguments$gam)) arguments$gam <- .25
  if(is.null(arguments$dval)) arguments$dval <- 2
  if(is.null(arguments$rule.reg)) arguments$rule.reg <- 'AND'
  if(is.null(arguments$method)) arguments$method <- "glm"
  if(is.null(arguments$missings)) arguments$missings <- 'error'
  if(is.null(arguments$weights)) arguments$weights <- NA
  if(is.null(arguments$ret.warn)) arguments$ret.warn <- FALSE
  if(is.null(arguments$binary.sign)) arguments$binary.sign <- FALSE
  
  # Override some arguments
  arguments$pbar <- FALSE

  # get some info from the data
  n <- nrow(arguments$data)
  p <- ncol(arguments$data)
  
  # If no N is specified, use N = nrow(data)
  if(is.null(N)) N <- nrow(arguments$data)
  
  # logical check subsampling
  if(replace==FALSE & N > n) stop('For subsampling (without replacement) one cannot sample more than n observations.')
  
  
  # --------- Do the Bootstrap ----------
  
  b_list <- list() # storage for bootstrap models
  ind_list <- list() # storage for indices defining bootstrap samples
  a_wadj <- array(dim = c(p, p, B))  # separate storage for weighted adjacency matrix (to make it easier to access for the user)
  
  # Progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=B, initial=0, char="-", style = 3)
  
  for(b in 1:B) {
    
    # Different subsampling for MGM and VAR (for Var -1 case, for lag1 VAR)
    if(arguments$VAR) {
      ind_list[[b]] <- sample(1:(n-1), size = N, replace=replace)
    } else {
      ind_list[[b]] <- sample(1:n, size = N, replace=replace) 
    }

    b_list[[b]] <- mgmfit_core(data = arguments$data,
                         type = arguments$type,
                         lev = arguments$lev,
                         lambda.sel = arguments$lambda.sel,
                         folds = arguments$folds,
                         gam = arguments$gam,
                         d = arguments$dval,
                         rule.reg = arguments$rule.reg,
                         pbar = arguments$pbar,
                         method = arguments$method,
                         missings = arguments$missings,
                         weights = arguments$weights,
                         ret.warn = arguments$ret.warn,
                         binary.sign = arguments$binary.sign,
                         VAR = arguments$VAR,
                         rs_indicator = ind_list[[b]])

    a_wadj[ , , b] <- b_list[[b]]$wadj * b_list[[b]]$signs_recov # add edge signs
    a_wadj[ , , b][is.na(a_wadj[ , , b])] <- 0
    
    if(pbar==TRUE) { setTxtProgressBar(pb, b) }
  }
  
  
  # --------- Preprocess output ----------
  
  # --- Preprocessing for MGM ---
  if(!arguments$VAR) {
    
    # Get Column names, if not present, assign them
    names <- colnames(b_list$call$data)
    if(is.null(names)) names <- paste0('Var ', 1:p)
    
    # Create Storage
    l_edgeNames <- list()
    l_edgeIDs <- list()
    l_edgeWeights <- list()
    l_nonzero <- list()
    
    # list all edges
    counter <- 1
    for(i in 1:p) {
      for(j in i:p) {
        if(j!=i) {
          l_edgeNames[[counter]] <- paste(names[i], '-', names[j]) # character vector listing all edges
          l_edgeIDs[[counter]] <- c(i,j) # 2-column matrix listing all edges
          l_edgeWeights[[counter]] <- ew <- a_wadj[i, j , ] # number(edges)-list with B edge-weights estimated on B resamples
          l_nonzero[[counter]] <- sum(ew) > 0 # indicator whether any of the B edge-weights are nonzero (used below)
          counter <- counter + 1
        }
      }
    }

  # --- Preprocessing for VAR ---
  } else {
    
    # Get Column names, if not present, assign them
    names <- colnames(b_list$call$data)
    if(is.null(names)) names <- paste0('Var ', 1:p)
    
    # Create Storage
    l_edgeNames <- list()
    l_edgeIDs <- list()
    l_edgeWeights <- list()
    l_nonzero <- list()
    
    # list all edges
    counter <- 1
    for(i in 1:p) {
      for(j in 1:p) {
          l_edgeNames[[counter]] <- paste(names[i], '->', names[j]) # character vector listing all edges
          l_edgeIDs[[counter]] <- c(i,j) # 2-column matrix listing all edges
          l_edgeWeights[[counter]] <- ew <- a_wadj[i, j , ] # number(edges)-list with B edge-weights estimated on B resamples
          l_nonzero[[counter]] <- sum(ew) > 0 # indicator whether any of the B edge-weights are nonzero (used below)
          counter <- counter + 1
      }
    }
    
  }
  

  
  # --------- Return Output list ----------
  
  outlist <- list('models' = b_list, 
                  'B_wadj' = a_wadj,
                  'B_samples' = ind_list,
                  'call'= arguments,
                  'call_rs'= list('B'=B, 'N'=N, 'replace'=replace),
                  'edgeNames' = l_edgeNames,
                  'edgeIDs' = l_edgeIDs, 
                  'edgeWeights' = l_edgeWeights, 
                  'edgeNonZero'= l_nonzero)
  
  if(arguments$VAR) { class(outlist) <- c('rs', 'var') } else { class(outlist) <- c('rs', 'mgm') }
  
  return(outlist)
  
  
} # end of Function

















