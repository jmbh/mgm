


boot.mgmfit <- function(..., B = 100, pbar=TRUE) {
  
  
  # --------- Checking & Preparing Input ----------
  
  arguments <- list(...) # get mgmfit arguments from (...) construct
  
  # Check whether minimum required input is provided
  if(is.null(arguments$data)) stop('Specify the data as for mgmfit(), see ?mgmfit')
  if(is.null(arguments$type)) stop('Specify the type of the variabales as for mgmfit(), see ?mgmfit')
  if(is.null(arguments$lev)) stop('Specify the levels of the variables as for mgmfit(), see ?mgmfit')

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
  # Override some arguments
  arguments$pbar <- FALSE

  
  # get some info from the data
  n <- nrow(arguments$data)
  p <- ncol(arguments$data)
  
  
  # --------- Do the Bootstrap ----------
  
  b_list <- list() # storage for bootstrap models
  ind_list <- list() # storage for indices defining bootstrap samples
  a_wadj <- array(dim = c(p, p, B))  # separate storage for weighted adjacency matrix (to make it easier to access for the user)
  
  # Progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=B, initial=0, char="-", style = 3)
  
  for(b in 1:B) {
    
    ind_list[[b]] <- sample(1:n, size = n, replace=TRUE)

    b_list[[b]] <- mgmfit(data = arguments$data[ind_list[[b]],],
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
                         ret.warn = arguments$ret.warn)

    a_wadj[ , , b] <- b_list[[b]]$wadj * b_list[[b]]$signs_recov # add edge signs
    a_wadj[ , , b][is.na(a_wadj[ , , b])] <- 0
    
    if(pbar==TRUE) { setTxtProgressBar(pb, b) }
  }
  
  
  # --------- Preprocess output ----------
  
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
        l_edgeNames[[counter]] <- paste(names[i], '-', names[j])
        l_edgeIDs[[counter]] <- c(i,j)
        l_edgeWeights[[counter]] <- ew <- a_wadj[i, j , ]
        l_nonzero[[counter]] <- sum(ew) > 0
        counter <- counter + 1
      }
    }
  }
  
  length(l_edgeIDs)
  length(l_edgeWeights)
  do.call(rbind, l_edgeWeights)
  
  # --------- Return Output list ----------
  
  outlist <- list('models' = b_list, 
                  'B_wadj' = a_wadj,
                  'B_samples' = ind_list,
                  'call'= arguments, 
                  'edgeNames' = l_edgeNames,
                  'edgeIDs' = l_edgeIDs, 
                  'edgeWeights' = l_edgeWeights, 
                  'edgeNonZero'= l_nonzero)
  
  class(outlist) <- c('boot', 'mgm')
  
  return(outlist)
  
  
} # end of Function

















