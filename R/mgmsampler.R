

mgmsampler <- function(factors,
                       interactions,
                       thresholds,
                       sds,
                       type,
                       level,
                       N,
                       nIter = 250,
                       pbar = TRUE,
                       divWarning = 10^3,
                       returnChains = FALSE)
  
  
{
  
  
  # ----- Compute basic aux vars -----
  
  p <- length(type)
  n_order <- length(lapply(interactions, function(x) 1))
  
  
  # ---------- Input Checks ----------
  
  # Check whether level are correctly specified according to type
  if(any(level[type!="c"] != 1)) stop("The levels of all non-categorical variable has to be specified as c(1).")
  
  # Check whether all interactions are specified in the correct dimension
  for(ord in 1:n_order) {
    
    rows_ord <- nrow(factors[[ord]])
    
    # Only check if there are actually factors of a given dimension
    if(!is.null(rows_ord)) {
      for(row in 1:rows_ord) {
        
        if(!all(level[factors[[ord]][row, ]] == dim(interactions[[ord]][[row]]))) {
          stop(paste0("The dimensions of the interaction ",
                      paste0(factors[[ord]][row, ], collapse = "-"), 
                      " is specified incorrectly. Please correct the 'interactions' argument."
          ))
        }
        
      } # end for: rows (interactions)      
    } # end if: 0> rows?
    
  } # end for: ord
  
  # browser()
  
  # Are all parameters specified as matrices?
  for(i in 1:n_order) {
    n_ints <- length(interactions[[i]])
    if(n_ints>0) for(row in 1:n_ints) if(!(class( interactions[[i]][[row]]) %in% c("array","matrix"))) stop("The parameters of each interaction have to be provided as k-dimensional array.")
  }
  
  
  # ----- Input Checks -----
  
  # 1) Do matrices in factors match the corresponding order?
  for(i in 1:n_order) if(!is.null(factors[[i]])) if(ncol(factors[[i]]) != i+1) stop(paste0('Matrix specifying ',i+1,'-way interctions has incorrect dimensions.' ))
  
  # 2) Are parameter of right dimension provided for all interactions
  for(i in 1:n_order) {
    
    n_ints <- nrow(factors[[i]])
    if(!is.null(n_ints)) {
      
      for(row in 1:n_ints) {
        
        check <- all.equal(level[factors[[i]][row,]],
                           dim(interactions[[i]][[row]]))
        
        if(any(check == FALSE)) stop(paste0('Incorrect dimensions of parameter array for interaction: ', paste(level[factors[[i]][row,]], collapse = ' ')))
      }
    }
  }
  
  # 3) Are interactions defined twice?
  for(i in 1:n_order) {
    
    if(!is.null(factors[[i]])) {
      
      symflag <- FlagSymmetric(factors[[i]])
      dup <- duplicated(symflag)
      if(any(dup)) stop(paste0('Interaction ', paste(factors[[i]][min(which(dup)),], collapse = ' '), ' is specified twice.' ))
      
    }
  }
  
  
  # 4) Are variances for conditional Gaussians specified properly?
  if(!('g' %in% type)) {
    if(missing(sds)) sds <- NULL
  } else {
    if(any(is.na(sds[type=='g']))) stop('Missing values in standard deviations')
    if(any(sds[type=='g'] < 0)) stop('Specified standard deviations need to be positive')
    if(any(!is.finite(sds[type=='g']))) stop('Specified standard deviations need to be finite')
  }
  
  # 5) Other basics
  if(length(type) != length(level)) stop('Length of type does not match length of level.')
  if(length(thresholds) != length(type)) stop('Length of threshold list does notmatch length of type vector')
  if(!all.equal(unlist(lapply(thresholds, length)), level)) stop('Specified thresholds do not match the number of categories in level')
  
  
  # Transforms all entries in factors and interactions into matrices to avoid trouble with nrow()
  for(i in 1:n_order) {
    
    if(!is.null(factors[[i]])) {
      factors[[i]] <- as.matrix(factors[[i]], ncol = ncol(factors[[i]]), nrow = nrow(factors[[i]]))
      
    }
  }
  
  
  # Create Output Object
  mgmsamp_obj <- list('call' = NULL,
                      'data' = NULL, 
                      'chains' = NULL)
  
  # Copy the call
  mgmsamp_obj$call <- list('factors' = factors,
                           'interactions' = interactions,
                           'thresholds' = thresholds,
                           'sds' = sds,
                           'type' = type,
                           'N' = N,
                           'level' = level,
                           'nIter' = nIter,
                           'pbar' = pbar,
                           'divWarning' = divWarning)
  
  
  # ---------- Optional: Enforce positive definite Gaussian subgraph ----------
  
  # browser()
  
  # ---------- Call C++ Gibbs Sampler (for now: R version) ----------
  
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=N, initial=0, char="-", style = 3)
  
  # Return chains of Gibbs sampler
  if(returnChains) chains <- array(NA, dim = c(N, nIter, p))
  
  data <- matrix(NA, nrow = N, ncol = p)
  
  for(case in 1:N) {
    
    sampling_seq <- matrix(NA, nrow = nIter, ncol = p)
    
    
    # Random Start
    for(v in 1:p) {
      if(type[v] == 'g') sampling_seq[1, v] <- rnorm(1)
      if(type[v] == 'p') sampling_seq[1, v] <- rpois(1, lambda = .5)
      if(type[v] == 'c') sampling_seq[1, v] <- sample(1:level[v], size = 1)
    }
    
    # Loop over iterations of the Gibbs sampler
    for(iter in 2:nIter) {
      
      # Looping over variables
      for(v in 1:p) {
        
        # A) ----- Continuous conditional -----
        if(type[v] != 'c') {
          
          l_natPar <- list() # list collecting the interaction terms separate for each order
          l_natPar[[1]] <- thresholds[[v]] # first entry = order 'zero', the intercept
          
          # Loop over order of interactions
          for(ord in 1:n_order) {
            
            # Loop over factors for fixed ord
            n_rows <- nrow(factors[[ord]])
            l_row_terms <- list() # collects terms for interactions of given order
            row_counter <- 1
            
            if(!is.null(n_rows)) { # only loop over factors if there are factors for given ord
              for(row in 1:n_rows) {
                
                f_ind <- factors[[ord]][row, ]
                
                # is the current variable in that row?
                if(v %in% f_ind) {
                  
                  # find out: which is current var, where is it, fill in 1 for continuous
                  # (hack; this will be slightly more complicated for categorical response)
                  fill_in <- rep(NA, ord + 1)
                  get_cont <- rep(NA, ord + 1) # get values of continuous variables
                  k_counter <- 1
                  
                  for(k in f_ind) { # loop over variables in factor
                    
                    # Gibbs algorithm: take always the fresh data
                    if(k >= v) gibbs_shift <- iter - 1 else gibbs_shift <- iter # equal case doesn't matter as excluded below
                    
                    if(level[k]==1) fill_in[k_counter] <- 1 else fill_in[k_counter] <- sampling_seq[gibbs_shift, k]
                    if(level[k]==1) get_cont[k_counter] <- sampling_seq[gibbs_shift, k] else get_cont[k_counter] <- 1
                    k_counter <- k_counter + 1
                    
                  }
                  
                  # Convergence Check
                  if(any(abs(sampling_seq[!is.na(sampling_seq)]) > divWarning)) warning('Gibbs Sampler diverged; adapt parameterization of continuous variables.')
                  
                  
                  # prepare fill in matrix:
                  m_fill_in <- matrix(fill_in, ncol=length(fill_in))
                  
                  # get the parameter for ord-interaction out the the ord-order parameter array
                  row_int <- interactions[[ord]][[row]][m_fill_in]
                  
                  # multiply it by continuous variables if in interaction and save in list
                  l_row_terms[[row_counter]] <- prod(get_cont[-which(f_ind == v)]) * row_int
                  
                  # update row counter (those rows (factors) in which v is contained)
                  row_counter <- row_counter + 1
                  
                } # end if: variable in given row?
                
              } # end for: rows in factors of given ord
              
            } # end if: isnul n_rows
            
            # fill in all terms of given order ord
            l_natPar[[ord + 1]] <- unlist(l_row_terms)
            
          } # end loop: ord
          
          
          # Compute Natural Parameter (collapse and sum all interaction terms)
          natPar <- sum(unlist(l_natPar))
          
          # Take a draw
          if(type[v] == 'g') sampling_seq[iter, v] <- rnorm(1, mean = natPar, sd = sds[v])
          
          if(type[v] == 'p') {
            if(natPar <= 0) stop(paste0('Lambda parameter of Poisson variable ', v, ' is nonpositive. Sampler returns NA.'))
            sampling_seq[iter, v] <- rpois(1, lambda = natPar)
          }
          
          
        } #end if: v == continuous?
        
        
        # A) ------ Categorical conditional ------
        
        if(type[v] == 'c') {
          
          v_Potential <- rep(NA, level[v]) # Storage to collect potentials
          
          for(cat in 1:level[v]) {
            
            l_natPar <- list() # list collecting the interaction terms separate for each order
            l_natPar[[1]] <- thresholds[[v]][cat] # first entry = order 'zero', the intercept for category cat
            
            # Loop over order of interactions
            for(ord in 1:n_order) {
              
              # Loop over factors for fixed ord
              n_rows <- nrow(factors[[ord]])
              l_row_terms <- list() # collects terms for interactions of given order
              row_counter <- 1
              
              
              if(!is.null(n_rows)) { # only loop over factors if there are factors for given ord
                for(row in 1:n_rows) {
                  
                  f_ind <- factors[[ord]][row, ]
                  
                  # is the current (v) variable in that row? If not the interaction doesn't matter for (v)
                  if(v %in% f_ind) {
                    
                    # find out: which is current var, where is it, fill in 1 for continuous
                    # (hack; this will be slightly more complicated for categorical response)
                    fill_in <- rep(NA, ord + 1) # here we collect the value in the previous iteration of the respective variables in f_ind 
                    get_cont <- rep(NA, ord + 1) # get values of continuous variables
                    k_counter <- 1
                    
                    for(k in f_ind) { # loop over variables in factor
                      
                      # Gibbs algorithm: take always the fresh data
                      if(k >= v) gibbs_shift <- iter - 1 else gibbs_shift <- iter # equal case doesn't matter as excluded below
                      
                      # level[k]==1 indicates a continuous varibale
                      if(level[k]==1) fill_in[k_counter] <- 1 else fill_in[k_counter] <- sampling_seq[gibbs_shift, k]
                      if(level[k]==1) get_cont[k_counter] <- sampling_seq[gibbs_shift, k] else get_cont[k_counter] <- 1 # numeric value for cat is always one for the present category (indicator function)
                      k_counter <- k_counter + 1
                      
                    }
                    
                    # Convergence Check
                    if(any(abs(sampling_seq[!is.na(sampling_seq)]) > divWarning)) warning('Gibbs Sampler diverged; adapt parameterization of continuous varibales.')
                    
                    # fill in current category
                    fill_in[which(v == f_ind)] <- cat # we set the category of the variable at hand to the cat we loop over; we fix the category of v to cat in 1:level, since we compute all level potentials
                    
                    # prepare fill in matrix:
                    m_fill_in <- matrix(fill_in, ncol=length(fill_in))
                    
                    # get the parameter for ord-interaction out the the ord-order parameter array
                    row_int <- interactions[[ord]][[row]][m_fill_in]
                    
                    # multiply it by continuous variables if in interaction and save in list
                    l_row_terms[[row_counter]] <- prod(get_cont[-which(f_ind == v)]) * row_int
                    
                    # update row counter (those rows (factors) in which v is contained)
                    row_counter <- row_counter + 1
                    
                  } # end if: variable in given row?
                  
                } # end for: rows in factors of given ord
                
              } # end if: isnul n_rows
              
              l_natPar[[ord + 1]] <- unlist(l_row_terms)
              
            } # end loop: ord
            
            v_Potential[cat] <- sum(unlist(l_natPar))
            
          } # end: for cat
          
          
          v_probabilities <- exp(v_Potential) / sum(exp(v_Potential)) # calc probability for each category
          sampling_seq[iter, v] <- sample(1:level[v], size = 1, prob = v_probabilities) # sample from multinomial distribution
          
        }
        
      } # end for: v
      
    } # end for: iter
    
    
    if(returnChains) chains[case, , ] <- sampling_seq
    
    data[case, ] <- sampling_seq[nIter, ] # Fill in last sampling iteration
    
    if(pbar==TRUE) setTxtProgressBar(pb, case)
    
  } # end for: case
  
  
  
  # ---------- Export ----------
  
  mgmsamp_obj$data <- data
  if(returnChains) mgmsamp_obj$chains <- chains
  
  return(mgmsamp_obj)
  
  
} # eoF
