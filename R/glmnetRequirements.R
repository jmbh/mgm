

glmnetRequirements <- function(data,
                               type,
                               weights, 
                               bootstrap = FALSE,
                               b = NULL, 
                               seed_b = NULL, 
                               silent = FALSE) {
  
  # if silent = FALSE; the function returns error messages
  # if silent = TRUE: the function returns a logical with TRUE=at least one error; 
  #   this is used in resampling schemes to discard bootstrap samples that do not allow to fit the model
  
  var_names <- colnames(data)
  n <- nrow(data)
  
  # 1) Nonzero variance
  var_check <- apply(data, 2, var)
  ind_zero_var <- which(var_check == 0)
  
  check_var_1 <- length(ind_zero_var) > 0
  
 if(!silent)  if(check_var_1) {
    if(bootstrap) cat(paste0("In boostrap sample ", b, " with seed ", seed_b, " the following error occured:\n"))
    stop(paste0('Please only provide variables with nonzero variance. Variable(s) with zero variance: ', paste(var_names[ind_zero_var], collapse = ', ')))
  } 
  
  
  # 2) > 1 events per category
  
  check_var_2 <- FALSE # for the case of: (1) no categorical variables and (2) silent = TRUE (when called in resample() )
  
  if('c' %in% type) {
    ind_cat <- which(type == 'c')
    l_frqu <- list()
    for(i in 1:length(ind_cat)) l_frqu[[i]] <- table(data[,ind_cat[i]]) # this does not catch the case where one category is not present at all; but this is catched by comparing specified levels and real levels
    v_check <- unlist(lapply(l_frqu, function(x) {
      frq_norm <- x / sum(x)
      ind_min <- which.min(frq_norm)
      check1 <- !(x[ind_min] > 1)
    }))
    # Error Msg Check 1:
    
    check_var_2 <- sum(v_check) > 0
    
    if(!silent) if(check_var_2) {
      ind_check1 <- ind_cat[v_check == TRUE]
      stop(paste0('At least 2 events required for each category. Requirement not met for variable(s): ',paste(var_names[ind_cat[ind_check1]], collapse = ', ')))
    }
  }
  
  
  # 3) For each category: p(K=l) > 10^-5
  
  check_var_3 <- FALSE # for the case of: (1) no categorical variables and (2) silent = TRUE (when called in resample() )
  
  if('c' %in% type) {
    
    # Function to compute weighted table:
    wtable <- function(x, weights) {
      n_level <- length(unique(x))
      v_level <- unique(x)
      n_obs <- length(x)
      v_wfrq <- rep(NA, n_level)
      for(i in 1:n_level) v_wfrq[i] <- sum(rep(1, sum(x == v_level[i])) * weights[x == v_level[i]]) / n
      return(v_wfrq)
    }
    
    
    ind_cat <- which(type == 'c')
    check2 <- rep(NA, length(ind_cat))
    for(i in 1:length(ind_cat)) check2[i] <- min(wtable(data[,ind_cat[i]], weights)) < 10^-5 # smaller = error
    
    check_var_3 <- sum(check2)>0
    
    if(!silent) if(check_var_3) stop(paste0('Each category has to have probability > 10^-5. Requirement not met for variable(s): ',paste(var_names[ind_cat[check2]], collapse = ', ')))
    
  }
  

  # If any of the three checks fails, return TRUE
  if(silent) return(any(check_var_1, check_var_2, check_var_3))
  
  
} # eoF





