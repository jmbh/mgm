# jonashaslbeck@gmail.com; July 2019

# Input:
# - the nodemodel of node i
# - matrix of fixed values
# Output: 
# - conditioned nodemodel object


condition_core <- function(i = i,
                           model_i, 
                           m_fixed_values) {
  
  
  names_i <- rownames(model_i)
  n_terms <- nrow(model_i)
  nCond <- nrow(m_fixed_values)
  
  # --- Get a nicer object for main effects / interactions --
  
  effects <- matrix(NA, nrow = n_terms-1, ncol=7)
  colnames(effects) <- c("Variable1", "Variable2", "Fixed1", "Fixed2", "Parameter", "Type1", "Type2")
  
  names_aux1 <- strsplit(names_i[-1], ":")
  names_aux2 <- lapply(names_aux1, function(x) gsub("V", "", x))
  
  names_aux3 <- lapply(names_aux2, function(x) {
    x_out <- rep(NA, length(x))
    for(v in 1:length(x)) if(substr(x[v], start = nchar(x[v]), nchar(x[v])) == ".") x_out[v] <- gsub("\\.", "", x[v]) else x_out[v] <- x[v]
    return(x_out)
  })
  
  n_var_i <- unlist(lapply(names_aux3, length))
  for(q in 1:length(n_var_i)) effects[q, 1:n_var_i[q]] <- as.numeric(unlist(names_aux3[[q]]))
  
  ## Fill in fixed values
  
  # type of predictor (cat vs con)?
  for(q in 1:(n_terms-1)) effects[q, 6] <- ifelse(effects[q, 1] == round(effects[q, 1]), 1, 0) # 1 = continuous, 0 = categorical
  for(q in which(n_var_i==2)) effects[q, 7] <- ifelse(effects[q, 2] == round(effects[q, 1]), 1, 0) # 1 = continuous, 0 = categorical
  
  # Fill in continuous predictors
  for(q in 1:(n_terms-1)) for(f in 1:nCond) if(effects[q, 1] == m_fixed_values[f, 1]) effects[q, 3] <- m_fixed_values[f, 2]
  for(q in which(n_var_i==2)) for(f in 1:nCond) if(effects[q, 2] == m_fixed_values[f, 1]) effects[q, 4] <- m_fixed_values[f, 2]
  
  # Fill in categorical predictors
  for(q in 1:(n_terms-1)) {
    if(effects[q, 6] == 0) {
      var_cat <- strsplit(as.character(effects[q, 1]), "\\.")[[1]]
      for(f in 1:nCond) if(as.numeric(var_cat[1]) == m_fixed_values[f, 1]) if(as.numeric(var_cat[2]) == m_fixed_values[f, 2]) effects[q, 3] <- 1 else effects[q, 3] <- 0  
    }
  }
  for(q in which(n_var_i==2)) {
    if(effects[q, 7] == 0) {
      var_cat <- strsplit(as.character(effects[q, 2]), "\\.")[[1]]
      for(f in 1:nCond) if(as.numeric(var_cat[1]) == m_fixed_values[f, 1]) if(as.numeric(var_cat[2]) == m_fixed_values[f, 2]) effects[q, 4] <- 1 else effects[q, 4] <- 0  
    }
  }
  
  
  # browser()
  
  # Fill in parameter values
  effects[, 5] <- model_i[-1, 1]
  
  # --- Fill (new) conditioned model i ---
  
  l_cPars <- vector("list", length = n_terms)
  l_cPars <- lapply(l_cPars, function(x) list() ) # list structure, since we don't "know" in advance how many terms we'll have
  
  # Only calculate new parameters for regressions on variables that are not fixed
  if(!i %in% m_fixed_values[, 1]) {
    
    # Copy intercept
    l_cPars[[1]][[1]] <- model_i[1, 1]
    
    for(q in 1:(n_terms-1)) {
      
      # main effects
      if(n_var_i[q] == 1) {
        
        if(is.na(effects[q, 3])) {
          # I) no fixed value: just copy main effect
          l_cPars[[q+1]][[length(l_cPars[[q+1]])+1]] <- effects[q, 5]
        } else {
          # II) fixed value: multiply times fixed valye & copy to intercept
          l_cPars[[1]][[length(l_cPars[[1]])+1]] <- effects[q, 5] * effects[q, 3]
        }
        
        
      } # end if: main effects
      
      # interaction effects
      if(n_var_i[q] == 2) {
        
        ind_spec <- sum(c(is.na(effects[q, 3]), is.na(effects[q, 4])))
        
        # I) nothing: just copy interaction effect
        if(ind_spec == 2) 
          l_cPars[[q+1]][[length(l_cPars[[q+1]])+1]] <- effects[q, 5]
        
        # II) one of the two: add to respective main effect
        if(ind_spec==1) {
          ind_specified <- !is.na(c(effects[q, 3], effects[q, 4]))
          
          ind_leftover_mainE <- which(effects[n_var_i==1, 1]==effects[q, 1:2][!ind_specified]) # indicates the row of the main effect to which we add the present moderation effect
          
          l_cPars[[ind_leftover_mainE+1]][[length(l_cPars[[ind_leftover_mainE+1]])+1]] <- effects[q, 5] * effects[q, 3:4][ind_specified]
        } 
        
        # III) both: add to intercept
        if(ind_spec == 0) l_cPars[[1]][[length(l_cPars[[1]])+1]] <- effects[q, 5] * effects[q, 3] * effects[q, 4]
        
      } # end if: interaction effects

      
    } # end for: loop parameters
    
  } # end if: variable still random (not fixed)?
  
  # Collapse lists into new model object
  model_i_new <- matrix(NA, nrow=n_terms, ncol=1)
  rownames(model_i_new) <- names_i
  for(q in 1:n_terms) model_i_new[q, 1] <- sum(unlist(l_cPars[[q]]))
  
  return(model_i_new)
  
} # eoF




