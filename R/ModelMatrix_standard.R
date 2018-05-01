# jonashaslbeck@gmail.com; May 2018

ModelMatrix_standard <- function(data, 
                                 d, # maximum neighborhood size in model = k-1
                                 v, # node on which current nodewise regression is performed
                                 moderators # moderators, if specified
                                 
)
  
  
{
  
  # -------- Input checks --------
  
  # ...
  
  # -------- Case I: No Moderation --------
  
  p <- ncol(data)
  
  if(is.null(moderators)) {
    
    if(d > (p - 1)) {
      stop("Order of interactions cannot be larger than the number of predictors.")
    } else if (d == 1){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
    } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^",d)) }
    
    # Construct standard design matrix (to get number of parameters for tau threshold)
    X_standard <- model.matrix(form, data = data)[, -1] # delete intercept (added by glmnet later)
    
  } else {
    
    
    # -------- Case II: Moderation --------
    
    # Simple predictors
    pred_simple <- paste(colnames(data)[-v], collapse = " + ")
    
    # Moderation terms
    l_mods <- list()
    n_mods <- length(moderators)
    n_l_mods <- 0
    
    
    if(v %in% moderators) {
      other_comb_pairw <- t(combn((1:p)[-v], 2))
      l_mods[[1]] <- paste0("V", other_comb_pairw[, 1], ".",  " * ", "V", other_comb_pairw[, 2], ".", collapse = " + ")
    } else {      
      for(i in 1:n_mods) {  # loop over moderators
        l_mods[[i]] <- paste0(colnames(data)[-c(v, moderators[i])], "*", "V", moderators[i], ".", collapse = " + ")
        n_l_mods <- n_l_mods + 1
      }
    } # end if: v = moderator?
    
    if(n_l_mods > 1) for(i in 1:(n_l_mods - 1)) l_mods[[i]] <- paste0(l_mods[[i]], " + ") # add plus between terms but not in the end
    v_mods <- do.call(paste0, l_mods)
    
    # Stitch together model formula
    if(length(v_mods) == 0) { # if there is no moderation for that variable (the case if v is the only moderator)
      form <- paste(colnames(data)[v], "~",
                    pred_simple)
      
    } else {
      form <- paste(colnames(data)[v], "~",
                    pred_simple,
                    " + ",
                    v_mods)
    }
    
    form <- as.formula(form)
    
    # Create design matrix
    X_standard <- model.matrix(form, 
                               data = data)[, -1] # delete intercept (added by glmnet later)
    
  } # end if: moderation?
  
  
  # -------- Output --------
  
  return(X_standard)
  
  
} # end of function



