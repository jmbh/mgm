# jonashaslbeck@gmail.com; May 2018

ModelMatrix_standard <- function(data, 
                                 type,
                                 d, # maximum neighborhood size in model = k-1
                                 v, # node on which current nodewise regression is performed
                                 moderators # moderators, if specified
                                 
)

  
{
  
  # -------- Input checks --------
  
  # ...
  
  # -------- Case I: No Moderation --------
  
  p <- ncol(data)
  for(i in 1:p) if(type[i] == "c") data[, i] <- as.factor(data[, i])
  
  if(is.null(moderators)) {
    
    if(d > (p - 1)) {
      stop("Order of interactions cannot be larger than the number of predictors.")
    } else if (d == 1){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
    } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^", d)) }
    
    # Construct standard design matrix (to get number of parameters for tau threshold)
    X_standard <- model.matrix(form, data = data)[, -1] # delete intercept (added by glmnet later)
    
  } else {
    
    
    # -------- Case II: Moderation --------
    
    # mSpec <- ifelse(class(moderators) %in% c("integer","numeric"), "vector", "matrix")
    mSpec <- ifelse(any(class(moderators) %in% c("integer","numeric")), "vector", "matrix")
    # browser()
    
    
    # Terms for interactions/moderation
    n_mods <- ifelse(mSpec == "vector", length(moderators), nrow(moderators))
    
    # IIa) Moderator specification: Vector
    if(mSpec == "vector") {
      
      n_l_mods <- 0
      l_mods <- list()
      
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
      
    } # end if: vector specification
    
    
    # IIb) Moderator specification: Matrix
    if(mSpec == "matrix") {
      
      nrow_mods <- nrow(moderators)
      ind_v_inMod <- as.logical(apply(moderators, 1, function(x) v %in% x  ))
      
      if(sum(ind_v_inMod)>0) { # if variable v is involved in at least one interaction
      # get interaction terms for node v
      int_terms <- t(apply(matrix(moderators[ind_v_inMod], ncol=3), 1, function(x) x[x!=v]))
      n_terms <- length(int_terms)
      v_mods <- paste0("V", int_terms[, 1], ".",  " * ", "V", int_terms[, 2], ".", collapse = " + ")
      } else {
        v_mods <- NULL
      }
      
      
    } # end if: matrix specification
    
    
    
    # -------- Final Model formula: main effects + interactions --------
    
    # Terms for main effects
    pred_simple <- paste(colnames(data)[-v], collapse = " + ")
    
    # Put together
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



