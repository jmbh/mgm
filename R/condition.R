# jonashaslbeck@gmail.com; July 2019

# input:
# - any mgm object with k <= 3
# - a list of variables and values at which to fix those variables

# output:
# - the conditional model (possibly still including 3-way interactions)

condition <- function(object, 
                      values) {
  
  # ---------- Basic Info -----------
  
  type <- object$call$type
  level <- object$call$level
  p <- length(type) # Number of Variables
  nCond <- length(values)
  
  # values: turn list into matrix
  m_fixed_values <- matrix(NA, nrow=nCond, ncol=2)
  m_fixed_values[, 1] <- as.numeric(names(values))
  m_fixed_values[, 2] <- unlist(values)
  
  
  # ---------- Create Output Object-----------
  
  object_new <- object
  object_new$call$condition <- values
  
  
  # ---------- Input Checks -----------
  
  if(object$call$k>3) stop("This function is only implemented for first-order moderation (3-way interactions).")
  if(! ("core" %in% class(object)) ) stop("condition() is currently only implemented for mgm() objects.")
  
  # Categorical variables: only condition on categories that exist
  for(cat in 1:nCond) {
    if(type[m_fixed_values[cat, 1]] == "c") {
      if(!(m_fixed_values[cat, 2] %in% object$call$unique_cats[[m_fixed_values[cat, 1]]])) stop("Fixed category does not exist in the data.")
    }  
  }
  
  # Continuous variables: give warning if one conditions outside 99% quantiles
  
  
  
  # ---------- Loop over response variables -----------
  
  
  for(i in 1:p) {
    
    # ----- Case I) Gaussian response -----
    
    if(type[i] == "g") {
      
      # Access node model
      model_i <- object$nodemodels[[i]]$model
      
      
      # Apply tau-thresholding & AND rule
      model_i <- applyTauAND(i = i,
                             object = object, 
                             model_i = model_i)
      
      
      # Condition / fix values
      model_i_new <- condition_core(i = i,
                                    model_i = model_i, 
                                    m_fixed_values = m_fixed_values)
      
      # Overwrite model object    
      object_new$nodemodels[[i]]$model <- model_i_new
      
      
    } # end if: response gaussian?
    
    
    
    # ----- Case II: Categorical response -----
    
    if(type[i] == "c") {
      
      # Retrieve nodemodel i
      model_i <- object$nodemodels[[i]]$model
      n_resp <- length(model_i)
      
      # Loop over response categories
      for(cat in 1:n_resp) {
        
        model_i_cat <- model_i[[cat]]
        
        # Apply tau-thresholding & AND rule
        model_i_cat <- applyTauAND(i = i,
                                   object = object, 
                                   model_i = model_i_cat)
        
        # Condition / fix values
        model_i_new <- condition_core(i = i,
                                      model_i = model_i_cat, 
                                      m_fixed_values = m_fixed_values)
        
        # Overwrite model object    
        object_new$nodemodels[[i]]$model[[cat]] <- model_i_new
        
      } # end for: response cats
      
    } # end if: response categorical?
    
  } # end for: response variables
  
  
  
  
  # ---------- Aggregation across regressions -----------
  
  object_new2 <- Reg2Graph(object_new)
  
  
  # ---------- Prepare output & return -----------
  
  return(object_new2)
  
  
} # eoF