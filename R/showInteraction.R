

showInteraction <- function(object, 
                            int) {
  
  
  # ----- Input Checks -----
  
  class_obj <- class(object)[2]
  
  # Check whether object is mgm model object
  if(!(class_obj %in% c("core", "tvmgm", "mvar", "tvmvar") )) stop("Please provide an mgm model object (mgm, tvmgm, mvar, tvmvar) as input.")
  
  # Check whether interactions are specified as integers
  if(!all(int == round(int))) stop("Interactions have to be specified as sets of integers in 1:p")
  if(any(int < 1)) stop("Interactions have to be specified as sets of integers in 1:p")
  
  # Current limitation of function:
  if(class_obj %in% c("core", "tvmgm")) if(object$call$k > 2) stop("showInteraction() currently only supports pairwise interactions.")
  if(class_obj %in% c("tvmgm", "mvar", "tvmvar")) stop("showInteraction() currently only supports mgm() objects.")
  
  # ----- Compute Aux variables -----
  
  int <- sort(int) # sort, because interactions are saved sorted 1:p in "object$interactions$weights[[n_order-1]][[int_row]]"
  levelNames <- list()
  levelNames[[1]] <- object$call$levelNames[[int[1]]]
  levelNames[[2]] <- object$call$levelNames[[int[2]]]
  

  # ------------------------------------------------------------------------------------------------
  # ---------- a) mgm-objects ----------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------
  
  if(class_obj == "core") {
    
    n_order <- length(int)
    type_int <- object$call$type[int]
    level_int <- object$call$level[int]
    
    
    # Input checks
    # if(n_order > object$call$k) stop("The order of requested interactions has to match the order of interactions in the model.")
    
    # Create outlist
    outlist <- list("variables" = int, 
                    "order" = n_order, 
                    "type" = type_int, 
                    "level" = level_int, 
                    "edgeweight" = NULL, 
                    "sign" = NULL,
                    "parameters" = list())
    
  
    # Get row of specified interaction
    int_row <- which(apply(matrix(object$interactions$indicator[[n_order-1]], ncol=n_order), 1, function(x) all(x %in% int))) # get row of interaction in "int" in interaction list
    
    outlist$edgeweight <- object$interactions$weightsAgg[[n_order-1]][[int_row]]
    outlist$sign <- object$interactions$sign[[n_order-1]][int_row]
    
    if(n_order > 2) {
      
      outlist$parameters = NULL
      
    } else {
      
      # ------- Collect & label parameter estimates -------
      
      if(length(int_row) == 0) stop("The specified interaction has been estimated to be absent.")
      
      for(i in 1:n_order) {
        
        # Get parameters of regression on i
        int_i  <- object$interactions$weights[[n_order-1]][[int_row]][[i]]
        
        # Create empty array with correct dimensions
        par_mat <- matrix(NA, level_int[i], level_int[-i])
        
        # Fill array
        if(type_int[i] == "c") { # if response = categorical
          
          if(type_int[-i] == "c") { # if predictor = categorical (only works for k=2 order MGMs ..)
            for(i_resp in 1:length(int_i)) {
              if(object$call$overparameterize) par_mat[i_resp, ] <- int_i[[i_resp]] else par_mat[i_resp, -1] <- int_i[[i_resp]]
            } 
          } else {
            par_mat[, 1] <- int_i
          }
          
        } else { # if response = continuous
          
          if(type_int[-i] == "c") { # if predictor = categorical
            par_mat[1, -1] <- int_i
          } else {
            par_mat[1, 1] <- int_i
          }
        }
        
        # browser()
        # Set dimension names in array
        if(type_int[i] == "c") row.names(par_mat) <- paste0(int[i], ".", levelNames[[i]]) else row.names(par_mat) <- int[i]
        if(type_int[-i] == "c") colnames(par_mat) <- paste0(int[-i], ".", levelNames[[abs(i-3)]]) else colnames(par_mat) <- int[-i]
        
        # Save array in output list
        outlist$parameters[[paste0("Predict_", int[i])]] <- par_mat
        
      } # end for: n_order
      
    } # if: order == 2
    
    
  } # end if: core?
  
  
  # ------------------------------------------------------------------------------------------------
  # ---------- b) tvmgm-objects ----------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------
  

  # ------------------------------------------------------------------------------------------------
  # ---------- c) mvar-objects ----------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------

  
  # ------------------------------------------------------------------------------------------------
  # ---------- d) tvmvar-objects ----------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------


  # ----- Return -----
  
  class(outlist) <- "int"
  
  return(outlist)
  
} # eoF