
# For now only for pairwise mgms and tvmgms

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
  
  int <- sort(int) # sort, because interactions are saved sorted 1:p in "object$rawfactor$weights[[n_order-1]][[int_row]]"
  levelNames <- list()
  levelNames[[1]] <- object$call$levelNames[[int[1]]]
  levelNames[[2]] <- object$call$levelNames[[int[2]]]
  
  # ----- Give summary for interaction -----
  
  # ----- a) mgm-objects -----
  
  if(class_obj == "core") {
    
    n_order <- length(int)
    
    # Input checks
    if(n_order > object$call$k) stop("The order of requested interactions has to match the order of interactions in the model.")
    
    # Create outlist
    outlist <- list()
    outlist$variables <- int
    outlist$type <- type <- object$call$type[int]
    outlist$level <- level <- object$call$level[int]
    outlist$parameters <- list()
    
    # Collect parameter estimates

    int_row <- which(apply(object$rawfactor$indicator[[n_order-1]], 1, function(x) all(x %in% int))) # get row of interaction in "int" in interaction list
    
    if(length(int_row) == 0) stop("The specified interaction has been estimated to be absent.")

    for(i in 1:n_order) {

    # Get parameters of regression on i
    int_i  <- object$rawfactor$weights[[n_order-1]][[int_row]][[i]]
    
    # Create empty array with correct dimensions
    par_mat <- matrix(NA, level[i], level[-i])
    
    # Fill array
    if(type[i] == "c") { # if response = categorical
      
      if(type[-i] == "c") { # if predictor = categorical
        for(i_resp in 1:length(int_i)) par_mat[i_resp, -1] <- int_i[[i_resp]]
      } else {
        par_mat[, 1] <- int_i
      }

      } else { # if response = continuous

        if(type[-i] == "c") { # if predictor = categorical
          par_mat[1, -1] <- int_i
        } else {
          par_mat[1, 1] <- int_i
        }
      }
    
    
    # browser()
    # Set dimension names in array
    if(type[i] == "c") row.names(par_mat) <- paste0(int[i], ".", levelNames[[i]]) else row.names(par_mat) <- int[i]
    if(type[-i] == "c") colnames(par_mat) <- paste0(int[-i], ".", levelNames[[abs(i-3)]]) else colnames(par_mat) <- int[-i]
    
    # Save array in output list
    outlist$parameters[[paste0("Predict_", int[i])]] <- par_mat
      
    } # end for: n_order
    
    
  } # end if: core?
  

  
  # ----- b) tvmgm-objects -----
  
  # ----- c) mvar-objects -----
  
  # ----- d) tvmvar-objects -----
  
  
  
  
  
  # ----- Prepare output -----
  
  
  return(outlist)
  
} # eoF