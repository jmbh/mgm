# jonashaslbeck@gmail.com; July 2019

# Input:
# - fitted mgm model object (mgm() output)
# - node i at hand
# Output:
# - Nodemodels with tau-threshold and AND rule applied

# This is used in condition() (for now) and in predict() (planned later)


applyTauAND <- function(i, 
                        object, 
                        model_i) {
  
  
  # ----- Apply: tau-threshold -----
  
  model_i[abs(model_i) < object$nodemodels[[i]]$tau] <- 0
  
  
  # ----- Apply AND-rule ------
  
  # if(i==12) browser()
  
  if(object$call$ruleReg == "AND") {
    
    # print(i)
    
    ## This creates an object with the same dimensionality as the glmnet output matrix
    ## Below we fill into this matrix whether the corresponding parameter has been set to zero by the AND rule in the fitted model
    ## We take that information from the indicator matrix
    
    # Aux matrix
    aux_m_AND <- matrix(NA, nrow(model_i)-1, 5) # columns: parameter value, variable 1, variable 2, order of interaction, parameter value after AND rule
    aux_m_AND[, 1] <- model_i[-1, 1] # copy parameters, without intercept
    
    # get variables out:
    names_i <- rownames(model_i)[-1]
    names_aux1 <- strsplit(names_i, ":")
    names_aux2 <- lapply(names_aux1, function(x) gsub("V", "", x))
    names_aux3 <- lapply(names_aux2, function(x) {
      out <- strsplit(x, "\\.")
      if(length(out) == 1) c(out[[1]][1], NA) else c(out[[1]][1], out[[2]][1])
    })
    names_aux4 <- do.call(rbind, names_aux3)
    names_aux4 <- apply(names_aux4, 2, as.numeric)
    
    aux_m_AND[, 2:3] <- names_aux4
    aux_m_AND[, 4] <- apply(aux_m_AND[, 2:3], 1, function(x) sum(!is.na(x))) # order of interaction
    n_2way <- sum(aux_m_AND[, 4] == 1)
    n_3way <- sum(aux_m_AND[, 4] == 2)
    
    
    # -- Loop over indicator matrix --
    
    # 2-way interactions
    ind_2way <- object$interactions$indicator[[1]]
    ind_2way <- matrix(ind_2way, ncol=2) # for the case of: only 1 e-way interaction
    out <- apply(ind_2way, 1, function(x) {
      if(i %in% x) {
        TRUE
      } else { 
        FALSE
      }
    })
    
    if(sum(out)>0) { # all of this needed only if at least one pairwise interaction is estimated to be nonzero
      
      ind_2way_i <- matrix(ind_2way[out, ], ncol=2) # subset nonzero estimated 2-way interactions that contain variable i in indicator matrix
      preds_2way <- aux_m_AND[aux_m_AND[, 4] == 1, 2] # subset 2-way interactions in model object
      
      for(v in 1:n_2way) {
        out_v <- apply(ind_2way_i, 1, function(x) preds_2way[v] %in% x)
        if(sum(out_v) > 0) aux_m_AND[v, 5] <- 1 else aux_m_AND[v, 5] <- 0
      }
      
    } else {
      
      aux_m_AND[1:n_2way, 5] <- 0
      
    }
    
    # 3-way interactions
    
    # only execute if there are any 3-way interactions _specified_ involving variable i
    # (it can happen that no 3-way interaction is specified if the moderators are specified by matrix input)
    if(n_3way > 0) {
      
      ind_3way <- object$interactions$indicator[[2]]
      ind_3way <- matrix(ind_3way, ncol=3) # for the case of: only one 3-way interaction
      out <- apply(ind_3way, 1, function(x) {
        if(i %in% x) {
          TRUE
        } else { 
          FALSE
        }
      })
      
      if(sum(out)>0) { # all of this needed only if at least one 3-way interaction is estimated to be nonzero
        
        ind_3way_i <- matrix(ind_3way[out, ], ncol=3)
        preds_3way <- matrix(aux_m_AND[aux_m_AND[, 4] == 2, 2:3], ncol=2)
        
        for(v in 1:n_3way) {
          out_v <- apply(ind_3way_i, 1, function(x) preds_3way[v, 1] %in% x & preds_3way[v, 2] %in% x)
          if(sum(out_v) > 0) aux_m_AND[v+n_2way, 5] <- 1 else aux_m_AND[v+n_2way, 5] <- 0
        }
        
      } else {
        
        aux_m_AND[(n_2way+1):(n_2way+n_3way), 5] <- 0
        
      }
      
    }
    
    
    # Finally, threshold
    model_i[-1, 1][aux_m_AND[, 5] == 0] <- 0
    
    
  } # end if: ruleReg = "AND" ?
  
  return(model_i)
  
} # eoF
