# jonashaslbeck@gmail.com; July 2019

# Input: 
#  - moderated MGM (max k = 3)
#  - a value for each of the included moderators

# Output:
#  - an pairwise mgm object, conditioned on the provided values of the moderators


FixModerator <- function(object, 
                         values) {
  
  
  
  # ---------- Get Basic Info ---------
  
  p <- length(object$call$type) # Number of Variables
  moderators <- object$call$moderators
  nMods <- length(moderators)
  type <- object$call$type
  level <- object$call$level
  
  # ---------- Input Checks ---------
  
  # Does number of values match number of moderators?
  if(nMods != length(values)) stop("Number of specified values does not match number of moderators.")
  if(length(values)!=nMods) stop("Please specify a value for each moderator in the model.")
  
  if(object$call$k>3) stop("This function is only implemented for first-order moderation (3-way interactions).")
  
  # ADD: only allow stationary models for now
  
  # ADD: not supported for poisson nodes
  
  
  # ---------- Create empty MGM object ---------
  
  object_new <- list("pairwise"=list("wadj"=NULL,
                                     "signs"=NULL,
                                     "edgecolor"=NULL),
                     # "interactions"=list("indicator"=NULL,
                     #                     "weightsAgg"=NULL,
                     #                     "weights"=NULL,
                     #                     "signs"=NULL),
                     "nodemodels"=vector("list", length=p))
  
  
  
  # -------- New approach: follow the nodewise estimation --------
  
  # All variables
  for(i in 1:p) {
    
    # ----- I) Gaussian response -----
    
    if(type[i] == "g") {
      
      # Load nodemodel of moderated MGM
      nodemodel_i <- object$nodemodels[[i]]$model
      names_i <- rownames(nodemodel_i)
      
      # Storage: new model object
      nodemodel_new_i <- matrix(NA, nrow=1+(p-1), ncol=1)
      nodemodel_new_i[1, 1] <- nodemodel_i[1, 1] # copy intercept
      
      # Loop over predictor variables
      for(j in 1:(p-1)) {
        
        j_var <- ((1:p)[-i])[j] # all variables except j
        
        
        # -- Get moderation effects for variable i --
        v_modeff_j_m <- rep(0, nMods)
        
        # loop over moderators
        for(m in 1:nMods) {
          
          m_fix <- moderators[m]
          if(j_var==m_fix) next # we don't have quadratic effects in the model
          
          # get the moderators that contain predictor j and moderator m
          ind_j_m <- grepl(paste0("V", j_var), names_i) & grepl(paste0("V", m_fix), names_i)
          
          # Get fixed moderator value
          mod_value <- values[[which(names(values) == m_fix)]]
          
          if(sum(ind_j_m) > 0) {
            
            # Moderator: Continuous
            if(type[m_fix] == "g") {
              
              v_modeff_j_m[m] <- mod_value * nodemodel_i[ind_j_m, 1] # fixed value of moderator x moderation effect
              
            } else {
              # Moderator: Categorical
              # here we have to match the fixed value to the appropriate parameter
              
              # get all characters after second . (bit complicated, because should also work for >9 categories)
              fpoint <-  sub("^[^.]*", "", names_i[ind_j_m])
              for(q in 1:length(fpoint)) fpoint[q] <- substr(fpoint[q], 2, nchar(fpoint[q]))
              aux_ind <- sub("^[^.]*", "", fpoint)
              for(q in 1:length(aux_ind)) aux_ind[q] <- substr(aux_ind[q], 2, nchar(aux_ind[q]))
              ind_par_mod_val <- as.numeric(aux_ind)
              
              # two cases: reference class (no parameter associated) vs. others
              if(mod_value == 0) {
                v_modeff_j_m[m] <- 0 # nothing happens if we are in the reference category
              } else {
                v_modeff_j_m[m] <- 1 * nodemodel_i[ind_j_m, 1][aux_ind == mod_value]
              }
              
            } # end if: moderator cont/cat
            
          } # end if: any moderation effects on variable i?
          
        } # end for: moderators
        
        # put it all together
        nodemodel_new_i[j+1, 1] <- nodemodel_i[j+1, 1] + sum(v_modeff_j_m)
        
      } # end for: j
      
      object_new$nodemodels[[i]] <- nodemodel_new_i # variable doesn't predict itself
      
    } # end if: gaussian response
    
    
    
    # ----- II) Categorical response -----
    
    if(type[i] == "c") {
      
      browser()
      
      # Load nodemodel of moderated MGM
      nodemodel_i <- object$nodemodels[[i]]$model
      names_i <- rownames(nodemodel_i)
      
      # Storage: new model object
      nodemodel_new_i <- matrix(NA, nrow=1+(p-1), ncol=1)
      nodemodel_new_i[1, 1] <- nodemodel_i[1, 1] # copy intercept
      
      
      
      # Loop over predictor variables
      
      
      # loop over moderators
      
      
      
      
      
    } # end if: categorical response
    
  } # end for: i
  
  
  
  # browser()
  
  
  
  # ---------- Compute aggregated (pairwise) edge weights ---------
  
  # Put all in one adj matrix
  m_wadj <- m_signs <- matrix(NA, p, p)
  for(i in 1:p) {
    m_wadj[i, -i] <- object_new$nodemodels[[i]][-1]
  }
  m_wadj <- (m_wadj + t(m_wadj)) / 2
  diag(m_wadj) <- 0
  
  # Create sign matrix
  m_signs <- sign(m_wadj)
  for(i in 1:p) if(type[i]=="c") m_signs[i, ] <- 0
  m_signs[m_wadj==0] <- NA
  
  # Create edgecolor matrix
  edgecolor <- matrix("black", p, p)
  edgecolor[m_signs==1] <- "darkgreen"
  edgecolor[m_signs==-1] <- "red"
  edgecolor[m_signs==0] <- "grey"
  
  
  # ---------- Prepare Output ---------
  
  object_new$pairwise$wadj <- m_wadj
  object_new$pairwise$signs <- m_signs
  object_new$pairwise$edgecolor <- edgecolor
  
  return(object_new)
  
  
} # eoF






