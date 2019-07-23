# jonashaslbeck@gmail.com; July 2019

# Input: 
# - glmnet output for each response variables (nodemodel object)
# - output object (so far only filled with function call)
# Output: 
# - output object with processed glmnet output


Reg2Graph <- function(mgmobj) {
  
  
  # -------------------- Processing glmnet Output -------------------
  
  # Basic info from model object
  p <- length(mgmobj$call$type) # number of variables
  n <- mgmobj$call$n
  d <- mgmobj$call$k - 1
  moderators <- mgmobj$call$moderators
  type <- mgmobj$call$type
  level <- mgmobj$call$level
  threshold <- mgmobj$call$threshold
  npar_standard <- mgmobj$call$npar
  binarySign <- mgmobj$call$binarySign
  ruleReg <- mgmobj$call$ruleReg
  overparameterize <- mgmobj$call$overparameterize
  ind_cat <- mgmobj$call$ind_cat
  ind_binary <- mgmobj$call$ind_binary
  
  mSpec <- ifelse(class(mgmobj$callmoderators) == "numeric", "vector", "matrix")
  
  
  # Storage
  Pars_ind <- list() # storage for interaction-indicators
  Pars_values <- list() # storage for actual parameters associate for the interactions indexed in Pars_ind
  mgmobj$intercepts <- vector('list', length = p)
  
  for(v in 1:p) {
    
    model_obj <- mgmobj$nodemodels[[v]]$model
    
    
    
    # ----- Create empty Storage for parameters -----
    
    predictor_set <- (1:p)[-v] # set of predictors for variable v (all others)
    
    # A) Indicator Object for all possible interactions: List contains Matrix
    v_Pars_ind <- vector('list', length = d) #  1 = pairwise, 2= threeway, etc.
    
    ## Create Indicator list for Present interactions (separate for k-order MGM and moderated MGM)
    if(is.null(moderators)) {
      
      # If there are no moderators: All interactions of specified degree (k)
      for(ord in 1:d) v_Pars_ind[[ord]] <- t(combn(predictor_set, ord))
      
    } else {
      
      v_Pars_ind[[1]] <- matrix(predictor_set, ncol=1) # main effects
      
      # Moderation effects
      if(mSpec == "vector") {
        
        d <- 2
        
        if(v %in% moderators) {
          ind_mods <- t(combn((1:p)[-v], 2)) # if moderator, all combinations of other variables
        } else {
          ind_mods <- expand.grid((1:p)[-v], moderators[moderators!=v]) # if not moderator, all variables times 
        }
        
        ind_mods <- ind_mods[!apply(ind_mods, 1, function(x) x[1] == x[2]), ] # remove rows with equal entries
        
      }
      
      if(mSpec == "matrix") {
        
        ind_v_inMod <- as.logical(apply(moderators, 1, function(x) v %in% x  ))
        ind_mods <- t(apply(matrix(moderators[ind_v_inMod], ncol=3), 1, function(x) x[x!=v]))
        
      }

      v_Pars_ind[[2]] <- ind_mods
      # no interactions k > 2 allowed, if moderators are specified
      
    } # end if: moderators?
    
    # Make sure all entries of "v_Pars_ind" are matrices
    for(j in 1:d) v_Pars_ind[[j]] <- matrix(as.matrix(v_Pars_ind[[j]]), ncol=j)
    
    no_interactions <- unlist(lapply(v_Pars_ind, nrow))
    
    # B) Parameter Object: Same structure as (A), but now with a list entry for each matrix row
    v_Pars_values <- vector('list', length = d)
    for(ord in 1:d) v_Pars_values[[ord]] <- vector('list', length = no_interactions[ord])
    
    # ----- Fill (B) with parameter estimates -----
    
    # separate for categorical/ non-categorical response node, because in former case we predict K categories
    # and hence need a different data structure
    
    # Categorical Case
    if(type[v] == 'c') {
      
      n_cat <- level[v]
      
      for(cat in 1:n_cat) {
        
        model_obj_i <- as.numeric(model_obj[[cat]]) # select parameters for predicting category i
        interaction_names <- rownames(model_obj[[cat]])[-1] # -1 because we remove the intercept
        interaction_order <- str_count(interaction_names, ":") + 1 # get order of each interaction parameter; +1 to put on same scale as ord, so ord=2 = pairwise interaction
        
        mgmobj$intercepts [[v]][[cat]] <- model_obj_i[1]
        
        # Thresholding:
        # p = number of covariances, as it should be
        model_obj_i_ni <- model_obj_i[-1] # remove intercept, this is no covariate
        if(threshold == 'LW') tau <- sqrt(d) * sqrt(sum(model_obj_i_ni^2)) * sqrt(log(npar_standard[v]) / n)
        if(threshold == 'HW') tau <- d * sqrt(log(npar_standard[v]) / n)
        if(threshold == 'none') tau <- 0
        model_obj_i[abs(model_obj_i) < tau] <- 0 # set all parameter estimates below threshold to zero
        mgmobj$nodemodels[[v]]$tau <- tau # Save tau
        
        
        for(ord in 1:d) {
          
          part_ord <- model_obj_i[-1][ord == interaction_order] # parameters for interaction of order = ord+1
          gotThemall <- rep(TRUE, length(part_ord)) # sanity check: did I copy all parameters
          int_names_ord <- interaction_names[ord == interaction_order]
          
          if(no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Internal Error: Error in parameter extraction type 1')
          
          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          
          for(t in 1:no_interactions[ord]) {
            
            # indicates location of parameters for given interaction
            
            for(cc in 1:ord) find_int_dummy[, cc] <- grepl(paste0('V', v_Pars_ind[[ord]][t, cc], '.'), int_names_ord,
                                                           int_names_ord,
                                                           fixed = TRUE) # not only single chars have to be contained, but exact order)
            select_int <- rowSums(find_int_dummy) == ord # because threeway interaction has 2 predictors; ord = order of interactions in joint distribution
            
            # fill in paramters + rownames into list
            parmat <- matrix(part_ord[select_int], ncol = 1)
            rownames(parmat) <- int_names_ord[select_int]
            v_Pars_values[[ord]][[t]][[cat]] <- parmat
            
            gotThemall[select_int == TRUE] <- FALSE
            
          }
          
          if(sum(gotThemall) > 0) stop('Internal Error: Error in parameter extraction type 2')
          
        } # end for: ord
        
      } # end for: categories
      
      # Continuous Case
    } else {
      
      model_obj_i <- as.numeric(model_obj) # select parameters for predicting category i
      interaction_names <- rownames(model_obj)[-1] # -1 because we remove the intercept
      interaction_order <- str_count(interaction_names, ":") + 1 # on same scale as ord, so ord=2 = pairwise interaction
      
      mgmobj$intercepts [[v]] <- model_obj_i[1]
      
      # Thresholding:
      # p = number of covariances, as it should be
      model_obj_i_ni <- model_obj_i[-1] # remove intercept, this is no covariate
      if(threshold == 'LW') tau <- sqrt(d) * sqrt(sum(model_obj_i_ni^2)) * sqrt(log(npar_standard[v]) / n)
      if(threshold == 'HW') tau <- d * sqrt(log(npar_standard[v]) / n)
      if(threshold == 'none') tau <- 0
      model_obj_i[abs(model_obj_i) < tau] <- 0 # set all parameter estimates below threshold to zero
      mgmobj$nodemodels[[v]]$tau <- tau # Save tau
      
      for(ord in 1:d) {
        
        if(no_interactions[ord] > 0) { # can be zero for moderated MGMs
          
          part_ord <- model_obj_i[-1][ord == interaction_order] # parameters for interaction of order = ord
          gotThemall <- rep(TRUE, length(part_ord)) # sanity check: did I copy all parameters
          int_names_ord <- interaction_names[ord == interaction_order]
          
          if(no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Fuckup in parameter extraction 1')
          
          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          
          for(t in 1:no_interactions[ord]) {
            
            # indicates location of parameters for given interaction
            for(cc in 1:ord) find_int_dummy[, cc] <- grepl(paste0('V', v_Pars_ind[[ord]][t, cc], '.'),
                                                           int_names_ord,
                                                           fixed = TRUE) # not only single chars have to be contained, but exact order
            select_int <- rowSums(find_int_dummy) == (ord) # because threeway interaction has 2 predictors; ord = order of interactions in joint distribution
            
            # fill in paramters + rownames into list
            parmat <- matrix(part_ord[select_int], ncol = 1)
            rownames(parmat) <- int_names_ord[select_int]
            v_Pars_values[[ord]][[t]] <- parmat
            
            gotThemall[select_int == TRUE] <- FALSE
            
          }
          
          if(sum(gotThemall) > 0) stop('Fuckup in parameter extraction 2')
          
        } # end if: no_interactions[ord] > 0
        
      } # end for: ord
      
    }
    
    Pars_ind[[v]] <- v_Pars_ind
    Pars_values[[v]] <- v_Pars_values
    
  } # end for: v
  
  
  
  # --------------------------------------------------------------------------------------------
  # -------------------- Postprocess Regression Estimates into (Factor) Graph Structure --------
  # --------------------------------------------------------------------------------------------
  
  
  # ----- Reduce to 1 parameter per Factor, apply AND rule and get sign -----
  
  # Combine interactions from all neighborhood regressions in one structure
  
  # We turn around the nesting to be able to collapse across over v
  # In addition we append the estimated node to 'complete' the interaction
  
  # Storage for indicator list
  Pars_ind_flip <- vector('list', length = d)
  dummy_list_flip <- vector('list', length = p)
  Pars_ind_flip <- lapply(Pars_ind_flip, function(x) dummy_list_flip)
  
  # Storage for value list
  Pars_values_flip <- vector('list', length = d)
  Pars_values_flip <- lapply(Pars_values_flip, function(x) dummy_list_flip)
  
  
  # Reordering so I can use do.call() below on inner list level
  for(v in 1:p) {
    for(ord in 1:d) {
      
      # Reordering indicator list
      Pars_ind_part <- Pars_ind[[v]][[ord]]
      colnames(Pars_ind_part) <- NULL
      Pars_ind_part <- as.matrix(Pars_ind_part)
      
      Pars_ind[[v]][[ord]] <- cbind(rep(v, nrow(Pars_ind_part)), Pars_ind_part)
      Pars_ind_flip[[ord]][[v]] <- Pars_ind[[v]][[ord]]
      
      # Reordering value list
      Pars_values_flip[[ord]][[v]] <- Pars_values[[v]][[ord]]
      
    }
  }
  
  # Collapse indicator list across nodes
  Pars_ind_flip_red <- lapply(Pars_ind_flip, function(x) do.call(rbind, x)  )
  # Collapse value list across nodes
  Pars_values_flip_red <- lapply(Pars_values_flip, function(x) do.call(c, x))
  
  # browser()
  
  # 1) Select each interaction
  
  # Compute number of interactions for each order
  #     Note that we could take this information also from the design matrices; however, we compute it theoretically as a sanity check
  n_terms_d <- rep(NA, d)
  if(is.null(moderators)) {
    for(ord in 1:d) n_terms_d[ord] <- choose(p, ord+1)
  } else {
    n_terms_d[1] <- choose(p, 1+1) # all pairwise interactions
    
    # Select moderation terms: Standard specification
    if(mSpec == "vector") {
      mod_terms <- expand.grid((1:p), (1:p), moderators) 
      id_uni <- FlagSymmetricFast(mod_terms)
      mod_terms2 <- mod_terms[!duplicated(id_uni), ]
      ind_diff <- as.numeric(apply(mod_terms2, 1, function(x) !any(duplicated(x))))
      mod_terms3 <- mod_terms2[ind_diff == 1, ]  
    }
    
    # Select moderation terms: Custom specification
    if(mSpec == "matrix") {
      mod_terms3 <- mgmobj$call$moderators
    }
    
    n_terms_d[2] <- nrow(mod_terms3) # ok to have interactions double; will remove them below using FlagSymmetricFast()
  }
  
  # Set up target data structure
  l_factors <- list() # saves all unique interavtions
  for(ord in 1:d) l_factors[[ord]] <- matrix(NA, nrow = n_terms_d[ord], ncol = ord+1)
  l_factor_par <- list() # saves parameters associated with all unique interactions
  for(ord in 1:d) l_factor_par[[ord]] <- vector('list', length = n_terms_d[ord])
  l_sign_par <- list() # saves sign (if defined) of all unique interactions
  for(ord in 1:d) l_sign_par[[ord]] <- rep(NA, n_terms_d[ord])
  
  l_factor_par_full <- l_factor_par_AggNodewise <- l_factor_par_SignNodewise <- l_factor_par # for un-aggregated parameter esimates
  
  # Define set of continous and binary variables: for interactions between these we can assign a sign
  # Depends on binarySign
  if(binarySign) {
    set_signdefined <- c(which(type == 'p'), which(type == 'g'), ind_cat[ind_binary])
  } else {
    set_signdefined <- c(which(type == 'p'), which(type == 'g'))
  }
  
  
  
  # Loop over: order of interactions (ord = 1 = pairwise)
  for(ord in 1:d) {
    
    set_int_ord <- Pars_ind_flip_red[[ord]]
    set_par_ord <- Pars_values_flip_red[[ord]]
    row.names(set_int_ord) <- NULL
    
    ids <- FlagSymmetricFast(x = set_int_ord) # BOTTLE NECK, now better with native solution, but still problematic for huge p
    
    # Get set of unique interactions
    unique_set_int_ord <- cbind(set_int_ord, ids)[!duplicated(ids), ]
    unique_set_int_ord <- matrix(unique_set_int_ord, ncol = ord+1+1)
    n_unique <- nrow(unique_set_int_ord)
    
    # loop over: unique interaction of order = ord
    for(i in 1:n_unique) {
      
      l_w_ind <- list()
      l_w_par <- list()
      ind_inter <- which(ids == i)
      
      # loop over the k estimates for a given k-order interaction
      for(j in 1:(ord+1)) {
        l_w_ind[[j]] <- set_int_ord[ind_inter[j], ]
        l_w_par[[j]] <- set_par_ord[[ind_inter[j]]]
      }
      
      # Mapping: Regression parameters -> Edge parameter (mean of absolute values of parameters)
      m_par_seq <- unlist(lapply(l_w_par, function(x) mean(abs(unlist(x)))))
      m_sign_seq <- unlist(lapply(l_w_par, function(x) {
        x <- unlist(x)
        if(length(x)>1) 0 else sign(x)
      } ))
      
      # Apply AND / OR rule
      if(ruleReg == 'AND') parcompr <- mean(m_par_seq) * !(0 %in% m_par_seq)
      if(ruleReg == 'OR') parcompr <- mean(m_par_seq)
      
      # Compute Sign if defined
      if(mean(m_par_seq) != 0) { # only relevant for nonzero parameters
        
        pair <- l_w_ind[[1]] # order doesn't matter, could take any but the first entry is always filled independent of "ord", so first
        
        if(sum(!(pair %in% set_signdefined)) == 0) { # check of all involved varibales are g, p, or binary
          
          # Computes combined sign (if defined) over k terms for same interaction
          sign_object <- getSign(l_w_ind, 
                                 l_w_par,
                                 type,
                                 set_signdefined,
                                 overparameterize,
                                 ord)
          
          int_sign <- sign_object$voteSign
          
          
        } else {
          int_sign <- 0 # if not defined (set_signdefined): 0
        }
      } else {
        int_sign <- NA # if no edge present: NA (defined but zero parameter estimate, so no sign available)
      }
      
      ## Get sign
      l_sign_par[[ord]][i] <- int_sign
      
      # Save indicator
      l_factors[[ord]][i, ] <- l_w_ind[[1]] # just choose first one, doesn't matter
      
      # Save edge weight
      
      for(i_ord in 1:(ord+1)) {
        l_factor_par_AggNodewise[[ord]][[i]][[i_ord]] <- m_par_seq[i_ord]
        l_factor_par_SignNodewise[[ord]][[i]][[i_ord]] <- m_sign_seq[i_ord]
      }
      
      l_factor_par[[ord]][[i]] <- parcompr
      l_factor_par_full[[ord]][[i]] <- l_w_par
      
    } # end for: i (unique interactions)
    
  } # end for: ord
  
  
  # -------------------- Compute Weighted Adjacency matrix (pairwise Interactions) -------------------
  
  # We copy the objects from above, and delete rows in l_factors if l_factor_par is zero
  l_factors_nz <- l_factors
  l_factor_par_nz <- l_factor_par
  l_factor_par_full_nz <- l_factor_par_full
  l_sign_par_nz <- l_sign_par
  
  for(ord in 1:d) {
    zero_indicator <-  which(unlist(lapply(l_factor_par[[ord]], function(x) x == 0 )))
    
    # Delete rows in l_factors
    if(length(zero_indicator) == 0) {
      l_factors_nz[[ord]] <- l_factors[[ord]]
    } else {
      l_factors_nz[[ord]] <- l_factors[[ord]][-zero_indicator,]
    }
    
    # Delete corresponding entries in l_factor_par
    l_factor_par_nz[[ord]] <- l_factor_par_full_nz[[ord]] <- list()
    counter <- 1
    for(k in 1:length(l_factor_par[[ord]])) {
      if(!(k %in% zero_indicator)) {
        l_factor_par_nz[[ord]][[counter]] <- l_factor_par[[ord]][[k]]
        l_factor_par_full_nz[[ord]][[counter]] <- l_factor_par_full[[ord]][[k]]
        counter <- counter + 1
      }
    }
    
    # Delete entries in l_sign_par_nz
    if(length(zero_indicator) == 0) {
      l_sign_par_nz[[ord]] <- l_sign_par[[ord]]
    } else {
      l_sign_par_nz[[ord]] <- l_sign_par[[ord]][-zero_indicator]
    }
    
  }
  
  # Save in output
  mgmobj$interactions$indicator <- l_factors_nz
  mgmobj$interactions$weightsAgg <- l_factor_par_nz
  mgmobj$interactions$weights <- l_factor_par_full_nz
  mgmobj$interactions$signs <- l_sign_par_nz
  
  # browser()
  
  
  # ---------- Fill into p x p Matrix ---------
  
  m_signs <- matrix(NA, p, p)
  wadj <-  matrix(0, p, p)
  edges <- matrix(l_factors_nz[[1]], ncol=2) # to avoid error if only 1 row and list entry = numeric
  n_edges <- nrow(edges)
  
  if(n_edges > 0) {
    for(i in 1:n_edges) {
      wadj[edges[i,1], edges[i,2]] <- wadj[edges[i,2], edges[i,1]] <- l_factor_par_nz[[1]][[i]]
      m_signs[edges[i,1], edges[i,2]] <- m_signs[edges[i,2], edges[i,1]] <- l_sign_par_nz[[1]][[i]]
    }
  }
  
  # Create sign color matrix
  sign_colors <- matrix('darkgrey', p, p)
  sign_colors[m_signs == 1] <- 'darkgreen'
  sign_colors[m_signs == -1] <- 'red'
  
  # Save in output
  mgmobj$pairwise$wadj <- wadj
  mgmobj$pairwise$signs <- m_signs
  mgmobj$pairwise$edgecolor <- sign_colors
  
  
  # ---------- Fill into p x p Nodewise Matrix ---------
  
  m_wadj <-  m_signs <- matrix(0, p, p)
  
  # get table of unique pairwise interactions (copied from above)
  ord <- 1
  set_int_ord <- Pars_ind_flip_red[[ord]]
  row.names(set_int_ord) <- NULL
  ids <- FlagSymmetricFast(x = set_int_ord) # BOTTLE NECK, now better with native solution, but still problematic for huge p
  unique_set_int_ord <- cbind(set_int_ord, ids)[!duplicated(ids), ]
  unique_set_int_ord <- matrix(unique_set_int_ord, ncol = ord+1+1)
  
  n_edges <- nrow(unique_set_int_ord)
  ED <- unique_set_int_ord
  
  for(i in 1:n_edges) {
    m_wadj[ED[i,1], ED[i,2]] <- l_factor_par_AggNodewise[[1]][[i]][[2]]
    m_wadj[ED[i,2], ED[i,1]] <- l_factor_par_AggNodewise[[1]][[i]][[1]]
    
    m_signs[ED[i,1], ED[i,2]] <- l_factor_par_SignNodewise[[1]][[i]][[2]]
    m_signs[ED[i,2], ED[i,1]] <- l_factor_par_SignNodewise[[1]][[i]][[1]]
  }
  
  
  # Create sign color matrix
  sign_colors <- matrix('darkgrey', p, p)
  sign_colors[m_signs == 1] <- 'darkgreen'
  sign_colors[m_signs == -1] <- 'red'
  
  # Save in output
  mgmobj$pairwise$wadjNodewise <- m_wadj
  mgmobj$pairwise$signsNodewise <- m_signs
  mgmobj$pairwise$edgecolorNodewise <- sign_colors
  
  
  return(mgmobj)
  
  
} # end of function

