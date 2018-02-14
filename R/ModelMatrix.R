



ModelMatrix <- function(data,  # matrix
                        type,  # type vector (I think not needed, level should be sufficient)
                        level, # level vector
                        labels,
                        d      , # largest neighborhood size
                        moderators,
                        v,
                        allCats = FALSE # if true, the model matrix does not use all unique categories, but the categories specified in level, this exists because I use this function within the sampling function
)
  
  
{
  
  
  # ---------- Calculate Auxilliary variables ----------
  
  p <- ncol(data) # note that we have only the predictors here!!
  n <- nrow(data)
  
  # ---------- Input Checks ----------
  
  if(p != length(type)) stop('Length of type has to match the number of columns in data.')
  if(p != length(level)) stop('Length of level has to match the number of columns in data.')
  if(!(class(level) == 'integer' | class(level) == 'numeric')) stop('level has to be an integer vector.')
  
  
  # ---------- Calculate Indicator Functions for all Variables ----------
  
  l_ind_datasets <- list()
  for(j in 1:p) {
    
    if(type[j] != 'c') {
      
      l_ind_datasets[[j]] <- as.matrix(data[, j])
      colnames(l_ind_datasets[[j]]) <- paste0(labels[j])
      
    } else {
      
      if(allCats==FALSE) {
        
        unique_labels <- unique(data[, j])
        unique_labels_sorted <- sort(unique_labels)
        n_labels <- length(unique_labels_sorted)
        ind_matrix <- matrix(NA, nrow = n, ncol=n_labels)
        for(s in 1:n_labels) ind_matrix[, s] <- data[, j] == unique_labels_sorted[s]
        
      } else {
        
        # This is used in mvarsampler() and avoids that the design matrix is too small cases early in the time series, where not all categories are seen yet
        unique_labels <- 1:level[j]
        n_labels <- level[j]
        ind_matrix <- matrix(NA, nrow = n, ncol=n_labels)
        for(s in 1:n_labels) ind_matrix[, s] <- data[, j] == unique_labels[s]
        
      }
      
      # Add Varnames
      cn <- paste0(labels[j], 1:n_labels)
      colnames(ind_matrix) <- cn
      l_ind_datasets[[j]] <- ind_matrix
      
    }
    
  }
  
  # ---------- Collect d = 1 interaction Terms ----------
  
  # In this case we just use the indicator functions computed above
  
  l_ind_datasets_nV <- l_ind_datasets
  Xd1 <- do.call(cbind, l_ind_datasets_nV)
  
  # ---------- Compute all d>1 Interaction Terms ----------
  
  # Here we loop over all orders, in each order over all interactions, and in each interaction over all combinations of the levels of all variables
  
  
  if(!is.null(moderators)) d <- 2 # If first-order moderation (three way interactions) is specified, we fix k = 3, d = k - 1 = 2
  
  if(d > 1) {
    
    ## List all possible Interactions
    l_interactions <- vector('list', length = d)
    
    # Case I: No Moderation
    if(is.null(moderators)) {
      for(ord in 1:d) l_interactions[[ord]] <- combn((1:p), ord, simplify = FALSE) # note that the numbers here refer to the all variables minus variable v
    } else {

      l_interactions[[1]] <- list()
      for(i in 1:p) l_interactions[[1]][[i]] <- i
      
      # Case II: No Moderation
    # if(v == 2)  browser()
      
      if(v %in% moderators) {
        
        l_interactions[[2]] <- combn(1:p, 2, simplify = FALSE) # all combinations of the remaining (here denoted by 1:p) variables

      } else {
        
        ind_mod_MM <- (1:(p+1) %in% moderators)[-v] # indicator(moderator?) for 1:p predictors
        n_mods <- sum(ind_mod_MM)
        which_mod <- which(ind_mod_MM)
        
        l_mods <- list()
        for(i in 1:n_mods) l_mods[[i]] <- expand.grid((1:p)[-which_mod[i]], which_mod[i]) # loop over moderators
        m_mods <- do.call(rbind, l_mods)
        m_mods <- as.matrix(m_mods)
        
        l_mods_combn <- list()
        for(i in 1:nrow(m_mods)) l_mods_combn[[i]] <- m_mods[i, ] # turn each row into list entry
        
        l_interactions[[2]] <- l_mods_combn
        
      } # end: v = moderator?

    } # end if: moderators?
     
    
    
    # storage for all interactions of all orders
    l_collect_terms <- list()
    
    # Loop over order of interactions;
    for(ord in 2:d) {
      
      # if(ord==3) browser()
      
      n_terms <- length(l_interactions[[ord]])
      
      # storage for all interactions of fixed order
      l_ord_terms <- list()
      
      # Loop over the interactions of a given order
      for(it in 1:n_terms) {
        
        # Storage: collect here all interaction terms of one interaction
        l_it_terms <- list()
        
        # For fixed interaction: which variables are involved?
        inter_it <- l_interactions[[ord]][[it]] # select interactions one by one
        
        # Compute amount of levels of each variable (continuous=1)
        l_indicator_it <- list()
        for(i in 1:ord) l_indicator_it[[i]] <- 1:level[inter_it[i]]
        
        # List all combination of levels of variables in the interaction it
        all_combs <- expand.grid(l_indicator_it)
        n_combs <- nrow(all_combs)
        
        # Loop over all combinations of levels of variables in interaction and collect in matrix
        for(comb in 1:n_combs) {
          tarmat <- matrix(NA, nrow = n, ncol = ord)
          for(i in 1:ord) tarmat[, i] <- as.matrix(l_ind_datasets[[inter_it[i]]])[,all_combs[comb, i]]
          l_it_terms[[comb]] <- apply(tarmat, 1, prod)
        }
        
        # combine all level combinations for interaction it
        it_data <- do.call(cbind, l_it_terms)
        
        # Create Names for all (combinations of) interactions
        all_combs_char <- apply(all_combs, 2, function(x) {
          if(length(unique(x))==1) {
            out <- rep("", length(x))
          } else {
            out <- x
          }
        })
        cn <- rep(NA, n_combs)
        
        # DEV: check whether this still works:
        
        for(comb in 1:n_combs) cn[comb] <- paste0(labels[inter_it], matrix(all_combs_char, ncol=ord)[comb,], collapse=':')
        # the point above added to deliminate category flag
        
        colnames(it_data) <- cn
        l_ord_terms[[it]] <- it_data
        
      }
      
      # Collapse over interactions of fixed order
      l_collect_terms[[ord]] <- do.call(cbind, l_ord_terms)
      
    }
    
    # Collapse across order of interactions
    # l_collect_terms[[1]] <- NULL
    all_HOI_terms <- do.call(cbind, l_collect_terms)
    
  }
  
  # Combine with d=1 size neighborhoods (singletons)
  
  if(d > 1) {
    X <- cbind(Xd1, all_HOI_terms)
  } else {
    X <- Xd1
  }
  
  
  # ---------- Output----------
  
  return(X)
  
  
} # end of function


