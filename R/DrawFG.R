
# ---------- 1) Standard version for mgm() ----------

DrawFG <- function(object,
                   PairwiseAsEdge = TRUE)
  
  
{
  
  
  # Take from model object
  list_ind <- object$factors$indicator
  list_weights <- object$factors$weightsAgg
  list_signs <- object$factors$signs
  
  p <- length(object$call$type)
  
  # Create empty matrix of node entities
  n_factors <- length(list_ind)
  size_factors <- rep(NA, n_factors)
  for(or in 1:n_factors) size_factors[or] <- length(list_ind[[or]]) / (or+1) # can't use apply, because list_ind might have only 1 row
  
  pF <- sum(size_factors) # how many factors in total
  
  # Subset of orders of interactions with at least 1 estimated interaction
  nz_ind <- which(size_factors > 0)
  
  # Make all list entries matrices to avoid trouble below
  one_ind <- which(size_factors == 1)
  for(i in one_ind) list_ind[[i]] <- matrix(list_ind[[i]], nrow = 1)
  
  # Create empty factor graph
  G <- Gw <- matrix(0, p+pF, p+pF)
  Gsign <- matrix(NA, p+pF, p+pF)
  
  counter <- p+1
  
  # loop over order of interactions
  for(or in nz_ind) {
    
    # Loop over rows in fixed order
    Nro <- nrow(list_ind[[or]])  #  now nrow() possible in case there is only 1 interaction, then error because no 2d object
    
    if(length(Nro) > 0) {
      
      for(r in 1:Nro) {
        
        # Loop over k connections of k-order factor
        for(or2 in 1:(or+1)) {
          
          G[counter, list_ind[[or]][r, or2]] <- G[list_ind[[or]][r, or2], counter] <- 1
          Gw[counter, list_ind[[or]][r, or2]] <- Gw[list_ind[[or]][r, or2], counter] <- list_weights[[or]][[r]]
          Gsign[counter, list_ind[[or]][r, or2]] <- Gsign[list_ind[[or]][r, or2], counter] <- list_signs[[or]][[r]]
          
        }
        
        counter <- counter + 1
        
      }
      
    }
  }
  
  
  # Define shape
  nodetype <- c(rep(0, p), rep(1, pF))
  
  
  # Caclulate order indicator
  n_order <- length(size_factors)
  l_oi <- list()
  for(i in 1:n_order) l_oi[[i]] <- rep(i, size_factors[i])
  
  order_ind <- c(rep(0,p), unlist(l_oi))
  
  # Calculate Color Matrix
  Gcol <- matrix('darkgrey', p+pF, p+pF)
  Gcol[Gsign == 1] <- 'darkgreen'
  Gcol[Gsign == -1] <- 'red'
  

  
  
  # Simplify Factor Graph by dropping factors for pairwise (2way) interactions, and putting edges back in  
  if(PairwiseAsEdge) {
    
    # browser()
    
    # Edges
    Gw <- Gw[order_ind!=1, order_ind!=1]
    Gw[1:p, 1:p] <- object$pairwise$wadj
    
    Gsign <- Gsign[order_ind!=1, order_ind!=1]
    Gsign[1:p, 1:p] <- object$pairwise$signs
    
    Gcol <- Gcol[order_ind!=1, order_ind!=1]
    Gcol[1:p, 1:p] <- object$pairwise$edgecolor
    
    nodetype <- nodetype[order_ind != 1]
    order_ind <- order_ind[order_ind != 1]
    
  }
  
  # Export
  outlist <- list("graph" = G,
                  'weightedgraph' = Gw,
                  'signs' = Gsign,
                  'signcolor' = Gcol,
                  'nodetype' = nodetype,
                  'order' = order_ind)
  
  
  return(outlist)
  
} # eOF











# ---------- 2) Time-varying version for mgm() ----------

# Includes all factors that are at least once present in the time series
# (Important to investigate factor weight across time)



DrawFGtv <- function(object, # list of all interactions that are estimated nonzer at at least one estimation point
                     PairwiseAsEdge = TRUE, 
                     estpoint # number of (initial) variables
                     
)
  
  
{
  
  
  list_ind <- object$interactions$indicator
  list_weights <- object$interactions$weightsAgg
  list_signs <- object$interactions$signs
  
  # ----------- Obtain set of all interactions that are present at least once ----------
  
  d <- object$call$k - 1
  p <- length(object$call$type)
  
  no_estpoints <- length(object$call$estpoints)
  
  # Obtain set of all interactions that occur at least once in the time series
  l_all_indicators <- list()
  l_unique_indicators <- vector('list', length = d)
  l_unique_indicators <- lapply(l_unique_indicators, function(x) vector('list', length=no_estpoints))
  for(i in 1:no_estpoints)  {
    l_all_indicators[[i]] <- object$tvmodels[[i]]$interactions$indicator
    for(ii in 1:d) l_unique_indicators[[ii]][[i]] <- l_all_indicators[[i]][[ii]]
  }
  
  # get unique set of interactions for each order
  l_unique_indicators <- lapply(l_unique_indicators, function(x) {
    x_comb <- do.call(rbind, x)
    x_comb[!duplicated(x_comb),]
  })
  
  
  # Select output of selected estpoint:
  list_ind <- object$tvmodels[[estpoint]]$interactions$indicator
  list_weights <- object$tvmodels[[estpoint]]$interactions$weightsAgg
  list_signs <- object$tvmodels[[estpoint]]$interactions$signs
  list_ind_all <- l_unique_indicators
  
  
  # ----------- Loop over Time points ----------
  

  # Create empty matrix of node entities
  n_factors <- length(list_ind_all)
  # all factor list
  size_factors <- rep(NA, n_factors)
  for(or in 1:n_factors) size_factors[or] <- length(list_ind_all[[or]]) / (or+1)
  # factor at hand list
  size_factors_fixed <- rep(NA, n_factors)
  for(or in 1:n_factors) size_factors_fixed[or] <- length(list_ind[[or]]) / (or+1)
  
  pF <- sum(size_factors) # how many factors in total
  
  # Subset of orders of interactions with at least 1 estimated interaction
  nz_ind <- which(size_factors_fixed > 0)
  one_ind <- which(size_factors_fixed == 1)
  
  # Make all list entries matrices to avoid trouble below
  for(i in one_ind) list_ind[[i]] <- matrix(list_ind[[i]], nrow = 1)
  
  # Create empty factor graph
  Gw <- matrix(0, p+pF, p+pF)
  Gsign <- matrix(NA, p+pF, p+pF)
  
  counter <- p+1
  
  
  
  # loop over order of interactions
  for(or in nz_ind) {
    # Loop over rows in fixed order
    
    Nro <- nrow(list_ind_all[[or]])  #  now nrow() possible in case there is only 1 interaction, then error because no 2d object
    
    if(length(Nro) > 0) {
      
      for(r in 1:Nro) {
        
        # is the present row in list_ind_all contained in list_ind?
        ind_where <- apply(list_ind[[or]], 1, function(x, want) isTRUE(all.equal(x, want)), list_ind_all[[or]][r, ])
        
        # if yes:
        if(any(ind_where)) {
          
          # Loop over k connections of k-order factor
          for(or2 in 1:(or+1)) {
            
            # browser()
            Gw[counter, list_ind_all[[or]][r, or2]] <- Gw[list_ind_all[[or]][r, or2], counter] <- list_weights[[or]][[which(ind_where)]]
            Gsign[counter, list_ind_all[[or]][r, or2]] <- Gsign[list_ind_all[[or]][r, or2], counter] <- list_signs[[or]][[which(ind_where)]]
          }
          
        } # if no: just leave the zero
        
        counter <- counter + 1
        
      }
      
    } # end if: at least 1 row present?
    
  }
  
  
  # Define shape
  nodetype <- c(rep(0, p), rep(1, pF))
  
  # Caclulate order indicator
  n_order <- length(size_factors)
  l_oi <- list()
  for(i in 1:n_order) l_oi[[i]] <- rep(i, size_factors[i])
  
  order_ind <- c(rep(0,p), unlist(l_oi))
  
  # Calculate Color Matrix
  Gcol <- matrix('darkgrey', p+pF, p+pF)
  Gcol[Gsign == 1] <- 'darkgreen'
  Gcol[Gsign == -1] <- 'red'
  
  # browser()
  
  # Simplify Factor Graph by dropping factors for pairwise (2way) interactions, and putting edges back in  
  if(PairwiseAsEdge) {
    
    # browser()
    
    # Edges
    Gw <- Gw[order_ind!=1, order_ind!=1]
    Gw[1:p, 1:p] <- object$pairwise$wadj[, , estpoint]
    
    Gsign <- Gsign[order_ind!=1, order_ind!=1]
    Gsign[1:p, 1:p] <- object$pairwise$signs[, , estpoint]
    
    Gcol <- Gcol[order_ind!=1, order_ind!=1]
    Gcol[1:p, 1:p] <- object$pairwise$edgecolor[, , estpoint]
    
    nodetype <- nodetype[order_ind != 1]
    order_ind <- order_ind[order_ind != 1]
    
  }
  
  # Export
  outlist <- list('weightedgraph' = Gw,
                  'signs' = Gsign,
                  'signcolor' = Gcol,
                  'nodetype' = nodetype,
                  'order' = order_ind)
  
  return(outlist)
  
} # eOF





