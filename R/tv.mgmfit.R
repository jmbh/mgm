


tv.mgmfit <- function(data, # input tv mgm
                      type, 
                      lev,
                      timepoints = NA,
                      estpoints = 'default', # either default, integer, or vector indicating where should be estimated
                      bandwidth,
                      gam = .25,                       
                      d = 2, 
                      rule.reg = 'AND', 
                      pbar = TRUE,
                      method ='glm',
                      missings = 'error',
                      ret.warn = FALSE) 

{
  
  
  # ---------- Compute Aux Variables ----------
  
  n <- nrow(data)
  nNode <- ncol(data)
  
  
  # ---------- Assign each data point a value on [0,1] normed time interval ----------
  
  if(!is.na(timepoints[1])) {
    if(length(timepoints) != nrow(data)) stop('Length of timepoints vector must match the number of rows in data.')
    
    # Scale to [0,1]
    timepoints_sc <- timepoints - min(timepoints)
    timepoints_sc <- timepoints_sc / max(timepoints_sc)
    
    # Input check: monotone increasing?
    mon_inc_check <- (n-1)==sum((timepoints_sc - c(0, timepoints_sc[-length(timepoints_sc)]))[-1]>0) # all differences have to be positive
    if(!mon_inc_check) stop('The vector timepoints must contain a positive and strictly increasing sequence.')
    
    tpoints <- timepoints_sc
    
  } else {
    
    tpoints <- seq(0, 1, length = n) # n sequence before and after for folded gaussian weight
    
  }
  
  
  # ---------- Determine Location where model is estimated ----------
  
 
  # We select a reasonable number of equally time-spaced estimation locations
  if(estpoints[1]=='default') {
    
    est_points_norm <- seq(0, 1, length = round(2/bandwidth)) # where do we estimate?
  
  # We use a user-specified number of equally time-spaced estimation locations  
  } else if(round(estpoints[1])==estpoints[1] & length(estpoints)==1) {
    
    # Input Check
    if(estpoints < 3) stop('Please choose a number of estimated time points > 2')
    est_points_norm <- seq(0, 1, length = estpoints)
  
  # We use a user-specified locations (which are not necessarily equally time-spaced)  
  } else if(length(estpoints) > 1) {
    
    # Input Check
    if(sum(estpoints>0) < length(estpoints)) stop('Time points at which the model should be estimated have to be positive')
    
    estpoints_sc <- estpoints - min(estpoints)
    estpoints_sc <- estpoints_sc / max(estpoints_sc)
   
    est_points_norm <- estpoints_sc
     
  }

  
  Nest <- length(est_points_norm) # number of estimated Locations
  
  # ---------- Fit weighted Model at each time point ----------
  
  if(pbar==TRUE) {pb <- txtProgressBar(min = 0, max = Nest, initial=0, char="-", style = 3) } # Set up: Progress bar
  
  # Create Storage
  a_wadj <- a_signs <- a_edgecol <- array(NA, dim=c(nNode, nNode, Nest))
  l_mpar.matrix <- list()
  tv_list <- list() 
  
  for(t in 1:Nest) {
    
    kernel_mean <- est_points_norm[t]
    
    # Define weights with folded gaussian kernel    
    weights <- dnorm(tpoints, 
                     mean = kernel_mean, 
                     sd = bandwidth) # we use unit time interval, so we choose the bw relative to that
    weights <- weights / max(weights) # doesn't really matter, but nice to have it on the same scale for all models
    
    # Estimate weighted models at Nest time locations
    tv_list[[t]] <- mgmfit_core(data = data, 
                                type = type, 
                                lev = lev, 
                                lambda.sel = "EBIC",  # force EBIC, CV makes no sense in nonstationary data
                                folds = 10, 
                                gam = gam, 
                                d = d, 
                                rule.reg = rule.reg, 
                                pbar = FALSE, 
                                method = method, 
                                missings = missings, 
                                weights = weights, 
                                ret.warn = FALSE)
    
    # From list to array
    a_wadj[,,t] <- tv_list[[t]]$wadj
    a_signs[,,t] <- tv_list[[t]]$signs
    a_edgecol[,,t] <- tv_list[[t]]$edgecolor
    l_mpar.matrix[[t]] <- tv_list[[t]]$mpar.matrix
    
    if(pbar==TRUE) setTxtProgressBar(pb, t)  # Update progress bar
    
  }
  
  # write model parameter matrix in array
  a_mpar.matrix <- array(NA, dim=c(dim(l_mpar.matrix[[t]])[1], dim(l_mpar.matrix[[t]])[2],Nest))
  
  # get weighted N for each time step; this will later also be used as a missing data signal
  l_wN <- unlist(lapply(tv_list, function(x) sum(x$call$weights)))
  
  
  # ---------- Output ----------
  
  l_call <- list('type' = type, 
                 'lev' = lev, 
                 'estpoints' = estpoints,
                 'est_points_norm' = est_points_norm,
                 'bandwidth' = bandwidth, 
                 'lambda.sel' = "EBIC",
                 'gam' = gam, 
                 'd' = d, 
                 'rule.reg' = rule.reg, 
                 'method' = method)
  
  outlist <- list('call' = l_call,
                  'wadj' = a_wadj, 
                  'mpar.matrix' = l_mpar.matrix, 
                  'sign' = a_signs, 
                  'edgecolor' = a_edgecol, 
                  't.models' = tv_list, 
                  'Nt' = l_wN)
  
  class(outlist) <- c('mgm', 'tv.mgm')
  
  # Return estimation messages:
  estimation_msg('tv.mgmfit') # note about where signs are stored
  
  return(outlist) 
  
}






