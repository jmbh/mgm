

tv.mgmfit <- function(data, # input tv mgm
                   type, 
                   lev,
                   timepoints = NA,
                   tsteps = 'default', 
                   bandwidth,
                   gam = .25,                       
                   d = 2, 
                   rule.reg = 'AND', 
                   pbar = TRUE,
                   method ='glm',
                   missings = 'error',
                   ret.warn = FALSE) {
  
  n <- nrow(data)
  nNode <- ncol(data)
  
  # create time sequence
  if(!is.na(timepoints[1])) {
    if(length(timepoints) != nrow(data)) stop('Length of timepoints vector must match the number of rows in data.')
    # scale to [0,1]
    timepoints_sc <- timepoints - min(timepoints)
    timepoints_sc <- timepoints_sc / max(timepoints_sc)
    
    # check whether monotone increasing
    mon_inc_check <- (n-1)==sum((timepoints_sc - c(0, timepoints_sc[-length(timepoints_sc)]))[-1]>0) # all differences have to be positive
    if(!mon_inc_check) stop('The vector timepoints must contain a positive and strictly increasing sequence.')
    tpoints <- timepoints_sc
  } else {
    tpoints <- seq(0,1, length=n) # n sequence before and after for folded gaussian weight
  }

  if(tsteps=='default') {tsteps <- round(2/bandwidth)} # reasonable coveres the data for any choice of bandwidth
  a_wadj <- a_signs <- a_edgecol <-  array(NA, dim=c(nNode,nNode,tsteps))
  l_mpar.matrix <- list()
  est_steps <- seq(0, 1, length = tsteps) # times at which estimation takes places
  
  if(pbar==TRUE) {pb <- txtProgressBar(min = 0, max=tsteps, initial=0, char="-", style = 3) } #progress bar
  
  # fit weighted mgm model
  tv_list <- list() # storage
  for(t in 1:tsteps) {
    
    kernel_mean <- est_steps[t]
    
    # define weights with folded gaussian kernel    
    weights <- dnorm(tpoints, mean = kernel_mean, sd=bandwidth) # we use unit time interval
    weights <- weights / max(weights) 

    # estimate using mgm
    tv_list[[t]] <- mgmfit_core(data, type, lev, lambda.sel="EBIC", folds = 10, gam = gam,  # force EBIC, CV makes no sense in nonstationary data
         d=d, rule.reg = rule.reg, pbar = FALSE, method = method, missings = missings, 
         weights=weights, ret.warn=FALSE)
    
    a_wadj[,,t] <- tv_list[[t]]$wadj
    a_signs[,,t] <- tv_list[[t]]$signs
    a_edgecol[,,t] <- tv_list[[t]]$edgecolor
    l_mpar.matrix[[t]] <- tv_list[[t]]$mpar.matrix
    
    if(pbar==TRUE) { setTxtProgressBar(pb, t) } # update progress bar
  
  }
  
  # write model parameter matrix in array
  a_mpar.matrix <- array(NA, dim=c(dim(l_mpar.matrix[[t]])[1], dim(l_mpar.matrix[[t]])[2],tsteps))
  
  # get weighted N for each time step; this will later also be used as a missing data signal
  l_wN <- unlist(lapply(tv_list, function(x) sum(x$call$weights)))
  
  # outputlist
  l_call <- list('type'=type, 'lev'=lev, 'tsteps'=tsteps, 'bandwidth'=bandwidth, 'lambda.sel'="EBIC",
                 'gam'=gam, 'd'=d, 'rule.reg'=rule.reg, 'method'=method)
  outlist <- list('call'=l_call, 'wadj'=a_wadj, 'mpar.matrix'=l_mpar.matrix, 
                  'sign'=a_signs, 'edgecolor'=a_edgecol, 't.models'=tv_list, 'Nt'=l_wN)
  
  class(outlist) <- c('mgm', 'tv.mgm')
  
  return(outlist) 
  
}






