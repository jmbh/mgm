

tv_var.mgm <- function(data, # input tv mgm
                      type, 
                      lev,
                      lags = 1, 
                      tsteps = 'default', 
                      bandwidth,
                      gam = .25,                       
                      d = 2, 
                      pbar = TRUE,
                      method ='glm',
                      missings = 'error',
                      ret.warn = FALSE, 
                      ... )

{
  
  
  # ---------- Compute Aux Variables ----------
  
  n <- nrow(data)
  nNode <- ncol(data)
  
  arg_pass <- list(...)
  if(is.null(arg_pass$signWarn)) arg_pass$signWarn <- TRUE
  
  
  # ---------- Assign each data point a value on [0,1] normed time interval ----------
  
  # Note: this part is much easier than in the t-v MGM model, because the VAR model assumes equally spaced measurements
  
  tpoints <- seq(0, 1, length = (n - 1)) # n sequence before and after for folded gaussian weight

  if(tsteps=='default') {tsteps <- round(2/bandwidth)} # reasonable coveres the data for any choice of bandwidth
  est_steps <- seq(0, 1, length = tsteps) # times at which estimation takes places
  
  # Storage
  a_wadj <- a_signs <- a_edgecol <-  array(NA, dim=c(nNode,nNode,tsteps))
  l_mpar.matrix <- list()
  tv_list <- list() 
  
  # Progress bar
  if(pbar==TRUE) {pb <- txtProgressBar(min = 0, max=tsteps, initial=0, char="-", style = 3) }
  
  
  # ---------- Estimate Time Points at specified locations ----------
  
  for(t in 1:tsteps) {
    
    kernel_mean <- est_steps[t]
    
    # define weights with folded gaussian kernel    
    weights <- dnorm(tpoints, mean = kernel_mean, sd = bandwidth) # we use unit time interval
    weights <- weights / max(weights) 
    
    # estimate using mgm
    tv_list[[t]] <- mgmfit_core(data = data, 
                                type = type, 
                                lev = lev, 
                                lambda.sel = "EBIC", 
                                folds = 10, 
                                gam = gam,  # force EBIC, CV makes no sense in nonstationary data
                                d = d, 
                                rule.reg = 'AND', 
                                pbar = FALSE, 
                                method = method, 
                                missings = missings, 
                                weights = weights, 
                                ret.warn = FALSE, 
                                VAR = TRUE)
    
    
    a_wadj[,,t] <- tv_list[[t]]$wadj
    a_signs[,,t] <- tv_list[[t]]$signs
    a_edgecol[,,t] <- tv_list[[t]]$edgecolor
    l_mpar.matrix[[t]] <- tv_list[[t]]$mpar.matrix
    
    if(pbar==TRUE) { setTxtProgressBar(pb, t) } # update progress bar
    
  }
  
  
  # ---------- Prepare output  ----------
  
  # write model parameter matrix in array
  a_mpar.matrix <- array(NA, dim=c(dim(l_mpar.matrix[[t]])[1], dim(l_mpar.matrix[[t]])[2],tsteps))
  
  # get weighted N for each time step; this will later also be used as a missing data signal
  l_wN <- unlist(lapply(tv_list, function(x) sum(x$call$weights)))
  
  # outputlist
  l_call <- list('type' = type, 
                 'lev' = lev, 
                 'tsteps' = tsteps, 
                 'bandwidth' = bandwidth, 
                 'lambda.sel' = "EBIC",
                 'gam' = gam, 
                 'd' = d, 
                 'method'=method)
  
  
  outlist <- list('call' = l_call, 
                  'wadj' = a_wadj, 
                  'mpar.matrix' = l_mpar.matrix, 
                  'sign' = a_signs, 
                  'edgecolor' = a_edgecol, 
                  't.models' = tv_list, 
                  'Nt' = l_wN)
  
  class(outlist) <- c('mgm', 'var', 'tv.var')
  
  
  # Return estimation messages:
  if(arg_pass$signWarn) estimation_msg('tv_var.mgm') # note about where signs are stored
  
  return(outlist) 
  
}
