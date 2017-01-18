

bwSelect <- function(data, # Data
                     test, # indicator vector: training cases
                     bwSeq, # sequence of bandwidth values to be tested
                     VAR = FALSE, # indicator value: VAR-1 model or 'contemporaneous' model 
                     ... # arguments passed to tv.mgmfit / tv_var.mgm and other stuff
                     )
  
{
  
  
  # ----- Input Checks -----
  
  arg_pass <- list(...)

  # Necessary Arguments  
  if(is.null(arg_pass$lev)) stop('Please specify the level vector. See tv.mgmfit()')
  if(is.null(arg_pass$type)) stop('Please specify the type vector. See tv.mgmfit()')

  # Fill in defaults if not specified
  if(is.null(arg_pass$timepoints)) arg_pass$timepoints <- NA
  if(is.null(arg_pass$gam)) arg_pass$gam <- .25
  if(is.null(arg_pass$d)) arg_pass$d <- 1
  if(is.null(arg_pass$rule.reg)) arg_pass$rule.reg <- 'AND'
  if(is.null(arg_pass$pbar)) arg_pass$pbar <- TRUE
  if(is.null(arg_pass$method)) arg_pass$method <- 'glm'
  
  
  # Actual Imput Checks
  if(arg_pass$method == 'linear') stop('bwSelect is only supported for method == \'glm\' ') 

  # Later: Build in: Error in case of missing data! this would fuck up the train/test set division
  
  # Fix some inputs
  arg_pass$ret.warn <- FALSE

  
  # ----- Compute Auxilliary Variables -----  

  # estpoints has different input depending on timepoints input
  if(is.na(arg_pass$timepoints)) {
    estpoints <- test  
  } else {
    estpoints <- arg_pass$timepoints[test]
  }


  # ----- Search bw-path -----  
  
  # Storage
  n_bw <- length(bwSeq)
  l_bwModels <- l_bwPreds <- list()
  
  data_train <- data
  data_train[test,] <- NA # knock out test samples (will be weighted to zero by NA-handling)
  
  # Set Pbar
  if(arg_pass$pbar) pb <- txtProgressBar(min = 0, max=n_bw, initial=0, char="-", style = 3)
  
  for(bw in 1:n_bw) {
    
    # add pbar 
    
    
    
    # Contemporanous Models
    if(!VAR) {
      
      l_bwModels[[bw]] <- tv.mgmfit(data = data_train,
                                    type = arg_pass$type, 
                                    lev = arg_pass$lev, 
                                    timepoints = arg_pass$timepoints,
                                    estpoints = estpoints, 
                                    bandwidth = bwSeq[bw], 
                                    gam = arg_pass$gam,  
                                    d = arg_pass$d, 
                                    rule.reg = arg_pass$rule.reg, 
                                    missings = 'casewise.zw', 
                                    pbar = FALSE, 
                                    ret.warn = arg_pass$ret.warn,
                                    signWarn = FALSE)
      
      l_bwPreds[[bw]] <- tvPredict(data = data,
                                   type = type,
                                   test = test, 
                                   model = l_bwModels[[bw]])
      
    
        
    # VAR Model
    } else {
      
      
      
    }
    
    
    if(arg_pass$pbar)  setTxtProgressBar(pb, bw) 
    
  } # end of for: bw
  
  
  
  # ----- Evaluate -----  
  
  
  loss_bw <- unlist(lapply(l_bwPreds, function(x) x$perMean))
  
 # plot(bwSeq, loss_bw)

  
  
  
  # ----- Output -----  
  
  outlist <- list('PredAll' = l_bwPreds, 
                  'PredBw' = loss_bw,
                  'bwSeq' = bwSeq,
                  'Models' = l_bwModels)
  
  return(outlist)
  
  
}

