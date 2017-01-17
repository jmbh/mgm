

bwSelect <- function(data, # Data
                     train, # indicator vector: training cases
                     bwSeq, # sequence of bandwidth values to be tested
                     VAR = FALSE, # indicator value: VAR-1 model or 'contemporaneous' model 
                     ... # arguments passed to tv.mgmfit / tv_var.mgm
                     )
  
{
  
  
  # ----- Input Checks -----
  
  # Actual Imput Checks
  if(arg_pass$method == 'linear') stop('bwSelect is only supported for method == \'glm\' ') 
  
  arg_pass <- list(...)
  arg_pass <- list()

  # Necessary Arguments  
  if(is.null(arg_pass$lev)) stop('Please specify the level vector. See tv.mgmfit()')
  if(is.null(arg_pass$type)) stop('Please specify the type vector. See tv.mgmfit()')

  # Fill in defaults
  if(is.null(arg_pass$timepoints)) arg_pass$timepoints <- NA
  if(is.null(arg_pass$gam)) arg_pass$gam <- .25
  if(is.null(arg_pass$d)) arg_pass$d <- 1
  if(is.null(arg_pass$rule.reg)) arg_pass$rule.reg <- 'AND'
  if(is.null(arg_pass$missings)) arg_pass$missings <- 'error'
  
  # Fix some inputs
  arg_pass$pbar <- FALSE
  arg_pass$ret.warn <- FALSE
  
  
  
  # ----- Compute Auxilliary Variables -----  

  # 
  # # Compute the times of
  # if(is.na(arg_pass$timepoints)) {
  #   train_points <- train
  # } else {
  #   timepoints_norm <- timepoints - min(timepoints) 
  #   timepoints_norm <- timepoints_norm / max(timepoints_norm)
  #   train_points <- timepoints_norm[train]
  # }
  # 
  # 
  # 

  # ----- Search bw-path -----  
  
  bwSeq <- seq(0.001, .4, length = 10)
  n_bw <- length(bwSeq)
  
  l_bwModels <- l_bwPreds <- list() # Storage

  
  for(bw in 1:n_bw) {
    
    
    # Contemporanous Models
    if(!VAR) {
      
      l_bwModels[[bw]] <- tv.mgmfit(data = data,
                                    type = arg_pass$type, 
                                    lev = arg_pass$lev, 
                                    timepoints = )
      
    
        
    # VAR Model
    } else {
      
      
      
    }
    
    
    
    
    
    
  }
  
  
  
  
  
  # ----- Output -----  
  
  
  # loss for each node, for each training time point, for each bandwidth
  # and summaries thereof
  
  
  
  
}

