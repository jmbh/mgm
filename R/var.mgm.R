

var.mgm <- function(
  data, # data matrix, col=variables
  type, # data type for col 1:ncol; c=categorical, g=gaussian, p=poisson, e=exponential
  lev, # number of categories of categorical variables, continuous variables have level=1
  lags = 1, # currently limited to 1 
  lambda.sel = "EBIC", # method for penalization parameter (lambda) -selection 
  folds = 10, # folds in case CV is used for lambda selection
  gam = .25, # tuning parameter for EBIC, in case EBIC is used for lambda selection
  d = 2, # maximal degree of the true graph
  pbar = TRUE, # shows a progress bar if TRUE
  method = 'glm',  # which method should be used for each nodewise regression?
  missings = 'error', # handling of missing data
  weights = NA, # weights for observations 
  ret.warn = TRUE, # TRUE returns warnings, makes sense to switch off for time varying wrapper
  binary.sign = FALSE, # see help file
  ...
)

{
  
  # ---------- VAR-specific input checks ----------
  
  if((nrow(data) - max(lags)) != length(weights)) stop('Provide nrow(data) - lags weights. If there are k lags, the first k rows are removed.')
  
  
  # ---------- Call mgmfit core ----------  
  
  outlist <- mgmfit_core(data = data, 
                         type = type, 
                         lev = lev, 
                         lambda.sel = lambda.sel, 
                         folds = folds, 
                         gam = gam, 
                         d = d, 
                         rule.reg = "AND", 
                         pbar = pbar, 
                         method = method, 
                         missings = missings, 
                         weights = weights, 
                         ret.warn = ret.warn, 
                         binary.sign = binary.sign,
                         VAR = TRUE) # use standard mgm.fit; no AR model
  
  
  # ---------- Export ----------
  
  # Return estimation messages:
  estimation_msg('var.mgm') # note about where signs are stored
  
  return(outlist)
} 

