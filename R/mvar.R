


mvar <- function(data,         # n x p data matrix
                 type,         # p vector indicating the type of each variable
                 level,        # p vector indivating the levels of each variable
                 lambdaSeq,    # sequence of considered lambda values (default to glmnet default)
                 lambdaSel,    # way of selecting lambda: CV vs. EBIC
                 lambdaFolds,  # number of folds if lambdaSel = 'CV'
                 lambdaGam,    # EBIC hyperparameter gamma, if lambdaSel = 'EBIC'
                 alphaSeq,     # sequence of considered alpha values (elastic net), default = 1 = lasso
                 alphaSel,     # way of selecting lambda: CV vs. EBIC
                 alphaFolds,   # number of folds if alphaSel = 'CV'
                 alphaGam,     # EBIC hyperparameter gamma, if alphaSel = 'EBIC',
                 lags,
                 consec,
                 weights,      # p vector of observation weights for weighted regression
                 threshold,    # defaults to 'LW', see helpfile
                 method,       # glm vs. 'linear'; for now only implement glm
                 binarySign,   # should a sign be computed for binary models (defaults to NO)
                 scale,        # If TRUE, scales Gaussians
                 verbatim,     # turns off all notifications
                 pbar,         #
                 warnings,     #
                 saveModels,    # defaults to TRUE, saves all estimated models
                 saveData,     # defaults to FALSE, saves the data, =TRUE makes sense for easier prediction routine in predict.mgm()
                 overparameterize,
                 signInfo,
                 ...
)



{
  
  
  # -------------------- Input Checks Global -------------------
  
  # ----- Compute Aux Variables -----
  
  # browser()
  
  n <- nrow(data)
  p <- ncol(data)
  n_var <- n - max(lags)
  n_lags <- length(lags)
  
  args <- list(...)
  
  # ----- Fill in Defaults -----
  
  if(missing(lambdaSeq)) lambdaSeq <- NULL
  if(missing(lambdaSel)) lambdaSel <- 'CV'
  if(missing(lambdaFolds)) lambdaFolds <- 10
  if(missing(lambdaGam)) lambdaGam <- .25
  if(missing(alphaSeq)) alphaSeq <- 1
  if(missing(alphaSel)) alphaSel <- 'CV'
  if(missing(alphaFolds)) alphaFolds <- 10
  if(missing(alphaGam)) alphaGam <- .25
  if(missing(lags)) lags <- 1
  if(missing(consec)) consec <- NULL
  # no AND rule necessary
  if(missing(weights)) weights <- rep(1, n-max(lags))
  if(missing(threshold)) threshold <- 'LW'
  if(missing(method)) method <- 'glm'
  if(missing(binarySign)) binarySign <- FALSE
  if(missing(verbatim)) verbatim <- FALSE
  if(missing(pbar)) pbar <- TRUE
  if(missing(warnings)) warnings <- TRUE
  if(missing(saveModels)) saveModels <- TRUE
  if(missing(saveData)) saveData <- FALSE
  if(missing(overparameterize)) overparameterize <- FALSE
  if(missing(signInfo)) signInfo <- TRUE
  
  if(missing(scale)) scale <- TRUE
  
  if(verbatim) pbar <- FALSE
  if(verbatim) warnings <- FALSE
  
  weights_initial <- weights
  
  
  # ----- Compute Auxilliary Variables II -----
  
  # Empirical Levels of each variable
  emp_lev <- rep(NA, p)
  ind_cat <- which(type=='c')
  if(length(ind_cat)>0) for(i in 1:length(ind_cat)) emp_lev[ind_cat][i] <-  length(unique(data[,ind_cat[i]])) # no apply() because of case of 1 categorical
  emp_lev[which(type!='c')] <- 1
  
  
  if(!missing(level)) {
    # Check whether provided levels are equal to levels in the data
    level_check <- level != emp_lev
    if(sum(level_check) > 0) stop(paste0('Provided levels not equal to levels in data for variables ',paste((1:p)[level_check], collapse = ', ')))
    # if not provided, do nothing, because the argument is not actually necessary
  }
  level <- emp_lev
  
  # Normalize weights (necessary to ensure that nadj makes sense)
  if(!missing(weights)) weights <- weights / max(weights)
  nadj <- sum(weights) # calc adjusted n
  
  # Scale Gaussians
  ind_Gauss <- which(type == 'g')
  if(scale) for(i in ind_Gauss) data[, i] <- scale(data[, i])
  
  # ----- Basic Input Checks -----
  
  if(!is.matrix(data)) stop('The data has to be provided as a n x p matrix (no data.frame)')
  
  if(!(threshold %in% c('none', 'LW', 'HW'))) stop('Please select one of the three threshold options "HW", "LW" and "none" ')
  
  if(nrow(data) < 2) ('The data matrix has to have at least 2 rows.')
  
  if(is.null(lags)) stop('No lags specified')
  if(any(duplicated(lags))) stop('No duplicates allowed in specified lags.')
  
  if(any(lags<1)) stop('Specified lags have to be in {1, 2, ..., n -1}')
  if(any(round(lags)!=lags)) stop('Specified lags have to be positive integers')
  if(missing(data)) stop('No data provided')
  if(ncol(data)<2) stop('At least 2 variables required')
  if(missing(type)) stop('No type vector provided.')
  
  if(sum(!(type %in% c('g', 'c', 'p')))>0) stop("Only Gaussian 'g', Poisson 'p' or categorical 'c' variables permitted.")
  if(any(is.na(data))) stop('No missing values permitted.')
  if(any(!is.finite(as.matrix(data)))) stop('No infinite values permitted.')
  if(any(!(apply(data, 2, class) %in% c('numeric', 'integer')))) stop('Only integer and numeric values permitted.')
  
  if(ncol(data) != length(type)) stop('Number of variables is not equal to length of type vector.')
  if(!missing(level)) if(ncol(data) != length(level)) stop('Number of variables is not equal to length of level vector.')
  if((nrow(data)-max(lags)) != length(weights)) stop('Weights vector has to be equal to the number of observations n - max(lags)')
  
  # Are Poisson variables integers?
  if('p' %in% type) {
    ind_Pois <- which(type == 'p')
    nPois <- length(ind_Pois)
    v_PoisCheck <- rep(NA, length=nPois)
    for(i in 1:nPois) v_PoisCheck[i] <- sum(data[, ind_Pois[i]] != round(data[, ind_Pois[i]])) > 0
    if(sum(v_PoisCheck) > 0) stop('Only integers permitted for Poisson variables.')
  }
  
  # ----- Binary Sign => values have to be in {0,1} -----
  # (compute anyway, because used later for sign extraction)
  
  # Find the binary variables
  ind_cat <- which(type == 'c')
  ind_binary <- rep(NA, length(ind_cat))
  ind_binary <- as.logical(ind_binary)
  if(length(ind_cat)>0) {
    for(i in 1:length(ind_cat)) ind_binary[i] <- length(unique(data[, ind_cat[i]])) == 2
  }
  
  # Check if they are coded in {0,1}
  if(sum(ind_binary)>0){
    check_binary <- rep(NA, sum(ind_binary))
    for(i in 1:sum(ind_binary)) check_binary[i] <- sum(!(unique(data[, ind_cat[ind_binary][i]]) %in% c(0,1)))
  }
  if(binarySign) {
    if(sum(check_binary)>0) stop(paste0('If binarySign = TRUE, all binary variables have to be coded {0,1}. Not satisfied in variable(s) ',paste(ind_cat[ind_binary][check_binary>0], collapse = ', ')))
  }
  
  
  
  # -------------------- Create VAR data structure -------------------
  
  # ----- Give Names to Variables & create data.frame -----
  
  colnames(data)[1:p] <- paste("V", 1:p, '.', sep = "")
  data <- as.data.frame(data)
  
  # Categoricals into factors (Needed to use formula to construct design matrix)
  for(i in which(type=='c')) data[, i] <- as.factor(data[, i])
  
  
  # ----- Split up predictor Sets by lags -----
  
  # Divide Data in several parts: one response set, and one set of each lag
  data_lagged <- lagData(data = data, 
                         lags = lags, 
                         consec = consec)
  
  data_response <- data_lagged$data_response
  l_data_lags <- data_lagged$l_data_lags
  
  n_design <- nrow(data_response)
  
  # # set weights for cases with not enough preceding complete measurements to zero
  # if(!is.null(consec)) {
  #   weights[data_lagged$included == FALSE] <- 0
  #   nadj <- sum(weights) # calc adjusted n again for case of non-consecutive measurements
  # }
  
  # do different, to make below bootstrap scheme simpler;
  # Subset instead:
  ind_included_wo_begin <- data_lagged$included[-c(1:n_lags)] # the weights-vector has alread length nrow-max(lags)
  
  data_response <- data_response[ind_included_wo_begin, ]
  l_data_lags <- lapply(l_data_lags, function(x) x[ind_included_wo_begin, ])
  weights <- weights[ind_included_wo_begin]
  nadj <- sum(weights)
  
  # browser()
  
  # ----- Use bootstrap instead of original data (called from resample()) -----

  if(!is.null(args$bootstrap)) {
    if(args$bootstrap) {

      # overwrite data with bootstrap sample, passed on from resample()
      data_response <- data_response[args$boot_ind, ]
      l_data_lags <- lapply(l_data_lags, function(x) x[args$boot_ind, ])

      # overwrite weights
      weights <- weights[args$boot_ind]
      
    }
  }

  
  # -------------------- Input Checks Local (for each set of predictors) -------------------
  
  # This checks glmnet requirements wrt categorical predictors, for each lag set
  
  n_lags <- length(lags)
  
  for(lag in 1:n_lags) {
    
    glmnetRequirements(data = l_data_lags[[lag]],
                       type = type,
                       weights = weights)
    
  }
  
  
  # ----- Storage: Create empty mgm object -----
  
  mvarobj <- list('call' = NULL, # fill in later
                  'wadj' = NULL,
                  'signs' = NULL,
                  'edgecolor' = NULL,
                  'rawlags' = list(),
                  'intercepts' = NULL,
                  'nodemodels' = NULL)
  
  
  
  # ----- Save the Call -----
  
  mvarobj$call <- list('data' = NULL,
                       'data_lagged' = NULL,
                       'type' = type,
                       'level' = level,
                       'lambdaSeq' = lambdaSeq,
                       'lambdaSel' = lambdaSel,
                       'lambdaFolds' = lambdaFolds,
                       'lambdaGam' = lambdaGam,
                       'alphaSeq' = alphaSeq,
                       'alphaSel' = alphaSel,
                       'alphaFolds' = alphaFolds,
                       'alphaGam' = alphaGam,
                       'lags' = lags,
                       'weights' = weights_initial,
                       'weights_design' = weights,
                       'threshold' = threshold,
                       'method' = method,
                       'binarySign' = binarySign,
                       'scale' = scale,
                       'verbatim' = verbatim,
                       'pbar' = pbar,
                       'warnings' = warnings,
                       'saveModels' = saveModels,
                       'saveData' = saveData,
                       'overparameterize' = overparameterize,
                       'signInfo' = signInfo)
  
  if(saveData)
  {
    mvarobj$call$data <- data
    mvarobj$call$data_lagged <- data_lagged
  }
  
  
  
  
  # -------------------- Estimate -------------------
  
  # Progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=p, initial=0, char="-", style = 3)
  
  
  # Save number of parameters of standard (non-overparameterized) design matrices for tau threshold
  npar_standard <- rep(NA, p)
  
  for(v in 1:p) {
    
    
    # ----- Create VAR Design Matrix -----
    
    # append response with predictors
    y <- data_response[,v] # response variable v
    data_v <- cbind(y, do.call(cbind, l_data_lags)) # combine
    data <- data_v[, -1]
    
    # Dummy coding
    form <- as.formula('y ~ (.)')
    
    # Construct standard design matrix (to get number of parameters for tau threshold)
    X_standard <- model.matrix(form, data = data_v)[, -1] # delete intercept (added by glmnet later)
    npar_standard[v] <- ncol(X_standard)
    
    if(overparameterize) {
      
      # Compute augmented type and level vectors
      type_aug <- rep(type, n_lags)
      level_aug <- rep(level, n_lags)
      
      # Construct over-parameterized design matrix
      X_over <- ModelMatrix(data = data,
                            type = type_aug,
                            level = level_aug,
                            labels = colnames(data),
                            d = 1)
      X <- X_over
      
    } else {
      
      X <- X_standard
      
    }
    

    
    
    # ----- Fit each neighborhood using the same procedure as in mgm -----
    
    # ----- Tuning Parameter Selection (lambda and alpha) -----
    
    n_alpha <- length(alphaSeq) # length of alpha sequence
    
    # alpha Section via CV
    
    if(alphaSel == 'CV') {
      
      l_alphaModels <- list() # Storage
      ind <- sample(1:alphaFolds, size = n_design, replace = TRUE) # fold-indicators, use same for each alpha
      
      v_mean_OOS_deviance <- rep(NA, n_alpha)
      
      if(n_alpha>1) {
        
        
        # For: alpha
        for(a in 1:n_alpha) {
          
          l_foldmodels <- list()
          v_OOS_deviance <- rep(NA, alphaFolds)
          
          for(fold in 1:alphaFolds) {
            
            # Select training and test sets
            train_X <- X[ind != fold, ]
            train_y <- y[ind != fold]
            test_X <- X[ind == fold, ]
            test_y <- y[ind == fold]
            
            # Recompute variables for training set
            n_train <- nrow(train_X)
            nadj_train <- sum(weights[ind != fold])
            
            l_foldmodels[[fold]] <- nodeEst(y = train_y,
                                            X = train_X,
                                            lambdaSeq = lambdaSeq,
                                            lambdaSel = lambdaSel,
                                            lambdaFolds = lambdaFolds,
                                            lambdaGam = lambdaGam,
                                            alpha = alphaSeq[a],
                                            weights = weights[ind != fold],
                                            n = n_train,
                                            nadj = nadj_train,
                                            v = v,
                                            type = type,
                                            level = level,
                                            emp_lev = emp_lev,
                                            overparameterize = overparameterize)
            
            # Calculte Out-of-sample deviance for current fold
            LL_model <- calcLL(X = test_X,
                               y = test_y,
                               fit = l_foldmodels[[fold]]$fitobj,
                               type = type,
                               level = level,
                               v = v,
                               weights = weights[ind == fold],
                               lambda = l_foldmodels[[fold]]$lambda,
                               LLtype = 'model')
            
            LL_saturated <- calcLL(X = test_X,
                                   y = test_y,
                                   fit = l_foldmodels[[fold]]$fitobj,
                                   type = type,
                                   level = level,
                                   v = v,
                                   weights = weights[ind == fold],
                                   lambda = l_foldmodels[[fold]]$lambda,
                                   LLtype = 'saturated')
            
            v_OOS_deviance[fold] <- 2 * (LL_saturated - LL_model)
            
          }
          
          v_mean_OOS_deviance[a] <- mean(v_OOS_deviance)
          
        }
        
        alpha_select <- alphaSeq[which.min(v_mean_OOS_deviance)]
        
      } else {
        
        alpha_select <- alphaSeq # in case alpha is just specified
        
      }
      
      # Refit Model on whole data, with selected alpha
      
      # browser()
      
      model <- nodeEst(y = y,
                       X = X,
                       lambdaSeq = lambdaSeq,
                       lambdaSel = lambdaSel,
                       lambdaFolds = lambdaFolds,
                       lambdaGam = lambdaGam,
                       alpha = alpha_select,
                       weights = weights,
                       n = n_design,
                       nadj = nadj,
                       v = v,
                       type = type,
                       level = level,
                       emp_lev = emp_lev,
                       overparameterize = overparameterize)
      
      mvarobj$nodemodels[[v]] <- model
      
    }
    
    
    
    # alpha Section via EBIC
    
    if(alphaSel == 'EBIC') {
      
      l_alphaModels <- list()
      EBIC_Seq <- rep(NA, n_alpha)
      
      # For: alpha
      for(a in 1:n_alpha) {
        
        
        l_alphaModels[[a]] <- nodeEst(y = y,
                                      X = X,
                                      lambdaSeq = lambdaSeq,
                                      lambdaSel = lambdaSel,
                                      lambdaFolds = lambdaFolds,
                                      lambdaGam = lambdaGam,
                                      alpha = alphaSeq[a],
                                      weights = weights,
                                      n = n,
                                      nadj = nadj,
                                      v = v,
                                      type = type,
                                      level = level,
                                      emp_lev = emp_lev,
                                      overparameterize = overparameterize)
        
        EBIC_Seq[a] <- l_alphaModels[[a]]$EBIC
        
      }
      
      ind_minEBIC_model <- which.min(EBIC_Seq)
      mvarobj$nodemodels[[v]] <- l_alphaModels[[ind_minEBIC_model]]
      
      
    } # end if: alpha EBIC?
    
    
    # Update Progress Bar
    if(pbar==TRUE) setTxtProgressBar(pb, v)
    
    
  } # end for: p
  
  # -------------------- Processing glmnet Output -------------------
  
  # ----- For now: only extract lag 1 -----
  # (but the structure in lists and variable names allows easy extraction of higher lags)
  
  
  # Storage
  no_lags <- length(lags)
  l_Par <- vector('list', length = p)
  l_preds_dummy <- vector('list', p*no_lags) # number of predictors (on variable level, not design matrix)
  l_Par <- lapply(l_Par, function(x) l_preds_dummy)
  l_intercepts <- vector('list', length = p) # collect intercepts
  
  pred_names <- unlist(lapply(l_data_lags, colnames))
  
  # Sanity
  if(length(pred_names) != (no_lags * p)) stop('Length of predictor variable names does not match calculated number of predictor varibales')
  
  
  
  
  # Fill parameter list
  
  for(v in 1:p) {
    for(v2 in 1:(p*no_lags)) {
      
      if(type[v] == 'c') {
        
        n_cats <- level[v]
        for(cat in 1:n_cats) {
          
          # all model parameter without intercept
          par_part <- mvarobj$nodemodels[[v]]$model[[cat]]
          par_part_ni <- par_part[-1]
          l_intercepts[[v]][[cat]] <- par_part[1]
          
          # compute threshold
          d <- 1 # at the moment no higher order interactions implemented, so max neighborhood degree = 1
          if(threshold == 'LW') tau <- sqrt(d) * sqrt(sum(par_part_ni^2)) * sqrt(log(npar_standard[v]) / n)
          if(threshold == 'HW') tau <- d * sqrt(log(npar_standard[v]) / n)
          if(threshold == 'none') tau <- 0
          
          # threshold
          par_part[abs(par_part) < tau] <- 0
          mvarobj$nodemodels[[v]]$tau <- tau # Save tau threshold
          
          ind_v2 <- grepl(pred_names[v2], rownames(par_part)[-1], fixed = TRUE) # indicator: where is that parameter; minus 1 to remove intercept
          
          l_Par[[v]][[v2]][[cat]] <- par_part[-1,][ind_v2] # minus 1 to remove intercept
          
        }
        
        
      } else {
        
        # Select set of parameters of category cat
        par_part <- mvarobj$nodemodels[[v]]$model
        par_part_ni <- par_part[-1]
        l_intercepts[[v]] <- par_part[1]
        
        # compute threshold
        d <- 1 # at the moment no higher order interactions implemented, so max neighborhood degree = 1
        if(threshold == 'LW') tau <- sqrt(d) * sqrt(sum(par_part_ni^2)) * sqrt(log(npar_standard[v]) / n_var)
        if(threshold == 'HW') tau <- d * sqrt(log(npar_standard[v]) / n_var)
        if(threshold == 'none') tau <- 0
        
        # threshold
        par_part[abs(par_part) < tau] <- 0
        mvarobj$nodemodels[[v]]$tau <- tau # Save tau threshold
        
        ind_v2 <- grepl(pred_names[v2], rownames(par_part)[-1], fixed = TRUE) # indicator: where is that parameter
        
        l_Par[[v]][[v2]] <- par_part[-1,][ind_v2] # minus 1 to remove intercept
        
      }
      
    } # end for: v2
  } # end for: v
  
  
  # ----- Recover signs where defined and reduce -----
  
  # Storage
  l_Par_red <- vector('list', length = p)
  l_signs <- vector('list', length = p)
  l_preds_dummy <- vector('list', p*no_lags)
  l_Par_red <- lapply(l_Par_red, function(x) l_preds_dummy)
  l_signs <- lapply(l_signs, function(x) l_preds_dummy)
  
  # Define set of variables which allow sign
  type_long <- rep(type, times=no_lags)
  ind_binary_long <- rep(1:p %in% ind_cat[ind_binary], times=no_lags)
  ind_cat_long <- rep(ind_cat, times=no_lags)
  
  # Define set of continous and binary variables: for interactions between these we can assign a sign
  # Depends on binarySign
  if(binarySign) {
    set_signdefined <- c(which(type_long == 'p'), which(type_long == 'g'), ind_cat_long[ind_binary_long])
  } else {
    set_signdefined <- c(which(type_long == 'p'), which(type_long == 'g'))
  }
  
  
  # Loop again:
  for(v in 1:p) {
    for(v2 in 1:(p*no_lags)) {
      
      # Reduce
      l_Par_red[[v]][[v2]] <- mean(abs(unlist(l_Par[[v]][[v2]])))
      
      # Extract Sign
      if(l_Par_red[[v]][[v2]] != 0) { # if zero parameter, set to NA
        if(sum(!(c(v, v2) %in% set_signdefined))==0) { # if both variables binary/continuous, consider sign, otherwise set to 0
          
          
          # Much easier as in mgm, because we don't have to integrate across two estimates
          if(overparameterize) {
            
            if(type[v]=='c') {
              
              if(type[v2]=='c') {
                if(l_Par[[v]][[v2]][[1]][1] != 0) sign_sel <- sign(l_Par[[v]][[v2]][[1]][1])
                if(l_Par[[v]][[v2]][[1]][2] != 0) sign_sel <- - sign(l_Par[[v]][[v2]][[1]][2])
              } else {
                sign_sel <- sign(l_Par[[v]][[v2]][2])
              }
              
            } else {
              sign_sel <- sign(l_Par[[v]][[v2]])
            }
            
            
          } else {
            
            # Sign extraction for standard parameterization
            if(type[v]=='c') {
              sign_sel <- sign(l_Par[[v]][[v2]][[2]]) # there will be two entries, second one is for category = 1 (other=0)
            } else {
              sign_sel <- sign(l_Par[[v]][[v2]])
            }
            
          }
          
          l_signs[[v]][[v2]] <- sign_sel
          
        } else {
          
          l_signs[[v]][[v2]] <- 0
          
        }
        
      } else {
        
        l_signs[[v]][[v2]] <- NA
        
      } # end if: nonzero coefficient?
      
    } # end for: v2
  } # end for: v
  
  
  
  
  # ----- Split up coefficients of each lag and save in all sorts of ways -----
  
  # Create lag indicator
  n_lags <- length(lags)
  lag_indicator <- rep(1:n_lags, each = p)
  
  # Storage
  list_wadj <- list_signs <- list()
  
  a_wadj <- a_signs <- array(dim = c(p, p, n_lags))
  a_edgecolor <- array('darkgrey', dim = c(p, p, n_lags))
  mvarobj$rawlags <- vector('list', length = n_lags)
  
  for(lag in 1:n_lags) {
    
    # weighted version
    lag_seq <- (1:(p*n_lags))[lag_indicator==lag]
    
    # Storage for raw Parameter output list
    l_lag <- vector('list', length = p)
    l_dum <- vector('list', length = length(lag_seq))
    l_lag <- lapply(l_lag, function(x) l_dum)
    
    wadj <- m_signs <- matrix(NA, p, p)
    for(i in 1:p) {
      k <- 1
      for(j in lag_seq) {
        a_wadj[i, k, lag] <- l_Par_red[[i]][[j]]
        l_lag[[i]][[j]] <- l_Par[[i]][[j]]
        a_signs[i, k, lag] <- l_signs[[i]][[j]]
        
        # edge color array
        if(!is.na(a_signs[i, k, lag])) {
          if(a_signs[i, k, lag] == -1) a_edgecolor[i, k, lag] <- 'red'
          if(a_signs[i, k, lag] == 1) a_edgecolor[i, k, lag] <- 'darkgreen'
        }
        
        k <- k + 1
        
      }
    }
    
    # save in list
    mvarobj$wadj <- a_wadj
    mvarobj$signs <- a_signs
    mvarobj$edgecolor <- a_edgecolor
    mvarobj$rawlags[[lag]] <- l_lag
    
  } # end for: lag
  
  # browser()
  
  
  # -------------------- Output -------------------
  
  # Save Node Models and extracted raw factors?
  if(!saveModels) {
    mvarobj$nodemodels <- NULL
    mvarobj$rawfactor <- NULL
  }
  
  # Save intercepts
  mvarobj$intercepts <- l_intercepts
  
  if(pbar) {
    if(signInfo) cat('\nNote that the sign of parameter estimates is stored separately; see ?mvar')    
  } else {
    if(signInfo) cat('Note that the sign of parameter estimates is stored separately; see ?mvar')    
  }
  
  
  class(mvarobj) <- c('mgm', 'mvar')
  
  
  return(mvarobj)
  
}
