

mgm <- function(data,         # n x p data matrix
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
                k,            # order of modeled interactions, 1 = pairwise
                ruleReg,      # rule to combine d+1 neighborhood estimates (defaults to 'AND')
                weights,      # p vector of observation weights for weighted regression
                threshold,    # defaults to 'LW', see helpfile
                method,       # glm vs. 'linear'; for now only implement glm
                binarySign,   # should a sign be computed for binary models (defaults to NO)
                scale,        # should gaussian variables be scaled? defaults to TRYE
                verbatim,     # turns off all notifications
                pbar,         #
                warnings,     #
                saveModels,   # defaults to TRUE, saves all estimated models
                saveData,     # defaults to FALSE, saves the data, =TRUE makes sense for easier prediction routine in predict.mgm()
                overparameterize, # if TRUE, uses the over-parameterized version,
                signInfo,
                ...
)


{
  
  # -------------------- Input Checks -------------------
  
  # ----- Compute Auxilliary Variables I -----
  
  p <- ncol(data)
  n <- nrow(data)
  
  # Give Names to variables (Needed to use formula to construct design matrix and to give informative error messages)
  colnames(data)[1:p] <- paste("V", 1:p, '.', sep = "")
  # colnames(data)[type == 'c'] <- paste("V", (1:p)[type == 'c'], '.', sep = "") # dot deliminator for categorical variables; needed to identify parameters assiciated with some k-order interaction below
  
  # Catch other passed on arguments
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
  if(missing(k)) k <- 2
  if(missing(ruleReg)) ruleReg <- 'AND'
  if(missing(weights)) weights <- rep(1, n)
  if(missing(threshold)) threshold <- 'LW'
  if(missing(method)) method <- 'glm'
  if(missing(binarySign)) {
    if(!is.null(args$binary.sign)) binarySign <- args$binary.sign else binarySign <- FALSE
  }
  
  if(!is.null(args$binary.sign)) {
    warning("The argument 'binary.sign' is deprecated Use 'binarySign' instead.")
  }
  

  if(missing(scale)) scale <- TRUE
  if(missing(verbatim)) verbatim <- FALSE
  if(missing(pbar)) pbar <- TRUE
  if(missing(warnings)) warnings <- TRUE
  if(missing(saveModels)) saveModels <- TRUE
  if(missing(saveData)) saveData <- FALSE
  if(missing(overparameterize)) overparameterize <- FALSE
  if(missing(signInfo)) signInfo <- TRUE
  
  
  if(verbatim) pbar <- FALSE
  if(verbatim) warnings <- FALSE
  
  # Switch all warnings off
  if(!warnings) {
    oldw <- getOption("warn")
    options(warn = -1)
  }
  
  
  
  # ----- Compute Auxilliary Variables II -----
  
  
  # From k to d
  d <- k - 1 # k = largest order of interaction in joint model; d = largest neighborhood size
  
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
  
  
  # ----- Basic Checks -----
  
  if(!is.matrix(data)) stop('The data has to be provided as a n x p matrix (no data.frame)')
  
  if(!(threshold %in% c('none', 'LW', 'HW'))) stop('Please select one of the three threshold options "HW", "LW" and "none" ')
  
  if(nrow(data) < 2) ('The data matrix has to have at least 2 rows.')
  
  if(missing(data)) stop('No data provided.')
  if(k<2) stop('The order of interactions should be at least k = 2 (pairwise interactions)')
  if(ncol(data)<3) stop('At least 3 variables required')
  if(missing(type)) stop('No type vector provided.')
  
  if(sum(!(type %in% c('g', 'c', 'p')))>0) stop("Only Gaussian 'g', Poisson 'p' or categorical 'c' variables permitted.")
  if(any(is.na(data))) stop('No missing values permitted.')
  if(any(!is.finite(as.matrix(data)))) stop('No infinite values permitted.')
  
  
  if(any(!(apply(data, 2, class) %in% c('numeric', 'integer')))) stop('Only integer and numeric values permitted.')
  
  if(ncol(data) != length(type)) stop('Number of variables is not equal to length of type vector.')
  if(!missing(level)) if(ncol(data) != length(level)) stop('Number of variables is not equal to length of level vector.')
  if(nrow(data) != length(weights)) stop('Number of observations is not equal to length of weights vector.')
  
  # Are Poisson variables integers?
  if('p' %in% type) {
    ind_Pois <- which(type == 'p')
    nPois <- length(ind_Pois)
    v_PoisCheck <- rep(NA, length=nPois)
    for(i in 1:nPois) v_PoisCheck[i] <- sum(data[, ind_Pois[i]] != round(data[, ind_Pois[i]])) > 0
    if(sum(v_PoisCheck) > 0) stop('Only integers permitted for Poisson variables.')
  }
  
  # ----- Checking glmnet minimum Variance requirements -----
  
  glmnetRequirements(data = data,
                     type = type,
                     weights = weights)
  
  
  
  # ----- Binary Sign => values have to be in {0,1} -----
  
  # compute anyway, because used later for sign extraction
  
  # browser()
  
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
    
    if(binarySign) {
      if(sum(check_binary)>0) stop(paste0('If binarySign = TRUE, all binary variables have to be coded {0,1}. Not satisfied in variable(s) ',paste(ind_cat[ind_binary][check_binary>0], collapse = ', ')))
    }
  }
  
  # ----- Storage: Create empty mgm object -----
  
  mgmobj <- list('call' = NULL,
                 'pairwise' = list('wadj' = NULL,
                                   'signs' = NULL,
                                   'edgecolor'= NULL),
                 'factorgraph' = list('graph' = NULL,
                                      'signs' = NULL,
                                      'edgecolor' = NULL,
                                      'order' = NULL),
                 'rawfactor' = list('indicator' = NULL,
                                    'weightsAgg' = NULL,
                                    'weights' = NULL,
                                    'signs' = NULL),
                 'intercepts' = NULL,
                 'nodemodels' = list())
  
  
  # ----- Save the Call -----
  
  
  mgmobj$call <- list('data' = NULL,
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
                      'k' = k,
                      'ruleReg' = ruleReg,
                      'weights' = weights,
                      'threshold' = threshold,
                      'method' = method,
                      'binarySign' = binarySign,
                      'scale' = scale,
                      'verbatim' = verbatim,
                      'pbar' = pbar,
                      'warnings' = warnings,
                      'saveModels' = saveModels,
                      'saveData' = saveData,
                      'overparameterize' = overparameterize)
  
  if(saveData) mgmobj$call$data <- data
  
  
  # ----- Some more variable Transforms -----
  
  data <- as.data.frame(data)
  
  # Categoricals into factors (Needed to use formula to construct design matrix)
  for(i in which(type=='c')) data[, i] <- as.factor(data[, i])
  
  # -------------------- Estimation -------------------
  
  # Progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=p, initial=0, char="-", style = 3)
  
  # Save number of parameters of standard (non-overparameterized) design matrices for tau threshold
  npar_standard <- rep(NA, p)
  
  for(v in 1:p) {
    
    # ----- Construct Design Matrix -----
    
    if(d > (p - 1)) {
      stop("Order of interactions cannot be larger than the number of predictors.")
    } else if (d == 1){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
    } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^",d)) }
    
    # Construct standard design matrix (to get number of parameters for tau threshold)
    X_standard <- model.matrix(form, data = data)[, -1] # delete intercept (added by glmnet later)
    npar_standard[v] <- ncol(X_standard)
    
    
    if(overparameterize) {
      
      # Construct over-parameterized design matrix
      X_over <- ModelMatrix(data = data[, -v],
                            type = type[-v],
                            level = level[-v],
                            labels = colnames(data)[-v],
                            d = d)
      
      X <- X_over
    } else {
      X <- X_standard
    }
    
    y <- as.numeric(data[, v])
    
    
    
    # ----- Tuning Parameter Selection (lambda and alpha) -----
    
    n_alpha <- length(alphaSeq) # length of alpha sequence
    
    # alpha Section via CV
    
    if(alphaSel == 'CV') {
      
      l_alphaModels <- list() # Storage
      ind <- sample(1:alphaFolds, size = n, replace = TRUE) # fold-indicators, use same for each alpha
      
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
      
      model <- nodeEst(y = y,
                       X = X,
                       lambdaSeq = lambdaSeq,
                       lambdaSel = lambdaSel,
                       lambdaFolds = lambdaFolds,
                       lambdaGam = lambdaGam,
                       alpha = alpha_select,
                       weights = weights,
                       n = n,
                       nadj = nadj,
                       v = v,
                       type = type,
                       level= level,
                       emp_lev = emp_lev,
                       overparameterize = overparameterize)
      
      mgmobj$nodemodels[[v]] <- model
      
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
      mgmobj$nodemodels[[v]] <- l_alphaModels[[ind_minEBIC_model]]
      
      
    } # end if: alpha EBIC?
    
    
    # Update Progress Bar
    if(pbar==TRUE) setTxtProgressBar(pb, v)
    
  } # end for: p
  
  
  
  # -------------------- Processing glmnet Output -------------------
  
  # Input: p regression vectors
  # Putput: hoo-object, wadj matrix (latter contained in former)
  
  Pars_ind <- list()
  Pars_values <- list()
  list_Intercepts <- vector('list', length = p)
  
  for(v in 1:p) {
    
    model_obj <- mgmobj$nodemodels[[v]]$model
    
    # ----- Create empty Storage for parameters -----
    
    predictor_set <- (1:p)[-v] # set of predictors for variable v (all others)
    
    # A) Indicator Object for all possible interactions: List contains Matrix
    v_Pars_ind <- vector('list', length = d) #  1 = pairwise, 2= threeway, etc.
    
    # Fill each with empty entries for all possible interactions
    for(ord in 1:d) v_Pars_ind[[ord]] <- t(combn(predictor_set, ord))
    
    no_interactions <- unlist(lapply(v_Pars_ind, nrow))
    
    # B) Parameter Object: Same structure as (A), but now with a list entry for each matrix row
    v_Pars_values <- vector('list', length = d)
    for(ord in 1:d) v_Pars_values[[ord]] <- vector('list', length = no_interactions[ord])
    
    # ----- Fill (B) with parameter estimates -----
    
    # separate for categorical/ non-categorical response node, because in former case we predict K categories
    # and hence need a different data structure
    
    # browser()
    
    
    # Categorical Case
    if(type[v] == 'c') {
      
      n_cat <- level[v]
      
      for(cat in 1:n_cat) {
        
        model_obj_i <- as.numeric(model_obj[[cat]]) # select parameters for predicting category i
        interaction_names <- rownames(model_obj[[cat]])[-1] # -1 because we remove the intercept
        interaction_order <- str_count(interaction_names, ":") + 1 # get order of each interaction parameter; +1 to put on same scale as ord, so ord=2 = pairwise interaction
        
        list_Intercepts[[v]][[cat]] <- model_obj_i[1]
        
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
          
          if(no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Fuckup in parameter extraction 1')
          
          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          
          for(t in 1:no_interactions[ord]) {
            
            # indicates location of parameters for given interaction
            
            for(cc in 1:ncol(v_Pars_ind[[ord]])) find_int_dummy[, cc] <- grepl(paste0('V', v_Pars_ind[[ord]][t, cc], '.'), int_names_ord,
                                                                               int_names_ord,
                                                                               fixed = TRUE) # not only single chars have to be contained, but exact order)
            select_int <- rowSums(find_int_dummy) == ord # because threeway interaction has 2 predictors; ord = order of interactions in joint distribution
            
            # fill in paramters + rownames into list
            parmat <- matrix(part_ord[select_int], ncol = 1)
            rownames(parmat) <- int_names_ord[select_int]
            v_Pars_values[[ord]][[t]][[cat]] <- parmat
            
            gotThemall[select_int == TRUE] <- FALSE
            
          }
          
          if(sum(gotThemall) > 0) stop('Fuckup in parameter extraction 2')
          
          
          
        } # end for: ord
        
      } # end for: categories
      
      # Continuous Case
    } else {
      
      model_obj_i <- as.numeric(model_obj) # select parameters for predicting category i
      interaction_names <- rownames(model_obj)[-1] # -1 because we remove the intercept
      interaction_order <- str_count(interaction_names, ":") + 1 # on same scale as ord, so ord=2 = pairwise interaction
      
      list_Intercepts[[v]] <- model_obj_i[1]
      
      # Thresholding:
      # p = number of covariances, as it should be
      model_obj_i_ni <- model_obj_i[-1] # remove intercept, this is no covariate
      if(threshold == 'LW') tau <- sqrt(d) * sqrt(sum(model_obj_i_ni^2)) * sqrt(log(npar_standard[v]) / n)
      if(threshold == 'HW') tau <- d * sqrt(log(npar_standard[v]) / n)
      if(threshold == 'none') tau <- 0
      model_obj_i[abs(model_obj_i) < tau] <- 0 # set all parameter estimates below threshold to zero
      mgmobj$nodemodels[[v]]$tau <- tau # Save tau
      
      for(ord in 1:d) {
        
        part_ord <- model_obj_i[-1][ord == interaction_order] # parameters for interaction of order = ord
        gotThemall <- rep(TRUE, length(part_ord)) # sanity check: did I copy all parameters
        int_names_ord <- interaction_names[ord == interaction_order]
        
        if(no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Fuckup in parameter extraction 1')
        
        find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
        
        for(t in 1:no_interactions[ord]) {
          
          # indicates location of parameters for given interaction
          for(cc in 1:ncol(v_Pars_ind[[ord]])) find_int_dummy[, cc] <- grepl(paste0('V', v_Pars_ind[[ord]][t, cc], '.'),
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
        
        
      } # end for: ord
      
    }
    
    Pars_ind[[v]] <- v_Pars_ind
    Pars_values[[v]] <- v_Pars_values
    
    
  } # end for: v
  
  
  
  
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
      Pars_ind_part <-  Pars_ind[[v]][[ord]]
      Pars_ind[[v]][[ord]] <- cbind(rep(v, nrow(Pars_ind_part)), Pars_ind[[v]][[ord]])
      Pars_ind_flip[[ord]][[v]] <- Pars_ind[[v]][[ord]]
      
      # Reordering value list
      Pars_values_flip[[ord]][[v]] <- Pars_values[[v]][[ord]]
      
    }
  }
  
  # Collapse indicator list across nodes
  Pars_ind_flip_red <- lapply(Pars_ind_flip, function(x) do.call(rbind, x)  )
  
  # Collapse value list across nodes
  Pars_values_flip_red <- lapply(Pars_values_flip, function(x) do.call(c, x))
  
  
  # 1) Select each interaction
  
  # Set up target data structure
  l_factors <- list() # saves all unique interavtions
  for(i in 1:d) l_factors[[i]] <- matrix(NA, nrow = choose(p, i+1), ncol = i+1)
  l_factor_par <- list() # saves parameters associated with all unique interactions
  for(i in 1:d) l_factor_par[[i]] <- vector('list', length = choose(p, i+1))
  l_sign_par <- list() # saves sign (if defined) of all unique interactions
  for(i in 1:d) l_sign_par[[i]] <- rep(NA, choose(p, i+1))
  
  l_factor_par_full <- l_factor_par # for un-aggregated parameter esimates
  
  # Define set of continous and binary variables: for interactions between these we can assign a sign
  # Depends on binarySign
  if(binarySign) {
    set_signdefined <- c(which(type == 'p'), which(type == 'g'), ind_cat[ind_binary])
  } else {
    set_signdefined <- c(which(type == 'p'), which(type == 'g'))
  }
  
  
  # loop over: order of interactions (ord = 1 = pairwise)
  for(ord in 1:d) {
    
    set_int_ord <- Pars_ind_flip_red[[ord]]
    set_par_ord <- Pars_values_flip_red[[ord]]
    
    ids <- FlagSymmetricFast(set_int_ord) # BOTTLE NECK, now better with native solution, but still problematic for huge p
    
    unique_ids <- unique(ids)
    
    # loop over: unique interaction of order = ord
    for(i in 1:length(unique_ids)) {
      
      l_w_ind <- list()
      l_w_par <- list()
      ind_inter <- which(ids == i)
      
      # loop over the k estimates for a given k-order interaction
      for(j in 1:(ord+1)) {
        
        l_w_ind[[j]] <- set_int_ord[ind_inter[j], ]
        l_w_par[[j]] <- set_par_ord[[ind_inter[j]]]
        
      }
      
      # Mapping: Regression parameters -> Edge parameter
      m_par_seq <- unlist(lapply(l_w_par, function(x) mean(abs(unlist(x)))))
      
      # Apply AND / OR rule
      if(ruleReg == 'AND') parcompr <- mean(m_par_seq) * !(0 %in% m_par_seq)
      if(ruleReg == 'OR') parcompr <- mean(m_par_seq)
      
      # Compute Sign if defined
      if(mean(m_par_seq) != 0) { # only relevant for nonzero parameters
        
        pair <- l_w_ind[[1]]
        
        if(sum(!(pair %in% set_signdefined)) == 0) { # check of all involved varibales are g, p, or binary
          
          # if(sum(pair == c(1,3))==2) browser()
          
          
          int_sign <- getSign(l_w_ind,
                              l_w_par,
                              type,
                              set_signdefined,
                              overparameterize,
                              ord)
          

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
      l_factor_par[[ord]][[i]] <- parcompr
      l_factor_par_full[[ord]][[i]] <- l_w_par
      
      
      
    } # end for: i (unique interactions)
    
  } # end for: ord
  
  
  
  # ----- Compute (weighted) adjacency matrix -----
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
  mgmobj$rawfactor$indicator <- l_factors_nz
  mgmobj$rawfactor$weightsAgg <- l_factor_par_nz
  mgmobj$rawfactor$weights <- l_factor_par_full_nz
  mgmobj$rawfactor$signs <- l_sign_par_nz
  
  
  # -------------------- Compute Weighted Adjacency matrix (pairwise Interactions) -------------------
  
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
  
  
  
  # -------------------- Compute Factor Graph -------------------
  
  # browser()
  
  FG <- DrawFG(list_ind = l_factors_nz,
               list_weights = l_factor_par_nz,
               list_signs = l_sign_par_nz,
               p = p)
  
  
  # Save in output
  mgmobj$factorgraph$graph <- FG$weightedgraph
  mgmobj$factorgraph$nodetype <- FG$nodetype
  mgmobj$factorgraph$order <- FG$order
  
  mgmobj$factorgraph$signs <- FG$signs
  mgmobj$factorgraph$edgecolor <- FG$signcolor
  
  
  
  # -------------------- Output -------------------
  
  # Save Node Models and extracted raw factors?
  if(!saveModels) {
    mgmobj$nodemodels <- NULL
    mgmobj$rawfactor <- NULL
  }
  
  # Save intercepts
  mgmobj$intercepts <- list_Intercepts
  
  # Switch warings back on
  if(!warnings) {
    options(warn = oldw)
  }
  
  if(pbar) {
    if(signInfo) cat('\nNote that the sign of parameter estimates is stored separately; see ?mgm')    
  } else {
    if(signInfo) cat('Note that the sign of parameter estimates is stored separately; see ?mgm')    
  }
  
  
  # Assign class
  class(mgmobj) <- c('mgm', 'core')
  
  return(mgmobj)
  
  
} # eoF




