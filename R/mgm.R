

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
                moderators,   # Vector specifying first-order moderators to be included in the model
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
                thresholdCat, # TRUE if overparameterize=FALSE; FALSE if overparamterize=TRUE; this argument overwrites this
                signInfo,
                ...
)


{
  
  # -------------------- Input Checks -------------------
  
  # ----- Compute Auxilliary Variables I -----
  
  p <- ncol(data)
  n <- nrow(data)
  data <- as.matrix(data)
  
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
  if(missing(moderators)) moderators <- NULL
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
  if(missing(thresholdCat)) if(overparameterize) thresholdCat <- TRUE else thresholdCat <- TRUE # always better
  
  if(verbatim) pbar <- FALSE
  if(verbatim) warnings <- FALSE
  
  # Switch all warnings off
  if(!warnings) {
    oldw <- getOption("warn")
    options(warn = -1)
  }
  
  # ----- Basic Checks I -----
  
  # Checks on data
  if(nrow(data) < 2) ('The data matrix has to have at least 2 rows.')
  if(any(!(apply(data, 2, class) %in% c('numeric', 'integer')))) stop('Only integer and numeric values permitted.')
  if(missing(data)) stop('No data provided.')
  if(any(!is.finite(as.matrix(data)))) stop('No infinite values permitted.')
  if(any(is.na(data))) stop('No missing values permitted.')
  if(!all(moderators %in% 1:p)) stop("Specified moderators are larger than number of variables in the data.")
  
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
  
  
  
  # ----- Basic Checks II -----
  
  # Checks on other arguments
  if(!(threshold %in% c('none', 'LW', 'HW'))) stop('Please select one of the three threshold options "HW", "LW" and "none" ')
  if(k<2) stop('The order of interactions should be at least k = 2 (pairwise interactions)')
  if(ncol(data)<3) stop('At least 3 variables required')
  if(missing(type)) stop('No type vector provided.')
  if(sum(!(type %in% c('g', 'c', 'p')))>0) stop("Only Gaussian 'g', Poisson 'p' or categorical 'c' variables permitted.")
  if(ncol(data) != length(type)) stop('Number of variables is not equal to length of type vector.')
  if(!missing(level)) if(ncol(data) != length(level)) stop('Number of variables is not equal to length of level vector.')
  if(nrow(data) != length(weights)) stop('Number of observations is not equal to length of weights vector.')
  if(!is.null(moderators) & k > 2) stop("Moderator specification is only implemented for first-order moderators. See ?mgm")
  
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
                 'interactions' = list('indicator' = NULL,
                                       'weightsAgg' = NULL,
                                       'weights' = NULL,
                                       'signs' = NULL),
                 'intercepts' = NULL,
                 'nodemodels' = list())
  
  
  # ----- Save the Call -----
  
  
  mgmobj$call <- list('data' = NULL,
                      'type' = type,
                      'level' = level,
                      "levelNames" = NULL,
                      'lambdaSeq' = lambdaSeq,
                      'lambdaSel' = lambdaSel,
                      'lambdaFolds' = lambdaFolds,
                      'lambdaGam' = lambdaGam,
                      'alphaSeq' = alphaSeq,
                      'alphaSel' = alphaSel,
                      'alphaFolds' = alphaFolds,
                      'alphaGam' = alphaGam,
                      'k' = k,
                      "moderators" = moderators, 
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
                      'overparameterize' = overparameterize,
                      "thresholdCat" = thresholdCat,
                      "signInfo" = signInfo)
  
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
  l_mods_ind <- list() # collect moderator terms for later output-processing
  
  for(v in 1:p) {
    
    # ----- Construct Design Matrix -----
    
    # Case I: No Moderation
    
    if(is.null(moderators)) {
      
      # Case I.a: Standard parameterization
      if(d > (p - 1)) {
        stop("Order of interactions cannot be larger than the number of predictors.")
      } else if (d == 1){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
      } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^",d)) }
      
      # Construct standard design matrix (to get number of parameters for tau threshold)
      X_standard <- model.matrix(form, data = data)[, -1] # delete intercept (added by glmnet later)
      npar_standard[v] <- ncol(X_standard)
      
      # Case I.b: Overparameterization (no intercept, all categories in design matrix)
      if(overparameterize) {
        
        # Construct over-parameterized design matrix
        X_over <- ModelMatrix(data = data[, -v],
                              type = type[-v],
                              level = level[-v],
                              labels = colnames(data)[-v],
                              d = d, 
                              moderators = NULL,
                              v = v)
        
        X <- X_over
      } else {
        X <- X_standard
      }
      
      
      # Case II: Moderation
    } else {
      
      # Case II.a: Standard parameterization
      
      # Simple predictors
      pred_simple <- paste(colnames(data)[-v], collapse = " + ")
      
      # Moderation terms
      l_mods <- list()
      n_mods <- length(moderators)
      n_l_mods <- 0
      
      
      if(v %in% moderators) {
        other_comb_pairw <- t(combn((1:p)[-v], 2))
        l_mods[[1]] <- paste0("V", other_comb_pairw[, 1], ".",  " * ", "V", other_comb_pairw[, 2], ".", collapse = " + ")
      } else {      
        for(i in 1:n_mods) {  # loop over moderators
          l_mods[[i]] <- paste0(colnames(data)[-c(v, moderators[i])], "*", "V", moderators[i], ".", collapse = " + ")
          n_l_mods <- n_l_mods + 1
        }
      } # end if: v = moderator?
      
      if(n_l_mods > 1) for(i in 1:(n_l_mods - 1)) l_mods[[i]] <- paste0(l_mods[[i]], " + ") # add plus between terms but not in the end
      v_mods <- do.call(paste0, l_mods)
      
      # if(v == 7) browser()
      # print(v)
      
      # Stitch together model formula
      if(length(v_mods) == 0) { # if there is no moderation for that variable (the case if v is the only moderator)
        form <- paste(colnames(data)[v], "~",
                      pred_simple)
        
      } else {
        form <- paste(colnames(data)[v], "~",
                      pred_simple,
                      " + ",
                      v_mods)
      }
      
      form <- as.formula(form)
      
      # Create design matrix
      X_standard <- model.matrix(form, data = data)[, -1] # delete intercept (added by glmnet later)
      npar_standard[v] <- ncol(X_standard)
      
      # Case II.b: Overparameterization (all category-configuration in design matrix)
      
      if(overparameterize) {
        
        X_over <- ModelMatrix(data = data[, -v],
                              type = type[-v],
                              level = level[-v],
                              labels = colnames(data)[-v],
                              d = d,
                              moderators = moderators, 
                              v = v)
        
        X <- X_over
        
      } else {
        
        X <- X_standard
        
        
      } # end if: overparameterize
      
      
    } # end if: Moderation
    
    
    ## Scale Gaussian variables after computing design matrix
    # Compute vector that tell us which interactions are purely consisting of continuous variables?
    if(scale) {
      # browser()
      if(any(type == "g")) {
        cn <- colnames(X)
        l_ints_split <- strsplit(cn, ":")
        ind_allc <- unlist(lapply(l_ints_split, function(x) {
          x2 <- sub("V", "", x)
          vars <- as.numeric(unlist(lapply(strsplit(x2, "[.]"), function(y) y[1] )))
          all(type[vars] == "g")
        }))
        
        X[, ind_allc] <- apply(matrix(X[, ind_allc], ncol = sum(ind_allc)), 2, scale)
      }
    }
    
    # Response 
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
                                            overparameterize = overparameterize,
                                            thresholdCat = thresholdCat)
            
            
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
        
        
        # If there is no search for alpha, goes continue and use default alpha
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
                       overparameterize = overparameterize,
                       thresholdCat = thresholdCat)
      
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
                                      overparameterize = overparameterize,
                                      thresholdCat = thresholdCat)
        
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
  # Putput: hoo-object, wadj matrix (latter defined by former)
  
  Pars_ind <- list() # storage for interaction-indicators
  Pars_values <- list() # storage for actual parameters associate for the interactions indexed in Pars_ind
  list_Intercepts <- vector('list', length = p)
  
  for(v in 1:p) {
    
    model_obj <- mgmobj$nodemodels[[v]]$model
    
    # ----- Create empty Storage for parameters -----
    
    predictor_set <- (1:p)[-v] # set of predictors for variable v (all others)
    
    # A) Indicator Object for all possible interactions: List contains Matrix
    v_Pars_ind <- vector('list', length = d) #  1 = pairwise, 2= threeway, etc.
    
    ## Create Indicator list for Present interactions (separate for k-order MGM and moderated MGM)
    if(is.null(moderators)) {
      
      # If there are no moderators: All interactions of specified degree (k)
      for(ord in 1:d) v_Pars_ind[[ord]] <- t(combn(predictor_set, ord))
      
    } else {
      
      d <- 2
      
      v_Pars_ind[[1]] <- matrix(predictor_set, ncol=1) # standard terms
      
      if(v %in% moderators) {
        ind_mods <- t(combn((1:p)[-v], 2)) # if moderator, all combinations of other variables
      } else {
        ind_mods <- expand.grid((1:p)[-v], moderators[moderators!=v]) # if not moderator, all variables times 
      }
      
      ind_mods <- ind_mods[!apply(ind_mods, 1, function(x) x[1] == x[2]), ] # remove rows with equal entries
      v_Pars_ind[[2]] <- ind_mods
      
      # no interactions k > 2 allowed, if moderators are specified
      
      
    } # end if: moderators?
    
    # Make sure all entries of "v_Pars_ind" are matrices
    for(j in 1:d) v_Pars_ind[[j]] <- matrix(v_Pars_ind[[j]], ncol=j)

    no_interactions <- unlist(lapply(v_Pars_ind, nrow))
    
    # B) Parameter Object: Same structure as (A), but now with a list entry for each matrix row
    v_Pars_values <- vector('list', length = d)
    for(ord in 1:d) v_Pars_values[[ord]] <- vector('list', length = no_interactions[ord])
    
    # browser()
    
    # ----- Fill (B) with parameter estimates -----
    
    # separate for categorical/ non-categorical response node, because in former case we predict K categories
    # and hence need a different data structure
    
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
        
        # browser()
        
        for(ord in 1:d) {
          
          part_ord <- model_obj_i[-1][ord == interaction_order] # parameters for interaction of order = ord+1
          gotThemall <- rep(TRUE, length(part_ord)) # sanity check: did I copy all parameters
          int_names_ord <- interaction_names[ord == interaction_order]
          
          if(no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Internal Error: Error in parameter extraction type 1')
          
          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          
          for(t in 1:no_interactions[ord]) {
            
            # indicates location of parameters for given interaction
            if(ord==2) browser()
            
            for(cc in 1:ord) find_int_dummy[, cc] <- grepl(paste0('V', v_Pars_ind[[ord]][t, cc], '.'), int_names_ord,
                                                                               int_names_ord,
                                                                               fixed = TRUE) # not only single chars have to be contained, but exact order)
            select_int <- rowSums(find_int_dummy) == ord # because threeway interaction has 2 predictors; ord = order of interactions in joint distribution
            
            # fill in paramters + rownames into list
            parmat <- matrix(part_ord[select_int], ncol = 1)
            rownames(parmat) <- int_names_ord[select_int]
            v_Pars_values[[ord]][[t]][[cat]] <- parmat
            
            gotThemall[select_int == TRUE] <- FALSE
            
          }
          
          if(sum(gotThemall) > 0) stop('Internal Error: Error in parameter extraction type 2')
          
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
        
        if(no_interactions[ord] > 0) { # can be zero for moderated MGMs
          
          part_ord <- model_obj_i[-1][ord == interaction_order] # parameters for interaction of order = ord
          gotThemall <- rep(TRUE, length(part_ord)) # sanity check: did I copy all parameters
          int_names_ord <- interaction_names[ord == interaction_order]
          
          if(no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Fuckup in parameter extraction 1')
          
          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          
          for(t in 1:no_interactions[ord]) {
            
            # if(ord == 2) browser()
            
            # indicates location of parameters for given interaction
            for(cc in 1:ord) find_int_dummy[, cc] <- grepl(paste0('V', v_Pars_ind[[ord]][t, cc], '.'),
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
          
        } # end if: no_interactions[ord] > 0
        
      } # end for: ord
      
    }
    
    # if(v == 7) browser()
    
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
      Pars_ind_part <- Pars_ind[[v]][[ord]]
      colnames(Pars_ind_part) <- NULL
      Pars_ind_part <- as.matrix(Pars_ind_part)
      
      Pars_ind[[v]][[ord]] <- cbind(rep(v, nrow(Pars_ind_part)), Pars_ind_part)
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
  
  # Compute number of interactions for each order
  #     Note that we could take this information also from the design matrices; however, we compute it theoretically as a sanity check
  n_terms_d <- rep(NA, d)
  if(is.null(moderators)) {
    for(ord in 1:d) n_terms_d[ord] <- choose(p, ord+1)
  } else {
    n_terms_d[1] <- choose(p, 1+1) # all pairwise interactions
    mod_terms <- expand.grid((1:p), (1:p), moderators) 
    id_uni <- FlagSymmetricFast(mod_terms)
    mod_terms2 <- mod_terms[!duplicated(id_uni), ]
    ind_diff <- as.numeric(apply(mod_terms2, 1, function(x) !any(duplicated(x))))
    mod_terms3 <- mod_terms2[ind_diff == 1, ]
    n_terms_d[2] <- nrow(mod_terms3) # ok to have interactions double; will remove them below using FlagSymmetricFast()
  }
  
  # Set up target data structure
  l_factors <- list() # saves all unique interavtions
  for(ord in 1:d) l_factors[[ord]] <- matrix(NA, nrow = n_terms_d[ord], ncol = ord+1)
  l_factor_par <- list() # saves parameters associated with all unique interactions
  for(ord in 1:d) l_factor_par[[ord]] <- vector('list', length = n_terms_d[ord])
  l_sign_par <- list() # saves sign (if defined) of all unique interactions
  for(ord in 1:d) l_sign_par[[ord]] <- rep(NA, n_terms_d[ord])
  
  l_factor_par_full <- l_factor_par # for un-aggregated parameter esimates
  
  # Define set of continous and binary variables: for interactions between these we can assign a sign
  # Depends on binarySign
  if(binarySign) {
    set_signdefined <- c(which(type == 'p'), which(type == 'g'), ind_cat[ind_binary])
  } else {
    set_signdefined <- c(which(type == 'p'), which(type == 'g'))
  }
  
  
  # browser()
  
  
  # Loop over: order of interactions (ord = 1 = pairwise)
  for(ord in 1:d) {
    
    set_int_ord <- Pars_ind_flip_red[[ord]]
    set_par_ord <- Pars_values_flip_red[[ord]]
    row.names(set_int_ord) <- NULL
    
    # browser()
    
    ids <- FlagSymmetricFast(x = set_int_ord) # BOTTLE NECK, now better with native solution, but still problematic for huge p
    
    # Get set of unique interactions
    unique_set_int_ord <- cbind(set_int_ord, ids)[!duplicated(ids), ]
    unique_set_int_ord <- matrix(unique_set_int_ord, ncol = ord+1+1)
    n_unique <- nrow(unique_set_int_ord)
    
    # print(v)
    # browser()
    
    # loop over: unique interaction of order = ord
    for(i in 1:n_unique) {
      
      l_w_ind <- list()
      l_w_par <- list()
      ind_inter <- which(ids == i)
      
      # loop over the k estimates for a given k-order interaction
      for(j in 1:(ord+1)) {
        l_w_ind[[j]] <- set_int_ord[ind_inter[j], ]
        l_w_par[[j]] <- set_par_ord[[ind_inter[j]]]
      }
      
      # Mapping: Regression parameters -> Edge parameter (mean of absolute values of parameters)
      m_par_seq <- unlist(lapply(l_w_par, function(x) mean(abs(unlist(x)))))
      
      # Apply AND / OR rule
      if(ruleReg == 'AND') parcompr <- mean(m_par_seq) * !(0 %in% m_par_seq)
      if(ruleReg == 'OR') parcompr <- mean(m_par_seq)
      
      # Compute Sign if defined
      if(mean(m_par_seq) != 0) { # only relevant for nonzero parameters
        
        pair <- l_w_ind[[1]] # order doesn't matter, could take any but first entry is always filled independent of "ord", so first
        
        if(sum(!(pair %in% set_signdefined)) == 0) { # check of all involved varibales are g, p, or binary
          
          # computes combined sign (if defined) over k terms for same interaction
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
      
      # if(v == 7 & ord ==1) browser()
      
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
  
  # browser()
  
  # Save in output
  mgmobj$interactions$indicator <- l_factors_nz
  mgmobj$interactions$weightsAgg <- l_factor_par_nz
  mgmobj$interactions$weights <- l_factor_par_full_nz
  mgmobj$interactions$signs <- l_sign_par_nz
  
  
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
  
  
  # -------------------- Output -------------------
  
  # Save Node Models and extracted raw factors?
  if(!saveModels) {
    mgmobj$nodemodels <- NULL
    mgmobj$factors <- NULL
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
  
  # Return level names (used in showInteraction())
  levelNames <- list()
  for(i in 1:p) {
    if(type[i] == "c") levelNames[[i]] <- sort(as.numeric(as.character(unique(data[, i])))) else levelNames[[i]] <- NA
  }
  
  mgmobj$call$levelNames <- levelNames
  
  # Assign class
  class(mgmobj) <- c('mgm', 'core')
  
  return(mgmobj)
  
  
} # eoF




