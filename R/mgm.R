

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
  
  # --------------------------------------------------------------------------------------------
  # -------------------- Input Checks ----------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
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
  
  # browser()
  
  # Checks on moderators
  if(!is.null(moderators)) {
    if(!all(moderators == round(moderators))) stop("Moderators have to be specified as integers mapping to the column numbers of variables in the data set.")
    if(!all(moderators %in% 1:p)) stop("Specified moderators are larger than number of variables in the data.")
    if(class(moderators) == "matrix") if(ncol(moderators) != 3) stop("Custom moderators have to be specified in a M x 3 matrix, for M moderators.")
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
  
  
  # Get unique values for all categorical variables (used in condition.R)
  unique_cats <- list()
  for(i in 1:p) if(type[i]=="c") unique_cats[[i]] <- unique(data[, i])

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
  if(!is.null(moderators) & k > 2) stop("Moderator specification is only implemented for first-order moderation (3-way interactions). See ?mgm")
  
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
  
  
  # --------------------------------------------------------------------------------------------
  # -------------------- Create Output Objects -------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  
  # ----- Storage: Create empty mgm object -----
  
  mgmobj <- list('call' = NULL,
                 'pairwise' = list('wadj' = NULL,
                                   'signs' = NULL,
                                   'edgecolor'= NULL, 
                                   "wadjNodewise" = NULL,
                                   "signsNodewise" = NULL,
                                   "edgecolorNodewise" = NULL),
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
                      "signInfo" = signInfo, 
                      "npar"=NULL,
                      "n"=nrow(data), 
                      "ind_cat" = ind_cat, 
                      "ind_binary" = ind_binary, 
                      "unique_cats" = unique_cats)
  
  if(saveData) mgmobj$call$data <- data
  
  
  # ----- Some more variable Transforms -----
  
  data <- as.data.frame(data)
  
  # Categoricals into factors (Needed to use formula to construct design matrix)
  for(i in which(type=='c')) data[, i] <- as.factor(data[, i])
  
  
  # --------------------------------------------------------------------------------------------
  # -------------------- Nodewise Estimation ---------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  
  # Progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=p, initial=0, char="-", style = 3)
  
  # Save number of parameters of standard (non-overparameterized) design matrices for tau threshold
  npar_standard <- rep(NA, p)
  
  l_mods_ind <- list() # collect moderator terms for later output-processing
  
  for(v in 1:p) {
    
    # ----- Construct Design Matrix -----
    
    # browser()
    
    X_standard <- X <- ModelMatrix_standard(data = data,
                                            type = type,
                                            d = d, 
                                            v = v, 
                                            moderators = moderators)
    
    npar_standard[v] <- ncol(X_standard)
    
    
    if(overparameterize) {
      
      X_over <- ModelMatrix(data = data, # fix that input, that's stupid
                            type = type,
                            level = level,
                            labels = colnames(data),
                            d = d, 
                            moderators = moderators,
                            v = v)
      
      X <- X_over
      
    } # end if: overparameterize?
    

    ## Scale Gaussian variables AFTER computing design matrix
    # Compute vector that tell us which interactions are purely consisting of continuous variables?
    if(scale) {
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
  
  mgmobj$call$npar <- npar_standard
  
  
  # --------------------------------------------------------------------------------------------
  # -------------------- Processing glmnet output ----------------------------------------------
  # --------------------------------------------------------------------------------------------  
  
  mgmobj <- Reg2Graph(mgmobj = mgmobj)

  # --------------------------------------------------------------------------------------------
  # -------------------- Output ----------------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  
  # Save Node Models and extracted raw factors?
  if(!saveModels) {
    mgmobj$nodemodels <- NULL
    mgmobj$factors <- NULL
  }
  
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




