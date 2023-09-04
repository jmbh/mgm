

nodeEst <- function(y,
                    X,
                    fam,
                    lambdaSeq,
                    lambdaSel,
                    lambdaFolds,
                    lambdaGam,
                    alpha,
                    weights,
                    n,
                    nadj,
                    v,
                    type,
                    level,
                    emp_lev,
                    overparameterize, 
                    thresholdCat)


{
  
  # ---------- Calc Aux Variables ----------
  
  # Define Exponential Family of Node at Hand
  if(type[v] == 'c') fam = 'multinomial'
  if(type[v] == 'g') fam = 'gaussian'
  if(type[v] == 'p') fam = 'poisson'
  
  
  # Set threshold (intercept) parameter to zero?
  if(type[v] == 'c') {
    intercept <- thresholdCat
  } else {
    intercept <- TRUE # for continuous variables always estimated
  }
  
  
  # ---------- Lambda selection via EBIC ----------
  
  if(lambdaSel == 'EBIC') {
    
    
    # if(v==3) browser()FA
    
    # ----- Fit Model -----
    
    fit <- glmnet(x = X,
                  y = y,
                  family = fam,
                  alpha = alpha,
                  weights = weights,
                  lambda = lambdaSeq,
                  intercept = intercept)
    
    
    
    n_lambdas <- length(fit$lambda) # length of fitted lambda sequence

    # ----- Calc EBIC of model: Fast Alternative by using -----
    
    # Calculate LL of Null Model
    LL_null <- calcLL(X = X,
                      y = y,
                      fit = fit,
                      type = type,
                      level = level,
                      v = v,
                      weights = weights,
                      lambda = fit$lambda[1], # any is fine, lambda has no influence on null (intercept) model
                      LLtype = 'nullmodel')
    
    LL_sat <- 1/2 * fit$nulldev + LL_null # calculate LL of saturated model
    deviance <- (1 - fit$dev.ratio) * fit$nulldev # note: dangerous in glmnet: fit$dev = fit$dev.ratio
    LL_lambda_models <- - 1/2 * deviance + LL_sat # length = length of lambda sequence
    
    
    n_neighbors <- rep(NA, n_lambdas)
    for(i in 1:n_lambdas)  n_neighbors[i] <- calcNeighbors(fit = fit,
                                                           lambda = fit$lambda[i],
                                                           type = type,
                                                           level = level,
                                                           v = v)
    # Note: n_neighbors = fit$df
    EBIC_lambda <- - 2 * LL_lambda_models + n_neighbors * log(nadj) + 2 * lambdaGam * n_neighbors * log(ncol(X))
    
    EBIC_min <- min(EBIC_lambda)
    
    ind_lambda_min <- which.min(EBIC_lambda)
    lambda_min <- fit$lambda[ind_lambda_min]
    lambad_min_model <- coef(fit, s = lambda_min)
    
    # ----- Output -----
    
    outlist <- list('EBIC' = EBIC_min,
                    'deviance' = deviance[which.min(EBIC_lambda)],
                    'lambda' = lambda_min,
                    'alpha' = alpha,
                    'model' = lambad_min_model,
                    'fitobj' = fit)
    
  } # end if: EBIC
  
  
  # ---------- Lambda selection via CV ----------
  
  if(lambdaSel == 'CV') {
    
    # ----- Fit Model -----
    
    fit <- cv.glmnet(x = X,
                     y = y,
                     family = fam,
                     alpha = alpha,
                     weights = weights,
                     nfolds = lambdaFolds,
                     type.measure = "deviance",
                     lambda = lambdaSeq,
                     intercept = intercept)
    
    lambda_min <-  fit$lambda.min
    lambad_min_model <- coef(fit, s = lambda_min)
    
    # ----- Calc Deviance of model -----
    # (used in alpha selection via EBIC)
    
    LL_model <- calcLL(X = X,
                       y = y,
                       fit = fit,
                       type = type,
                       level = level,
                       v = v,
                       weights = weights,
                       lambda = lambda_min,
                       LLtype = 'model')
    
    n_neighbors <- calcNeighbors(fit = fit,
                                 lambda = lambda_min,
                                 type = type,
                                 level = level,
                                 v = v)
    
    EBIC <- - 2 * LL_model + n_neighbors * log(nadj) + 2 * lambdaGam * log(ncol(X))
    
    
    # ----- Output -----
    
    outlist <- list('EBIC' = EBIC,
                    'deviance' = NULL,
                    'lambda' = lambda_min,
                    'alpha' = alpha,
                    'model' = lambad_min_model,
                    'fitobj' = fit)
    
  } # end if: CV
  
  
  # ----- Return -----
  
  return(outlist)
  
  
  
} # eoF
