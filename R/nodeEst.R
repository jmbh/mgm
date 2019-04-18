

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

  # browser()
  
# if(v == 1) browser()
  
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


    # ----- Fit Model -----
    
    fit <- glmnet(x = X,
                  y = y,
                  family = fam,
                  alpha = alpha,
                  weights = weights,
                  lambda = lambdaSeq,
                  intercept = intercept)


    
    n_lambdas <- length(fit$lambda) # length of fitted lambda sequence

    # coef(fit, s = fit$lambda[61])
    #
    # if(v==5)browser()

    # ----- Calc EBIC of model directly (to check what I cooked up below) -----

    # l_LL_model <- list()
    # l_n_neighbors <- list()
    # for(i in 1:n_lambdas) {
    #
    #   l_LL_model[[i]] <- calcLL(X = X,
    #                             y = y,
    #                             fit = fit,
    #                             type = type,
    #                             v = v,
    #                             weights = weights,
    #                             lambda = fit$lambda[i])
    #
    #   l_n_neighbors[[i]] <- calcNeighbors(fit = fit,
    #                                lambda = fit$lambda[i],
    #                                type = type,
    #                                v = v)
    # }
    #
    #
    # v_EBIC <- - 2 * unlist(l_LL_model) + unlist(l_n_neighbors) * log(nadj) + 2 * lambdaGam * log(ncol(X))
    #
    # EBIC_min <- min(v_EBIC)
    #
    # ind_lambda_min <- which.min(v_EBIC)
    # lambda_min <- fit$lambda[ind_lambda_min]
    # lambad_min_model <- coef(fit, s = lambda_min)


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

    # Calc LL_NULL manually here

    # coef(fit, s = fit$lambda[1])
    #
    # LL_null_parts <- dnorm(y, mean = 0.005466318, sd = 1, log = TRUE)
    # sum(LL_null_parts)
    #
    # hist(y)


    LL_sat <- 1/2 * fit$nulldev + LL_null # calculate LL of saturated model
    deviance <- (1 - fit$dev.ratio) * fit$nulldev # note: dangerous in glmnet: fit$dev = fit$dev.ratio
    LL_lambda_models <- - 1/2 * deviance + LL_sat # length = length of lambda sequence


    n_neighbors <- rep(NA, n_lambdas)
    for(i in 1:n_lambdas)  n_neighbors[i] <- calcNeighbors(fit = fit,
                                                           lambda = fit$lambda[i],
                                                           type = type,
                                                           level = level,
                                                           v = v)
    # btw: n_neighbors = fit$df


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

    return(outlist)

  }


  # ---------- Lambda selection via CV ----------

  if(lambdaSel == 'CV') {

    # ----- Fit Model -----
    
    # browser()
    
    # browser()

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

    return(outlist)


  }





} # eoF
