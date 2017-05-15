
bwSelPredict <- function(data,
                         type,
                         level,
                         obj,
                         test,
                         modeltype,
                         ...) # is lags passed on?



{

  # -------------------- Input checks & calc Aux Vars -------------------

  p <- ncol(data)

  # Fetch arguments
  args <- list(...)

  args$d <- args$k - 1
  n_test <- length(test)

  # Scale Gaussians
  ind_Gauss <- which(type == 'g')
  for(i in ind_Gauss) data[, i] <- scale(data[, i])

  data_df <- as.data.frame(data)

  # Categoricals into factors (Needed to use formula to construct design matrix)
  for(i in which(type=='c')) data_df[, i] <- as.factor(data_df[, i])


  # -------------------- Get VAR data structure -------------------

  if(modeltype == 'mvar') {

    # Glue all lagged variables together

    n_lags <- length(args$lags)

    data_lagged <- lagData(data_df, args$lags)

    data_response <- data_lagged$data_response
    l_data_lags <- data_lagged$l_data_lags

    data_response <- apply(data_response, 2, as.numeric) # to avoid confusion with labels of categories if there are factors

  }

  # -------------------- Compute Preditions -------------------

  m_pred <- matrix(NA, nrow = n_test, ncol = p)

  for(v in 1:p) {

    # ----- Create Design Matrices -----

    if(modeltype == 'mvar') {

      # append response with predictors
      y <- data_response[,v] # response variable v
      data_v <- cbind(y, do.call(cbind, l_data_lags)) # combine
      data <- data_v[, -1]

      # Dummy coding
      form <- as.formula('y ~ (.)')

      # Construct standard design matrix (to get number of parameters for tau threshold)
      X_standard <- model.matrix(form, data = data_v)[, -1] # delete intercept (added by glmnet later)

      if(args$overparameterize) {

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

    }

    if(modeltype == 'mgm') {

      # Create HOI design matrix (copied from mgm)

      if(args$d > (p - 1)) {
        stop("Order of interactions cannot be larger than the number of predictors.")
      } else if (args$d == 1){ form <- as.formula(paste(colnames(data_df)[v],"~ (.)"))
      } else { form <- as.formula(paste(colnames(data_df)[v],"~ (.)^", args$d)) }

      # Construct standard design matrix (to get number of parameters for tau threshold)
      X_standard <- model.matrix(form, data = data_df)[, -1] # delete intercept (added by glmnet later)

      if(args$overparameterize) {

        # Construct over-parameterized design matrix
        X_over <- ModelMatrix(data = data[, -v],
                              type = type[-v],
                              level = level[-v],
                              labels = colnames(data)[-v],
                              d = args$d)  # defined above from  k
        X <- X_over
      } else {
        X <- X_standard
      }

      y <- data[, v]

    } # end if: modeltype mgm





    # ----- Make Predictions for test points -----

    for(i in 1:n_test) {

      coefs <- obj$tvmodels[[i]]$nodemodels[[v]]$model
      row <- test[i]

      if(type[v] == 'c') {

        n_cat <- level[v]

        # Compute exp(potentials)
        Potentials <- rep(NA, n_cat)
        for(k in 1:n_cat)  Potentials[k] <- exp(coefs[[k]][1] + X[row,] %*% coefs[[k]][-1])

        # compute category-probabilities
        Probabilities <- Potentials / sum(Potentials)
        # Classify
        m_pred[i, v] <- sort(unique(y))[which.max(Probabilities)]


      } else {

        ## Prediction Continuous
        coefs <- as.numeric(coefs) # get coefficients
        m_pred[i, v] <- coefs[1] + X[row,] %*% coefs[-1] # predict

      }

    } # end for: i (over test locations)

  } # end for: v


  # -------------------- Compute Errors -------------------

  if(modeltype == 'mvar') m_true <- matrix(data_response[test, ], nrow = n_test) # enforce matrix in case n_test = 1
  if(modeltype == 'mgm') m_true <- matrix(data[test, ], nrow = n_test) # enforce matrix in case n_test = 1

  m_error <- matrix(NA, nrow = n_test, ncol = p)

  
  # 0/1 loss for categorical, RMSE for continuous (which here is just absolute error, because 1 element)
  for(i in 1:p) {
    if(type[v] == 'c') {
      m_error[, i] <- abs(m_true[, i] - m_pred[, i]) 
    } else {
      m_error[, i] <- abs(m_true[, i] - m_pred[, i])
    }
  }


  # -------------------- Output -------------------

  predOutlist <- list()
  predOutlist$errors <- m_error
  predOutlist$error_time <- rowMeans(m_error)
  predOutlist$error_mean <- mean(m_error)


  return(predOutlist)


}
