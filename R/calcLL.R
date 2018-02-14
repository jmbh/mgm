

calcLL <- function(X,
                   y,
                   fit, #glmnet fit object
                   type,
                   level,
                   v,
                   weights,
                   lambda,
                   LLtype = 'model')


  {

  if(missing(level)) stop('No levels passed to calcLL !')

  # This function calculates three different LL:
  # 1) LLtype = 'model': The LL of a given model via fit
  # 2) LLtype = 'nullmodel': The LL of the Null (Intercept) Model
  # 3) LLtype = 'saturated': The LL of the saturated model

  n <- nrow(X)

  if(LLtype == 'model') {

    if(type[v] == 'g') {
      beta_vector <- matrix(coef(fit, s = lambda), ncol = 1)
      predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
      LL_model <- dnorm(y, mean = predicted_mean, sd = 1, log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }

    if(type[v] == 'p') {
      beta_vector <- matrix(coef(fit, s = lambda), ncol = 1)
      predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }

    if(type[v] == 'c') {

      n_cats <- level[v] # number of levels

      ## Compute LL (see http://www.stanford.edu/~hastie/Papers/glmnet.pdf, equation 22)
      m_respdum <- matrix(NA, n, n_cats) # dummy for data
      m_coefs <- matrix(NA, n, n_cats) # dummy for coefficients
      cats <- unique(y)

      LL_n <- rep(NA, n) # Storage
      m_LL_parts <- matrix(NA, nrow = n, ncol=n_cats+1)

      for(catIter in 1:n_cats) {
        m_respdum[,catIter] <- (y==cats[catIter])*1 # dummy matrix for categories
        m_coefs[,catIter] <- cbind(rep(1, n), X) %*% matrix(coef(fit, s = lambda)[[catIter]], ncol = 1)
        m_LL_parts[,catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
      }

      m_LL_parts[, n_cats+1] <- - log(rowSums(exp(m_coefs))) # the log part, see eq (22)
      LL_n <- rowSums(m_LL_parts) # sum up n_cat + 1 parts
      mean_LL_model <- sum(LL_n * weights) # apply weighting

    }

  }

  if(LLtype == 'nullmodel') {

    if(type[v] == 'g') {
      beta_vector <- matrix(coef(fit, s = 1)[1], ncol = 1) # only intercept here
      predicted_mean <- rep(1, n) * as.vector(beta_vector)
      LL_model <- dnorm(y, mean = predicted_mean, sd = 1, log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }

    if(type[v] == 'p') {
      beta_vector <- matrix(coef(fit, s = 1)[1], ncol = 1)
      predicted_mean <- rep(1, n) * as.vector(beta_vector) # log mean actually
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }

    if(type[v] == 'c') {

      n_cats <- level[v] # number of levels

      ## Compute LL (see http://www.stanford.edu/~hastie/Papers/glmnet.pdf, equation 22)
      m_respdum <- matrix(NA, n, n_cats) # dummy for data
      m_coefs <- matrix(NA, n, n_cats) # dummy for coefficients
      cats <- unique(y)

      LL_n <- rep(NA, n) # Storage
      m_LL_parts <- matrix(NA, nrow = n, ncol=n_cats+1)

      for(catIter in 1:n_cats) {
        m_respdum[,catIter] <- (y==cats[catIter])*1 # dummy matrix for categories
        m_coefs[,catIter] <- cbind(rep(1, n), X) %*% matrix(coef(fit, s = 1)[[catIter]], ncol = 1)
        m_LL_parts[,catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
      }

      m_LL_parts[, n_cats+1] <- - log(rowSums(exp(m_coefs))) # the log part, see eq (22)
      LL_n <- rowSums(m_LL_parts) # sum up n_cat + 1 parts
      mean_LL_model <- sum(LL_n * weights) # apply weighting

    }

  }


  if(LLtype == 'saturated') {

    if(type[v] == 'g') {
      predicted_mean <- y
      LL_model <- dnorm(y, mean = predicted_mean, sd = 1, log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }

    if(type[v] == 'p') {
      predicted_mean <- y
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }

    if(type[v] == 'c') {

      mean_LL_model <- 0

      # For discrete RVs,the saturated model has Likelihood = 1 and LL = log(1) = 0
      # e.g. http://stats.stackexchange.com/questions/114073/logistic-regression-how-to-obtain-a-saturated-model

    }

  }



  return(mean_LL_model)

}
