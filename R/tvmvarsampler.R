

tvmvarsampler <- function(coefarray, # v x v2 x cat(v) x cat(v2) x lag array x N specifying lagged parameters
                          lags, # vector specifiying the lags
                          thresholds, # N x p matrix containint vector with a threshold for each category
                          sds, # N x p matrix of sds
                          type,
                          level,
                          pbar)




{

  # ---------- Input Checks ----------

  # Does coefarray have the right dimensions?
  if(length(dim(coefarray))!=6) stop('coefarray has to have 5 dimensions: p x p x max(level) x max(level) x n_lags x N.')
  if(dim(coefarray)[5] != length(lags)) stop('Dimensions for lags in coefarray have to match the number of lags specified in lags.')
  if(dim(coefarray)[1] != length(thresholds)) stop('Dimensions for variables in coefarray have to match the number of specified thresholds.')
  if(dim(coefarray)[1] != ncol(sds)) stop('Dimensions for variables in coefarray have to match the number of specified standard deviations.')
  if(dim(coefarray)[1] != length(type)) stop('Dimensions for variables in coefarray have to match the number of specified types.')
  if(dim(coefarray)[1] != length(level)) stop('Dimensions for variables in coefarray have to match the number of specified levels.')
  if(dim(coefarray)[1] != dim(coefarray)[2]) stop('The first two dimensions specifying the cross-lagged effects in coefarray have to have the same dimensionality.')

  if(class(thresholds) != 'list') stop('The thresholds have to be provided via a list with p entries. See ?tvmvarsampler')
  if(class(sds) != 'matrix') stop('The standard deviations of Gaussian variables have to be provided via a p x n matrix.')
  

  # ---------- Fill in Defaults ----------
  
  if(missing(pbar)) pbar <- TRUE
  

  # ---------- Create Output Object & copy the call ----------

  # Create Output Object
  mvarsamp_obj <- list('call' = NULL,
                       'data' = NULL)

  # Copy the call
  mvarsamp_obj$call <- list('coefarray' = coefarray,
                            'lags' = lags,
                            'thresholds' = thresholds,
                            'sds' = sds,
                            'type' = type,
                            'level' = level,
                            'pbar' = pbar)



  # ---------- Compute Auxilliary Variables ----------

  p <- dim(coefarray)[1]
  N <- dim(coefarray)[6]
  n_lags <- length(lags)
  labels <- paste0('V', 1:p, '.')
  max_lag <- max(lags)



  # ---------- Sample Data ----------

  data <- matrix(NA, ncol = p, nrow = N + max_lag)

  # Create random starting values
  for(r in 1:max_lag) {
    for(v in 1:p) {
      if(type[v] == 'c') data[r, v] <- sample(1:level[v], size = 1)
      if(type[v] == 'g') data[r, v] <- rnorm(1)
      if(type[v] == 'p') data[r, v] <- rpois(1, 1)
    }
  }


  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=N, initial=0, char="-", style = 3)


  for(r in 1:N) {
    r2 <- ((max_lag+1):(N+max_lag))[r]

    for(v in 1:p) {

      if(type[v] != 'c') {

        # Get design matrix
        design_mat <- ModelMatrix(data = data, 
                                  type = type, 
                                  level = level, 
                                  labels = labels, 
                                  d = 1, 
                                  v=v,
                                  allCats=TRUE)

        potential_lag <- list()
        for(lag in 1:n_lags) {

          # Get parameters in order
          l_parms <- list()
          for(par in 1:p) l_parms[[par]] <- coefarray[v, par, 1:level[v], 1:level[par], lag, r]
          v_parms <- unlist(l_parms)

          # Compute part of potential for given lag
          potential_lag[[lag]] <- design_mat[r2 - lags[lag], ] * v_parms

        } # end for: lag

        mu <- thresholds[[v]][r,1] + sum(unlist(potential_lag))

        if(type[v] == 'g') data[r2, v] <- rnorm(1, mean = mu, sd = sds[r, v])
        if(type[v] == 'p') data[r2, v] <- rpois(1, lambda = mu)

      } # end if: !='c'



      if(type[v] == 'c') {

        # Get design matrix
        design_mat <- ModelMatrix(data = data, 
                                  type = type, 
                                  level = level, 
                                  labels = labels, 
                                  d = 1, 
                                  v=v,
                                  allCats=TRUE)

        # Loop over categories
        l_potentials <- list()
        for(cat in 1:level[v]) {

          potential_lag <- list()
          for(lag in 1:n_lags) {

            # Get parameters in order
            l_parms <- list()
            for(par in 1:p) l_parms[[par]] <- coefarray[v, par, (1:level[v])[cat], 1:level[par], lag, r]
            v_parms <- unlist(l_parms)

            # Compute part of potential for given lag
            potential_lag[[lag]] <- design_mat[r2 - lags[lag], ] * v_parms

          } # end for: lag

          # put potential-parts from different lags together
          l_potentials[[cat]] <- thresholds[[v]][r, cat] + sum(unlist(potential_lag))

        } # end for: cat

        # if(r==2 & v == 1) browser()

        # compute probabilities
        v_potentials <- exp(unlist(l_potentials))
        probabilities <- v_potentials / sum(v_potentials)

        data[r2, v] <- sample(1:level[v], prob = probabilities, size = 1)

      } # end if: =='c'?


    } # end for: v variables

    if(pbar==TRUE) setTxtProgressBar(pb, r)

  } # end for: r rows


  # ---------- Output ----------


  mvarsamp_obj$data <- data[(max_lag+1):(N+max_lag), ]

  class(mvarsamp_obj) <- 'tvmvar'


  return(mvarsamp_obj)




} # eoF
