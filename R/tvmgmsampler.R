

tvmgmsampler <- function(factors,
                         interactions,
                         thresholds,
                         sds,
                         type,
                         level,
                         nIter = 250,
                         pbar = TRUE,
                         ...)


{

  # ---------- Input Checks ----------


  # ----- Calc Aux Variables -----

  n_timepoints <- nrow(thresholds[[1]])
  n_order <- length(factors)
  p <- length(level)

  data <- matrix(NA,
                 nrow = n_timepoints,
                 ncol = p)

  # ----- Input Checks -----

  # - same nrow for all thresholds
  # - match time dimension for factors, interaction, thresholds and sds
  # - check time vector: is it a integer sequence 1:N
  # - only checks on sd if at least one gaussian present!

  if(missing(sds)) sds <- NULL


  # ----- Create Outout object -----

  tvmgmsamp_obj <- list('call' = NULL,
                        'data' = NULL)

  # Copy the call
  tvmgmsamp_obj$call <- list('factors' = factors,
                             'interactions' = interactions,
                             'thresholds' = thresholds,
                             'sds' = sds,
                             'type' = type,
                             'level' = level,
                             'nIter' = nIter,
                             'pbar' = pbar)


  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max = n_timepoints, initial=0, char="-", style = 3)

  for(ts in 1:n_timepoints) {

    # ---------- Select Time dimension ----------


    # a) Take out time point from: interactions

    factors_ts <- list()
    interactions_ts <- list()

    for(ord in 1:n_order) {

      n_row <- dim(factors[[ord]])[1]

      if(is.null(n_row)) {

        interactions_ts[[ord]] <- NULL

      } else {

        # factors_ts[[ord]] <- matrix(NA, nrow = n_row, ncol = ord+1)
        interactions_ts[[ord]] <- list()

        for(row in 1:n_row) {

          # factors_ts[[ord]][row, ] <- fr <- factors[[ord]][row, , ts]
          fr <- factors[[ord]][row, ]

          # create list we use to access dynamic array: interactions
          l_query <- list()
          for(i2 in 1:(ord+1)) l_query[[i2]] <- 1:level[fr[i2]]
          l_query[[(ord+2)]] <- ts

          A <- interactions[[ord]][[row]]
          query_out <- do.call(function(...)A[...], l_query)
          dim_ts_fix_array <- unlist(lapply(l_query, length))[-(ord+2)]
          interactions_ts[[ord]][[row]] <- array(query_out, dim = dim_ts_fix_array)

        }

      } # end if else: at least 1 row?

    } # end for: ord


    # b) Take out time point of: thresholds and sds

    thresholds_ts <- list()
    for(v in 1:p) thresholds_ts[[v]] <- as.numeric(thresholds[[v]][ts, ])
    sds_ts <- sds[ts, ]

    # c) factors stays the same
    factors_ts <- factors


    # ---------- Call mgmsampler() ----------

    one_case <- mgmsampler(factors = factors_ts,
                             interactions = interactions_ts,
                             thresholds = thresholds_ts,
                             sds = sds_ts,
                             type = type,
                             level = level,
                             N = 1,
                             nIter = nIter,
                             pbar = FALSE)

    data[ts, ] <- one_case$data



    if(pbar==TRUE) setTxtProgressBar(pb, ts)

  } # end for: time steps

  # ---------- Output ----------

  tvmgmsamp_obj$data <- data

  class(tvmgmsamp_obj) <- 'tvmgm'

  return(tvmgmsamp_obj)


} # eoF














