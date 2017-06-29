
lagData <- function(data, lags, consec = NULL) {


  # ---------- Compute Aux Variables ----------

  max_lag <- max(lags) # maximum lag
  lags_ext <- 1:max(lags) # makes it easier to delete right columns

  n <- nrow(data)
  p <- ncol(data)
  n_var <- nrow(data) - max(lags)
  n_lags <- length(lags_ext)

  data_response <- data[-c(1:max_lag), ]

  if(!is.null(consec)) m_consec <- matrix(NA, nrow = n_var, ncol = n_lags)

  # browser()

  # ---------- Lag Variables ----------

  lag_pos <- 1 # to make sure that the list is filled successively, if not a full sequence (e.g. lags = c(1,5)); otherwise this leads to problems later in the code

  l_data_lags <- list()
  for(lag in lags) {

    front <- max_lag - lag
    end <- max_lag - front

    if(lag == max_lag) {

      l_data_lags[[lag_pos]] <- data[-((n-end+1):n), ]
      colnames(l_data_lags[[lag_pos]]) <- paste("V", 1:p, '.lag', lag, '.', sep = "")

    } else {

      l_data_lags[[lag_pos]] <- data[-c((1:front),((n-end+1):n)), ]
      colnames(l_data_lags[[lag_pos]]) <- paste("V", 1:p, '.lag', lag, '.', sep = "")

    }

    lag_pos <- lag_pos + 1
  }


  # ---------- Knock Out if not consecutive ----------

  if(!is.null(consec)) {

    for(lag in lags_ext) {

      front <- max_lag - lag
      end <- max_lag - front

      if(lag == max_lag) {
        m_consec[,lag] <- consec[-((n-end+1):n)]
      } else {
        m_consec[,lag] <- consec[-c((1:front),((n-end+1):n))]
      }
    }

    # calculate which cases to knock out
    m_consec_check <- cbind(consec[-c(1:max_lag)], m_consec)

    v_check <- apply(m_consec_check, 1, function(x) {

      check_row <- x[1] - x[-1] == 1:length(x[-1]) # check for extended lags 1:max(lags)
      check_row_relevant <- check_row[lags_ext %in% lags] # but then compute check only over the lags that are actually specified
      if(any(check_row_relevant == FALSE)) FALSE else TRUE # and return test result: any required previous measurement missing? => FALSE

    })

  }


  # ---------- Output ----------
  
  outlist <- list()
  outlist$data_response <- data_response
  outlist$l_data_lags <- l_data_lags

  if(!is.null(consec)) {
    outlist$included <- v_check
  } else {
    outlist$included <- rep(TRUE, n)
    outlist$included[lags_ext] <- FALSE
  }


  return(outlist)

}
