
lagData <- function(data, 
                    lags, 
                    consec = NULL) {


  # ---------- Compute Aux Variables ----------

  data <- as.matrix(data) # turn into matrix
  
  max_lag <- max(lags) # maximum lag
  lags_ext <- 1:max(lags) # makes it easier to delete right columns

  n <- nrow(data)
  p <- ncol(data)
  n_var <- nrow(data) - max(lags)
  n_lags <- length(lags_ext)

  data_response <- data #[-c(1:max_lag), ]

  if(!is.null(consec)) m_consec <- matrix(NA, 
                                          nrow = n, 
                                          ncol = n_lags)

  # ---------- Lag Variables ----------

  # Storage
  l_data_lags <- list()

  # Loop through lags
  lag_pos <- 1 # to make sure that the list is filled successively, if not a full sequence (e.g. lags = c(1,5)); otherwise this leads to problems later in the code  
  for(lag in lags) {
    
    lagged_data <- matrix(NA, nrow = n, ncol=p)
    lagged_data[(lag+1):n, ] <- data[-((n-lag+1) : n), ]
    lagged_data <- matrix(lagged_data, 
                          ncol = p, 
                          nrow = n) 
    colnames(lagged_data) <- paste("V", 1:p, '.lag', lag, '.', sep = "")
    
    l_data_lags[[lag_pos]] <- lagged_data

    lag_pos <- lag_pos + 1
  }

  # ---------- Knock Out if not consecutive ----------

  
  # browser()
  
  if(!is.null(consec)) {

    for(lag in lags_ext)  m_consec[(lag+1):n, lag] <- consec[-((n-lag+1) : n)]

    # Calculate which cases to knock out
    m_consec_check <- cbind(consec, m_consec)
    
    v_check <- apply(m_consec_check, 1, function(x) {
      
      if(any(is.na(x))) {
        FALSE
      } else {
        check_row <- x[1] - x[-1] == 1:length(x[-1]) # check for extended lags 1:max(lags)
        check_row_relevant <- check_row[lags_ext %in% lags] # but then compute check only over the lags that are actually specified
        if(any(check_row_relevant == FALSE)) FALSE else TRUE # and return test result: any required previous measurement missing? => FALSE
      }

    })

  } else {
    
    v_check <- rep(TRUE, n)
    v_check[1:n_lags] <- FALSE
    
  }


  # ---------- Output ----------
  
  outlist <- list()
  outlist$data_response <- data_response
  outlist$l_data_lags <- l_data_lags
  outlist$included <- v_check


  return(outlist)

}
