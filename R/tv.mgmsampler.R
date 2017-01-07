

tv.mgmsampler <- function(type, # p time vector
                          lev, # p level vector
                          graphs, # p x p x timestep array
                          threshs, # n list with each p threshold entries (same structure as in mgmsampler)
                          parmatrices = NA,
                          nIter = 250, 
                          varadj = .2,
                          exportGraph = FALSE # if TRUE exports list with data, graph with (transformed) gaussian sub graph and the parameter matrix, if FALSE just exports the data 
                          ) {
  
  # ----- Input Type: Edge weight or parameter matrix? -----
  if(!is.na(parmatrices)) {
    graphs <- NA
    # basic input info
    n <- dim(parmatrices)[3]
    p <- length(type)
  } else {
    # basic input info
    n <- dim(graphs)[3]
    p <- dim(graphs)[2]
    # input checks
    if(sum(!apply(graphs, 3, matrixcalc::is.symmetric.matrix))>0) stop('The weight matrix must be symmetric for every time step.')
    if(length(threshs)!=n) stop('A threshold for each time step has to be specified.')

  }
  
  # ----- Sampling -----
  data <- matrix(NA, n, p)
  l_dlist <- vector('list', n)
  
  for(timestep in 1:n) {
    
    if(exportGraph) {
      
      
      if(!is.na(parmatrices)) {
        l_dlist[[timestep]] <- mgmsampler(n=1, type=type, lev=lev, graph=NA, 
                                      thresh = threshs[[timestep]], parmatrix = parmatrices[[timestep]], 
                                      nIter = nIter, varadj = varadj, exportGraph = exportGraph)
      } else {
        l_dlist[[timestep]] <- mgmsampler(n=1, type=type, lev=lev, graph=graphs[,,timestep], 
                                      thresh = threshs[[timestep]], parmatrix = NA, 
                                      nIter = nIter, varadj = varadj, exportGraph = exportGraph)
      }
      
      
      
    } else {
      
      
      if(!is.na(parmatrices)) {
        data[timestep,] <- mgmsampler(n=1, type=type, lev=lev, graph=NA, 
                                      thresh = threshs[[timestep]], parmatrix = parmatrices[[timestep]], 
                                      nIter = nIter, varadj = varadj, exportGraph = exportGraph)
      } else {
        data[timestep,] <- mgmsampler(n=1, type=type, lev=lev, graph=graphs[,,timestep], 
                                      thresh = threshs[[timestep]], parmatrix = NA, 
                                      nIter = nIter, varadj = varadj, exportGraph = exportGraph)
      }
      
    }
    
  } # end time for-loop
  
  
  # ----- Data Export -----
  
  if(exportGraph) { 
    
    for(timestep in 1:n) data[timestep,] <- l_dlist[[timestep]]$Data
    l_graphs <- l_parmats <- vector('list', n)
    for(timestep in 1:n) l_graphs[[timestep]] <- l_dlist[[timestep]]$Graph
    for(timestep in 1:n) l_parmats[[timestep]] <- l_dlist[[timestep]]$ParMatrix
    
    outlist <- list('Data' = data,
                    'Graphs' = l_graphs,
                    'ParMatrices' = l_parmats)
    
    return(outlist)
    
  } else {
    
    return(data)
    
  }
  
  
}






