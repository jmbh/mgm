



mgmsampler <- function(n, #number of samples
                       type, #type of data/from which distribution comes the data
                       lev, # Number of levels
                       graph, #graph structure
                       thresh, #thresholds, for every node (& category)
                       parmatrix = NA, #possibility to provide costum function to create model parameter matrix
                       nIter = 250, #number of samples for each node
                       varadj = .2, # additive constant to conditional variances before normalization; avoids partial correlations close to 1,
                       exportGraph = FALSE # if TRUE exports list with data, graph with (transformed) gaussian sub graph and the parameter matrix, if FALSE just exports the data 
){

  # --------- Input Checks -------
    
  lev <- as.numeric(lev)

  stopifnot(length(type)==length(lev))
  if(is.na(parmatrix)==TRUE){
    stopifnot(isSymmetric(graph))
    stopifnot(length(type)==nrow(graph))
    stopifnot(sum(diag(graph)!=0)==0)
  }
  
  
  # --------- Ensure that the Gaussian Submatrix is Positive Definite -------
  
  # in case we have more than 1 gaussian variable in the graph, 
  # we have to rescale the conditional variances to 1
  if(sum("g"==type)>1) {
    
    # get the gaussian subgraph
    graph.g <- graph[type=="g", type=="g"]
    
    # all parameters have to be negative
    graph.g <- -abs(graph.g)
    
    # define diagonal in inverse covariance matrix
    diag(graph.g) <- - (colSums(graph.g) - rep(varadj, sum(type=="g"))) 
    
    #rescale to unit conditional variance
    graph.g.re <- round(cov2cor(graph.g),10) 
    
    # check whether Gaussian submatrix is now positive definite
    Eig <- eigen(graph.g.re)
    if(sum(Eig$values>0)!=ncol(graph.g)) stop('Gaussian Submatrix not positive definite.')
    
    graph.g.re.e <- graph.g.re; diag(graph.g.re.e) <- 0 #zero diagonal to put it back in the graph
    graph[type=="g", type=="g"] <- graph.g.re.e #put back in the graph
    
  } 
  
  # --------- Create Parameter Matrix Used in Gibbs Sampler -------
  
  if(is.na(parmatrix)==TRUE) {
    graphe <- potts_parameter(graph, type, lev, thresh) #create model.parameter.matrix
  } else {
    stopifnot(isSymmetric(parmatrix))
    graphe <- parmatrix
  }
  
  nNodes <- ncol(graph)  # number of nodes
  Data <- matrix(0, n, nNodes)  # create empty data matrix
  inde <- as.numeric(colnames(graphe))
  colnames(graphe) <- rownames(graphe) <- NULL
  
  #transform thresh into a matrix (C doesnt take lists)
  thresh_m <- matrix(0, nrow=length(thresh), ncol=max(lev))
  for(t in 1:nNodes) {
    thresh_m[t,1:length(thresh[[t]])] <- thresh[[t]]
  }
  
  #transform types into integers
  type_c <- numeric(nNodes)
  type_c[type=="c"] <- 1
  type_c[type=="g"] <- 2
  type_c[type=="p"] <- 3
  type_c[type=="e"] <- 4
  
  # --------- C-implementation of Gibbs Sampler -------

  c_out <- mMRFCsampler(Data, n, nNodes, type_c, lev, nIter=nIter, thresh_m, graphe, inde)
  
  
  
  # --------- Different Output Options -------
  
  if(exportGraph) {
    
    outlist <- list('Data' = c_out,
                    'Graph' = graph,
                    'ParMatrix' = graphe)
    
    return(outlist)
    
  } else {
    
    return(c_out)
    
  }
  
  
} # end of function
