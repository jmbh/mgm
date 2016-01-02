
mgmsampler <- function(
  n, #number of samples
  type, #type of data/from which distribution comes the data
  lev, # Number of levels
  graph, #graph structure
  thresh, #thresholds, for every node (& category)
  parmatrix = NA, #possibility to provide costum function to create model parameter matrix
  nIter = 250 #number of samples for each node
){
  
  lev <- as.numeric(lev)
  # some checks on the input
  stopifnot(length(type)==length(lev))
  if(is.na(parmatrix)==TRUE){
    stopifnot(isSymmetric(graph))
    stopifnot(length(type)==nrow(graph))
    stopifnot(sum(diag(graph)!=0)==0)
  }
  
  # check: is the within-gaussian covariance matrix positive definite?
  g_cov <- graph[type=="g",type=="g"] #select within-gaussian cov matrix
  
  if(length(g_cov)>1)
  {
    g_covm <- matrix(g_cov, nrow(g_cov), ncol(g_cov))
    diag(g_covm) <- 1 # we are working with unit variance
    stopifnot(is.positive.definite(g_covm))
  }
  
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
  
  #CALL C CORE
  c_out <- mMRFCsampler(Data, n, nNodes, type_c, lev, nIter=nIter, thresh_m, graphe, inde)
  
  return(c_out)
  
} # end of function
