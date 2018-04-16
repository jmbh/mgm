# jonashaslbeck@gmail.com; March 2016

FactorGraph <- function(object, 
                        colors, 
                        labels,
                        PairwiseAsEdge = TRUE, # box
                        DoNotPlot = FALSE, 
                        FactorLabels = TRUE, 
                        shapes = c('circle', 'square'), 
                        shapeSizes = c(8, 4), 
                        estpoint = NULL,
                        ...)
  
{
  
  
  
  # --------- Compute Aux Variables ---------
  
  p <- length(object$call$level)
  n_estpoints <- length(object$call$estpoints)
  
  # --------- Input Checks ---------
  
  if(!missing(labels)) if(length(labels) != p) stop("Number of provided labels has to match the number of variables.")
  
  # Checks for time-varying FactorGraph  
  if("tvmgm" %in% class(object)) {
    if(missing(estpoint)) stop("Specify the estimation point for which the factor graph should be visualized.")
    if(estpoint > n_estpoints) stop(paste0("The provided fit object has only ", n_estpoints, " estimation points."))
  } 
  
  
  # --------- Create FractorGraph object ---------
  
  call <- list("object" = object)
  
  FG_object <- list("call" = call,
                    "graph" = NULL,
                    "nodetype" = NULL,
                    "order" = NULL,
                    "signs" = NULL,
                    "edgecolor" = NULL, 
                    "qgraph" = NULL)
  
  # --------- Fill in defaults ---------
  
  if(missing(labels)) labels <- 1:p
  if(missing(colors)) colors <- c("white", "tomato", "lightgrey")
  layout <- "spring"
  cut <- 0
  
  
  # --------- Compute Factor Graph ----------
  
  # Call different DrawFG() version for stationary/time-varying
  if("tvmgm" %in% class(object)) {
    
    # Time-varying
    FG <- DrawFGtv(object = object,
                   PairwiseAsEdge = PairwiseAsEdge, 
                   estpoint = estpoint)
    
  } else {
    
    # Stationary
    FG <- DrawFG(object = object,
                 PairwiseAsEdge = PairwiseAsEdge)
    
  }
  
  
  # Save into FG_object
  FG_object$graph <- FG$weightedgraph
  FG_object$nodetype <- FG$nodetype
  FG_object$order <- FG$order
  FG_object$signs <- FG$signs
  FG_object$edgecolor <- edge.color <- FG$signcolor
  
  
  
  # Allow overwriting ...
  args <- list(...)
  if(!is.null(args$cut)) cut <- args$cut
  if(!is.null(args$layout)) layout <- args$layout
  if(!is.null(args$edge.color)) edge.color <- args$edge.color
  
  
  # --------- Plot & Return ---------
  
  if(!DoNotPlot){
    
    # ----- Compute stuff necessary for plotting -----
    
    # Create labels for factors (label = order of factor/interaction)
    ifelse(PairwiseAsEdge, ek <- 1, ek <- 0)
    if(FactorLabels) {
      tb <- table(FG_object$order)[-1]
      l_lf <- list()
      for(k in 1:length(tb)) l_lf[[k]] <- rep(k+1+ek, tb[k])
      labels_ex <- c(labels, unlist(l_lf))
    } else {
      labels_ex <- c(labels, rep('', sum(FG_object$nodetype)))
    }
    
    # ----- Call qgraph -----
    
    # browser()
    
    qgraph_object <- qgraph(FG_object$graph,
                            color = colors[FG_object$order + 1],
                            edge.color = edge.color,
                            layout = layout,
                            labels =  labels_ex,
                            shape = shapes[FG_object$nodetype + 1],
                            vsize = shapeSizes[FG_object$nodetype + 1], 
                            cut = cut,
                            ...)
    
    FG_object$qgraph <- qgraph_object
    
    
    invisible(FG_object) # return output object invisible
    
  } else {
    return(FG_object)
  }
  
  
  
} # eoF