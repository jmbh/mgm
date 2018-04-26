# jonashaslbeck@gmail.com; March 2016

FactorGraph <- function(object, 
                        colors, 
                        labels,
                        PairwiseAsEdge = TRUE, # box
                        Nodewise = FALSE,
                        DoNotPlot = FALSE, 
                        FactorLabels = TRUE, 
                        shapes = c('circle', 'square'), 
                        shapeSizes = c(8, 4), 
                        estpoint = NULL,
                        ...)
  
{
  
  # --------- Compute Aux Variables ---------
  
  if(Nodewise) PairwiseAsEdge <- FALSE
  
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
                    "nonzero" = NULL,
                    "qgraph" = NULL)
  
  # --------- Fill in defaults ---------
  
  if(missing(labels)) labels <- 1:p
  if(missing(colors)) colors <- c("white", "tomato", "lightblue")
  layout <- "spring"
  cut <- 0
  
  
  # --------- Compute Factor Graph ----------
  
  # Call different DrawFG() version for stationary/time-varying
  if("tvmgm" %in% class(object)) {
    
    # Time-varying
    FG <- DrawFGtv(object = object,
                   PairwiseAsEdge = PairwiseAsEdge, 
                   Nodewise = Nodewise,
                   estpoint = estpoint)
    
  } else {
    
    # Stationary
    FG <- DrawFG(object = object,
                 PairwiseAsEdge = PairwiseAsEdge, 
                 Nodewise = Nodewise)
    
  }
  
  # Save into FG_object
  FG_object$graph <- FG$weightedgraph
  FG_object$nodetype <- FG$nodetype
  FG_object$order <- FG$order
  FG_object$signs <- FG$signs
  FG_object$edgecolor <- edge.color <- FG$signcolor
  FG_object$nonzero <- FG$nonzero
  
  
  
  # Allow overwriting ...
  args <- list(...)
  if(!is.null(args$cut)) cut <- args$cut
  if(!is.null(args$layout)) layout <- args$layout
  if(!is.null(args$edge.color)) edge.color <- args$edge.color
  
  
  # browser()
  
  # Adapt edge labels for zero edges in Nodewise=TRUE
  if(!is.null(args$edge.labels)) { # if specified, otherwise set to FALSE
    if(is.logical(args$edge.labels)) { # if specified and logical, then adapt for nonzero or FALSE
      if(args$edge.labels) {
        edge.labels <- FG_object$graph
        edge.labels[FG_object$nonzero == 2] <- 0
        edge.labels <- round(edge.labels, 2)
      } else {
        edge.labels = FALSE
      }
    } else {
      # if not logical, take the input
      edge.labels <- args$edge.labels
    }
  } else {
    edge.labels = FALSE
  }
  
  
  # --------- Plot & Return ---------
  
  if(!DoNotPlot){
    
    # ----- Compute stuff necessary for plotting -----
    
    # Create labels for factors (label = order of factor/interaction)
    ifelse(PairwiseAsEdge, ek <- 1, ek <- 0)
    if(FactorLabels) {
      tb <- table(FG_object$order)[-1]
      
      if(length(tb)==0) { # For the case PairwiseAsEdge=FALSE and no 3-way interactions
        FL <- NULL
      } else {
        l_lf <- list()
        for(k in 1:length(tb)) l_lf[[k]] <- rep(k+1+ek, tb[k])
        FL <- unlist(l_lf)
      }
      
      labels_ex <- c(labels, FL)
    } else {
      labels_ex <- c(labels, rep('', sum(FG_object$nodetype)))
    }
    
    # ----- Call qgraph -----
    
    qgraph_object <- qgraph(FG_object$graph,
                            color = colors[FG_object$order + 1],
                            edge.color = edge.color,
                            lty = FG_object$nonzero,
                            layout = layout,
                            labels =  labels_ex,
                            shape = shapes[FG_object$nodetype + 1],
                            vsize = shapeSizes[FG_object$nodetype + 1], 
                            edge.labels = edge.labels,
                            cut = cut,
                            ...)
    
    FG_object$qgraph <- qgraph_object
    
    
    invisible(FG_object) # return output object invisible
    
  } else {
    return(FG_object)
  }
  
  
  
} # eoF