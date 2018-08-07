\name{FactorGraph}
\alias{FactorGraph}

\title{
Draws a factor graph of a (time-varying) MGM
}
\description{
Wrapper function around qgraph() that draws factor graphs for (time-varying) MGMs
}
\usage{
FactorGraph(object, colors, labels, PairwiseAsEdge = FALSE, 
            Nodewise = FALSE, DoNotPlot = FALSE, 
            FactorLabels = TRUE, shapes = c("circle", "square"), 
            shapeSizes = c(8, 4), estpoint = NULL, 
            negDashed = FALSE, ...)
}


\arguments{
  \item{object}{
The output object of \code{mgm()} or \code{tvmgm()}.
}
  \item{colors}{
A character vector of colors for nodes and factors. The first color is for nodes, the second for 2-way interactions, the third for 3-way interactions, etc. Defaults to \code{colors = c("white", "tomato", "lightblue")}.
}
  \item{labels}{
A character vector of (variable) node labels.
}
  \item{PairwiseAsEdge}{
If \code{TRUE}, pairwise interactions are not displayed as factors but as simple edges between nodes. Defaults to \code{PairwiseAsEdge = FALSE}.
}
  \item{Nodewise}{
If \code{TRUE}, the estimates from the individual nodewise regressions are displayed as a directed edge towards the node on which the respective nodewise regression was performed. This is useful to identify model misspecification (e.g. moderation effects / interaction parameters with largely different values across nodewise regressions). Defaults to \code{Nodewise = FALSE}.
}
  \item{DoNotPlot}{
If \code{DoNotPlot = TRUE} no factorgraph is plotted. This way the computed factor graph can be obtained without plotting. Defaults to \code{DoNotPlot = FALSE}.
}
  \item{FactorLabels}{
If \code{FactorLabels = TRUE} the factors are labeled by their order. If \code{FactorLabels = FALSE} no label is shown. Defaults to \code{FactorLabels = TRUE}.
}
  \item{shapes}{
A character vector of length two indicating the shapes for nodes and factors. Defaults to \code{shapes = c("circle", "square")}, which gives circles to nodes and squares to factors. See \code{?qgraph} for available shapes.
}
  \item{shapeSizes}{
A numeric vector of length two indicating the size of shapes for nodes and factors. Defaults to \code{shapeSizes = c(8, 4)}.
}
  \item{estpoint}{
An integer indicating the estimation point to display if the output object of a time-varying MGM is provided. 
}
  \item{negDashed}{
If \code{negDashed = TRUE}, edges with negative sign are dashed.
}
  \item{\dots}{
Arguments passed to qgraph.
}
}
\details{
\code{FactorGraph()} is a wrapper around \code{qgraph()} from the qgraph package. Therefore all arguments of \code{qgraph()} are available and can be provided as additional arguments.

To make time-varying factor graphs comparable across estimation points, the factor graph of each estimation point includes all factors that are estimated nonzero at least at one estimation point.
 
}
\value{
Plots the factor graph and returns a list including the arguments used to plot the factor graph using qgraph().

Specifically, a list is returned including:  \code{graph} contains a weighted adjacency matrix of a (bipartide) factor graph. If p is the number of variables and E the number of interactions (factors) in the model, this matrix has dimensions (p+E) x (p+E). The factor graph is furter specified by the following objects: \code{signs} is a matrix of the same dimensions as \code{graph} that indicates the sign of each interaction, if defined (see \code{pairwise} above). \code{edgecolor} is a matrix with the same dimension as \code{graph} that provides edge colors depending on the sign as above. \code{order} is a (p+E) vector indicating the order of interaction. The first p entries are set to zero. \code{qgraph} contains the qgraph object created while plotting.

}
%%\references{
%% ~put references to the literature/web site here ~
%%}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\seealso{
\code{mgm()}, \code{tvmgm()}, \code{qgraph()}
}


\examples{


\dontrun{

# Fit MGM with pairwise & threeway interactions to Autism Dataset
fit_k3 <- mgm(data = autism_data$data,
              type = autism_data$type,
              level = autism_data$lev,
              k = 3, 
              overparameterize = TRUE,
              lambdaSel = "EBIC", 
              lambdaGam = .5) 

# List of estimated interactions
fit_k3$interactions$indicator

# Plot Factor Graph 
nodeColors <- c("white", "orange", "lightblue")

FactorGraph(object = fit_k3, 
            PairwiseAsEdge = FALSE, 
            DoNotPlot = FALSE, 
            colors = nodeColors, 
            labels = 1:7, 
            layout="spring")

# For more examples see https://github.com/jmbh/mgmDocumentation
}

}


