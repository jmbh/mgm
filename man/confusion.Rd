\name{confusion}
\alias{confusion}
\title{
Confusion matrix}
\description{
Computes accuracy, sensitivity, specificity and precision from a true and an estimated adjacency matrix.}
\usage{
confusion(tg, eg, undirected = TRUE, diagonal = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{tg}{p x p adjacency matrix of the true graph}

  \item{eg}{p x p adjacency matrix of the estimated graph}

  \item{undirected}{Logical indicating whether the true and estimated graphs are undirected or directed.}

  \item{diagonal}{Logical indicating whether the diagonal should be taken into account. \code{diagonal} defaults to \code{FALSE} if \code{undirected = TRUE}.}

}

\value{A list containing accuracy, sensitivity, specificity and precision of the estimated graph with respect to the true graph.}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>}

\examples{

\dontrun{

# create some data from a mixed graphical model

n <- 100 # number of samples
type <- c("g", "c", "c") # one gaussian, two categorical
lev <- c(1, 3, 2) # the categorical variables have 3 and 2 categories, resp
graph <- matrix(0, now=3, ncol=3) 
graph[1, 2] <- graph[2, 1] <- .5 # edge weight .5
thresh <- list(c(0), c(0, 0, 0), c(0, 0))
data <- mgmsampler(n, type, lev, graph, thresh, parmatrix=NA, nIter=1000)

# fit a mixed graphical model
fit <- mgmfit(data, type, lev=lev, d=2)
grapht <- graph
grapht[grapht!=0] <- 1 # binarize true graph

# compute derivatives of confusion matrix
confusion(grapht, fit$adj)

}

}
