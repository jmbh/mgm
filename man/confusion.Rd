\name{confusion}
\alias{confusion}
\title{
Confusion matrix}
\description{
Computes accuracy, sensitivity, specificity and precision from a true and an estimated adjacency matrix.}
\usage{
confusion(tg, eg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tg}{
p x p adjacenty matrix of the true graph
}
  \item{eg}{
p x p adjacenty matrix of the estimated graph}
}
\value{
A list containing accuracy, sensitivity, specificity and precision of the estimated graph in reference to the true graph.}
\author{
Jonas Haslbeck}

\examples{

#create some data from a mixed graphical model
n <- 100 # number of samples
type <- c("g", "c", "c") # one gaussian, two categorical
lev <- c(1,3, 2) # the categorical variable have 3 and 2 categories, respectively
graph <- matrix(0,3,3) 
graph[1,2] <- graph[2,1] <- .5 # we have an edge with weight .5 between node 1 and 2
thresh <- list(c(0), c(0,0,0), c(0,0)) # all thresholds are zero (for the categorical variables each categoriy has a threshold)
data <- mgmsampler(n, type, lev, graph, thresh, parmatrix=NA, nIter=1000)

#fit a mixed graphical model
fit <- mgmfit(data, type, lev=lev, d=2)
grapht <- graph
grapht[grapht!=0] <- 1 #binarize true graph

confusion(grapht, fit$adj) #compute derivatives of confusion matrix


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
