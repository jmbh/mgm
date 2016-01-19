\name{mgm-package}
\alias{mgm-package}
\alias{mgm}
\docType{package}
\title{
Estimation of Mixed Graphical Models
}
\description{
Estimation of Mixed Graphical Models using L1-constrained neighborhood regression. 
}
\details{
\tabular{ll}{
Package: \tab mgm\cr
Type: \tab Package\cr
Version: \tab 1.1-4 \cr
Date: \tab 2016-01-06\cr
License: \tab ? \cr
}
}
\author{
Jonas Haslbeck

Maintainer: <jonashaslbeck@gmail.com>
}



\references{
Haslbeck, J., & Waldorp, L. J. (2015). Structure estimation for mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.05677.

Loh, P. L., & Wainwright, M. J. (2013). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. The Annals of Statistics, 41(6), 3022-3049.

Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).
}
\keyword{ package }
\examples{

\dontrun{
#create some data from a pairwise mixed graphical model
n <- 100 # number of samples
type <- c("g", "c", "c") # one gaussian, two categorical
lev <- c(1,3, 2) # the categorical variable have 3 and 2 categories, respectively
graph <- matrix(0,3,3) 
graph[1,2] <- graph[2,1] <- .5 # we have an edge with weight .5 between node 1 and 2
thresh <- list(c(0), c(0,0,0), c(0,0)) # all thresholds are zero 
data <- mgmsampler(n, type, lev, graph, thresh, parmatrix=NA, nIter=1000)

#fit a mixed graphical model
fit <- mgmfit(data, type, lev=lev, d=2)
grapht <- graph
grapht[grapht!=0] <- 1 #binarize true graph

confusion(grapht, fit$adj) #evaluate performance
}

}
