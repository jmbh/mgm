\name{mgmsampler}
\alias{mgmsampler}

\title{
Samples from a mixed Graphical Model}
\description{
Samples from a pairwise mixed Graphical Model}
\usage{
mgmsampler(n, type, lev, graph, thresh, parmatrix = NA, nIter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of samples
}
  \item{type}{
 x 1 vector specifying the type of distribution for each variable ("g" = Gaussian, "p" = Poisson, "e" = Exponential, "c" = Categorical)
}
  \item{lev}{
p x 1 vector specifying the number of levels for each variables (for continuous variables: 1)
}
  \item{graph}{
A p x p symmetric (weighted) adjacency matrix
}
  \item{thresh}{
A list in which each entry corresponds to a variable; categorical variables have as many thresholds as categories
}
  \item{parmatrix}{
Optional; A matrix specifying all parameters in the model. This provides a possibility to exactly specify all parameters instead of using the default mapping from the weighted adjacency matrix to the model parameter matrix as described below (Details). If provided, \code{graph} will be ignored.
}
  \item{nIter}{
Iterations in the Gibbs sampler
}
}
\details{
In case of interactions involving categorical variables with more than m=2 categories, there are more than one one parameter describing the interaction. As we use the overcomplete representation, in case of an interaction between two categorical variables, we have m^2 parameters, in case of an interaction between a cateogircal and a continuous variable we have m parameters. In case of Catogorical-Categorical interactions, we use the Potts-parameterization (see e.g. Wainwright & Jordan, 2009), where the value of all non-zero parameters is given bei the weighted adjacency matrix. In case of Categorical-Continous interactions, parameters for categories m > |m|/2 have a nonzero parameter with the value is given by the weighted adjacency matrix. All categories m >=  |m|/2 have zero parameters. Please note that this mapping is absolutely arbitrary. In the offical release of the package, the user will be able to specify the mapping herself.
}
\value{
n x p data matrix
}
\references{
Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).

Wainwright, M. J., & Jordan, M. I. (2008). Graphical models, exponential families, and variational inference. Foundations and Trends in Machine Learning, 1(1-2), 1-305.
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\examples{

n <- 10 # number of samples
type <- c("g", "c") # one Gaussian, one categorical
lev <- c(1,3) # the categorical variable has 3 categories
graph <- matrix(0,2,2) 
graph[1,2] <- graph[2,1] <- .5 # we have an edge with weight .5 between the two nodes
thresh <- list(c(0), c(0,0,0)) # all thresholds are zero 
data <- mgmsampler(n, type, lev, graph, thresh, parmatrix=NA, nIter=1000)
head(data)


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
