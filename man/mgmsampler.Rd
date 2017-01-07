\name{mgmsampler}
\alias{mgmsampler}

\title{
Samples from a mixed Graphical Model}
\description{
Samples from a pairwise Mixed Graphical Model}
\usage{
mgmsampler(n, type, lev, graph, thresh, parmatrix = NA, 
           nIter = 250, varadj = .2, exportGraph = FALSE)
}

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
The number of iterations in the Gibbs sampler
}
  \item{varadj}{
Additive constant to conditional Gaussian variances before normalization. This avoids partial correlations close to 1.
}
  \item{exportGraph}{
If \code{FALSE} exports only the data. If \code{TRUE} exports a list containing the data, the graph with a potentially transformed Gaussian submatrix and the explicit parameter matrix used in the Gibbs sampler.
}
}
\details{
For interactions involving categorical variables with more than m = 2 categories, an interaction is comprised by more than one parameter. Specifically, in the overcomplete representation used here, in case of an interaction between two categorical variables with each m categories, we have m^2 parameters, whereas in case of an interaction between a cateogircal and a continuous variable we have m parameters. 

We use the following (arbitrary) mapping from provided edge-weights to interaction parameters: for catogorical-categorical interactions, we use the parameterization of the Potts model (see e.g. Wainwright & Jordan, 2008) extended to the case of interactions between categorical variables with different number of categories. All non-zero elements in the parameter matrix are assigned the value specified in the weighted adjacency matrix. In case of Categorical-Continous interactions, parameters for categories m > |m|/2 have a nonzero parameter with the value is given by the weighted adjacency matrix. All categories m >=  |m|/2 have zero parameters. A customized parameterization of interactions including categorical variables can be specified using the \code{parmatrix} argument.

The specification of the mixed graphical model has to satisfy a number of contraints. For details see Yang et al. (2014) or Haslbeck and Waldorp (2015).

}
\value{
If \code{exportGraph = FALSE} a n x p data matrix. If \code{exportGraph = TRUE} in addition the graph with a potentially transformed Gaussian submatrix and the explicit parameter matrix used in the Gibbs sampler.
}
\references{
Haslbeck, J., & Waldorp, L. J. (2015). Structure estimation for mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.06871.

Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).

Wainwright, M. J., & Jordan, M. I. (2008). Graphical models, exponential families, and variational inference. Foundations and Trends in Machine Learning, 1(1-2), 1-305.
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\examples{

n <- 10 # number of samples
type <- c("g", "c") # one Gaussian, one categorical
lev <- c(1, 3) # the categorical variable has 3 categories
graph <- matrix(0, 2, 2) 
graph[1, 2] <- graph[2, 1] <- .5 # we have an edge with weight .5 between the two nodes
thresh <- list(c(0), c(0, 0, 0)) # all thresholds are zero 

data <- mgmsampler(n, type, lev, graph, thresh, parmatrix = NA, nIter = 1000)

head(data)


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
