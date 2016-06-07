\name{tv.mgmsampler}
\alias{tv.mgmsampler}

\title{
Sampling from a time-varying mixed Graphical Model}
\description{
Samples from a time-varying pairwise Mixed Graphical Model}
\usage{
tv.mgmsampler(type, lev, graphs, threshs, parmatrices = NA, nIter = 250, varadj = .2)
}

\arguments{
  \item{type}{
p (variables) x 1 vector specifying the type of distribution for each variable ("g" = Gaussian, "p" = Poisson, "e" = Exponential, "c" = Categorical)
}
  \item{lev}{
p x 1 vector specifying the number of levels for each variables (for continuous variables = 1)
}
  \item{graphs}{
A p x p x n (time steps) symmetric (weighted) adjacency matrix, where n are the number of time steps, i.e. a graph is specified for each time step.
}
  \item{threshs}{
A list in which each entry corresponds to a time step. For each time step there is a list in which each entry corresponds to one variable (as \code{thresh} in \code{mgmsampler}.)
}
  \item{parmatrices}{
Optional; Equivalent to \code{parmatrix} in \code{mgmsampler}, only with a parameter matrix provided for each time step \code{n}.  If provided, \code{graphs} will be ignored.
}
  \item{nIter}{
Number of iterations in the Gibbs sampler
}
  \item{varadj}{
Additive constant to conditional Gaussian variances before normalization. This avoids partial correlations close to 1.
}
}
%\details{
%
%}
\value{
n x p data matrix.
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

\dontrun{
# We sample from a graph with 400 time steps.
# In the first 200 steps the underlying graph is graph A, 
# and in the last 200 steps graph B

  p <- 4 # 4 variables
  n <- 400 # number of time steps
  graphA <- graphB <- matrix(0, p, p)
  graphA[2,1] <- graphA[1,2] <- 1 # graphA has edge (2,1)
  graphB[3,4] <- graphB[4,3] <- 1 # graphB has edge (3,4)
  graphs <- array(dim=c(p, p, n)) 
  graphs[,,1:(n/2)] <- graphA
  graphs[,,(n/2+1):n] <- graphB
  
  # specify type and thresholds
  type <- c('g', 'g', 'c', 'c') # two Gaussian, two Binary
  lev <- c(1, 1, 2, 2)
  # zero means (Gaussians) and thresholds (Binary)
  thresh <- list(0, 0 , c(0, 0), c(0, 0))
  # same means/thresholds at all time steps
  threshs <- list(); for(nIter in 1:n) threshs[[nIter]] <- thresh

  # sample
  set.seed(1)
  data <- tv.mgmsampler(type, lev, graphs, threshs)
  
  # look at sampled data
  dim(data)
  head(data)
}


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
