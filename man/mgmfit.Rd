\name{mgmfit}
\alias{mgmfit}

\title{
Estimation of Mixed Graphical Models
}
\description{
Estimation of Mixed Graphical Models using L1-constrained neighborhood regression.}
\usage{
mgmfit(data, type, lev, lambda.sel="EBIC", folds = 10, gam = .25, 
       d=1, rule.reg = "AND", rule.cat = "OR", pbar = "TRUE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x p data matrix
}
  \item{type}{
p-vector containing the types of variable ("g" Gaussian, "e" Exponential, "p" Poisson, "c" Categorical)
}
  \item{lev}{
p-vector containing the number of levels for each variable p (lev=1 for all continuous variables)
}

  \item{lambda.sel}{
Procedure to select the lambda-parameter for L1-penalized regressions. The two options are cross validation "CV" and "EBIC" (default). While cross-validation is used in the paper (see below), the Extended Bayesian Information Criterion (EBIC, see Barber et al., 2015) is a useful alternative in cases where the cross-validation requirement that all categories are present in each fold is not met.
}
  \item{folds}{
The number of folds to be used for cross-validation in case lambda.sel="CV". Defaults to folds=10.}
  \item{gam}{
Gamma hyperparameter when lambda.sel="EBIC". Defaults to gam=.25 (Barber et al., 2015)
}
  \item{d}{
Degrees of augmented interactions. The degree of augmented interactions reflects our belief about the maximal degree in the true graph. (see Loh & Wainwright, 2013)
}
  \item{rule.cat}{
Rule for combining parameters of interactions including categorical variables. The default, the "OR"-rule sets an edge to present when at least one parameter is non-zero. This corresponds to a conditional independence graph (Markov random field). Optionally, one can provide a matrix with the dimension of the parameter matrix to specify a costumized rule when to set an edge to present. Ones in this matrix indicate that the parameter has to be present to render the corresponding edge present and a zero indicates that the parameter does not have to be zero in order to render the corresponding edge present. A matrix of ones would indicate that for each interaction, all parameters have to be nonzero to render the corresponding edge present. A matrix of zeros is equivalent to the "OR"-rule.


The "OR"-rule determines conditional dependence if at least one parameter is non-zero. The "AND"-rule determines conditional dependence if all parameters are non-zero.
}
  \item{rule.reg}{
Rule for combining the two parameters obtained for each edge due to the neighborhood-regression approach. The "OR"-rule determines conditional dependence if at least one of the two parameters is non-zero, the "AND"-rule determines conditional dependence if both parameters are non-zero.}

  \item{pbar}{Shows a progress-bar if TRUE}
}



\value{
Returns a list containing:
\item{adj}{Adjacency matrix}
\item{wadj}{Weighted adjacency matrix; Please note that the parameters of interactions including categorical variables are not interpretable; for the actual interaction parameters inspect wpar.matrix}
\item{wpar.matrix}{Weighted model parameter matrix; This matrix is has the dimension of an overcomplete covariance matrix. By deleting empty columns/rows one obtains the minimal representation. It is given to enable the user to to inspect the exact interaction parameters involving categorical variables, which are (arbitrarily) aggregated in wadj by "rule.cat" (see above)}
\item{indvar.map}{List of length p, providing the mapping of indicator variables for categories and categorical variables in \code{wpar.matrix}. Can be used directly to plot indicator variables grouped by variables using the \code{qgraph} package.}
\item{indvar.map.labels}{List of length p, providing labels for the indicator variables in \code{indvar.map} of the form 1.1 1.2 ... 2 3.1 ...
}


\item{lambda}{p vector containing the L1-regularization parameter used for each nodewise-regression (selected with CV or EBIC)}
}
\references{
Barber, R. F., & Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. Electronic Journal of Statistics, 9, 567-607.

Loh, P. L., & Wainwright, M. J. (2013). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. The Annals of Statistics, 41(6), 3022-3049.

Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\seealso{
\code{confusion}, \code{mgmsampler}
}
\examples{

head(data_mixed2) #example data

type <- c("g", "p", "e", "c", "c") #Gaussian, Poission, ...
levs <- c(1,1,1,3,2) #Number of categories

fit <- mgmfit(data = data_mixed2, 
              type = type, 
              lev = levs, 
              lambda.sel = "CV", 
              folds = 10, 
              gam = .25, 
              d = 2, 
              rule.reg = "AND", 
              rule.cat = "OR")

}