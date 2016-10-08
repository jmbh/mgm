\name{var.mgm}
\alias{var.mgm}

\title{
Estimation of mixed Autoregressive Models
}
\description{
Estimation of mixed Autoregressive Models using L1-constrained regression.}
\usage{
var.mgm(data, type, lev, lags = 1, lambda.sel = "EBIC", 
        folds = 10, gam = .25, d=2, pbar = TRUE, 
        method = "glm", missings = 'error', 
        weights = NA, ret.warn = TRUE, binary.sign = FALSE)
}

% input
\arguments{
  \item{data}{
n x p data matrix
}
  \item{type}{
p-vector containing the types of variable ("g" Gaussian, "p" Poisson, "c" Categorical)
}
  \item{lev}{
p-vector containing the number of levels for each variable p (\code{lev=1} for all continuous variables)
}
  \item{lags}{Order of included lags. Currently limited to \code{lag = 1}.
}
  \item{lambda.sel}{
Procedure to select the lambda-parameter for L1-penalized regressions. The two options are cross validation "CV" and "EBIC" (default). While cross-validation (CV) is used in the paper (see below), the Extended Bayesian Information Criterion (EBIC, see Barber et al., 2015) is a useful alternative in cases where the CV requirement that all categories are present in each fold is not met or in cases where CV is computationally unfeasible.
}
  \item{folds}{
The number of folds to be used for cross-validation in case \code{lambda.sel="CV"}. Defaults to \code{folds=10}.}
  \item{gam}{
Gamma hyperparameter when \code{lambda.sel="EBIC"}. Defaults to \code{gam=.25} (Barber et al., 2015).
}
  \item{d}{
Degrees of augmented interactions. The degree of augmented interactions reflects our belief about the maximal degree in the true graph. (see Loh & Wainwright, 2013)
}

\item{pbar}{Shows a progress-bar if TRUE.}

\item{method}{For each neighborhood regression, \code{method = "glm"} uses the appropriate link function for each variable type. \code{method = "linear"} uses linear regression for each variable, no matter of which type it is (for categorical variables, this method predicts each indicator variable using linear regression).}

\item{missings}{Handling of missing values. The default \code{missings = 'error'} returns an error message in case there are missing data. \code{missings = 'casewise.zw'} sets the weight of missing cases to zero. For stationary graphs, this is the equivalent to casewise deletion. In time varying graphs, this avoids corruption of the time scale, which would be the consequence of casewise deletion. For details, please see Haslbeck and Waldorp (2015).}

\item{weights}{Weights for observations. Defaults to 1 for each observation.}

\item{ret.warn}{If \code{FALSE} warnings are surpressed.}

\item{binary.sign}{The default in mgm is to assign no sign to interactions involving categorical variables. This is because for the case where the number of categories m > 2, every interaction has to involve > 1 parameter and hence no sign of the interaction can be defined. However, for the case = 2 we have one parameter describing the interaction between the binary variable at hand and any other binary or continuous variable. If \code{binary.sign = TRUE} and if the binary variable is in (0,1) then the sign of this interaction will be recovered and returned in the sign matrix (see below).}

}

% output

\value{
Returns a list containing:
\item{call}{A list with all function inputs except the data}

\item{adj}{Directed adjacency matrix}

\item{wadj}{Weighted directed adjacency matrix; Please note that we collapsed the interaction parameters involving categorical variables into one value by taking the mean of the absolute values; for the actual interaction parameters inspect wpar.matrix}

\item{mpar.matrix}{This matrix includes all estimated parameters in one matrix. If there are only continuous variables in the dataset this matrix has the same dimensions as \code{wadj}. If there are categorical variables in the dataset, the matrix has the dimensions of the overcomplete representations. NA values indicate parameters that are not estimated (because they are redundant). For details see Friedman, Hastie and Tibshirani (2010).}

\item{signs}{A p x p matrix that gives the sign of the edge parameters, if defined. The sign is defined for edges between two continuous variables, where the edge weight is equal to the single parameter. It is not defined for interactions involving categorical variables, where dependency depends on more than one parameter. In case a sign is defined it is indicated as -1 or 1, if a sign is undefined there is a 0 entry. If an edge is absent, there is NA entry.}

\item{edgecolor}{A p x p matrix that is equivalent to the \code{signs} matrix, however, it contains color names. This is useful for plotting the graphical model (for instance with the qgraph package) and indicating the sign (negative = red, positive = green, undefined = grey) for each edge in the graph.}

\item{node.models}{A list that includes for each node: a list of estimated parameters, the selected regularization parameter, the threshold parameter and the EBIC of the model (if applicable).}

\item{par.labels}{A vector matching with the dimensions of \code{mpar.matrix} indicating which parameter belongs to the interaction between which variables.}

\item{warnings}{A list of all warnings produced during execution. Useful to have when warnings are surpressed in console output \code{(ret.warn = TRUE)}.}}



\references{
Barber, R. F., & Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. Electronic Journal of Statistics, 9, 567-607.

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of statistical software, 33(1), 1. Chicago.

Haslbeck, J., & Waldorp, L. J. (2015). Structure estimation for mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.05677.

Loh, P. L., & Wainwright, M. J. (2013). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. The Annals of Statistics, 41(6), 3022-3049.

Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\seealso{
\code{mgmfit}, \code{tv.mgmfit} \code{tv-var.mgm}
}
\examples{

\dontrun{

#set up graph
p <- 5 # number of variables
tsteps <- 200 # number of time steps
set.seed(3)
ParMat <- matrix(runif(p^2, .4, .8), p, p) * 
            matrix(sample(c(0,1), p^2, replace=TRUE, prob=c(.6,.4)), p, p)
ParMat
InnoVar <- colSums(ParMat) + 1 # innovation variance

# generate data from VAR model
set.seed(99)
data <- matrix(NA, tsteps, p)
data[1, ] <- rnorm(p, 0, InnoVar)
for(i in 1:tsteps) {
  for(v in 1:p) {
    data[i, v] <-  sum(data[i-1,] * ParMat[,v]) + rnorm(1,0,InnoVar[v])
  }
}

head(data) # first 6 time steps

# estimate VAR model
library(mgm)
var.fit <- var.mgm(data, type=rep('g', p), lev=rep(1,p), lags=1, gam=.25, d=2)

# compare true & estimated VAR model
library(qgraph)
par(mfrow=c(1,2))
Q <- qgraph(ParMat, minimum=0, maximum=range(ParMat)[2]) 
qgraph(var.fit$mpar.matrix, layout=Q$layout, minimum=0, maximum=range(ParMat)[2])

# observation:
# all parameters were recovered correctly; parameter weights are smaller due to regularization

}


}

