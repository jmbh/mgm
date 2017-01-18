\name{tv.mgmfit}
\alias{tv.mgmfit}

\title{
Estimation of time-varying Mixed Graphical Models
}
\description{
Estimation of time-varying Mixed Graphical Models using L1-constrained neighborhood regression.}
\usage{
tv.mgmfit(data, type, lev, timepoints = NA, estpoints = "default", 
          bandwidth, gam = .25, d = 2,  rule.reg = "AND", 
          pbar = TRUE, method = "glm", missings = 'error', 
          ret.warn = FALSE)
}
\arguments{
  \item{data}{
n (time step) x p (variable) data matrix
}
  \item{type}{
p-vector containing the types of variable ("g" Gaussian, "p" Poisson, "c" Categorical)
}
  \item{lev}{
p-vector containing the number of levels for each variable p (lev=1 for all continuous variables)
}
  \item{timepoints}{
n-vector containing the time points of the measurements. Has to be a non-negative, strictly increasing sequence. Default \code{timepoints = NA} assumes equally spaced time intervals.
}
  \item{estpoints}{There are several input options: if \code{estpoints = 'default'}, the model is estimated at 2/bandwidth equally spaced time points. If a positive integer is provided, then the model is estimated at that many equally spaced time points. If a vector of time points is provided, then the model is estimated at exactly these time points. The specified time points have to be on the same scale as the provided data: if a time-vector is provided via \code{timepoints}, \code{estpoints} have to be provided on the same time scale. If \code{timepoints = NA}, \code{estpoints} can be directly provided as the row-numbers.}
  
  \item{bandwidth}{Standard deviation of the Gaussian kernel function.}
  \item{gam}{
Gamma hyperparameter when \code{lambda.sel="EBIC"}. Defaults to gam=.25 (Barber et al., 2015)
}
  \item{d}{
Degrees of augmented interactions. The degree of augmented interactions reflects our belief about the maximal degree in the true graph. (see Loh & Wainwright, 2013)
}
  \item{rule.reg}{
Rule for combining the two parameters obtained for each edge due to the neighborhood-regression approach. The "OR"-rule determines conditional dependence if at least one of the two parameters is non-zero, the "AND"-rule determines conditional dependence if both parameters are non-zero.}
  \item{pbar}{Shows a progress-bar if TRUE}
\item{method}{For each neighborhood regression, \code{method = "glm"} uses the appropriate link function for each variable type. \code{method = "linear"} uses linear regression for each variable, no matter of which type it is (for categorical variables, we predict each indicator variable using linear regression).}
\item{missings}{Handling of missing values. The default \code{missings = 'error'} returns an error message in case there are missing data. \code{missings = 'casewise.zw'} sets the weight of missing cases to zero. This has the advantage over casewise (timestepwise) deletion that it does not corrupt the time scale. Using the weighting technique, the time scale remains intact and we use less data at a time step that is close to time steps with missing data. As a consequence, the algorithm is less sensitive at these time steps.}
\item{ret.warn}{\code{TRUE} prints all warning messsages from \code{mgmfit} in the console. Note that all warnings are saved in the \code{mgmfit} object.}
}

\details{
For time-varying graphs \code{lambda.sel} is fixed to EBIC, as it is not straight-forward to apply cross-validation (CV) to time-series data.
}


\value{
Returns a list containing:

\item{call}{A list with all function inputs except the data}

\item{wadj}{A p x p x n weighted adjacency array, where n are the number of estimated time steps; Please note that we collapsed the interaction parameters involving categorical variables into one value by taking the mean of the absolute values; for the actual interaction parameters inspect mpar.matrix in the node models for each estimated time step.}

\item{mpar.matrix}{A n list with par x par matrices, where n are the number of estimated time steps and par are the number of parameters in each row in the overcomplete representation at the estimated time point n. For details see \code{mpar.matrix} in \code{mgmfit}.}

\item{signs}{A p x p x n array that gives the sign of the edge parameters, if defined, for each estimated time point n. The sign is defined for edges between two continuous variables, where the edge weight is equal to the single parameter. It is not defined for interactions involving categorical variables, where dependency depends on more than one parameter. In case a sign is defined it is indicated as -1 or 1, if a sign is undefined there is a 0 entry. If an edge is absent, there is NA entry.}

\item{edgecolor}{A p x p x n array that is equivalent to the \code{signs} array, however, it contains color names. This is useful for plotting the graphical model (for instance with the qgraph package) and indicating the sign (negative = red, positive = green, undefined = grey) for each edge in the graph.}

\item{t.models}{A list with n entries containing the \code{mgmfit} output at each estimated time point. For details see the help file of \code{mgmfit}.}

\item{Nt}{A n vector, indicating the sum of weights used at each time steps. As weights are normalized to the range [0,1], this can be interpreted as the effectively used sample size for estimation at time step n. This is interesting to check in the presence of missing data in the time-series.}
}

\references{

Haslbeck, J., & Waldorp, L. J. (2016). mgm: Structure estimation for time-varying mixed graphical models in high-dimensional data. arXiv preprint http://arxiv.org/abs/1510.06871v2.

Barber, R. F., & Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. Electronic Journal of Statistics, 9, 567-607.

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of statistical software, 33(1), 1. Chicago.

Haslbeck, J., & Waldorp, L. J. (2015). Structure estimation for mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.05677.

Loh, P. L., & Wainwright, M. J. (2013). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. The Annals of Statistics, 41(6), 3022-3049.

Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).

Zhou, S., Lafferty, J., & Wasserman, L. (2010). Time varying undirected graphs. Machine Learning, 80(2-3), 295-319.

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\seealso{
 \code{mgmfit}, \code{tv.mgmsampler}, \code{mgmsampler}
}
\examples{

\dontrun{

# We generate samples from a time-varying graph 
# and try to recover the true model from the sampled data

# We sample from a graph with 400 time steps, in which one parameter smoothly changes from 1 to 0
# and another parameter from 0 to 1

# specify time-varying graph
p <- 4 # 4 nodes
n <- 600 # number of time steps
graph <- matrix(0, p, p)
k <- 15 # steepness of sigmoid function
sigm_in <- 1/ (1+ exp(-k*(seq(0,1,length=n)-.5))) # sigmoid curve
sigm_de <- 1/ (1+ exp(k*(seq(0,1,length=n)-.5)))
graphs <- array(dim=c(p, p, n)) 
graphs[,,1:n] <- graph
graphs[1,2,] <- graphs[2,1,] <- - sigm_in
graphs[3,4,] <- graphs[4,3,] <- - sigm_de

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

## 3.3 Estimation of time-varying Mixed Graphical Models

# estimate
tv.obj <- tv.mgmfit(data = data, 
                    type = type,
                    lev = lev, 
                    estpoints = 50, # estimate at 50 equally spaced time points
                    bandwidth = .8/n^(1/3),
                    gam = 0, 
                    d = 2)


# visual check

plot(tv.obj$wadj[2,1,], xlab = 'Estimated time points', ylab = 'Edge weight', 
     type = 'l', col = 'red', ylim=c(0,1), lwd=2, yaxt='n')
lines(tv.obj$wadj[3,4,], col='blue', lwd=2)
axis(side = 2, at = round(seq(0,1,length=5), 2), las=2)
# true parameters
estpoints <- 50
lines( 1/ (1+ exp(-k*(seq(0,1,length=estpoints)-.5))) , col='red', lty=2)
lines( 1/ (1+ exp(k*(seq(0,1,length=estpoints)-.5))), col='blue', lty=2)

legend(38,.55, c('(2,1) true', '(2,1) estimated', '(3,4) true', '(3,4) estimated'), lty=c(2,1,2,1), 
       col = c('red', 'red', 'blue', 'blue'), lwd=c(2,2, 2, 2))


}

}