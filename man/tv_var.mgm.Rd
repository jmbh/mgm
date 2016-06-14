\name{tv_var.mgm}
\alias{tv_var.mgm}

\title{
Estimation of time-varying mixed Autoregressive Models
}
\description{
Estimation of time-varying mixed Autoregressive Models using L1-constrained regression.}
\usage{
tv_var.mgm(data, type, lev, lags = 1, tsteps = "default", 
          bandwidth, gam = .25, d = 2, pbar = TRUE, 
          method = "glm", missings = 'error', ret.warn = FALSE)
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
\item{lags}{Order of included lags. Currently limited to \code{lag = 1}.
}
  \item{tsteps}{Number of equidistant time points at which the model is estimated. If \code{tsteps = 'default'} then the number of estimated time steps is 2/bandwidth. This provides a reasonable coverage of the time series for many choices of bandwidt.}
  \item{bandwidth}{Standard deviation of the Gaussian kernel function.}
  \item{gam}{
Gamma hyperparameter when lambda.sel="EBIC". Defaults to gam=.25 (Barber et al., 2015)
}
  \item{d}{
Degrees of augmented interactions. The degree of augmented interactions reflects our belief about the maximal degree in the true graph. (see Loh & Wainwright, 2013)
}

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

\item{wadj}{A p x p x n weighted adjacency matrix, where n are the number of estimated time steps; Please note that we collapsed the interaction parameters involving categorical variables into one value by taking the mean of the absolute values; for the actual interaction parameters inspect mpar.matrix in the node models for each estimated time step.}

\item{mpar.matrix}{A n list with par x par matrices, where n are the number of estimated time steps and par are the number of parameters in each row in the overcomplete representation at the estimated time point n. For details see \code{mpar.matrix} in \code{mgmfit}.}

\item{signs}{A p x p x n array that gives the sign of the edge parameters, if defined, for each estimated time point n. The sign is defined for edges between two continuous variables, where the edge weight is equal to the single parameter. It is not defined for interactions involving categorical variables, where dependency depends on more than one parameter. In case a sign is defined it is indicated as -1 or 1, if a sign is undefined there is a 0 entry. If an edge is absent, there is NA entry.}

\item{edgecolor}{A p x p x n array that is equivalent to the \code{signs} array, however, it contains color names. This is useful for plotting the graphical model (for instance with the qgraph package) and indicating the sign (negative = red, positive = green, undefined = grey) for each edge in the graph.}

\item{t.models}{A list with n entries containing the \code{mgmfit} output at each estimated time point. For details see the help file of \code{mgmfit}.}

\item{Nt}{A n vector, indicating the sum of weights used at each time steps. As weights are normalized to the range [0,1], this can be interpreted as the effectively used sample size for estimation at time step n. This is interesting to check in the presence of missing data in the time-series.}
}

\references{
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

# we set up a VAR model of 4 variables and 2 smoothly changing parameters:
# the autocorrelation 1-> smootly increases from 0 to .8
# the cross-correlation 2-> smootly decreases from .8 to 0

#set up graph
p <- 4 # 4 nodes
n <- 600 # number of time steps
graph <- matrix(0, p, p)
k <- 15 # steepness of sigmoid function
theta <- .8
sigm_in <- theta * 1/ (1+ exp(-k*(seq(0,1,length=n)-.5))) # sigmoid curve
sigm_de <- theta * 1/ (1+ exp(k*(seq(0,1,length=n)-.5)))
graphs <- array(dim=c(p, p, n)) 
graphs[,,1:n] <- graph
graphs[1,1,] <- sigm_in # autocorrelation 1->1
graphs[2,3,] <- sigm_de # cross correlation 2->3

InnoVar <- array(data=NA, dim=c(p,1,n)) #innocation variance
for(i in 1:n) InnoVar[,,i] <- colSums(graphs[,,i]) + 1 

# generate data from VAR model
set.seed(23)
data <- matrix(NA, n, p)
data[1, ] <- rnorm(p, 0, InnoVar)
for(i in 2:n) {
  for(v in 1:p) {
    data[i, v] <-  sum(data[i-1,] * graphs[,v,i]) + rnorm(1,0,InnoVar[v,,i])
  }
}

head(data) # first 6 time steps

# estimate

type <- c('g', 'g', 'g', 'g') # four Gaussians
lev <- c(1, 1, 1, 1)
tsteps <- 50

tv_var.obj <- tv_var.mgm(data = data, 
                     type = type,
                     lev = lev, 
                     tsteps = tsteps,
                     bandwidth = .8/n^(1/3),
                     gam = 0, 
                     d = 1)


# visual check
plot(tv_var.obj$wadj[1,1,], xlab = 'Estimated time points', ylab = 'Edge weight', 
     type = 'l', col = 'red', ylim=c(0,1), lwd=2, yaxt='n')
lines(tv_var.obj$wadj[2,3,], col='blue', lwd=2)
axis(side = 2, at = round(seq(0,1,length=5), 2), las=2)
# true parameters
tsteps <- 50
lines( theta* 1/ (1+ exp(-k*(seq(0,1,length=tsteps)-.5))) , col='red', lty=2)
lines( theta* 1/ (1+ exp(k*(seq(0,1,length=tsteps)-.5))), col='blue', lty=2)

legend(38,.55, c('1->1 true', '1->1 estimated', '2->3 true', '2->3 estimated'), lty=c(2,1,2,1), 
       col = c('red', 'red', 'blue', 'blue'), lwd=c(2,2, 2, 2))




}

}