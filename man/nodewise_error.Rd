\name{nodewise_error}
\alias{nodewise_error}

\title{
Computing node wise errors.
}

\description{
Computes the node wise error for all nodes in the graph.}
\usage{
nodewise_error(fitobj, data)
}

\arguments{
  \item{fitobj}{
The output object of \code{mgmfit}, \code{tv.mgmfit}, \code{var.mgm} or \code{tv_var.mgm}.
}
\item{data}{
A data matrix with the same properties as the data on which the model was fit. Note that different data can be used here to estimate the out of sample prediction error.
}

}

\details{
For time-varying graphs \code{lambda.sel} is fixed to EBIC, as it is not straight-forward to apply cross-validation (CV) to time-series data.
}


\value{
In the case of stationarty graphs the function returns a dataframe, in the case of time varying graphs the function returns a list of dataframes with the following columns:

\item{Variable}{The variable lable, in case it was provided.}

\item{Error}{The root mean squared error (RMSE) for continuous variables and the proportion of correct classification for categorical variables. For time-varying graphs, both types of errors are weighted by the same weights used during estimation.}

\item{ErrorType}{A character vector indicating which type of error is displayed for a given variable.}

}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\seealso{
 \code{mgmfit}, \code{tv.mgmfit}, \code{var.mgm},  \code{tv_var.mgm}, \code{summary.mgm}
}
\examples{

\dontrun{

set.seed(1)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100) + 4*x1
data <- cbind(x1, x2, x3)
fitobj <- mgmfit(data, type=c('g', 'g', 'g'), lev=rep(1,3), d=1, ret.warn = FALSE)

fitobj # We get a strong dependency between x1-x3

# and as expected we are doing well in predicting x3 by x1 and vice versa:
nodewise_error(fitobj, data)

}

}