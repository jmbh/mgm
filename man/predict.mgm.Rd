\name{predict.mgm}
\alias{predict.mgm}

\title{
Nodewise predictions and errors
}

\description{
Computes predictions and error for all nodes in the graph.
}
\usage{
\method{predict}{mgm}(object, data, variables='all', 
        error.continuous = 'RMSE', 
        error.categorical = 'CorrectClass', ...)
}

\arguments{
  \item{object}{
The output object of \code{mgmfit}, \code{tv.mgmfit}, \code{var.mgm} or \code{tv_var.mgm}.
}
\item{data}{
A data matrix with the same properties as the data on which the model was fit. Note that different data can be used here to estimate the out of sample prediction error.}
\item{variables}{A character vector of variable names or a numeric vector of column numbers for which the the nodewise predictions and nodewise error should be computed. The default \code{variables='all'} computes predictions and error for all variables.}

\item{error.continuous}{Error measure for continuous variables. Default is \code{error.continuous = 'RMSE'} Root Mean Squared Error (RMSE). The alternative is \code{error.continuous = 'VarExpl'} proportion of explained variance 1 - var(y - yhat)/var(y).}


\item{error.categorical}{Error measure for categorica variables. Default is the proportion of correct classification \code{error.categorical = 'CorrectClass'}. The alternative is a normalized version of correct classification \code{error.categorical = 'CorrectClassNorm'}, which takes the marginal distribution into account. Specifically, CorrectClassNorm = (CorrectClass - max(RelativeFrequencies)) / (1-max(RelativeFrequencies)).

\code{error.categorical = 'VarExpl'} proportion of explained variance 1 - var(y - yhat)/var(y).}

\item{...}{Other prediction arguments.}
}

\value{
In the case of stationarty graphs the function returns a list, in the case of time varying graphs the function returns a list of lists with the following entries:

\item{pred}{p x n matrix of predictions, where p are the selected variables and n are the number of observations in the provided dataset.}

\item{error}{A data frame with the columns: Variable, Error and ErrorType. As error we report the root mean squared error (RMSE) for continuous variables and the proportion of correct classification for categorical variables. For time-varying graphs, both types of errors are weighted by the same weights used during estimation.}


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

fitobj$wadj # We get a strong dependency between x1-x3

# We compute nodewise predictions and errors
predobj <- predict.mgm(fitobj, data)

# And look at nodewise errors
predobj$error

}

}