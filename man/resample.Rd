\name{resample}
\alias{resample}

\title{
Resampling scheme for mgm objects
}
\description{
Fits mgm model types (mgm, mvar, tvmgm, tvmvar) to a specified number of bootstrap samples
}
\usage{
resample(object, data, nB, seeds, blocks, quantiles, pbar, ...)
}

\arguments{
  \item{object}{
An mgm model object, the output of mgm(), tvmgm(), mvar(), tvmvar(). The model specifications for all fitted models are taken from this model object.
}
  \item{data}{
The n x p data matrix.
}
  \item{nB}{
The number of bootstrap samples.
}
  \item{seeds}{
A integer vector of length \code{nB}, providing random seeds for each bootstrap sample. Defaults to \code{seeds = 1:nB}.
}
  \item{blocks}{
The number of blocks for the block bootstrap used for time-varying models.
}

  \item{quantiles}{
A vector with values in [0, 1], specifying quantiles computed on the bootstrapped sampling distributions. Defaults to \code{quantiles = c(.05, .95)}
}


  \item{pbar}{
If \code{TRUE}, a progress bar is shown. Defaults to \code{pbar = TRUE}.
}
  \item{verbatim}{
If \code{TRUE}, the seed of the current bootstrap sample is printed in the console. Useful to exclude zero-variance bootstrap samples in datasets with low variance.
}
  \item{\dots}{
Additional arguments.
}
}

\details{

\code{resample()} fits a model specified via the \code{object} argument to \code{nB} bootstrap samples obtained from the orginial dataset. For stationary models (mgm() and mvar()) objects, we use the standard bootstrap. For time-varying models (tvmgm() and tvmvar()) we use the block bootstrap.

For mvar models, \code{bootParameters} is a p x p x nlags x nB array, where p is the number of variables, nlags is the number of specified lags, and nB is the number of bootstrap samples. Thus \code{bootParameters[7, 3, 2, ]} returns the bootstrapped sampling distribution of the lagged effect from variable 3 on 7 for the 2nd specified lag. See also \code{?mvar}.

For tvmar models, \code{bootParameters} is a p x p x nlags x nestpoints x nB array, analogously to mvar models. nestpoints is the number of specified estpoints. See also \code{?tvmvar}.

Resampling is currently only supported for pairwise MGMs (\code{k = 2}). For mgms, \code{bootParameters} is a p x p x nB array. For tvmgms, \code{bootParameters} is a p x p x nestpoint x nB array.

}
\value{
The output consists of a list with the entries

\item{call}{Contains the function call.}

\item{models}{A list with \code{nB} entries, containing the models fit to the bootstrapped sampels.}

\item{bootParameters}{Contains all the bootstrapped sampling distribution of all parameters. The dimension of this object depends on the type of model. Specifically, this object has the same dimension as the main  parameter output of the corresponding estimation function, with one dimension added for the bootstrap iterations. See "Details" above.
}

\item{bootQuantiles}{Contains the specified quantiles of the bootstrapped sampling distribution for each parameter. Has the same structure as \code{bootParameters}. See "Details" above. }

\item{Times}{Returns the running time for each bootstrap model in seconds.}

\item{totalTime}{Returns the running time for all bootstrap models together in seconds.}

}


\references{

Haslbeck, J., & Waldorp, L. J. (2017). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\examples{

\dontrun{


# 1) Fit mgm to example dataset (true edges: 1-4, 2-3, 1-2)

mgm_obj <- mgm(data = mgm_data$data, 
               type = c('g', 'c', 'c', 'g'),
               level = c(1, 2, 4, 1),
               k = 2, 
               lambdaSel = 'CV', 
               threshold = 'none')


# 2) Take 50 bootstrap samples using resample()

res_obj <- resample(object = mgm_obj, 
                    data = mgm_data$data,
                    nB = 50)


# 3) Plot histogram of bootstrapped sampling distribution of edge 1-4

hist(res_obj$bootParameters[1, 4, ], 
     main = "", 
     xlab = "Parameter Estimate")


}


}
