\name{resample}
\alias{resample}

\title{
Resampling scheme for mgm objects
}
\description{
Fits mgm model types (mgm, mvar, tvmgm, tvmvar) to a specified number of bootstrap samples.
}
\usage{
resample(object, data, nB, blocks, quantiles, 
         pbar, verbatim, ...)
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
  \item{blocks}{
The number of blocks for the block bootstrap used for time-varying models.
}

  \item{quantiles}{
A vector with two values in [0, 1], specifying quantiles computed on the bootstrapped sampling distributions. Defaults to \code{quantiles = c(.05, .95)}
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

The seeds for the bootstrap samples are randomly sampled. For MGMs, the seeds are resampled until there are nB bootstrap samples on which an MGM can be estimated. This resampling has been implemented, because especially for small data sets, one can obtain bootstrap samples in which one or several variables have zero variance. For the other model classes, an informative error is returned in case the respective model cannot be estimated on one or more of the bootstrap samples.

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

\item{seeds}{\code{nB} integers indicating the seeds used to sample the \code{nB} bootstrap samples.}

}


\references{

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. Journal of Statistical Software, 93(8), pp. 1-46. DOI: 10.18637/jss.v093.i08


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

# 4) Plot summary of all sampling distributions
plotRes(object = res_obj, 
        quantiles = c(0.05, .95))


# For more examples see https://github.com/jmbh/mgmDocumentation



}


}
