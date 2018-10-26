\name{tvmvarsampler}
\alias{tvmvarsampler}

\title{
Sampling from a time-varying mixed VAR model
}
\description{
Function to sample from a time-varying mixed VAR (mVAR) model
}
\usage{
tvmvarsampler(coefarray, lags, thresholds,
              sds, type, level, pbar)
}

\arguments{
  \item{coefarray}{
A p x p x max(level) x max(level) x n_lags x N array, where p are the number of variables, level is the input argument \code{level} and n_lags is the number of specified lags in \code{lags}, so n_lags = length(n_lags), and N is the number of time points in the time series. The first four dimensions specify the parameters involved in the cross-lagged effects of the lag specified in the 5th dimension. I.e. coefarray[5, 6, 1, 1, 3, 100] indicates the cross-lagged effect of variable 6 on variable 5 (if both are continuous), for the third lag specified in \code{lags} at time point 100. If variable 1 and 3 are categorical with m = 2 and  = 4 categories, respectively, then coefarray[1, 3, 1:2, 1:4, 1, 250] indicates the m*s=8 parameters specifying this interaction for the first lag specified in \code{lags} at time point 250. See the examples below for an illustration.
}
  \item{lags}{
A vector indicating the lags in the mVAR model. E.g. \code{lags = c(1, 4, 9)} specifies lags of order 1, 3, 9. The number of specified lags has to match the 5th dimension in \code{coefarray}.
}
  \item{thresholds}{
A list with p entries, each consisting of a matrix indicating a threshold for each category of the given variable (column) and time point (row). For continuous variable, the matrix has 1 column.
}
  \item{sds}{
A N x p matrix specifying the standard deviation of Gaussian variables (columns) at each time point (rows)If non-Gaussian variables are included in the mVAR model, the corresponding columns are ignored.
}
  \item{type}{
p vector indicating the type of variable for each column in \code{data}. "g" for Gaussian, "p" for Poisson, "c" for categorical.
}
  \item{level}{
p vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{pbar}{
If \code{pbar = TRUE}, a progress bar is shown.
}
}
\details{
We sample from the mVAR model by separately sampling from its corresponding p conditional distributions.
}
\value{
A list with two entries:

\item{call}{The function call}

\item{data}{The sampled n x p data matrix}

}
\references{

Haslbeck, J., & Waldorp, L. J. (2018). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}



\examples{


\dontrun{

## We specify a tvmvar model, sample from it and recover it

# a) Set up time-varying mvar model

p <- 6 # Six variables
type <- c("c", "c", "c", "c", "g", "g") # 4 categorical, 2 gaussians
level <- c(2, 2, 4, 4, 1, 1) # 2 categoricals with 2 categories, 2 with 5
max_level <- max(level)

lags <- c(1, 3, 9) # include lagged effects of order 1, 3, 9
n_lags <- length(lags)

N <- 5000

# Specify thresholds
thresholds <- list()
thresholds[[1]] <- matrix(0, ncol=2, nrow=N)
thresholds[[2]] <- matrix(0, ncol=2, nrow=N)
thresholds[[3]] <- matrix(0, ncol=4, nrow=N)
thresholds[[4]] <- matrix(0, ncol=4, nrow=N)
thresholds[[5]] <- matrix(0, ncol=1, nrow=N)
thresholds[[6]] <- matrix(0, ncol=1, nrow=N)

# Specify standard deviations for the Gaussians
sds <- matrix(NA, ncol=6, nrow=N)
sds[,5:6] <- 1

# Create coefficient array
coefarray <- array(0, dim=c(p, p, max_level, max_level, n_lags, N))

# a.1) interaction between continuous 5<-6, lag=3
coefarray[5, 6, 1, 1, 2, ] <- c(rep(.5, N/2), rep(0, N/2))
# a.2) interaction between 1<-3, lag=1
m1 <- matrix(0, nrow=level[2], ncol=level[4])
m1[1, 1:2] <- 1
m1[2, 3:4] <- 1
coefarray[1, 3, 1:level[2], 1:level[4], 1, ] <- m1
# a.3) interaction between 1<-5, lag=9
coefarray[1, 5, 1:level[1], 1:level[5], 3, ] <- c(0, 1)

dim(coefarray)


# b) Sample
set.seed(1)
dlist <- tvmvarsampler(coefarray = coefarray,
                       lags = lags,
                       thresholds = thresholds,
                       sds = sds,
                       type = type,
                       level = level,
                       pbar = TRUE)


# c) Recover: time-varying mVAR model
set.seed(1)
tvmvar_obj <- tvmvar(data = dlist$data,
                     type = type,
                     level = level,
                     lambdaSel = "CV",
                     lags = c(1, 3, 9),
                     estpoints = seq(0, 1, length=10),
                     bandwidth = .05)

tvmvar_obj$wadj[5, 6, 2, ] # parameter goes down, as specified
tvmvar_obj$wadj[1, 3, 1, ]
tvmvar_obj$wadj[1, 5, 3, ]

# For more examples see https://github.com/jmbh/mgmDocumentation

}

}


