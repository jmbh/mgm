\name{mvarsampler}
\alias{mvarsampler}

\title{
Sampling from a mixed VAR model
}
\description{
Function to sample from a mixed VAR (mVAR) model
}
\usage{
mvarsampler(coefarray, lags, thresholds, 
            sds, type, level, N, pbar)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{coefarray}{
A p x p x max(level) x max(level) x n_lags array, where p are the number of variables, level is the input argument \code{level} and n_lags is the number of specified lags in \code{lags}, so n_lags = length(n_lags). The first four dimensions specify the parameters involved in the cross-lagged effects of the lag specified in the 5th dimension. I.e. coefarray[5, 6, 1, 1, 3] indicates the cross-lagged effect of variable 6 on variable 5 (if both are continuous), for the third lag specified in \code{lags}. If variable 1 and 3 are categorical with m = 2 and  = 4 categories, respectively, then coefarray[1, 3, 1:2, 1:4, 1] indicates the m*s=8 parameters specifying this interaction for the first lag specified in \code{lags}. See the examples below for an illustration.
}
  \item{lags}{
A vector indicating the lags in the mVAR model. E.g. \code{lags = c(1, 4, 9)} specifies lags of order 1, 3, 9. The number of specified lags has to match the 5th dimension in \code{coefarray}.
}
  \item{thresholds}{
A list with p entries, each consisting of a vector indicating a threshold for each category of the given variable. For continuous variable, the vector has length 1.
}
  \item{sds}{
A vector of length p indicating the standard deviations of the included Gaussian nodes. If non-Gaussian variables are included in the mVAR model, the corresponding entries are ignored.
}
  \item{type}{
p vector indicating the type of variable for each column in \code{data}. "g" for Gaussian, "p" for Poisson, "c" for categorical.
}
  \item{level}{
p vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{N}{
The number of samples to be drawn from the specified mVAR model.
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

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. Journal of Statistical Software, 93(8), pp. 1-46. DOI: 10.18637/jss.v093.i08

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\examples{

\dontrun{


## Generate data from mixed VAR model using mvarsampler() and recover model using mvar()

# 1) Define mVAR model

p <- 6 # Six variables
type <- c("c", "c", "c", "c", "g", "g") # 4 categorical, 2 gaussians
level <- c(2, 2, 4, 4, 1, 1) # 2 categoricals with m=2, 2 categoricals with m=4, two continuous
max_level <- max(level)

lags <- c(1, 3, 9) # include lagged effects of order 1, 3, 9
n_lags <- length(lags)

# Specify thresholds
thresholds <- list()
thresholds[[1]] <- rep(0, level[1])
thresholds[[2]] <- rep(0, level[2])
thresholds[[3]] <- rep(0, level[3])
thresholds[[4]] <- rep(0, level[4])
thresholds[[5]] <- rep(0, level[5])
thresholds[[6]] <- rep(0, level[6])

# Specify standard deviations for the Gaussians
sds <- rep(NULL, p)
sds[5:6] <- 1

# Create coefficient array
coefarray <- array(0, dim=c(p, p, max_level, max_level, n_lags))

# a.1) interaction between continuous 5<-6, lag=3
coefarray[5, 6, 1, 1, 2] <- .4
# a.2) interaction between 1<-3, lag=1
m1 <- matrix(0, nrow=level[2], ncol=level[4])
m1[1,1:2] <- 1
m1[2,3:4] <- 1
coefarray[1, 3, 1:level[2], 1:level[4], 1] <- m1
# a.3) interaction between 1<-5, lag=9
coefarray[1, 5, 1:level[1], 1:level[5], 3] <- c(0, 1)


# 2) Sample
set.seed(1)
dlist <- mvarsampler(coefarray = coefarray,
                     lags = lags,
                     thresholds = thresholds,
                     sds = sds,
                     type = type,
                     level = level,
                     N = 200,
                     pbar = TRUE)

# 3) Recover
set.seed(1)
mvar_obj <- mvar(data = dlist$data,
                 type = type,
                 level = level,
                 lambdaSel = "CV",
                 lags = c(1, 3, 9),
                 signInfo = FALSE,
                 overparameterize = F)

# Did we recover the true parameters?
mvar_obj$wadj[5, 6, 2] # cross-lagged effect of 6 on 5 over lag lags[2]
mvar_obj$wadj[1, 3, 1] # cross-lagged effect of 3 on 1 over lag lags[1]
mvar_obj$wadj[1, 5, 3] # cross-lagged effect of 1 on 5 over lag lags[3]


# For more examples see https://github.com/jmbh/mgmDocumentation


}

}






