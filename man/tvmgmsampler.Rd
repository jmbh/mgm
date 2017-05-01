\name{tvmgmsampler}
\alias{tvmgmsampler}

\title{
Sample from time-varying k-order Mixed Graphical Model
}
\description{
Generates samples from a time-varying k-order Mixed Graphical Model
}
\usage{
tvmgmsampler(factors, interactions, thresholds, sds, type,
             level, nIter = 250, pbar = TRUE, ...)
}

\arguments{
  \item{factors}{
The same object as \code{factors} in \code{mgmsampler()}. An interaction is specified in \code{factors} if it should be nonzero at least at one time point in the time series. The values of each parameter at each time point is specified via \code{interactions}.
}
  \item{interactions}{
The same object as \code{factors} in \code{mgmsampler()}, except that each array indicating the parameters of an interaction has an additional (the last) dimension, indicating time. Corresponding to the time vector in \code{factors}, the time vector has to be a sequence of integers \{1, 2, ..., N\}. For an illustration see the examples below.
}
  \item{thresholds}{
A list with p entries for p variables, each of which contains a N x m matrix. The columns contain the m thresholds for m categories (for continuous variables m = 1 and the entry contains the threshold/intercept). The rows indicate how the thresholds change over time.
}
  \item{sds}{
N x p matrix indicating the standard deviations of Gaussians specified in \code{type} for \{1, ..., N\} time points. Entries not referring to Gaussians are ignored.
}
  \item{type}{
p character vector indicating the type of variable for each column in \code{data}. 'g' for Gaussian, 'p' for Poisson, 'c' of each variable.
}
  \item{level}{
p integer vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{nIter}{
Number of iterations in the Gibbs sampler until a sample is drawn.
}
  \item{pbar}{
  If \code{pbar = TRUE} a progress bar is shown. Defaults to \code{pbar = TRUE}.
}

  \item{...}{
  Additional arguments.
}

}
\details{

\code{tvmgmsampler} is a wrapper function around \code{mgmsampler}. Its input is the same as for \code{mgmsampler}, except that each object has an additional dimension for time. The number of time points is specified via entries in the additional time dimension.

}
\value{

A list containing:

\item{call}{
Contains all provided input arguments.
}

\item{data}{
The N x p data matrix of sampled values
}


}

\references{

Haslbeck, J., & Waldorp, L. J. (2016). mgm: Structure Estimation for time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

Yang, E., Baker, Y., Ravikumar, P., Allen, G. I., & Liu, Z. (2014, April). Mixed Graphical Models via Exponential Families. In AISTATS (Vol. 2012, pp. 1042-1050).

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\examples{


\dontrun{


# --------- Example 1: p = 4 dimensional Gaussian ---------

# ----- 1) Specify Model -----

# a) General Graph Info
type = c('g', 'g', 'g', 'g') # Four Gaussians
level = c(1, 1, 1, 1)
n_timepoints = 500 #  Number of time points

# b) Define Interaction
factors <- list()
factors[[1]] <- array(NA, dim=c(2, 2)) # two pairwise interactions
factors[[1]][1, 1:2] <- c(3,4)
factors[[1]][2, 1:2] <- c(1,2)

# Two parameters, one linearly increasing from 0 to 0.8, another one lin decreasing from 0.8 to 0
interactions <- list()
interactions[[1]] <- vector('list', length = 2)
interactions[[1]][[1]] <- array(0, dim = c(level[1], level[2], n_timepoints))
interactions[[1]][[1]][1,1, ] <- seq(.8, 0, length = n_timepoints)
interactions[[1]][[2]] <- array(0, dim = c(level[1], level[2], n_timepoints))
interactions[[1]][[2]][1,1, ] <- seq(0, .8, length = n_timepoints)

# c) Define Thresholds
thresholds <- vector('list', length = 4)
thresholds <- lapply(thresholds, function(x) matrix(0, ncol = level[1], nrow = n_timepoints))

# d) Define Standard deviations
sds <- matrix(1, ncol = length(type), nrow = n_timepoints) # constant across variables and time


# ----- 2) Sample cases -----

set.seed(1)
dlist <- tvmgmsampler(factors = factors,
                      interactions = interactions,
                      thresholds = thresholds,
                      sds = sds,
                      type = type,
                      level = level,
                      nIter = 75,
                      pbar = TRUE)


# ----- 3) Recover model from sampled cases -----

set.seed(1)
tvmgm_obj <- tvmgm(data = dlist$data,
                   type = type,
                   level = level,
                   estpoints = seq(1, n_timepoints, length = 15),
                   bandwidth = .2,
                   k = 3,
                   lambdaSel = 'CV',
                   ruleReg = 'AND')

# How well did we recover those two time-varying parameters?
plot(tvmgm_obj$pairwise$wadj[3,4,], type='l', ylim=c(0,.8))
lines(tvmgm_obj$pairwise$wadj[1,2,], type='l', col='red')
# Looks quite good


# --------- Example 2: p = 5 binary; one 3-way interaction ---------


# ----- 1) Specify Model -----

# a) General Graph Info
p <- 5 # number of variables
type = rep('c', p) # all categorical
level = rep(2, p) # all binary
n_timepoints <- 1000

# b) Define Interaction
factors <- list()
factors[[1]] <- NULL # no pairwise interactions
factors[[2]] <- array(NA, dim = c(1,3)) # one 3-way interaction
factors[[2]][1, 1:3] <- c(1, 2, 3)

interactions <- list()
interactions[[1]] <- NULL # no pairwise interactions
interactions[[2]] <- vector('list', length = 1)  # one 3-way interaction
# 3-way interaction no1
interactions[[2]][[1]] <- array(0, dim = c(level[1], level[2], level[3], n_timepoints))
theta <- 2
# fill in time-varying parameter: linear increase from 0 to 2
interactions[[2]][[1]][1, 1, 1, ] <- seq(0, 2, length = n_timepoints)

# c) Define Thresholds
thresholds <- list()
for(i in 1:p) thresholds[[i]] <- matrix(0, nrow = n_timepoints, ncol = level[i])


# ----- 2) Sample cases -----

set.seed(1)
dlist <- tvmgmsampler(factors = factors,
                      interactions = interactions,
                      thresholds = thresholds,
                      type = type,
                      level = level,
                      nIter = 150,
                      pbar = TRUE)


# ----- 3) Check Marginals -----

dat <- dlist$data[1:round(n_timepoints/2),]
table(dat[,1], dat[,2], dat[,3])

dat <- dlist$data[round(n_timepoints/2):n_timepoints,]
table(dat[,1], dat[,2], dat[,3])

# Observation: much stronger effect in second hald of the time-series,
# which is what we expect


# ----- 4) Recover model from sampled cases -----

set.seed(1)
tvmgm_obj <- tvmgm(data = dlist$data,
                   type = type,
                   level = level,
                   estpoints = seq(1, n_timepoints, length = 15),
                   bandwidth = .2,
                   k = 3,
                   lambdaSel = 'CV',
                   ruleReg = 'AND')

tvmgm_obj$rawfactor$indicator
# Seems very difficult to recover this time-varying 3-way binary interaction
# See also the corresponding difficulties in the stationary version (examples of ?mgmsampler)




}


}


