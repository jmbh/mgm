\name{mgmsampler}
\alias{mgmsampler}

\title{
Sample from k-order Mixed Graphical Model
}
\description{
Generates samples from a k-order Mixed Graphical Model
}
\usage{
mgmsampler(factors, interactions, thresholds, sds, type,
           level, N, nIter = 250, pbar = TRUE, 
           divWarning = 10^3, returnChains = FALSE)
}


\arguments{
  \item{factors}{
This object indicates which interactions are present in the model.  It is a list in which the first entry corresponds to 2-way interactions, the second entry corresponds to 3-way interactions, etc. and the kth entry to the k+1-way interaction. Each entry contains a matrix with dimensions order x number of interaction of given order. Each row in the matrix indicates an interaction, e.g. (1, 3, 7, 9) in the matrix in list entry three indicates a 4-way interaction between the variables 1, 3, 7 and 9.
}
  \item{interactions}{
This object specifies the parameters associated to the interactions specified in \code{factors}. Corresponding to the structure in \code{factors}, this object is a list, where the kth entry corresponds to k+1-way interactions. Each list entry contains another list, with entries equal to the number of rows in the corresponding matrix in \code{factors}. Each of these list entries (for a fixed k) contains a k-dimensional array that specifies the parameters of the given k-order interaction. For instance, if we have a 3-way interaction (1, 2, 3) and all variables are binary, we have a 2 x 2 x 2 array specifiying the parameters for each of the 2^3 = 8 possible configurations. If all variables are continuous, we have a 1 x 1 x 1 array, so the interaction is specified by a single parameter. See the examples below for an illustration.
}
  \item{thresholds}{
A list with p entries corresponding to p variables in the model. Each entry contains a vector indicating the threshold for each category (for categorical variables) or a numeric value indicating the threshold/intercept (for continuous variables).
}
  \item{sds}{
A numeric vector with p entries, specifying the variances of Gaussian variables. If variables 6 and 13 are Gaussians, then the corresponding entries of \code{sds} have to contain the corresponding variances. Other entries are ignored.
}
  \item{type}{
p character vector indicating the type of variable for each column in \code{data}. "g" for Gaussian, "p" for Poisson, "c" of each variable.
}
  \item{level}{
p integer vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{N}{
Number of samples that should be drawn from the distribution.
}
  \item{nIter}{
Number of iterations in the Gibbs sampler until a sample is drawn.
}
  \item{pbar}{
  If \code{pbar = TRUE} a progress bar is shown. Defaults to \code{pbar = TRUE}.
}
  \item{divWarning}{
  \code{mgmsampler()} returns a warning message if the absolute value of a continuous variable the chain of the gibbs sampler is larger than \code{divWarning}. To our best knowledge there is no theory yet defining a parameter space that ensures a proper probability density and hence a converging chain. Defaults to \code{divWarning = 10^3}.
}

  \item{returnChains}{
  If \code{returnChains = TRUE}, the sampler provides the entire chain of the Gibbs sampler, for each sampled case. Can be used to check convergence of the Gibbs sampler. Defaults to \code{returnChains = FALSE}.}
}

\details{
We use a Gibbs sampler to sample from the join distribution introduced by Yang and colleageus (2014). Note that the contraints on the parameter space necessary to ensure that the joint distribution is normalizable are to our best knowledge unknown. Yang and colleagues (2014) give these constraints for a number of simple pairwise models. In practice, an "improper joint density" will lead to a sampling process that approaches infinity, and hence \code{mgmsampler()} will return \code{Inf} / \code{-Inf} values.
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
Haslbeck, J., & Waldorp, L. J. (2018). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

Yang, E., Baker, Y., Ravikumar, P., Allen, G. I., & Liu, Z. (2014, April). Mixed Graphical Models via Exponential Families. In AISTATS (Vol. 2012, pp. 1042-1050).

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\examples{

\dontrun{

# --------- Example 1: p = 10 dimensional Gaussian ---------

# ----- 1) Specify Model -----

# a) General Graph Info
p <- 10 # number of variables
type = rep("g", p) # type of variables
level = rep(1, 10) # number of categories for each variable (1 = convention for continuous)

# b) Define interactions
factors <- list()
factors[[1]] <- matrix(c(1,2,
                         1,3,
                         4,5,
                         7,8), ncol=2, byrow = T) # 4 pairwise interactions
interactions <- list()
interactions[[1]] <- vector("list", length = 4)

# all pairwise interactions have value .5
for(i in 1:4) interactions[[1]][[i]] <- array(.5, dim=c(1, 1))

# c) Define Thresholds
thresholds <- vector("list", length = p)
thresholds <- lapply(thresholds, function(x) 0 ) # all means are zero

# d) Define Variances
sds <- rep(1, p) # All variances equal to 1


# ----- 2) Sample cases -----

data <- mgmsampler(factors = factors,
                   interactions = interactions,
                   thresholds = thresholds,
                   sds = sds,
                   type = type,
                   level = level,
                   N = 500,
                   nIter = 100,
                   pbar = TRUE)


# ----- 3) Recover model from sampled cases -----

set.seed(1)
mgm_obj <- mgm(data = data$data,
               type = type,
               level = level,
               k = 2,
               lambdaSel = "EBIC",
               lambdaGam = 0.25)

mgm_obj$interactions$indicator # worked!



# --------- Example 2: p = 3 Binary model with one 3-way interaction ---------

# ----- 1) Specify Model -----

# a) General Graph Info
type = c("c", "c", "c")
level = c(2, 2, 2)

# b) Define Interaction
factors <- list()
factors[[1]] <- NULL # no pairwise interactions
factors[[2]] <- matrix(c(1,2,3), ncol=3, byrow = T) # one 3-way interaction

interactions <- list()
interactions[[1]] <- NULL
interactions[[2]] <- vector("list", length = 2)
# threeway interaction no1
interactions[[2]][[1]] <- array(0, dim = c(level[1], level[2], level[3]))
theta <- 2
interactions[[2]][[1]][1, 1, 1] <- theta # fill in nonzero entries
# thus: high probability for the case that x1 = x2 = x3 = 1

# c) Define Thresholds
thresholds <- list()
thresholds[[1]] <- rep(0, level[1])
thresholds[[2]] <- rep(0, level[2])
thresholds[[3]] <- rep(0, level[3])


# ----- 2) Sample cases -----

set.seed(1)
dlist <- mgmsampler(factors = factors,
                    interactions = interactions,
                    thresholds = thresholds,
                    type = type,
                    level = level,
                    N = 500,
                    nIter = 100,
                    pbar = TRUE)


# ----- 3) Check: Contingency Table -----

dat <- dlist$data
table(dat[,1], dat[,2], dat[,3]) # this is what we expected


# ----- 4) Recover model from sampled cases -----

mgm_obj <- mgm(data = dlist$data,
               type = type,
               level = level,
               k = 3,
               lambdaSel = "EBIC",
               lambdaGam = 0.25, 
               overparameterize = TRUE)

mgm_obj$interactions$indicator # recovered, plus small spurious pairwise 1-2


# --------- Example 3: p = 5 Mixed Graphical Model with two 3-way interaction ---------

# ----- 1) Specify Model -----

# a) General Graph Info
type = c("g", "c", "c", "g")
level = c(1, 3, 5, 1)
# b) Define Interaction
factors <- list()
factors[[1]] <- NULL # no pairwise interactions
factors[[2]] <- matrix(c(1,2,3,
                         2,3,4), ncol=3, byrow = T) # no pairwise interactions
interactions <- list()
interactions[[1]] <- NULL
interactions[[2]] <- vector("list", length = 2)
# 3-way interaction no1
interactions[[2]][[1]] <- array(0, dim = c(level[1], level[2], level[3]))
interactions[[2]][[1]][,,1:3] <- rep(.8, 3) # fill in nonzero entries
# 3-way interaction no2
interactions[[2]][[2]] <- array(0, dim = c(level[2], level[3], level[4]))
interactions[[2]][[2]][1,1,] <- .3
interactions[[2]][[2]][2,2,] <- .3
interactions[[2]][[2]][3,3,] <- .3
# c) Define Thresholds
thresholds <- list()
thresholds[[1]] <- 0
thresholds[[2]] <- rep(0, level[2])
thresholds[[3]] <- rep(0, level[3])
thresholds[[4]] <- 0
# d) Define Variances
sds <- rep(.1, length(type))


# ----- 2) Sample cases -----

set.seed(1)
data <- mgmsampler(factors = factors,
                   interactions = interactions,
                   thresholds = thresholds,
                   sds = sds,
                   type = type,
                   level = level,
                   N = 500,
                   nIter = 100,
                   pbar = TRUE)


# ----- 3) Check: Conditional Means -----

# We condition on the categorical variables and check whether
# the conditional means match what we expect from the model:

dat <- data$data

# Check interaction 1
mean(dat[dat[,2] == 1 & dat[,3] == 1, 1]) # (compare with interactions[[2]][[1]])
mean(dat[dat[,2] == 1 & dat[,3] == 5, 1])
# first mean higher, ok!

# Check interaction 2
mean(dat[dat[,2] == 1 & dat[,3] == 1, 4]) # (compare with interactions[[2]][[2]])
mean(dat[dat[,2] == 1 & dat[,3] == 2, 4])
# first mean higher, ok!


}

}

