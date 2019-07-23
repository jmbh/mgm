\name{condition}
\alias{condition}

\title{
Computes mgm object conditional on a set of variables
}
\description{
The function takes an mgm object and a set of variables fixed to given values as input and returns the conditional mgm object.
}
\usage{
condition(object, values)
}

\arguments{
  \item{object}{
An mgm object, the output of the \code{mgm()} function.
}
  \item{values}{
A list, where the entry name indicates the variable that should be fixed, and the entry value indicates the value to which the corresponding variable should be fixed. See below for an example.
}
}
\value{
The function returns an mgm object that is conditional on the provided values. The new mgm object can again be used as input in \code{predict()}, \code{print()}, \code{showInteraction()}, etc.. 
}

\details{
The new conditional object still contains the variables that were fixed, however, they are not related to any of the random variables anymore. We kept the variables in the object to avoid confusion with variable labels and plotting. Also note that \code{mgm()} by default scales all Gaussian variables to mean=0, sd=1. Thus, fixed values should be selected based on the scaled version of variables.
}

\references{
Haslbeck, J., & Waldorp, L. J. (2019). mgm: Estimating time-varying mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.06871.
}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}



\seealso{
\code{mgm}
}


\examples{

\dontrun{

# --- Create Mixture of two Gaussians ---

set.seed(1)
n <- 500
library(MASS)

# Component A
Sigma_a <- diag(2)
Sigma_a[1, 2] <- Sigma_a[2, 1] <- .5
Xa <- mvrnorm(n = n, mu = rep(0, 2), Sigma = Sigma_a)

# Component B
Sigma_b <- diag(2)
Sigma_b[1, 2] <- Sigma_b[2, 1] <- 0
Xb <- mvrnorm(n = n, mu = rep(0, 2), Sigma = Sigma_b)

data <- as.data.frame(cbind(rbind(Xa, Xb), c(rep(0, n), rep(1, n))))
colnames(data) <- c("x1", "x2", "x3")


# --- Fit MGM ---

# with mgm
mgm_obj <- mgm(data = data, 
               type = c("g","g","c"), 
               level = c(1, 1, 2), 
               moderator = c(3), 
               lambdaSel = "EBIC")

# --- Condition on / fix values of variable 3 ---

# Fix x3=0
mgm_obj_x3.0 <- condition(object = mgm_obj, 
                          values = list("3"=0))
mgm_obj_x3.0$pairwise$wadj

# Fix x3=1
mgm_obj_x3.1 <- condition(object = mgm_obj, 
                          values = list("3"=1))
mgm_obj_x3.1$pairwise$wadj

}

}












