
\name{bwSelect}
\alias{bwSelect}

\title{
Select optimal bandwidth for time-varying MGMs and mVAR Models
}

\description{
Selects the bandwidth parameter with lowest out of sample prediction error for MGMs and mVAR Models.
}

\usage{
bwSelect(data, type, level, bwSeq, bwFolds,
         bwFoldsize, modeltype, pbar, ...)
}

\arguments{
  \item{data}{
A n x p data matrix.
}
  \item{type}{
p vector indicating the type of variable for each column in \code{data}. "g" for Gaussian, "p" for Poisson, "c" for categorical.
}
  \item{level}{
p vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{bwSeq}{
A sequence with candidate bandwidth values (0, s] with s < Inf. Note that the bandwidth is applied relative to the unit time interval [0,1] and hence a banwidth of > 2 corresponds roughly to equal weights for all time points and hence gives similar estimates as the stationary model estimated via \code{mvar()}.
}
  \item{bwFolds}{
The number of folds (see details below).
}
  \item{bwFoldsize}{
The size of each fold (see details below).
}
  \item{modeltype}{
If \code{modeltype = "mvar"} model, the optimal bandwidth parameter for a \code{tvmvar()} model is selected. If \code{modeltype = "mgm"} model, the optimal bandwidth parameter for a \code{tvmgm()} model is selected. Additional arguments to \code{tvmvar()} or \code{tvmgm()} can be passed via the \code{\dots} argument.
}

  \item{pbar}{
If TRUE a progress bar is shown. Defaults to \code{pbar = "TRUE"}.
}

  \item{\dots}{
Arguments passed to \code{tvmgm} or \code{tvmvar}.
}
}

\details{

Performs a cross-validation scheme that is specified by \code{bwFolds} and \code{bwFoldsize}. In the first fold, the test set is defined by an equally spaced sequence between [1,  n - \code{bwFolds}] of length \code{bwFoldsize}. In the second fold, the test set is defined by an equally spaced sequence between [2, n - \code{bwFolds} + 1] of length \code{bwFoldsize}, etc. . Note that if \code{bwFoldsize} = n / \code{bwFolds}, this procedure is equal to \code{bwFolds}-fold cross valildation. However, full cross validation is computationally very expensive and a single split in test/training set by setting \code{bwFolds = 1} is sufficient in many situations. The procedure selects the bandwidth with the lowest prediction error, averaged over variables and time points in the test set.

\code{bwSelect} computes the absolute error (continuous) or 0/1-loss (categorical) for each time point in the test set defined by \code{bwFoldsize} as described in the previous paragraph for every fold specified in \code{bwFolds}, separately for each variable. The computed errors are returned in different levels of aggregation in the output list (see below). Note that continuous variables are scaled (centered and divided by their standard deviation), hence the absolute error and 0/1-loss are roughly on the scale scale.

Note that selecting the bandwidth with the EBIC is no alternative. This is because the EBIC always selects the intercept model with the lowest bandwidth.  The reason is that the unregularized intercept closely models the noise in the data and hence the penalty sets all other parameters to zero. This problem is solved by using out of sample prediction error in the cross validation scheme.


}
\value{

The function returns a list with the following entries:

  \item{call}{
Contains all provided input arguments. If \code{saveData = TRUE}, it also contains the data.
}

  \item{bwModels}{
Contains the models estimated at the time points in the tests set. For details see \code{tvmvar} or \code{tvmgm}.
}

  \item{fullErrorFolds}{
List with number of entries equal to the length of \code{bwSeq} entries. Each entry contains a list with \code{bwFolds} entries. Each of those entries contains a contains a \code{bwFoldsize} times p matrix of out of sample prediction errors.
}

  \item{fullError}{
The same as \code{fullErrorFolds} but pooled over folds.
}

  \item{meanError}{
List with number of entries equal to the length of \code{bwSeq} entries. Each entry contains the average prediction error over variables and time points in the test set.
}

  \item{testsets}{
List with \code{bwFolds} entries, which contain the rows of the test sample for each fold.
}

  \item{zeroweights}{
List with \code{bwFolds} entries, which contains the observation weights used to fit the model at the \code{bwFoldsize} time points.
}


}
\references{

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. In Advances in neural information processing systems (pp. 604-612).

Barber, R. F., & Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. Electronic Journal of Statistics, 9(1), 567-607.

Haslbeck, J., & Waldorp, L. J. (2018). mgm: Estimating for time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.


}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}




\examples{

\dontrun{


## A) bwSelect for tvmgm() 

# A.1) Generate noise data set
p <- 5
n <- 100
data_n <- matrix(rnorm(p*n), nrow=100)
head(data_n)

type <- c("c", "c", rep("g", 3))
level <- c(2, 2, 1, 1, 1)
x1 <- data_n[,1]
x2 <- data_n[,2]
data_n[x1>0,1] <- 1
data_n[x1<0,1] <- 0
data_n[x2>0,2] <- 1
data_n[x2<0,2] <- 0

head(data_n)

# A.2) Estimate optimal bandwidth parameter

bwobj_mgm <- bwSelect(data = data_n,
                      type = type,
                      level = level,
                      bwSeq = seq(0.05, 1, length=3),
                      bwFolds = 1,
                      bwFoldsize = 3,
                      modeltype = "mgm",
                      k = 3,
                      pbar = TRUE,
                      overparameterize = TRUE)


print.mgm(bwobj_mgm)



## B) bwSelect for tvmVar() 

# B.1) Generate noise data set

p <- 5
n <- 100
data_n <- matrix(rnorm(p*n), nrow=100)
head(data_n)

type <- c("c", "c", rep("g", 3))
level <- c(2, 2, 1, 1, 1)
x1 <- data_n[,1]
x2 <- data_n[,2]
data_n[x1>0,1] <- 1
data_n[x1<0,1] <- 0
data_n[x2>0,2] <- 1
data_n[x2<0,2] <- 0

head(data_n)

# B.2) Estimate optimal bandwidth parameter

bwobj_mvar <- bwSelect(data = data_n,
                       type = type,
                       level = level,
                       bwSeq = seq(0.05, 1, length=3),
                       bwFolds = 1,
                       bwFoldsize = 3,
                       modeltype = "mvar",
                       lags = 1:3,
                       pbar = TRUE,
                       overparameterize = TRUE)


print.mgm(bwobj_mvar)

# For more examples see https://github.com/jmbh/mgmDocumentation


}

}


