\name{tvmgm}
\alias{tvmgm}
\alias{tv.mgmfit}

\title{
Estimating time-varying Mixed Graphical Models
}

\description{
Estimates time-varying k-degree Mixed Graphical Models (MGMs) via elastic-net regularized kernel smoothed Generalized Linear Models
}

\usage{
tvmgm(data, type, level, timepoints, estpoints, bandwidth, ...)
}


\arguments{
  \item{data}{
A n x p data matrix.
}
  \item{type}{
p vector indicating the type of variable for each column in \code{data}. 'g' for Gaussian, 'p' for Poisson, 'c' for categorical.
}
  \item{level}{
p vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{timepoints}{
A strictly increasing numeric vector of length \code{nrow(data)} indicating time points for the measurements in \code{data}. If \code{timepoints} is not specified, it is assumed that the time points are equally spaced. For details, see Haslbeck and Waldorp (2016).
}
  \item{estpoints}{
Vector indicating estimation points on interval [0, n-max(lags)]. Lags are passed to via the \code{...} arguments. See \code{?mvar} for how to specify lags. n-max(lags) is the number of rows in the design matrix in a VAR model: the first max(lags) rows of the data matrix have to be omitted, otherwise we have no data to predict the first row.
}
  \item{bandwidth}{
We use a gaussian density on the unit time-interval [0,1] to determine the weights for each observation at each estimated time point. The bandwidth specifies the standard deviation the Gaussian density. To get some intuition, which bandwidth results in the combination of how many data close in time one can plot Gaussians on [0,1] for different bandwidths. The bandwidth can also be selected in a data driven way using the function (see \code{bwSelect}).

}
  \item{\dots}{
Arguments passed to \code{mgm}, specifying how each single model should be estimated. See \code{?mgm}.
}
}


\details{
Estimates a sequence of MGMs at the time points specified at the locations specified via \code{estpoints}. \code{tvmgm()} is a wrapper around \code{mgm()} and estimates a series of MGM with different weightings which are defined by the estimation locations in \code{estpoints} and the banwdith parameter specified in \code{bandwidth}. For details see Haslbeck and Waldorp (2016).


}
\value{
A list with the following entries:

\item{call}{
Contains all provided input arguments. If \code{saveData = TRUE}, it also contains the data.
}


\item{pairwise}{
Contains a list with all information about estimated pairwise interactions. \code{wadj} contains a p x p x estpoints array containing the weighted adjacency matrix for each estimation point specified in \code{estpoints}, if p is the number of variables in the network. \code{signs} has the same dimensions as \code{wadj} and contains the signs for the entries of \code{wadj}: 1 indicates a positive sign, -1 a negative sign and 0 an undefined sign. A sign is undefined if an edge is a function of more than one parameter. This is the case for interactions involving a categorical variable with more than 2 categories. \code{edgecolor} also has the same dimensions as \code{wadj} contains a color for each edge, depending on \code{signs}. It is provided for more convenient plotting. If only pairwise interactions are modeled (\code{d = 1}), \code{wadj} contains all conditional independence relations.
}

\item{factorgraph}{
Contains a list with all information about all estimated interactions (factors). \code{graph} contains an (p+E) x (p+E) x estpoints array containing weighted adjacency matrix of a (bipartide) factor graph for each estimation point specified in \code{estpoints}. If p is the number of variables and E the number of interactions (factors) in the model, each of these matrices (for a fixed estimation point) has dimensions (p+E) x (p+E). The factor graph is furter specified by the following objects: \code{signs} is a matrix of the same dimensions as \code{graph} that indicates the sign of each interaction, if defined (see \code{pairwise} above). \code{edgecolor} is a matrix with the same dimension as \code{graph} that provides edge colors depending on the sign as above. \code{order} contains a list with \code{length(estpoints)} entries; each element is a (p+E) vector indicating the order of interaction. The first p entries are set to zero and correspond to the variables. \code{nodetype} contains a list with \code{length(estpoints)} entries; each element is a (p+E) vector indicating whether a row/column corresponds to a variable (0) or a factor (1). The last two vectors are helpful to plot the factor graph (see examples below). Note that \code{factorgraph} contains the same information as \code{pairwise} if only pairwise interactions are estimated.
}

\item{rawfactor}{
Contains a list with entries equal to the number of estimated time points specified in \code{estpoints}; each entry is a list with three entries that relate each interaction in the model to all its parameters. This is different to the output provided in \code{factorgraph}, where one value is assigned to each interaction. \code{indicator} contains a list with d entries, one for each order of modeled interaction, which contain the estimated (nonzero) interactions. \code{weights} contains a list with d entries, which in turn contain R lists, where R is the number of interactions (and rows in the corresponding list entry in\code{indicator}) that were estimated (nonzero) in the given entry. \code{signs} has the same structure as \code{weights} and provides the sign of the interaction, if defined.
}

\item{intercepts}{
Contains a list with entries equal to the number of estimated time points specified in \code{estpoints}; each entry is a list with p entries, which contain the intercept/thresholds for each node in the network. In case a given node is categorical with m categories, there are m thresholds for this variable.
}


\item{tvmodels}{
Contains the MGM model estimated by \code{mgm()} at each time point specified via \code{estpoints}. See \code{?mgm} for a detailed description of this output.
}

}


\references{

Haslbeck, J., & Waldorp, L. J. (2016). mgm: Structure Estimation for time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

Haslbeck, J., & Waldorp, L. J. (2015). Structure estimation for mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.05677.

}

\author{

Jonas Haslbeck <jonashaslbeck@gmail.com>

}



\examples{

\dontrun{


## We specify a time-varying MGM and recover it using tvmgm()

# 1) Specify Model

# a) Define Graph
p <- 6
type = c('c', 'c', 'g', 'g', 'p', 'p')
level = c(2, 3, 1, 1, 1, 1)
n_timepoints <- 1000

# b) Define Interaction
factors <- list()
factors[[1]] <- matrix(c(1,2,
                         2,3,
                         3,4), ncol=2, byrow = T)  # no pairwise interactions
factors[[2]] <- matrix(c(1,2,3,
                         2,3,4), ncol=3, byrow = T) # one 3-way interaction

interactions <- list()
interactions[[1]] <- vector('list', length = 3)
interactions[[2]] <- vector('list', length = 2)
# 3 2-way interactions
interactions[[1]][[1]] <- array(0, dim = c(level[1], level[2], n_timepoints))
interactions[[1]][[2]] <- array(0, dim = c(level[2], level[3], n_timepoints))
interactions[[1]][[3]] <- array(0, dim = c(level[3], level[4], n_timepoints))
# 2 3-way interactions
interactions[[2]][[1]] <- array(0, dim = c(level[1], level[2], level[3], n_timepoints))
interactions[[2]][[2]] <- array(0, dim = c(level[2], level[3], level[4], n_timepoints))
theta <- .3
interactions[[1]][[1]][1, 1, ] <- theta
interactions[[1]][[2]][1, 1, ] <- theta
interactions[[1]][[3]][1, 1, ] <- seq(0, theta, length = n_timepoints)
interactions[[2]][[1]][1, 1, 1, ] <- theta
interactions[[2]][[2]][1, 1, 1, ] <- theta
# c) Define Thresholds
thresholds <- list()
thresholds[[1]] <- matrix(0, nrow = n_timepoints, ncol= level[1])
thresholds[[2]] <- matrix(0, nrow = n_timepoints, ncol= level[2])
thresholds[[3]] <- matrix(0, nrow = n_timepoints, ncol= level[3])
thresholds[[4]] <- matrix(0, nrow = n_timepoints, ncol= level[4])
thresholds[[5]] <- matrix(.1, nrow = n_timepoints, ncol= level[5])
thresholds[[6]] <- matrix(.1, nrow = n_timepoints, ncol= level[6])
# d) define sds
sds <- matrix(.2, ncol=p, nrow=n_timepoints)

# 2) Sample Data
iter <- 1
set.seed(iter)
d_iter <- tvmgmsampler(factors = factors,
                       interactions = interactions,
                       thresholds = thresholds,
                       sds = sds,
                       type = type,
                       level = level,
                       nIter = 100,
                       pbar = TRUE)

data <- d_iter$data
head(data)
# delete inf rows:
ind_finite <- apply(data, 1, function(x) if(all(is.finite(x))) TRUE else FALSE)
table(ind_finite) # all fine for this setup & seed
# in case of inf values (no theory on how to keep k-order MGM well-defined)
data <- data[ind_finite, ] 


# 3) Recover

mgm_c_cv <- tvmgm(data = data,
                  type = type,
                  level = level,
                  k = 3,
                  estpoints = seq(0, 1000, length=10),
                  bandwidth = .1,
                  lambdaSel = 'CV',
                  ruleReg = 'AND',
                  pbar = TRUE,
                  overparameterize = T,
                  signInfo = FALSE)

# Look at time-varying pairwise parameter 3-4
mgm_c_cv$pairwise$wadj[3,4,] # recovers increase


# 4) Predict values / compute nodewise Errors

pred_mgm_cv_w <- predict.mgm(mgm_c_cv,
                             data = data,
                             tvErrorType = 'weighted')
pred_mgm_cv_cM <- predict.mgm(mgm_c_cv,
                              data = data,
                              tvErrorType = 'closestModel')

pred_mgm_cv_w$errors
pred_mgm_cv_cM$errors



}

}



