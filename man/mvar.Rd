\name{mvar}
\alias{mvar}
\alias{var.mgm}

\title{
Estimating mixed Vector Autoregressive Model (mVAR)
}
\description{
Estimates Mixed Vector Autoregressive Model (mVAR) via elastic-net regularized Generalized Linear Models
}
\usage{
mvar(data, type, level, lambdaSeq, lambdaSel, lambdaFolds,
lambdaGam, alphaSeq, alphaSel, alphaFolds, alphaGam, lags,
consec, weights, threshold, method, binarySign, scale, verbatim,
pbar, warnings, saveModels, saveData, overparameterize, signInfo, ...)
}

\arguments{
  \item{data}{
n x p data matrix
}
  \item{type}{
p vector indicating the type of variable for each column in \code{data}. 'g' for Gaussian, 'p' for Poisson, 'c' for categorical.
}
  \item{level}{
p vector indicating the number of categories of each variable. For continuous variables set to 1.
}
  \item{lambdaSeq}{
A sequence of lambdas that should be searched (see also \code{lambdaSel}). Defaults to \code{NULL}, which lets \code{glmnet} select a reasonable lambda sequence (recommended).
}
  \item{lambdaSel}{
Procedure to select the tuning parameter controlling the Lq-penalization. The two options are cross validation 'CV' and the Extended Bayesian Information Criterion (EBIC) 'EBIC'. The EBIC performs well in selecting sparse graphs (see Barber and Drton, 2010 and Foygel and Drton, 2014). Note that when also searching the alpha parameter in the elastic net penalty, cross validation should be preferred, as the parameter vector will not necessarily be sparse anymore. The EBIC tends to be a bit more conservative than CV (see Haslbeck and Waldorp, 2016). CV can sometimes not be performed with categorical variables, because \code{glmnet} requires at least 2 events of each category of each categorical variable in each training-fold. Defaults to \code{lambdaSel = 'CV'}.
}
  \item{lambdaFolds}{
Number of folds in cross validation if \code{lambdaSel = 'CV'}
}
  \item{lambdaGam}{
Hyperparameter gamma in the EBIC if \code{lambdaSel = 'EBIC'}. Defaults to \code{lambdaGam = .25}.
}
  \item{alphaSeq}{
A sequence of alpha parameters for the elastic net penality in [0,1] that should be searched (see also \code{alphaSel}). Defaults to \code{alphaSeq = 1}, which means that the lasso is being used. \code{alphaSeq = 0} corresponds to an L2-penalty (Ridge regresion). For details see Friedman, Hastie and Tibshirani (2010).
}
  \item{alphaSel}{
Procedure to select the alpha parameter in the elastic net penalty. The two options are cross validation 'CV' and the Extended Bayesian Information Criterion (EBIC) 'EBIC'. The EBIC performs well in selecting sparse graphs (see Barber and Drton, 2010 and Foygel and Drton, 2014). Note that when also searching the alpha parameter in the elastic net penalty, cross validation should be preferred, as the parameter vector will not necessarily be sparse anymore. The EBIC tends to be a bit more conservative than CV (see Haslbeck and Waldorp, 2016). CV can sometimes not be performed with categorical variables, because \code{glmnet} requires at least 2 events of each category of each categorical variable in each training-fold. Defaults to \code{alphaSel = 'CV'}.

}
  \item{alphaFolds}{
Number of folds in cross validation if \code{alphaSel = 'CV'}
}
  \item{alphaGam}{
Hyperparameter gamma in the EBIC if \code{alphaSel = 'EBIC'}. Defaults to \code{alphaGam = .25}.
}
  \item{lags}{
Vector of positive integers indicating the lags included in the mVAR model (e.g. 1:3 or c(1,3,5))
}
  \item{consec}{
An integer vector of length n, indicating the measurement points of the rows in \code{data}. The time difference of 1 is what we model as a lag with size 1 in the VAR model. This vector is used to delete all rows from the VAR design matrix, for which no \code{max(lags)} complete previous measurements are available. This is useful in many applications where some time points are completely missing, or gaps in the time series have structural reasons (respondents only respond during the hours they are awake). The 'trimmed' dataset is returned in \code{call$data_lagged} if \code{saveData = TRUE}. Defaults to \code{consec = NULL}, which assumes that all measurements are consecutive, i.e. \code{consec = 1:n}. In this case only the first \code{max(lags)} lags are excluded to obtain the VAR design matrix.
}
  \item{weights}{
A vector with n - max(lags) entries, indicating the weight for each observation. The mVAR design matrix has with n - max(lags) rows, because the first row must be predictable by the highest lag. The weights have to be on the scale [0, n - max(lags) ].
}
  \item{threshold}{
  A threshold below which edge-weights are put to zero. This is done in order to guarantee a constant false-positive rate. \code{threshold = 'LW'} refers to the threshold in Loh and Wainwright (2013), which was used in all previous versions of mgm. \code{threshold = 'HW'} refers to the threshold in Haslbeck and Waldorp (2016). If \code{threshold = 'none'} no thresholding is applied.  Defaults to \code{threshold = 'LW'}.

}
  \item{method}{
Estimation method, currently only \code{method = 'glm'}.
}
  \item{binarySign}{
If \code{binarySign = TRUE} a sign for the interaction within binary nodes and between binary and continuous nodes is provided in the output. Note that in this case the two categories of the binary variables have to be coded in {0,1}. This is to ensure that the interpretation of the sign is unambigous: a positive sign of a parameter means that increasing the associated predictor results in a higher probability for category 1.
}


  \item{scale}{
If \code{scale = TRUE}, all Gaussian nodes (specified by \code{'g'} in \code{type}) are centered and divided by their standard deviation. Scaling is recommended, because otherwise the penalization of a parameter depends on the variance of the associated predictor.
}

  \item{verbatim}{
If \code{verbatim = TRUE} no warnings and no progress bar is shown. Defaults to \code{verbatim = FALSE}.
}
  \item{pbar}{
If \code{pbar = TRUE} a progress bar is shown. Defaults to \code{pbar = TRUE}.
}
  \item{warnings}{
If \code{warnings = TRUE} no warnigns are returned. Defaults to \code{warnings = FALSE}.
}
  \item{saveModels}{
If \code{saveModels = FALSE}, only information about the weighted adjacency matrix, and if d > 1 about the factor graph is provided in the output list. If \code{saveModels = TRUE}, all fitted parameters are additionally returned.
}
  \item{saveData}{
If \code{saveData = TRUE}, the data is saved in the output list. Defaults to \code{saveData = FALSE}.
}

  \item{overparameterize}{
If \code{overparameterize = TRUE}, \code{mgm()} uses over-parameterized design-matrices for each neighborhood regression; this means that a cross-lagged effect between two categorical variables with m and s categories is parameterized by m*s parameters. If \code{overparameterize = FALSE} the standard parameterization (in \code{glmnet}) with m*(s-1) parameters is used, where the first category of the predicting variable serves as reference category. If all variables are continuous both parameterizations are the same. The default is set to \code{overparameterize = FALSE}.
}


  \item{signInfo}{
If \code{signInfo = TRUE}, a message is shown in the console, indicating that the sign of estimates is stored separately. Defaults to \code{signInfo = TRUE}.
}

  \item{...}{
Additional arguments.
}

}



\details{
\code{mgm()} estimates an exponential mixed graphical model as introduced in Yang and colleagies (2014). Haslbeck and Waldorp (2015) show that the estimation procedure is consistent. See Haslbeck and Waldorp (2016) for details about how to rewrite the exponential MGM as a mVAR model.
}



\value{

The function returns a list with the following entries:

\item{call}{
Contains all provided input arguments. If \code{saveData = TRUE}, it also contains the data.
}


\item{wadj}{
A p x p x n_lags array, in which rows are predicted by columns, i.e. entry \code{wadj[1, 2, 4]} corresponds to the parameter(s) of variable 2 at time point t predicting variable 1 at time point t - z, where z is the fourth specified lag in \code{lags} and n_lags is the number of specified lags in \code{lags}. For interactions that involve more than two parameters (e.g. always for categorical variables with more than 2 categories), we take the arithmetic mean of the absolute value of all parameters. The full set of estimated parameters is saved in \code{rawlags} (see below).
}

\item{signs}{
A p x p x n_lags array, specifying the signs corresponding to the entries of \code{wadj} (if defined), where n_lags is the number of specified lags in \code{lags}. 1/-1 indicate positive and negative relationships, respectively. 0 indicates that no sign is defined, which is the case for interactions that involve a categorical variable where an interaction can have more than one parameter. If \code{binarySign = TRUE}, a sign is calculated for interactions between binary variables and binary and continuous variables, where the interaction is still defined by one parameter and hence a sign can be specified. \code{NA} indicates that the corresponding parameter in \code{wadj} is zero.
}

\item{edgecolor}{
A p x p x n_lags array of colors indicating the sign of each parameter. This array contains the same information is \code{signs} and is included for convenient plotting.
}

\item{rawlags}{
List with entries equal to the number of specified lags in \code{lags}. Each entry is a nested list, with each p entries: the first level indicates the predicted variable, the second level the predictor variable. In case of categorical variables, interactions have more than one parameter.
}

\item{intercepts}{
A list with p entries, which contain the intercept/thresholds for each node in the network. In case a given node is categorical with m categories, there are m thresholds for this variable.
}

\item{nodemodels}{
A list with p \code{glmnet()} models, from which all above output is computed. Also contains the coefficients \code{models} for the selected lambda and the applied tau threshold \code{tau}.
}


}




\references{


Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. In Advances in neural information processing systems (pp. 604-612).

Barber, R. F., & Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. Electronic Journal of Statistics, 9(1), 567-607.

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of statistical software, 33(1), 1.

Haslbeck, J., & Waldorp, L. J. (2015). Structure estimation for mixed graphical models in high-dimensional data. arXiv preprint arXiv:1510.05677.

Haslbeck, J., & Waldorp, L. J. (2016). mgm: Structure Estimation for time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

Loh, P. L., & Wainwright, M. J. (2012, December). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. In NIPS (pp. 2096-2104).

Yang, E., Baker, Y., Ravikumar, P., Allen, G. I., & Liu, Z. (2014, April). Mixed Graphical Models via Exponential Families. In AISTATS (Vol. 2012, pp. 1042-1050).



}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}



\examples{

\dontrun{


## We generate data from a mixed VAR model and then recover the model using mvar()

# 1) Define mVAR model

p <- 6 # Six variables
type <- c('c', 'c', 'c', 'c', 'g', 'g') # 4 categorical, 2 gaussians
level <- c(2, 2, 4, 4, 1, 1) # 2 categoricals with m=2, 2 categoricals with m=4, two continuous
max_level <- max(level)

lags <- c(1, 3, 9) # include lagged effects of order 1, 3, 9
n_lags <- length(lags)

# Specify thresholds
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
                 lambdaSel = 'CV',
                 lags = c(1, 3, 9),
                 signInfo = FALSE,
                 overparameterize = F)

# Did we recover the true parameters?
mvar_obj$wadj[5, 6, 2] # cross-lagged effect of 6 on 2 over lag lags[2]
mvar_obj$wadj[1, 3, 1] # cross-lagged effect of 3 on 1 over lag lags[1]
mvar_obj$wadj[1, 5, 3] # cross-lagged effect of 1 on 5 over lag lags[3]

# How to get the exact parameter estimates?
# Example: the full parameters for the crossed-lagged interaction of 2 on 1 over lag lags[1]
mvar_obj$rawlags[[1]][[1]][[2]] 



# 4) Predict / Compute nodewise Error

pred_mvar <- predict.mgm(mvar_obj, dlist$data)

head(pred_mvar$predicted) # predicted values

pred_mvar$errors # Nodewise errors

# For more examples see https://github.com/jmbh/mgmDocumentation

}



}




