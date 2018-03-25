\name{mgm}
\alias{mgm}
\alias{mgmfit}

\title{
Estimating Mixed Graphical Models
}
\description{
Function to estimate k-degree Mixed Graphical Models via neighborhood regression.
}
\usage{
mgm(data, type, level, lambdaSeq, lambdaSel, lambdaFolds,
    lambdaGam, alphaSeq, alphaSel, alphaFolds, alphaGam, 
    k, moderators, ruleReg, weights, threshold, method, 
    binarySign, scale, verbatim, pbar, warnings, saveModels, 
    saveData, overparameterize, thresholdCat, signInfo, ...)
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
  \item{k}{
Order up until including which interactions are included in the model. \code{k = 2} means that all pairwise interactions are included, \code{k = 3} means that all pairwise and all three-way interactions are included, etc. In previous versions of mgm the order of interactions was specified by the parameter d, the largest size or a neighborhood. Note that k = d + 1.
}

  \item{moderators}{
Integer vector with elements in \code{1:p}, specifying moderation effects to be included in the model. For instance, \code{moderators = c(4)} includes all linear moderation effects of variable 4. This is equivalent to including all 3-way interactions that include variable 4. Note that \code{moderators = 1:p} gives the same model as setting \code{k = 3} (see previous argument).
}

  \item{ruleReg}{
Rule to combine estimates from neighborhood regression. E.g. for pairwise interactions, two estimates (one from regressing A on B and one from B on A) have to combined in one edge parameter. \code{ruleReg = 'AND'} requires all estiamtes to be nonzero in order to set the edge to be present. \code{ruleReg = 'OR'} requires at least one estiamte to be nonzero in order to set the edge to be present. For higher order interactions, \code{k} estiamtes have to be combined with this rule.
}
  \item{weights}{
A n vector with weights for observations.
}
  \item{threshold}{
  A threshold below which edge-weights are put to zero. This is done in order to guarantee a constant false-positive rate. \code{threshold = 'LW'} refers to the threshold in Loh and Wainwright (2013), which was used in all previous versions of mgm. \code{threshold = 'HW'} refers to the threshold in Haslbeck and Waldorp (2016). If \code{threshold = 'none'} no thresholding is applied. Defaults to \code{threshold = 'LW'}.
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
If \code{saveModels = FALSE}, only information about the weighted adjacency matrix, and if k > 2 about the factor graph is provided in the output list. If \code{saveModels = TRUE}, all fitted parameters are additionally returned.
}
  \item{saveData}{
If \code{saveData = TRUE}, the data is saved in the output list. Defaults to \code{saveData = FALSE}.
}
  \item{overparameterize}{
If \code{overparameterize = TRUE}, \code{mgm()} uses over-parameterized design-matrices for each neighborhood regression; this means that an interaction between two categorical variables with m and s categories is parameterized by m*s parameters. If \code{overparameterize = FALSE} the standard parameterization (in \code{glmnet}) with m*(s-1) parameters is used, where the first category of the predicting variable serves as reference category. If all variables are continuous both parameterizations are the same. Note that the default is set to \code{overparameterize = FALSE}, to be consistent with the previous mgm versions. However note that in order to recover higher order interactions categorical variables the overparameterized parameterization is often advantageous. See the examples below for an illustration. Note that we can estimate the model despite the colinear columns in the design matrix because we use penalized regression.
}

  \item{thresholdCat}{
If \code{thresholdCat = FALSE}, the thresholds of categorical variables are set to zero. Defaults to \code{thresholdCat = TRUE} for which the thresholds are esimated.
}

  \item{signInfo}{
If \code{signInfo = TRUE}, a message is shown in the console, indicating that the sign of estimates is stored separately. Defaults to \code{signInfo = TRUE}.
}


  \item{...}{
Additional arguments.
}


}




\details{
\code{mgm()} estimates an exponential mixed graphical model as introduced in Yang and colleagies (2014). Note that MGMs are not normalizable for all parameter values. See Chen, Witten & Shojaie (2015) for an overview of when pairwise MGMs are normalizable. To our best knowledge, for MGMs with interactions of order > 2 that include non-categorical variables, the conditions for normalizablity are unknown.
}


\value{

The function returns a list with the following entries:

\item{call}{
Contains all provided input arguments. If \code{saveData = TRUE}, it also contains the data
}

\item{pairwise}{
Contains a list with all information about estimated pairwise interactions. \code{wadj} contains the p x p weighted adjacency matrix, if p is the number of variables in the network. \code{signs} has the same dimensions as \code{wadj} and contains the signs for the entries of \code{wadj}: 1 indicates a positive sign, -1 a negative sign and 0 an undefined sign. A sign is undefined if an edge is a function of more than one parameter. This is the case for interactions involving a categorical variable with more than 2 categories. \code{edgecolor} also has the same dimensions as \code{wadj} contains a color for each edge, depending on \code{signs}. It is provided for more convenient plotting. If only pairwise interactions are modeled (\code{d = 1}), \code{wadj} contains all conditional independence relations.
}

\item{interactions}{
A list with three entries that relate each interaction in the model to all its parameters. This is different to the output provided in \code{factorgraph}, where one value is assigned to each interaction. \code{indicator} contains a list with k-1 entries, one for each order of modeled interaction, which contain the estimated (nonzero) interactions. \code{weightsAgg} contains a list with k-1 entries, which in turn contain R lists, where R is the number of interactions (and rows in the corresponding list entry in\code{indicator}) that were estimated (nonzero) in the given entry. Each of these entries contains the mean of the absolute values of all parameters involved in this interaction. \code{weights} has the same structure as \code{weightsAgg}, but does contain all parameters involved in the interaction instead of the mean of their absolute values. \code{signs} has the same structure as \code{weightsAgg}/\code{weights} and provides the sign of the interaction, if defined.
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

Haslbeck, J., & Waldorp, L. J. (2017). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

Loh, P. L., & Wainwright, M. J. (2012, December). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. In NIPS (pp. 2096-2104).

Chen S, Witten DM & Shojaie (2015). Selection and estimation for mixed graphical models. Biometrika, 102(1), 47.

Yang, E., Baker, Y., Ravikumar, P., Allen, G. I., & Liu, Z. (2014, April). Mixed Graphical Models via Exponential Families. In AISTATS (Vol. 2012, pp. 1042-1050).

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\examples{

\dontrun{


## We fit a pairwise and 3-order MGM to the mixed Autism dataset (?autism_data)

# 1) Fit Pairwise MGM

# Call mgm()
fit_k2 <- mgm(data = autism_data$data,
              type = autism_data$type,
              level = autism_data$lev,
              k = 2) # ad most pairwise interacitons

# Weighted adjacency matrix
fit_k2$pairwise$wadj

# Visualize using qgraph()
library(qgraph)
qgraph(fit_k2$pairwise$wadj,
       edge.color = fit_k2$pairwise$edgecolor,
       layout = 'spring',
       labels =  autism_data$colnames)


# 2) Fit MGM with pairwise & threeway interactions
fit_k3 <- mgm(data = autism_data$data,
              type = autism_data$type,
              level = autism_data$lev,
              k = 3) # include all interactions up to including order 3

# List of estimated interactions
fit_k3$interactions$indicator

# 3) Plot Factor Graph 

nodeColors <- c("white", "orange", "lightblue")

FactorGraph(object = fit_k3, 
            PairwiseAsEdge = TRUE, 
            DoNotPlot = FALSE, 
            colors = nodeColors, 
            labels = 1:7, 
            layout="spring")

# 4) Predict values

pred_obj <- predict(fit_k3, autism_data$data)

head(pred_obj$predicted) # predicted values

pred_obj$errors # Nodewise errors



## Here we illustrate why we need to overparameterize the design matrix to 
## recover higher order interactions including categorical variables

# 1) Define Graph (one 3-way interaction between 3 binary variables)

# a) General Graph Info
type = c('c', 'c', 'c')
level = c(2, 2, 2)
# b) Define Interaction
factors <- list()
factors[[1]] <- NULL # no pairwise interactions
factors[[2]] <- matrix(c(1,2,3), ncol=3, byrow = T) # one 3-way interaction
interactions <- list()
interactions[[1]] <- NULL
interactions[[2]] <- vector('list', length = 1)
# threeway interaction no1
interactions[[2]][[1]] <- array(0, dim = c(level[1], level[2], level[3]))
theta <- .7
interactions[[2]][[1]][1, 1, 1] <- theta  #weight theta for conf (1,1,1), weight 0 for all others
# c) Define Thresholds
thresholds <- list()
thresholds[[1]] <- c(0, 0)
thresholds[[2]] <- c(0, 0)
thresholds[[3]] <- c(0, 0)


# 2) Sample from Graph

iter <- 1
set.seed(iter)
N <- 2000
d_iter <- mgmsampler(factors = factors,
                     interactions = interactions,
                     thresholds = thresholds,
                     type = type,
                     level = level,
                     N = N,
                     nIter = 50,
                     pbar = TRUE)


# 3.1) Estimate order 3 MGM usind standard parameterization

d_est_stand <- mgm(data = d_iter$data,
                   type = type,
                   level = level,
                   k = 3,
                   lambdaSel = 'CV',
                   ruleReg = 'AND',
                   pbar = TRUE, 
                   overparameterize = FALSE, 
                   signInfo = FALSE)

# We look at the nodewise regression for node 1 (same for all)
coefs_stand <- d_est_stand$nodemodels[[1]]$model
coefs_stand 
# we see that nonzero-zero pattern of parameter vector does not allow us to infer whether 
# interactions are present or not


# 3.2) Estimate order 3 MGM usind overparameterization

d_est_over <- mgm(data = d_iter$data,
                  type = type,
                  level = level,
                  k = 3,
                  lambdaSel = 'CV',
                  ruleReg = 'AND',
                  pbar = TRUE, 
                  overparameterize = TRUE, 
                  signInfo = FALSE)

# We look at the nodewise regression for node 1 (same for all)
coefs_over <- d_est_over$nodemodels[[1]]$model
coefs_over # recovers exactly the 3-way interaction


# For more examples see https://github.com/jmbh/mgmDocumentation


}


}
