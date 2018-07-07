\name{predict.mgm}
\alias{predict.mgm}



\title{
Compute predictions from mgm model objects
}
\description{
Computes predictions and prediction errors from a mgm model-object (\code{mgm}, \code{mvar}, \code{tvmgm} or \code{tvmvar}).

}
\usage{
\method{predict}{mgm}(object, data, errorCon, errorCat, 
        tvMethod, consec, beepvar, dayvar, \dots)
}

\arguments{
  \item{object}{
An mgm model object (the output of one of the functions \code{mgm()}, \code{mvar()}, \code{tvmgm()} or \code{tvmvar()})
}
  \item{data}{
A n x p data matrix with the same structure (number of variables p and types of variables) as the data used to fit the model.
}
  \item{errorCon}{
  Either a character vector specifying the types of nodewise errors that should be computed, where the two provided error functions for continuous varaibles are \code{errorCon = "RMSE"}, the Root Mean Squared Error, and \code{errorCon = "R2"}, the proportion of explained variance. The default is \code{errorCon = c("RMSE" "R2")}.

Alternatively, \code{errorCon} can be a list, where each list entry is a custom error function of the form \code{foo(true, pred)}, where \code{true} and \code{pred} are the arguments for the vectors of true and predicted values, respectively.  If predictions are made for a time-varying model and \code{tvMethod = "weighted"}, the weighted R2 or RMSE are computed. If a custom function is used, an additional argument for the weights has to be provided: \code{foo(true, pred, weights)}. Note that custom error functios can also be combined with the buildt-in functions, i.e. \code{errorCon = list("RMSE", "CustomError"=foo)}.}


\item{errorCat}{
Either a character vector specifying the types of nodewise errors that should be computed, where the two provided error functions for categorical variables are \code{errorCat = "CC"}, the proportion of correct classification (accuracy) and \code{errorCat = "nCC"}, the proportion of correct classification normalized by the marginal distribution of the variable at hand. Specifically, nCC = (CC - norm_constant) / (1 - norm_constant), where norm_constant is the highest relative frequency across categories. Another provided error is "CCmarg" which returns the accuracy of the intercept/marginal model. The default is to return all types of errors \code{errorCon = c("CC" "nCC", "CCmarg")}.

Alternatively, \code{errorCat} can be a list, where each list entry is a custom error function of the form \code{foo(true, pred)}, where \code{true} and \code{pred} are the arguments for the vectors of true and predicted values, respectively.  If predictions are made for a time-varying model and \code{tvMethod = "weighted"}, the weighted R2 or RMSE are computed. If a custom function is used, an additional argument for the weights has to be provided: \code{foo(true, pred, weights)}. Note that custom error functios can also be combined with the buildt-in functions, i.e. \code{errorCon = list("nCC", "CustomError"=foo)}.}

\item{tvMethod}{
Specifies how predictions and errors are computed for time-varying models: \code{tvMethod = "weighted"} computes errors by computing a weighted error over all cases in the time series at each estimation point specified in \code{estpoints} in \code{tvmgm()} or \code{tvmvar()}. The weighting corresponds to the weighting used for estimation (see \code{?tvmgm} or \code{?tvmvar}). \code{tvMethod = "closestModel"} determines for each time point the closest model and uses that model for prediction. See Details below for a more detailed explanation.
}

\item{consec}{Only relevant for (time-varying) mVAR models. An integer vector of length \code{nrow(data)}, indicating the sequence of measurement points in a time series. This is only relevant for mVAR models and time series with unequal time intervals. Defaults to \code{consec = NULL}, which assumes equal time intervals. \code{consec} is ignored if a \code{mgm} or \code{tvmgm} object is provided to \code{predict.mgm()}. For details see \code{?mvar}.}

\item{beepvar}{Together with the argument \code{dayvar}, this argument is an alternative to the \code{consec} argument (see above) to specify the consecutiveness of measurements. This is tailored to ecological momentary assessment (EMA) studies, where the consectutiveness is defined by the number of notification on a given day (\code{beepvar}) and the given day (\code{dayvar}).
}

\item{dayvar}{See \code{beepvar}.}

\item{...}{Additional arguments.}


}



\details{
Nodewise errors in time-varying models can be computed in two different ways: first, one computes the predicted value for each of the N cases in the time series for all models (estimated at different estimation points, see \code{?tvmgm} or \code{?tvmvar}). Then the error of each of the N cases for each of the models is weighted by the weight that has been used to estimate a given model at its estimation point. This means that the error of a data point close to the end of a time-series gets a high weight for models estimated in the end of the time-series and a low weight for models estimated in the beginning of the time series.

Second, we determine for each case in the time-series the closest estmation point, and use the model estimated at that estimation point to make predictions for that case.

Note that the error function normalized accuracy (nCC) is negative if the full model performs worse than the intercept model. This can happen if the model overfits the data.

}


\value{
A list with the following entries:

\item{call}{
Contains all provided input arguments.
}

\item{predicted}{
A n x p matrix with predicted values, matching the dimension of the true values in \code{true}.
}

\item{probabilities}{
A list with p entries corresponding to p nodes in the data. If a variable is categorical, the corresponding entry contains a n x k matrix with predicted probabilities, where k is the number of categories of the categorical variable. If a variable is continuous, the corresponding entry is empty.
}

\item{true}{
Contains the true values. For \code{mgm} and \code{tvmgm} objects these are equal to the data provided via \code{data}. For \code{mvar} and \code{tvmvar} objects, these are equal to the rows that can be predicted in a VAR model, depending on the largest specified lag and (if specified) the \code{consec} argument. 
}

\item{errors}{
A matrix containing the all types of errors specified via \code{errorCon} and \code{errorCat}, for each variable. If \code{tvMethod = "weighted"}, the matrix becomes an array, with an additional dimension for the estimation point.
}

\item{tverrors}{
If \code{tvMethod = "weighted"}, this list entry contains a list with errors of the format of \code{errors}, separately for each estimation point. The errors are computed from predictions of the model at the given estimation points and weighted by the weight-vector at that estimation point. If \code{tvMethod = "closestModel"}, this entry is empty.
}


}
\references{

Haslbeck, J., & Waldorp, L. J. (2018). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\examples{

\dontrun{
# See examples in ?mgm, ?tvmgm, ?mvar and ?tvmvar.
}

}

