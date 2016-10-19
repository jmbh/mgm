\name{summary.mgm}
\alias{summary.mgm}
\title{
Prints summary of an mgm object}
\description{
Print a selection of the output-list of the function \code{mgmfit}}
\usage{
\method{summary}{mgm}(object, data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The output object of \code{mgmfit}, \code{tv.mgmfit}, \code{var.mgm} or \code{tv_var.mgm}.}
  \item{data}{If a dataset matching the dataset used for estimation is provided, the summary function adds a column with node wise errors to the output.}
  \item{...}{Not sed.}
}

\value{
Outputs a dataframe (stationary models) or list of dataframes (time-varying models) with the columns Variable, Type, Degree, Lambda, Threshold, EBIC and if specified, nodewise errors}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>}

