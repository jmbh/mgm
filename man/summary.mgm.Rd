\name{summary.mgm}
\alias{summary.mgm}
\title{
Prints summary of an mgm object}
\description{
Print a selection of the output-list of the function \code{mgmfit}}
\usage{
\method{summary}{mgm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The output object of \code{mgmfit}, \code{tv.mgmfit}, \code{var.mgm} or \code{tv_var.mgm}.}
  \item{...}{Other summary arguments.}
}
\value{
Outputs a dataframe (stationary models) or list of dataframes (time-varying models) with the columns Variable, Type, Lambda, Threshold, EBIC.}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>}

