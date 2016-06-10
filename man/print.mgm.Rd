\name{print.mgm}
\alias{print.mgm}
\title{
Print method for mgm objects}

\description{
Returns the model class, number of nodes and (for time-varying models) the number of estimated time steps from a given mgm object.}

\usage{
\method{print}{mgm}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{mgm fit-object (output of functions \code{mgmfit}, \code{tv.mgmfit}, \code{var.mgm} or \code{tv_var.mgm}.)}
  \item{\dots}{additional print arguments}
}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>}

