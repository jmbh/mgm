\name{print.mgm}
\alias{print.mgm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print method for mgm objects
}
\description{
Returns basic information about fit objects, prediction objects and bandwidth-selection objects.
}
\usage{
\method{print}{mgm}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The output object of \code{mgm()}, \code{mvar()}, \code{tvmgm()}, \code{tvmvar()}, \code{predict.mgm()} or \code{bwSelect()}.
}
  \item{\dots}{
Additional arguments.
}
}
\value{
Writes basic information about the object in the console.
}


\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}
