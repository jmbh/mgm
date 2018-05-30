\name{plotRes}
\alias{plotRes}

\title{
Plot summary of resampled sampling distributions
}
\description{
Plots a summary of sampling distributions resampled with the resample() function
}
\usage{
plotRes(object, quantiles, labels = NULL, decreasing = TRUE, 
        cut = NULL, cex.label = 0.75, lwd.qtl = 2, 
        cex.mean = 0.5, cex.bg = 3.5, 
        axis.ticks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))
}

\arguments{
  \item{object}{
An output object from the \code{resample()} function.
}
  \item{quantiles}{
A numerical vector of length two specifying the desired lower/upper quantile, e.g. \code{c(.05, .95)}.
}
  \item{labels}{
A character vector of length p containing the label of each variable, where p is the number of variables.
}
  \item{decreasing}{
If \code{TRUE} (default), the edges are displayed by the mean of the sampling distribution in decreasing order. If \code{FALSE} they are ordered in increasing order.
}
  \item{cut}{
A sequence of integers, specifying which edges are represented. For instance, if \code{decreasing = TRUE} and \code{cut = 1:10}, summaries for the 10 edges with the largest parameter estimate are displayed.
}
  \item{cex.label}{
Text size of the labels.
}
  \item{lwd.qtl}{
Line width of line indicating the upper/lower quantiles.
}
  \item{cex.mean}{
Text size of the number indicating the proportion of the estimates whose absolute value is larger than zero.
}
  \item{cex.bg}{
Size of the white background of the number indicating the proportion of the estimates whose absolute value is larger than zero.
}
  \item{axis.ticks}{
A numeric vector indicating the axis ticks and labels for the x-axis.
}
}


\details{
Currently only supports summaries for resampled \code{mgm()} objects. 
}

\value{
Plots a figure that shows summaries of the resampled sampling distribution for (a set of) all edge parameters. These include the mean, a specified upper and lower quantile and the proportion of parameter estimates whose absolute value is larger than zero.
}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}


\seealso{
\code{resample()}, \code{mgm()}, \code{mvar()}, \code{tvmgm()}, \code{tvmar()}
}


\examples{

\dontrun{

# Fit initial model
fit_aut <- mgm(data = as.matrix(autism_data$data),
              type = autism_data$type,
              level = autism_data$lev,
              k = 2)


# Fit bootstrapped models
res_aut <- resample(object = fit_aut, 
                    data = as.matrix(autism_data$data), 
                    nB = 10)

# Plot Summary
plotRes(object = res_aut,
        qtl = c(.05, .95), 
        labels = NULL, 
        axis.ticks = c(-.25, 0, .25, .5, .75))

  }
}


