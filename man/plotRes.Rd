\name{plotRes}
\alias{plotRes}

\title{
Plot summary of resampled sampling distributions
}
\description{
Plots a summary of sampling distributions resampled with the resample() function
}
\usage{
plotRes(object, quantiles = c(.05, .95), labels = NULL, 
        decreasing = TRUE, cut = NULL, cex.label = 0.75, 
        lwd.qtl = 2, cex.mean = 0.55, cex.bg = 2.7, 
        axis.ticks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), 
        axis.ticks.mod = NULL, layout.width.labels = .2,
        layout.gap.pw.mod = .15)
}

\arguments{
  \item{object}{
An output object from the \code{resample()} function.
}
  \item{quantiles}{
A numerical vector of length two, specifying the desired lower/upper quantiles. Defaults to \code{quantiles = c(.05, .95)}.
}
  \item{labels}{
A character vector of length p, containing the label of each variable, where p is the number of variables.
}
  \item{decreasing}{
If \code{TRUE} (default), the edges are ordered by the arithmetic mean of the sampling distribution in decreasing order. If \code{FALSE} they are ordered in increasing order.
}
  \item{cut}{
A sequence of integers, specifying which edges are represented. For instance, if \code{decreasing = TRUE} and \code{cut = 1:10}, summaries for the 10 edges with the largest parameter estimate are displayed. The cut argument can also be used to present the boostrapped CIs in several figures.
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
  \item{axis.ticks.mod}{
A numeric vector indicating the axis ticks and labels for the x-axis for moderation effects. If \code{axis.ticks.mod=NULL}, the values from \code{axis.ticks} for pairwise interactions are used.
}
\item{layout.width.labels}{
A positive numeric value which specifies the width of the left-hand-side legend relative to the width of the data panel (or data panels, in case of a moderator model), which have width = 1. Defaults to \code{layout.width.labels = 0.2}.}
\item{layout.gap.pw.mod}{
A positive numeric value which specifies the width of the gap between the stability of pairwise effects and moderation effects. Defaults to \code{layout.gap.pw.mod = 0.15}.
}
}

\details{
Currently only supports summaries for resampled \code{mgm()} objects, and moderated MGMs with a single moderator.
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
                    nB = 10) # should be more in real applications 

# Plot Summary
plotRes(object = res_aut,
        quantiles = c(.05, .95), 
        labels = NULL, 
        axis.ticks = c(-.25, 0, .25, .5, .75))

  }
}


