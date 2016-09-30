\name{rs.mgmfit}
\alias{rs.mgmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stability analysis of Mixed Graphical Models using resampling schemes
}
\description{
Estimates Mixed Graphical Models on a number of resampled datasets. Can be used as input to \code{summary()} and \code{plot()} to investigate parameter-stability.
}
\usage{
rs.mgmfit(..., B = 100, N = NULL, replace = TRUE, pbar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Arguments as in \code{mgmfit}.
}
  \item{B}{
Number of resampling iterations.
}
  \item{N}{
Size of each resample.
}
  \item{replace}{
Subsampling with (bootstrap) or without replacement.
}
  \item{pbar}{
Shows a progress bar if \code{TRUE}}
}
\details{
Returns B mixed graphical models, each fitted on one bootstrap sample. Use the \code{summary()} and \code{plot()} function to get information about the stability of (a selection of) edges. If \code{N = nrow(data)} and \code{replace = TRUE} the resampling scheme is the nonparametric bootstrap.
}
\value{
Returns a list with
  \item{models}{A list with B entries containing the mgmfit() output for each of the B bootstrap samples.}
  \item{B_wadj}{A p x p x B array containing the estimated weighted adjacency matrix for each of the B bootstrap samples.}
  \item{B_samples}{A list with B entries containing n x 1 indicator vectors, defining the used bootstrap samples.}
  \item{call}{Returns the function call.}
}

\author{
Jonas Haslbeck
}

\examples{

\dontrun{

# Fit Mixed Graphical Model to Bootstrap Samples of the autism dataset
rs_obj <- rs.mgmfit(data = autism_data$data, 
                      type = autism_data$type,
                      lev = autism_data$lev,
                      B = 10,
                      N = nrow(autism_data$data), # for classical bootstrap
                      replace=TRUE)

# Summary of bootstrapped edge weights
summary(rs_obj)

# Plot bootstrapped edge weights using boxplots
plot(rs_obj)

# Pairwise scatter plots between bootstrapped edge weights to detect trade-offs
# ... in the case of correlated predictors (negative correlation between weights)

ind <- unlist(rs_obj$edgeNonZero) # only consider edges with nonzero estimates
m_bweights <- do.call(cbind, rs_obj$edgeWeights)[, ind] # B x k matrix, with k edges
colnames(m_bweights) <- unlist(rs_obj$edgeNames)[ind] # add edge names

library(psych) # psych package for pairs
pairs.panels(m_bweights)
# (recommended to print a PDF, because this can be a huge plot:
# pdf('filename.pdf', 20, 20))
# ... plotting
# dev.off()



}


}
