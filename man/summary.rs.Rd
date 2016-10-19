\name{summary.rs}
\alias{summary.rs}
\title{
Prints summary of an resampled mgm object}
\description{
Gives a summary of resampled MGMs / mixed VAR models.}
\usage{
\method{summary}{rs}(object, ZeroEdges = FALSE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The output object of \code{rs.mgm}}
  \item{ZeroEdges}{If \code{FALSE}, only edges which are at least once estimated as nonzero in B bootstrap samples are shown in the summary.}
  \item{...}{Not used.}
}

\value{
Gives a summary of edge estimates of models estiamted on B resampled datasets (see \code{?rs.mgm}).
}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>}

\examples{

\dontrun{

# Fit Mixed Graphical Model to Bootstrap Samples of the autism dataset
rs_obj <- rs.mgm(data = autism_data$data, 
                 type = autism_data$type,
                 lev = autism_data$lev,
                 VAR = FALSE, # fit mgm model
                 B = 10,
                 N = nrow(autism_data$data), # for classical bootstrap
                 replace=TRUE)

# Summary of bootstrapped edge weights
summary(rs_obj)

# Plot bootstrapped edge weights using boxplots
plot(rs_obj)

}

}