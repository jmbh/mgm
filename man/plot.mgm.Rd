\name{plot.mgm}
\alias{plot.mgm}

\title{
Visualization of bootstrapped edge weights
}
\description{
Visualizes the bootstrapped edge weights from boot.mgmfit() using boxplots
}
\usage{
\method{plot}{mgm}(x, y = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Output object from the \code{boot.mgmfit()} function.
}
  \item{y}{
Not used
}
  \item{\dots}{
Not used
}
}
\details{
Returns a boxplot of each edge-weight that was at least once nonzero in the B bootstrap samples.
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

# Plot bootstrapped edge weights using boxplots
plot(rs_obj)

# Summary of bootstrapped edge weights
summary(rs_obj)


}


}
