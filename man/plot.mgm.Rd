\name{plot.mgm}
\alias{plot.mgm}

\title{
Visualization of bootstrapped edge weights
}
\description{
Visualizes the bootstrapped edge weights from boot.mgmfit() using boxplots
}
\usage{
plot.mgm(x, y = NULL, ...)
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
boot_obj <- boot.mgmfit(data = autism_data$data, 
                   type = autism_data$type,
                   lev = autism_data$lev,
                   B = 10)

# Plot bootstrapped edge weights using boxplots
plot.mgm(boot_obj)

# Summary of bootstrapped edge weights
summary(boot_obj)


}


}
