\name{showInteraction}
\alias{showInteraction}

\title{
Retrieving details of interactions
}
\description{
Retrieves details of a specified interaction from mgm model objects.
}
\usage{
showInteraction(object, int)
}


\arguments{
  \item{object}{
The output of one of the estimation functions \code{mgm()}, \code{tvmgm()}, \code{mvar()}, \code{tvmvar()}.
}
  \item{int}{
An integer vector specifying the interaction. For mVAR models, this vector has length 2. For MGMs the vector can be larger to request details of interaction of order > 2.
}
}

\details{
Currently the function only returns details of pairwise interactions from output objects of \code{mgm()}.
}
\value{
  \item{variables}{Integer vector returning the variables specified via the argument \code{int}}
  \item{type}{Character vector returning the type of the specified variables variables}
  \item{level}{Integer vector returning the number of levels of the specified variables variables}
  \item{parameters}{A list of length equal to the order k of the specified interaction. The entries contain the set of parameters obtained from the nodewise regressions on the k variables. Depending on the type of the variables in the interaction, these sets can obtain one or several parameters. For details see \code{?mgm} or Haslbeck & Waldorp (2017).}
  


}
\references{
Haslbeck, J., & Waldorp, L. J. (2017). mgm: Estimating time-varying Mixed Graphical Models in high-dimensional Data. arXiv preprint arXiv:1510.06871.
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\seealso{
\code{mgm}, \code{tvmgm}, \code{mvar}, \code{tvmvar}
}


\examples{

\dontrun{


## We fit a pairwise and 3-order MGM to the mixed Autism dataset (?autism_data)

# 1) Fit Pairwise MGM

# Call mgm()
fit_d2 <- mgm(data = autism_data$data,
              type = autism_data$type,
              level = autism_data$lev,
              k = 2) # ad most pairwise interacitons

# Weighted adjacency matrix
fit_d2$pairwise$wadj # for instance, we see there is an interaction 1-2

# 2) Look at details of interaction 1-2
showInteraction(object = fit_d2, 
                int = c(1, 2))

# For more examples see https://github.com/jmbh/mgmDocumentation

}


}

