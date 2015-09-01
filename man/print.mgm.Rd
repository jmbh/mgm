\name{print.mgm}
\alias{print.mgm}
\title{
print a mgm object}
\description{
Print a selection of the output-list of the function \code{mgmfit}}
\usage{
\method{print}{mgm}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{fitted mgm object (output of function \code{mgmfit}.)}
  \item{\dots}{additional print arguments}
}
\value{
The adjacency matrix and the weighted adjacency matrix is printed}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>}

\examples{
head(data_mixed2) #example data

type <- c("g", "p", "e", "c", "c") #Gaussian, Poission, ...
levs <- c(1,1,1,3,2) #Number of categories

fit <- mgmfit(data = data_mixed2, 
              type = type, 
              lev = levs, 
              lambda.sel = "CV", 
              folds = 10, 
              gam = .25, 
              d = 2, 
              rule.reg = "AND", 
              rule.cat = "OR")

print(fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
