\name{glmnet-internal}
\title{Internal glmnet functions}
\alias{calc.EBIC, node.est, potts_parameter, VARreshape}

\description{Internal mgm functions.}

\author{Jonas Haslbeck}
\details{These are not intended for use by users. 

\code{calc.EBIC}
  calculates the EBIC from the sequence of lambda values in \code{glmnet} objects.
  
\code{node.est}
  Subroutine in \code{mgmfit} for neighborhood regression at each node.
  
\code{potts_parameter}
  Translate provided edge weights in \code{mgmsampler} into a paramter matrix. The details are in the help file for \code{mgmsampler}.

\code{VARreshape}
  Reshapes data in a lag-1 VAR format.
  
  }
\keyword{internal}
