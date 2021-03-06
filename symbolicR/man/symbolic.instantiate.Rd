% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00symbolic.000general.R
\name{symbolic.instantiate}
\alias{symbolic.instantiate}
\title{symbolic.instantiate}
\usage{
symbolic.instantiate(skeleton, dict)
}
\arguments{
\item{skeleton}{expression skeleton to be instantiated}

\item{dict}{dictionary with definition of instance variables}
}
\value{
the instantiated expression
}
\description{
Given a dictionary, instantiate a new expression which might have instance variables
The \code{':'(x)} is a example of instance variable.
If \code{x} equals to \code{2}, then the following expression
\code{sin(':'(x))} will instantiate to \code{sin(2)}
}
\seealso{
\code{\link{symbolic.simplify.gigo}}
}
\author{
jjxie
}
