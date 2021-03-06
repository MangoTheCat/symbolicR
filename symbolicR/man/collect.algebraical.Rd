% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00symbolic.000general.R
\name{collect.algebraical}
\alias{collect.algebraical}
\title{collect.algebraical}
\usage{
collect.algebraical(le)
}
\arguments{
\item{le}{list of expressions}
}
\value{
list of combined expressions
}
\description{
For a list of expression, match them to coeff*indeterminators, and then combine together the indeterminators\cr
\code{ [ A, B*C, 5*A ] -> [6*A, B*C] } \cr
Here coefficients are atomic by \code{'?a'}
}
\author{
jjxie
}
