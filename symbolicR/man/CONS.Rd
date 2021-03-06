% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00symbolic.000general.R
\name{CONS}
\alias{CONS}
\title{CONS : evaluate argument other than first one}
\usage{
CONS(x1, ...)
}
\arguments{
\item{x1}{: the function, can be a variable, or quoted symbol}

\item{...}{: other arguments}
}
\value{
the call
}
\description{
similar as QCONS, however, will evaluated the argument first, as following
}
\details{
\code{CONS('*', 1, 2+3)} \samp{->} \code{1 * 5}
}
\seealso{
\code{\link{QCONS}}
}
\author{
jjxie
}
