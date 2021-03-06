% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00symbolic.006Matrix.R
\name{tuple.to.list}
\alias{tuple.to.list}
\title{convert TUPLE to list which means a symbolic Vector
\code{ TUPLE( x1,x2,x3) => list(x1,x2,x3) }
\code{
TUPLE( 
   TUPLE(x1,x2,x3),
   TUPLE(x4,x5,x6)) => list( list(x1,x2,x3), list(x4,x5,x6)) }}
\usage{
tuple.to.list(e0)
}
\arguments{
\item{e0}{quoted expression started by TUPLE}
}
\value{
list
}
\description{
convert TUPLE to list which means a symbolic Vector
\code{ TUPLE( x1,x2,x3) => list(x1,x2,x3) }
\code{
TUPLE( 
   TUPLE(x1,x2,x3),
   TUPLE(x4,x5,x6)) => list( list(x1,x2,x3), list(x4,x5,x6)) }
}
