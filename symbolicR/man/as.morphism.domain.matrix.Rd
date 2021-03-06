% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15domain.080morphism.domain.R
\name{as.morphism.domain.matrix}
\alias{as.morphism.domain.matrix}
\title{setting functions for the domain of a morphism
assume \code{m} has 3 columns low, est, upper corresponding to NONMEM statements
return a domain description}
\usage{
as.morphism.domain.matrix(m, parameter.space = 1:NROW(m))
}
\arguments{
\item{m}{matrix}

\item{parameter.space}{the row number to use, or indices for the THETAs}
}
\value{
expression for domain
}
\description{
setting functions for the domain of a morphism
assume \code{m} has 3 columns low, est, upper corresponding to NONMEM statements
return a domain description
}
