% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00symbolic.003eqns.R
\name{analyse.connected.component.critical.points}
\alias{analyse.connected.component.critical.points}
\title{analyse.connected.component.first.critical.points}
\usage{
analyse.connected.component.critical.points(eqns,
  socp.filter = socp.filter.default, ...)
}
\arguments{
\item{...}{extra.rules}
}
\value{
found pk relations and left equations(might be for other use)
}
\description{
for a connected component, try to extract PK variables (as lhs) one by one
If no pk variables found, return eqns as is
actually, what we can get are not PK variables, but however the critical point varibles
}
