% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05graphic.000common.R
\name{topo.sort}
\alias{topo.sort}
\title{topo.sort
Topological sort according to dependency}
\usage{
topo.sort(G, rev = F)
}
\arguments{
\item{G}{graph}

\item{rev}{logical \code{TRUE} then reverse direction of all edges}
}
\value{
a list has two names \code{L} and \code{stage.number} \cr 
     where \code{L} is the sorted vertex and \code{N} are the corresponding stage number for the vertex in \code{L}
}
\description{
Assume \code{G} is directed graph as following: \cr
 \code{G=list(V=V,E=E)}
Do a topological sort of vertex of \code{G}
}
\author{
jjxie
}
