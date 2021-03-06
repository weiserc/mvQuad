% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvQuad.R
\name{size (size.NIGrid)}
\alias{dim.NIGrid}
\alias{size}
\alias{size (size.NIGrid)}
\alias{size.NIGrid}
\title{returns the size of an NIGrid-object}
\usage{
size(object, ...)

\method{size}{NIGrid}(object, ...)

\method{dim}{NIGrid}(x)
}
\arguments{
\item{object}{a grid of type \code{NIGrid}}

\item{...}{other arguments passed to the specific method}

\item{x}{object of type \code{NIGrid}}
}
\value{
Returns the grid size in terms of dimensions, number of grid points and used memory
}
\description{
Returns the size of an NIGrid-object
}
\examples{
myGrid <- createNIGrid(dim=2, type="GHe", level=5)
size(myGrid)
dim(myGrid)
}

