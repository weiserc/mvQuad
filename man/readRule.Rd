% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvQuad.R
\name{readRule}
\alias{readRule}
\title{reads a quadrature-rule from a text file}
\usage{
readRule(file = NULL)
}
\arguments{
\item{file}{file name of the text file containing the quadrature rule}
}
\value{
Returns an object of class 'customRule', which can be used for creating a 'NIGrid' (\code{\link{createNIGrid}})
}
\description{
\code{readRule} reads a quadrature-rule from a text file
}
\details{
The text file containing the quadrature rule has to be formatted in the following way:

The first line have to declare the domain \code{initial.domain a b}, where a and b denotes the lower and upper-bound for the integration domain.
This can be either a number or '-Inf'/'Inf' (for example \code{initial.domain 0 1} or \code{initial.domain 0 Inf})

Every following line contains one single node and weight belonging to one level of the rule (format: 'level' 'node' 'weight').
This example shows the use for the "midpoint-rule" (levels: 1 - 3).

> \code{initial.domain 0 1}

> \code{1 0.5 1}

> \code{2 0.25 0.5}

> \code{2 0.75 0.5}

> \code{3 0.166666666666667 0.333333333333333}

> \code{3 0.5 0.333333333333333}

> \code{3 0.833333333333333 0.333333333333333}
}
\examples{
\dontrun{myRule <- readRule(file="midpoint_rule.txt")}
\dontrun{nw <- createNIGrid(d=1, type = myRule.txt, level = 2)}
}
\seealso{
\code{\link{createNIGrid}}
}

