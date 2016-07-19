\name{mvQuad-package}
\alias{mvQuad-package}
\alias{mvQuad}
\docType{package}
\title{
Methods for multivariate Quadrature (numerical integration)
}
\description{
This package provides methods to construct multivariate grids, which can be
used for multivariate quadrature. This grids can be based on different quadrature
rules like Newton-Cotes formulas (trapezoidal-, Simpson-rule, ...) or Gauss-Quadrature
(Gauss-Hermite, Gauss-Legendre, ...). For the construction of the multidimensional
grid the product-rule or the combination-technique can be applied.
}
\details{
\tabular{ll}{
Package: \tab mvQuad\cr
Type: \tab Package\cr
Version: \tab 1.0-6\cr
Date: \tab 2016-07-19\cr
License: \tab GPL-3\cr
}
}
\author{
Constantin Weiser (HHU of Duesseldorf / Germany)
Maintainer: Constantin Weiser <weiserc@hhu.de>
}
\references{
Philip J. Davis, Philip Rabinowitz (1984): Methods of Numerical Integration

F. Heiss, V. Winschel (2008): Likelihood approximation by numerical integration on sparse grids, Journal of Econometrics

H.-J. Bungartz, M. Griebel (2004): Sparse grids, Acta Numerica

Peter Jaeckel (2005): A note on multivariate Gauss-Hermite quadrature
}
\keyword{ package }
\keyword{ math }
\examples{
myGrid <- createNIGrid(dim=2, type="GLe", level=5)
rescale(myGrid, domain=rbind(c(-1,1),c(-1,1)))
print(myGrid)
plot(myGrid, col="blue")
myFun <- function(x){
   1 - x[,1]^2 * x[,2]^2
}
quadrature(myFun, myGrid)
}
