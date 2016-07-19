   \name{QuadRules}
   \alias{QuadRules}
     \docType{data}
     \title{nodes and weights for 1D - Gauss-Quadrature}
     \description{
       This data set stores nodes an weights for Gauss-Quadrature.
       Syntax:

          QuadRules[['type']][['level']]


          \itemize{
            \item type="\bold{GLe}" Gauss-Legendre; interval [0,1]; max-level 45
            \item type="\bold{nLe}" nested-type Gauss-Legendre; interval [0,1]; max-level 25
            \item type="\bold{GKr}" Gauss-Kronrod; interval [0,1]; max-level 29
            \item type="\bold{GLa}" Gauss-Laguere; interval [0, Inf); max-level 30
            \item type="\bold{GHe}" Gauss-Hermite; interval (-Inf, Inf); max-level 45
            \item type="\bold{GHN}" Gauss-Hermite (as above, but pre-multiplied weights \eqn{\hat(w)_i = w_i * \phi(x_i)})
            \item type="\bold{nHe}" nested-type Gauss-Hermite; interval (-Inf, Inf) max-level 25
            \item type="\bold{nHN}" nested-type Gauss-Hermite (as above, but pre-multiplied weights \eqn{\hat(w)_i = w_i * \phi(x_i)})
            \item type="\bold{Leja}" Leja-points; interval [0,1]; max-level 141
          }
     }
     \examples{nw <- QuadRules[["GHe"]][[2]]}
     \format{list of nodes and weights (for organisation see "Syntax" in description section)}
     \source{ - http://keisan.casio.com/exec/system/1329114617 high precission computing (for G..-rules)

     - further information in \code{\link{createNIGrid}}}
     \keyword{datasets}
