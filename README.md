# mvQuad: Methods for Multivariate Quadrature

link to CRAN: https://cran.r-project.org/web/packages/mvQuad/

This package provides a collection of methods for (potentially) multivariate
quadrature in R. It's especially designed for use in _statistical problems_, 
characterized by integrals of the form $\int_a^b g(x)p(x) \ dx$, where $p(x)$
denotes a density-function. Furthermore the methods are also applicable to 
standard integration problems with finite, semi-finite and infinite intervals.

In general quadrature (also named: numerical integration, numerical
quadrature) work this way: The integral of interests is approximated by a 
weighted sum of function values.

$$ I = \int_a^b h(x) \ dx \approx A = \sum_{i=1}^n w_i \cdot h(x_i) $$

The so called nodes ($x_i$) and weights($w_i$) are defined by the chosen 
quadrature rule, which should be appropriate (better: optimal) for the 
integration problem in hand^[A rigorous introduction to numerical integration
can be found in Davis/Rabinowitz [-@DavisRabinowitz1984]].

This principle can be transferred directly to the multivariate case.

The methods provided in this package cover the following tasks:  

* creating an uni-/multivariate grid (grid: collection of nodes and weights) for a chosen quadrature rule     
* examining the created grid     
* rescaling the grid for appropriate/efficient use    
* computing the approximated integral   


## Quick Start
This section shows a typical workflow for quadrature with the `mvQuad`-package.
More details and additional features of this package are provided in the subsequent sections.
In this illustration the following two-dimensional integral will be approximated:
$$ I = \int_1^2 \int_1^2 x \cdot e^y \ dx dy $$


```{r }
  library(mvQuad)

  # create grid
  nw <- createNIGrid(dim=2, type="GLe", level=6)
  
  # rescale grid for desired domain
  rescale(nw, domain = matrix(c(1, 1, 2, 2), ncol=2))

  # define the integrand
  myFun2d <- function(x){
    x[,1]*exp(x[,2])
  }

  # compute the approximated value of the integral
  A <- quadrature(myFun2d, grid = nw)
```

**Explanation Step-by-Step**    

0. `mvQuad`-package is loaded    
1. with the `createNIGrid`-command a two-dimensional (`dim=2`) grid, based on Gauss-Legendre quadrature rule (`type="GLe"`) with a given accuracy level (`level=6`) is created and stored in
the variable `nw`    

The grid created above is designed for the domain $[0, 1]^2$ but we need one
for the domain $[1, 2]^2$     

2. the command `rescale` changes the domain-feature of the grid; the new domain
is passed in a matrix (`domain=...`)

3. the integrand is defined

4. the `quadrature`-command computes the weighted sum of function values as mentioned
in the introduction

## Supported Rules (build in)
The choice of quadrature rule is heavily related to the integration problem. Especially
the domain/support ($[a, b]$ (finite), $[a, \infty)$ (semi-finite), $(-\infty, \infty)$ (infinite)) determines the choice.

The `mvQuad`-packages provides the following quadrature rules.      

* __`cNC1, cNC2, ..., cNC6`__: closed Newton-Cotes Formulas of degree 1-6 (1=trapezoidal-rule; 2=Simpson's-rule; ...), finite domain

* __`oNC0, oNC1, ..., oNC3`__: open Newton-Cote Formula of degree 0-3 (0=midpoint-rule; ...),
 finite domain
* __`GLe, GKr`__:  Gauss-Legendre and Gauss-Kronrod rule, finite domain
* __`nLe`__: nested Gauss-Legendre rule (finite domain) [@Petras2003]
* __`Leja`__: Leja-Points (finite domain)
* __`GLa`__: Gauss-Laguerre rule (semi-finite domain)
* __`GHe`__: Gauss-Hermite rule (infinite domain)
* __`nHe`__: nested Gauss-Hermite rule (infinite domain)  [@GenzKeister1996]
* __`GHN`, `nHN`__: (nested) Gauss-Hermite rule as before but weights are pre-multiplied by the standard normal density ($\hat{w}_i = w_i * \phi(x_i)$).^[Those rules are computationally more efficent for integrands of the form: $\int_{-\infty}^{\infty} g(x)\phi(x)dx$ because the approximation reduces to $\sum \hat{w}_i \cdot g(x)$.]

For each rule grids can be created of different accuracy. The adjusting screw in
the `createNIGrid` is the `level`-option. In general, the higher `level` the more precise the approximation. For the Newton-Cotes rules an arbitrary level can be chosen. The other rules uses lookup-tables for the nodes and weights and are therefore restricted to a maximum level (see `help(QuadRules)`)
