# ‘BATSS’ R package

## Package description

BATSS - Bayesian Adaptive Trials Simulator Software - provides a
flexible structure for the fast simulation of Bayesian adaptive designs
for clinical trials.

This package can be used to define and evaluate the operating
characteristics of Bayesian adaptive designs for various different types
of primary outcomes (e.g., those that follow a normal, binary, Poisson
or negative binomial distribution or time-to-event outcomes) and can
incorporate the most common types of adaptations: stopping treatments
(or the entire trial) for efficacy or futility, and Bayesian response
adaptive randomisation - based on user-defined adaptation rules.

Other important features of this highly modular package include:
parallel processing, customisability, use on a cluster computer or
PC/Mac, and adjustment for covariates.

## Installation instructions

INLA needs to be installed prior to installing BATSS. For instructions,
visit <https://www.r-inla.org/download-install>

To install the stable version of the BATSS package from CRAN, run

``` r
install.packages("BATSS")
```

To install the development version of the BATSS package from GitHub, run

``` r
install.packages("devtools")
devtools::install_github("batss-dev/BATSS")
```

Note that a developement integer is added to the package version
sequence when BATSS is intalled via GitHub, so that the package version
breaks down into `major`.`minor`.`patch`.`dev`.
