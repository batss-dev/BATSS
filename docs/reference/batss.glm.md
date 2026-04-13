# Bayesian adaptive trial simulations for generalised linear models

Simulation of Bayesian adaptive trials with GLM endpoint using
Integrated Nested Laplace Approximation (INLA).

## Usage

``` r
batss.glm(
  model,
  var,
  var.control = NULL,
  family = "gaussian",
  link = "identity",
  beta,
  which,
  alternative = "greater",
  R = 10000,
  N,
  interim,
  prob0,
  delta.eff = 0,
  delta.fut = delta.eff,
  delta.RAR = 0,
  eff.arm,
  eff.arm.control = NULL,
  eff.trial = NULL,
  eff.trial.control = NULL,
  fut.arm,
  fut.arm.control = NULL,
  fut.trial = NULL,
  fut.trial.control = NULL,
  RAR = NULL,
  RAR.control = NULL,
  H0 = TRUE,
  computation = "parallel",
  mc.cores = getOption("mc.cores", 3L),
  extended = 0,
  ...
)
```

## Arguments

- model:

  an object of class '[formula](https://rdrr.io/r/stats/formula.html)'
  indicating a symbolic description of the model to be fitted (as in the
  [lm](https://rdrr.io/r/stats/lm.html) and
  [glm](https://rdrr.io/r/stats/glm.html) functions).

- var:

  A list. Each entry corresponds to a variable described under '`model`'
  and indicates the name of a function allowing to generate variates
  (like [rnorm](https://rdrr.io/r/stats/Normal.html) and
  [rexp](https://rdrr.io/r/stats/Exponential.html), for example). The
  list names have to match the variable names unded in '`model`' and its
  first element should correspond to the model outcome. The grouping
  variable corresponding to the target parameters has to be of class
  '[factor](https://rdrr.io/r/base/factor.html)' with levels
  corresponding to the names indicated in argument `prob0` (see below).

- var.control:

  An optional list of control parameters for the functions indicated in
  '`var`'. The names of the list items need to correspond to the names
  used in '`var`'. Each element is another list with names of the
  elements corresponding to the parameter names of the functions
  specified in '`var`'.

- family:

  A character string indicating the name of the conditional distribution
  as described in the package INLA (check
  [inla.list.models](https://rdrr.io/pkg/INLA/man/list-models.html)).
  Default set to '`gaussian`'.

- link:

  A character string describing the link function to be used in the
  model to relate the outcome to the set of predictors: 'identity',
  'log', 'logit', 'probit', 'robit', 'cauchit', 'loglog' and 'cloglog'
  are the currently available options. Default set to 'identity'.

- beta:

  A numerical vector of parameter values for the linear predictor. Its
  length has to match the number of column of the **X** matrix induced
  by the formula indicated under '`model`' (check
  [model.matrix](https://rdrr.io/r/stats/model.matrix.html)).

- which:

  A numerical vector indicating the position of the target `beta`
  parameters.

- alternative:

  A vector of strings providing the one-sided direction of the
  alternative hypothesis corresponding to each target parameter
  indicated under '`which`' (in the same order). Possibilities are
  'greater' (default) or 'less'. If the vector is of length 1, the same
  direction will be used for all target parameter tests.

- R:

  a vector of natural numbers to be used as seeds (check
  [set.seed](https://rdrr.io/r/base/Random.html)) for the different
  Monte Carlo trials (the vector length will thus correspond to the
  number of Monte Carlo trials). When `R` is a scalar, seeds `1` to `R`
  are used, where `R` corresponds to the number of Monte Carlo trials.

- N:

  A scalar indicating the maximum sample size.

- interim:

  A list of parameters related to interim analyses. Currently, only
  '`recruited`' is available. It consists in a vector of integers
  indicating the number of completed observations at each look, last
  excluded, in increasing order.

- prob0:

  A named vector with initial allocation probabilities. Names need to
  correspond to the levels of the grouping variable. If `RAR = NULL`,
  these probabilities/ratios will be used throughout (fixed allocation
  probabilities).

- delta.eff:

  A vector (of length equal to the number of looks (i.e., number of
  interims + 1)) of clinically meaningful treatment effect values (on
  the linear predictor scale) to be used to define the efficacy-related
  posterior probabilities for each target parameter at each look. If a
  scalar is provided, the same value is used at each look. The default
  is `delta.eff = 0`.

- delta.fut:

  A vector (of length equal to the number of looks (i.e., number of
  interims + 1)) of clinically meaningful treatment effect values (on
  the linear predictor scale) to be used to define the futility-related
  posterior probabilities for each target parameter at each look. If a
  scalar is provided, the same value is used at each look. The default
  is `delta.fut = delta.eff`.

- delta.RAR:

  A vector (of length equal to the number of looks (i.e., number of
  interims + 1)) of clinically meaningful treatment effect values (on
  the linear predictor scale) to be used to define the RAR-related
  posterior probabilities for each target parameter at each look. If a
  scalar is provided, the same value is used at each interim analysis.
  The default is `delta.RAR = 0`. Note that, when a vector is provided,
  its last value is ignored as no randomisation is made at the last
  look.

- eff.arm:

  A function defining if efficacy has been achieved at a given look
  given the information available at that stage a given target
  parameter. The output of this function must be a
  [logical](https://rdrr.io/r/base/logical.html) (of length 1).
  Arguments of this function will typically consider 'BATSS'
  ingredients. Check [eff.arm.simple](eff.arm.simple.md) and
  [eff.arm.infofract](eff.arm.infofract.md) for examples.

- eff.arm.control:

  An optional list of parameters for the function indicated in
  '`eff.arm`'.

- eff.trial:

  A function defining if the trial can be stopped for efficacy given the
  output of the function indicated in '`eff.arm`'. The output of this
  function must be a [logical](https://rdrr.io/r/base/logical.html) of
  length one. Arguments of this function will typically only consider
  the 'BATSS' ingredient `eff.target`. Check
  [eff.trial.all](eff.trial.all.md) and
  [eff.trial.any](eff.trial.any.md) for examples. When
  `eff.trial = NULL` (default), the trial stops for efficacy when *all*
  target parameters are found to be effective (like in
  [eff.trial.all](eff.trial.all.md)).

- eff.trial.control:

  An optional list of parameters for the function indicated in
  '`eff.trial`'.

- fut.arm:

  A function defining if futility has been achieved at a given look
  given the information available at that stage for each target
  parameter. The output of this function must be a
  [logical](https://rdrr.io/r/base/logical.html) (of length 1).
  Arguments of this function will typically consider 'BATSS'
  ingredients. Check [fut.arm.simple](fut.arm.simple.md) to see an
  example of such a function.

- fut.arm.control:

  An optional list of parameters for the function indicated in
  '`fut.arm`'.

- fut.trial:

  A function defining if the trial can be stopped for futility given the
  output of the function indicated in '`fut.arm`'. The output of this
  function must be a [logical](https://rdrr.io/r/base/logical.html) of
  length one. Arguments of this function will typically only consider
  the 'BATSS' ingredient `fut.target`. Check
  [fut.trial.all](fut.trial.all.md) for an example of such a function.
  When `fut.trial = NULL` (default), the trial stops for futility when
  *all* target parameters are found to be futile (like in
  [fut.trial.all](fut.trial.all.md)).

- fut.trial.control:

  An optional list of parameters for the function indicated in
  '`fut.trial`'.

- RAR:

  A function defining the response-adaptive randomisation probabilities
  of each group - reference group included - with the same group names
  and ordering as used in '`prob0`'. Arguments of this function will
  typically consider 'BATSS' ingredients. Check
  [RAR.trippa](RAR.trippa.md) and [RAR.optimal](RAR.optimal.md) for
  examples. If `RAR = NULL` (default), the probabilities/ratios
  indicated under `prob0` will be used throughout (fixed allocation
  probabilities).

- RAR.control:

  An optional list of control parameters for the function provided in
  '`RAR`'.

- H0:

  A logical indicating whether the simulation should also consider the
  case with all target parameters set to 0 to check the probability of
  rejecting the hypothesis that the target parameter value is equal to 0
  individually (pairwise type I error) or globally (family-wise error
  rate). Default set to `H0=TRUE`.

- computation:

  A character string indicating how the computation should be performed.
  Possibilities are 'parallel' or 'sequential' with default
  `computation="parallel"` meaning that the computation is split between
  `mc.cores`.

- mc.cores:

  An integer indicating the number of CPUs to be used when
  `computation="parallel"` (Default to 3 if no global '`mc.cores`'
  global option is available via
  [getOption](https://rdrr.io/r/base/options.html)).

- extended:

  an integer indicating the type of results to be returned. 0 (default)
  provides summary statistics, 1 adds the results of each Monte Carlo
  trial and 2 additionally returns each Monte Carlo dataset.
  [batss.combine](batss.combine.md) requires extended \> 0 as the
  function needs to merge results of different sets of seeds.

- ...:

  Additional arguments to control fitting in
  [inla](https://rdrr.io/pkg/INLA/man/inla.html).

## Value

The function batss.glm returns an S3 object of class 'batss' with
available print/summary/plot functions

- beta - A data frame providing information related to the beta
  parameter vector, like parameter names and values, for example.

- look - A data frame providing information related to looks, like
  sample size of a given interim (m) and cumulative sample size at a
  given interim (n), for example.

- par - A list providing different information, like the used seeds
  (seed) and the groups (group), for example.

- H1 - A list providing trial results under the alternative, like the
  estimates per target parameter when the corresponding arm was stopped
  (estimate), the efficacy and futility probabilites per target
  parameter and overall (target, efficacy and futility), the sample size
  per group and trial (sample), the probabilities associated to each
  combination of efficacy and futility per group (scenario), the
  detailed results per trial (trial), for example.

- H0 - A list providing trial results under the global null hypothesis
  (same structure as H1).

- call - The matched call.

- type - The type of 'BATSS' analysis (only 'glm' is currently
  available).

## See also

[summary.batss](summary.batss.md) and [plot.batss](plot.batss.md) for
detailed summaries and plots, and [batss.combine](batss.combine.md) to
combine different evaluations of batss.glm considering the same trial
design but different sets of seeds (useful for cluster computation).

## Examples

``` r
# \donttest{
# Example: 
# * Gaussian conditional distribution with sigma = 5
# * 3 groups with group means 'C' = 1 (ref), 'T1' = 2, 'T2' = 3,
#     where higher means correspond to better outcomes 
# * 5 interim analyses occurring when n = 100, 120, 140, 160, and 180
# * fixed and equal allocation probabilities per arm (i.e., no RAR)
# * max sample size = 200 
# * efficacy stop per arm when the prob of the corresponding parameter 
#     being greater than 0 is greater than 0.975 (?eff.arm.simple)
# * futility stop per arm when the prob of the corresponding parameter 
#     being greater than 0 is smaller than 0.05 (?fut.arm.simple) 
# * trial stop once all arms have stopped (?eff.trial.all and ?fut.trial.all)
#     or the max sample size was reached 

sim = batss.glm(model            = y ~ group,   
                var              = list(y     = rnorm,
                                        group = alloc.balanced),
                var.control      = list(y = list(sd = 5)),
                beta             = c(1, 1, 2),
                which            = c(2:3),
                alternative      = "greater",
                R                = 20,
                N                = 200,
                interim          = list(recruited = seq(100, 180, 20)),
                prob0            = c(C = 1/3, T1 = 1/3, T2 = 1/3),
                eff.arm          = eff.arm.simple,
                eff.arm.control  = list(b = 0.975),
                fut.arm          = fut.arm.simple,
                fut.arm.control  = list(b = 0.05),
                computation      = "parallel",
                H0               = TRUE,
                mc.cores         = 2)# better: parallel::detectCores()-1
#>     Initialisation
#>     Evaluation of H1
#>     Evaluation of H0
#>     Results
# }
```
