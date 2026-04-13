# Bayesian adaptive trial simulations for time-to-event data

Simulation of Bayesian adaptive trials with time-to-event endpoint using
Integrated Nested Laplace Approximation (INLA).

## Usage

``` r
batss.surv(
  model,
  family = "weibullsurv",
  surv = simsurv::simsurv,
  surv.control,
  var,
  var.control = NULL,
  fup = NULL,
  cens = NULL,
  cens.control = NULL,
  accr = runif,
  accr.control = list(min = 0, max = 3),
  accr.type = "fixed",
  hr,
  which,
  R = 10000,
  N,
  n.max = NULL,
  alternative = "less",
  RAR = NULL,
  RAR.control = NULL,
  interim,
  t.max = NULL,
  event.max = NULL,
  prob0,
  delta.eff = 0,
  delta.fut = delta.eff,
  delta.RAR = 0,
  eff.arm,
  eff.trial = NULL,
  eff.arm.control = NULL,
  eff.trial.control = NULL,
  fut.arm,
  fut.trial = NULL,
  fut.arm.control = NULL,
  fut.trial.control = NULL,
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
  indicating a symbolic description of the model to be fitted

- family:

  A character string indicating the log-survival likelihood as described
  in the package INLA. Default set to '`weibullsurv`', other options
  include '`exponential.surv`', '`loglogistic.surv`','`lognormal.surv`',
  '`coxph`'. See
  [inla.list.models](https://rdrr.io/pkg/INLA/man/list-models.html).

- surv:

  Time-to-event data simulator function. Currently,
  '[simsurv](https://rdrr.io/pkg/simsurv/man/simsurv.html)' is the only
  option and used by default.

- surv.control:

  A list of parameters for the call to
  '[simsurv](https://rdrr.io/pkg/simsurv/man/simsurv.html)'.

- var:

  A list. Each entry corresponds to a variable described under '`model`'
  and indicates the name of a function allowing to generate variates
  (like [rnorm](https://rdrr.io/r/stats/Normal.html) and
  [rexp](https://rdrr.io/r/stats/Exponential.html), for example). The
  list names have to match the variable names used in '`model`'. The
  grouping variable corresponding to the target parameters has to be of
  class '[factor](https://rdrr.io/r/base/factor.html)' with levels
  corresponding to the names indicated in argument `prob0` (see below).

- var.control:

  An optional list of control parameters for the functions indicated in
  '`var`'. The names of the list items need to correspond to the names
  used in '`var`'. Each element is another list with names of the
  elements corresponding to the parameter names of the functions
  specified in '`var`'.

- fup:

  A scalar specifying the follow-up time of the last patient recruited
  and hence the final analysis, this only applies if 'maxt=NULL' in
  'surv.control'.

- cens:

  Random number generating function for censoring, e.g.
  '[rexp](https://rdrr.io/r/stats/Exponential.html)'

- cens.control:

  Additional parameters for call to censoring function.

- accr:

  Random number generating function for entry times, e.g.
  '[rexp](https://rdrr.io/r/stats/Exponential.html)'

- accr.control:

  Additional parameters for call to entry times generating function.

- accr.type:

  One of 'random' or 'fixed', the former will draw entry times from the
  distribution defined in '`accr`', adding the resulting times results
  in a random study duration. Setting this parameter to 'fixed' will
  allow setting a fixed study length and '`accr`' will be used to
  describe the distribution within the interval from '`0`' to study
  duration.

- hr:

  A numerical vector of parameter values for the hazard ratios. Its
  length has to match the number of treatments (excluding control).

- which:

  A numerical vector indicating the position of the target parameters.

- R:

  a vector of natural numbers to be used as seeds (check
  [set.seed](https://rdrr.io/r/base/Random.html)) for the different
  Monte Carlo trials (the vector length will thus correspond to the
  number of Monte Carlo trials). When `R` is a scalar, seeds `1` to `R`
  are used, where `R` corresponds to the number of Monte Carlo trials.

- N:

  A scalar indicating the total sample size.

- n.max:

  A vector indicating the maximum sample size per arm, if a scalar is
  given the same size will be used for all arms.

- alternative:

  A vector of strings providing the one-sided direction of the
  alternative hypothesis corresponding to each target parameter
  indicated under '`which`' (in the same order). Possibilities are
  'greater' or 'less'. If the vector is of length 1, the same direction
  will be used for all target parameter tests.

- RAR:

  A character string indicating the how the response-adaptive randomised
  probabilities corresponding to each group, reference group included,
  are defined. If RAR=NULL, equal probabilities of being attributed to
  each (active) group are used. Check ?BATS::prob.fun for examples.

- RAR.control:

  An optional list of control parameters the function in 'RAR'.

- interim:

  A list of parameters related to interim analyses. Possible list items
  include, '`recruited`' a vector of integers indicating the number of
  recruited participants at each look, last excluded, in increasing
  order, '`time`' a vector of integers indicating the time points of
  interim analyses, and '`event`' a vector of integers indicating the
  number of events. '`event.type`' is a character string describing the
  groups that are counted towards the numbers in '`event`' with options
  options '`all`', counting all events,'`control`', counting events in
  the control group only, '`cplusone`', counting the control plus
  treatments individually (this will result in different time points for
  the interims of the treatment groups), and '`min`', counting the last
  group to reach the threshold. See details for explanations and
  combinations.

- t.max:

  A scalar of the maximum trial duration, safeguard if 'events' takes
  too long to reach

- event.max:

  A scalar indicating the maximum number of event.max in the control
  group, once reached the trial will stop recruiting

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

- eff.arm.control:

  An optional list of parameters for the function indicated in
  '`eff.arm`'.

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

- fut.arm.control:

  An optional list of parameters for the function indicated in
  '`fut.arm`'.

- fut.trial.control:

  An optional list of parameters for the function indicated in
  '`fut.trial`'.

- H0:

  A logical indicating whether the simulation should also consider the
  case with all target parameters set to 0 to check the probability of
  rejecting the hypothesis that the target parameter value is equal to 0
  individually (pairwise type I error) or globally (family-wise error
  rate). Default set to `H0=TRUE`.

- computation:

  A character string indicating how the computation should be performed.
  Possibilities are '`parallel`' or '`sequential`' with default
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

The function batss.surv returns an S3 object of class 'batss' with
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

- type - The type of 'BATSS' analysis, i.e. '`surv`'

## See also

[summary.batss](summary.batss.md) and [plot.batss](plot.batss.md) for
detailed summaries and plots, and [batss.combine](batss.combine.md) to
combine different evaluations of batss.surv considering the same trial
design but different sets of seeds (useful for cluster computation).

## Examples

``` r
# Example 1: TBC
```
