# Negative binomial endpoint

## Context

This case study is inspired by the *CHANGE-MS* study (Curtin et al.
2016), which was a multi-arm, phase II randomised controlled trial which
aimed to study the efficacy and safety of *temelimab* for the treatment
of relapsing-remitting multiple sclerosis (MS) and determine whether
*temelimab* slowed down or stopped the progression of MS.

The primary endpoint is the cumulative number of new active lesions that
were identified on brain MRI scans that were performed on 4 occasions
(once a month) from weeks 12-24 post-randomisation. The data were
overdispersed (Curtin et al. 2016; Hartung et al. 2022) and modelled
using a negative binomial GLM; we will use the same model here.

## Design

The chosen design has the following characteristics:

- **treatment arms**: Our case study has 4 treatment arms: placebo
  (control), low, medium and high doses of *temelimab* (intervention
  arms referred to as arms “A”, “B” and “C”, respectively).
- **alternative hypotheses**: This is a superiority trial where we are
  interested in demonstrating efficacy compared to the placebo arm and
  intervention arms will not be compared to one another.
- **interim analyses and maximum sample size**: We will begin the
  interim analyses at 100 participants completing follow-up, and perform
  interim analyses every 40 participants completing follow-up
  thereafter. Our maximum sample size is 260 participants, which is the
  same as the original study (Curtin et al. 2016).
- **power**: Similar to the original study, we aim to have 90% power to
  detect a 60% reduction (RR=0.4) in the mean cumulative number of
  lesions in the highest dose intervention arm (from a mean of 4
  (control) to 1.6 lesions).
- **endpoint conditional distribution**: the endpoint is assumed to
  follow a negative binomial distribution with nuisance parameter equal
  to $\phi = 0.5$ when considering the parametrisation leading to the
  following variance: $$\text{Var}(Y) = \mu + \frac{\mu^{2}}{\phi}.$$
- **group allocation**: we will consider equal allocation probabilities
  per group (no response-adaptive randomisation).
- **efficacy stopping rule**: Early stopping of intervention arms for
  efficacy may occur if there is a high posterior probability of (any)
  benefit at look j, i.e., when
  $$\text{Prob}\left( \beta_{k} < \delta_{\epsilon}|y,X \right) > 1 - b_{\epsilon}\left( \frac{\sum\limits_{k = 1}^{K}n_{k}}{N} \right)^{p_{\epsilon}}$$
  where
  - $\beta_{k}$ denotes the $k$th target parameter,
  - $\delta_{\epsilon} = 0$ denotes the (efficacy-related) clinically
    meaningful treatment effect value,
  - $N$ is the maximum sample size,
  - $n_{k}$ is the number of participants that have completed follow up
    in arm $k = 1,...,K$ at look $j$ so that
    $\left( \sum_{k = 1}^{K}n_{k} \right)/N$ corresponds to the
    information fraction,
  - $b_{\epsilon} = 0.009$ and $p_{\epsilon} = 3$ are tuning parameters
    that determine the shape of the efficacy function (Gotmaker et al.
    2019).
- **futility stopping rule**: Early stopping of intervention arms for
  futility may occur if there is a low posterior probability of seeing a
  minimum clinically important difference of a 20% reduction (i.e.,
  RR=0.8), i.e., when
  $$\text{Prob}\left( \beta_{k} < \delta_{f}|y,X \right) < b_{f}$$ where
  - $\delta_{f} = log(0.8)$ denotes the (futility-related) clinically
    meaningful treatment effect values,
  - $b_{f} = 0.2025$ is a the cut-off value used to declare futility.
- **trial stopping rule:** The trial will run until an efficacy or
  futility decision has been reached for each intervention arm or once
  the maximum sample size has been reached.

Note that, in the above,

- the parameter values of $b_{f}$, $b_{\epsilon}$ and $p_{\epsilon}$
  have been optimised via a grid search to lead to suitable operating
  characteristics,
- the parameters values of $\delta_{\epsilon}$ and $\delta_{f}$ have
  been defined by clinicians.

## Self-defined R functions

In the following, we define the group allocation, futility and efficacy
functions corresponding to the design described above

#### Group allocation

We need to generate a function that randomises the `m` participants of
the next look according to the allocation ratios `prob`, where `m` and
`prob` are ingredients described in the [Ingredients
section](Ingredients.md).

``` r
# function
treatalloc.fun = function(m,prob){
  prob = abs(prob)/sum(abs(prob)) 
  m0.g = floor(prob*m)
  m0   = sum(m0.g)
  factor(rep(names(prob),m0.g+rmultinom(1,m-m0,prob)),
         levels=names(prob))
}
# test on m = 40 patients and equal allocation per group
table(treatalloc.fun(m=40,prob=c(control=.25,A=.25,B=.25,C=.25)))
table(treatalloc.fun(m=40,prob=c(control=1,A=1,B=1,C=1)))
```

`treatalloc.fun` first allocates the largest possible number of units to
the different groups given their exact target probabilities and then
assigns randomly the remaining units to the different groups according
to multinomial draws.

Note that

- this corresponds to the function `alloc.balanced` available in the
  **BATSS** package.
- when no RAR rule is used, the allocation probabilities used throughout
  the trial by `batss.glm` are the ones indicated under `prob0`.

#### Arm efficacy stopping rule

We need to generate a function that leads to a logical output and takes
as input

- the ingredients
  - `posterior` for the posterior probability of the target parameter
    being smaller than `delta.eff = 0`,
  - `n` and `N`, respectively the sample size per arm at the look of
    interest and the max sample size to define the information fraction,
- the additional parameters (to be added to `eff.arm.control` in
  `batss.glm`)
  - $b_{\epsilon}$ that we will name `b.eff`,
  - $p_{\epsilon}$ that we will name `p.eff`.

``` r
# function
efficacy.arm.fun = function(posterior,n,N,b.eff,p.eff){
  posterior > (1-(b.eff*(sum(n)/N)^p.eff))
}
# test at interim 1 for a parameter with a posterior = 0.999
efficacy.arm.fun(0.999, n=c(control=25,A=25,B=25,C=25),N=260,
                 b.eff = 0.009,p.eff=3)
# test at interim 1 for a parameter with a posterior = 0.9995
efficacy.arm.fun(0.9995, n=c(control=25,A=25,B=25,C=25),N=260,
                 b.eff = 0.009,p.eff=3)
```

Note that this corresponds to the function `eff.arm.infofract` available
in the **BATSS** package.

#### Trial efficacy stopping rule

We need to generate a function that, based on the ingredient
`eff.target` (indicating if efficacy was reached for each target
parameter at the stage of interest or before) leads to a logical output
equal to `TRUE` if all target parameters reached efficacy and `FALSE`
otherwise.

``` r
# function
efficacy.trial.fun = function(eff.target){
  all(eff.target)
}

# test 
efficacy.trial.fun(c(A=TRUE,B=TRUE,C=TRUE))
efficacy.trial.fun(c(A=TRUE,B=TRUE,C=FALSE))
```

Note that this corresponds to the function `eff.trial.all` available in
the **BATSS** package and that this behaviour is the default one of the
function `batss.glm` when `eff.trial = NULL`.

#### Arm futility stopping rule

We need to generate a function that leads to a logical output and takes
as input

- the ingredients `posterior` for the posterior probability of the
  target parameter being smaller than `delta.fut = log(0.8)`,
- the additional parameters $b_{f}$ that we will name `b.fut` (and that
  needs to be added to `fut.arm.control` in `batss.glm`)

``` r
# function
futility.arm.fun = function(posterior,b.fut){
  posterior < b.fut
}
# test 
futility.arm.fun(0.9, b.fut=.2025)
futility.arm.fun(0.1, b.fut=.2025)
```

Note that this corresponds to the function `fut.arm.simple` available in
the **BATSS** package.

#### Trial futility stopping rule

We need to generate a function that, based on the ingredient
`fut.target` (indicating if futility was declared for each target
parameter at the stage of interest or before) leads to a logical output
equal to `TRUE` if all target parameters were declared futile and
`FALSE` otherwise.

``` r
# function
futility.trial.fun = function(fut.target){
  all(fut.target)
}

# test 
futility.trial.fun(c(A=TRUE,B=TRUE,C=TRUE))
futility.trial.fun(c(A=TRUE,B=TRUE,C=FALSE))
```

Note that this corresponds to the function `fut.trial.all` available in
the **BATSS** package and that this behaviour is the default one of the
function `batss.glm` when `fut.trial = NULL`.

## Monte Carlo Simulations

We consider three scenarios:

- **Scenario 1 = ‘global null’**: each arm has a cumulative mean number
  of lesions equal to 4,
- **Scenario 2 = ‘one treatment works’**: each arm has a cumulative mean
  number of lesions equal to 4 except the last one, “C”, that shows a
  60% reduction in the cumulative mean number of lesions compared to the
  other groups (RR = 0.4),
- **Scenario 3 = ‘better, best’**: the control arm has a cumulative mean
  number of lesions equal to 4 and the other arms respectively lead to a
  20%, 40% and 60% reduction in the cumulative mean number of lesions
  compared to the reference group (i.e., RR = 0.8 for “A”, RR = 0.6 for
  “B” and RR = 0.4 for “C”).

### Adaptive designs

We first show the code corresponding to the planned adaptive design for
each Scenario:

#### Scenario 1

``` r
R = 25

scenario1 = batss.glm(   
  model           = y~treatment,
  var             = list(y = rnbinom,          
                         treatment = treatalloc.fun),
  var.control     = list(y = list(size = 1/2)), 
  family          = "nbinomial",
  link            = "log",
  beta            = c(log(4),log(1),log(1),log(1)),
  which           = c(2:4),
  R               = R,
  alternative     = c("less"),
  RAR             = NULL,
  prob0           = c(control = .25, A = .25, B = .25, C = .25),
  N               = 260,
  interim         = list(recruited=c(100,140,180,220)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = 0, 
  eff.arm.control = list("b.eff"=0.009, "p.eff"=3),
  delta.fut       = log(0.8), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list("b.fut"=0.2025),
  computation     = "parallel",
  mc.cores        = 10,
  H0              = TRUE,
  extended        = 1) 
```

You can note that

- `y` is generated via the function `rnbinom` with nuisance parameter
  equal to 0.5 (`size = 1/2`) and expected values on the `log` scale of
  respectively
  - `log(4)` for the control group,
  - `log(4)+log(1) = log(4)` for all other groups (as contrasts of type
    treatment are used for the factor `treatment` in the self-defined
    function `treatalloc.fun`),
- the target parameters (corresponding to the shift in means of
  treatment arms “A”, “B” and “C” compared to the “control” on the log
  scale) are in position 2, 3 and 4 of the fitted coefficients obtained
  when using the formula `y~treatment`,
- `prob0` provides
  - the (equal) allocation probabilities throughout the trial (as
    `RAR = NULL`),
  - the names of the groups,
- `eff.arm` and `fut.arm` are set to the functions defined above (i.e.,
  `efficacy.arm.fun` and `futility.arm.fun`),
- `eff.trial` and `fut.trial` are not specified and therefore equal to
  `NULL` (default) which leads to the behaviour wished in this case,
- the values of the additional parameters of `efficacy.arm.fun` and
  `futility.arm.fun` (i.e., `b.eff`, `p.eff` and `b.fut`) are specified
  under `eff.arm.control` and `fut.arm.control`,
- `delta.eff` and `delta.fut` are respectively set to `0` and
  `log(0.8)`.

We chose here a low number of seeds/trials (`R=25`) to save time.

#### Scenario 2

``` r
R = 25

scenario2 = batss.glm(   
  model           = y~treatment,
  var             = list(y = rnbinom,          
                         treatment = treatalloc.fun),
  var.control     = list(y = list(size = 1/2)), 
  family          = "nbinomial",
  link            = "log",
  beta            = c(log(4),log(1),log(1),log(0.4)),
  which           = c(2:4),
  R               = R,
  alternative     = c("less"),
  RAR             = NULL,
  prob0           = c(control = .25, A = .25, B = .25, C = .25),
  N               = 260,
  interim         = list(recruited=c(100,140,180,220)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = 0, 
  eff.arm.control = list("b.eff"=0.009, "p.eff"=3),
  delta.fut       = log(0.8), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list("b.fut"=0.2025),
  computation     = "parallel",
  mc.cores        = 10,
  H0              = TRUE,
  extended        = 1) 
```

Same as above except for `beta` which now leads to a mean of

- `exp(log(4)) = 4` for the control group,
- `exp(log(4) + log(1)) = 4` (no effect) for the groups “A” and “B”,
- `exp(log(4) + log(0.4)) = 1.6` (target effect) for group “C” which
  corresponds to a RR of `1.6/4 = 0.4`.

#### Scenario 3

``` r
R = 25

scenario3 = batss.glm(   
  model           = y~treatment,
  var             = list(y = rnbinom,          
                         treatment = treatalloc.fun),
  var.control     = list(y = list(size = 1/2)), 
  family          = "nbinomial",
  link            = "log",
  beta            = c(log(4),log(0.8),log(0.6),log(0.4)),
  which           = c(2:4),
  R               = R,
  alternative     = c("less"),
  RAR             = NULL,
  prob0           = c(control = .25, A = .25, B = .25, C = .25),
  N               = 260,
  interim         = list(recruited=c(100,140,180,220)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = 0, 
  eff.arm.control = list("b.eff"=0.009, "p.eff"=3),
  delta.fut       = log(0.8), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list("b.fut"=0.2025),
  computation     = "parallel",
  mc.cores        = 10,
  H0              = TRUE,
  extended        = 1) 
```

Same as above except for `beta` which now leads to a mean of

- `exp(log(4)) = 4` for the control group,
- `exp(log(4) + log(0.8)) = 3.2` (small effect) for group “A” which
  corresponds to a RR of `3.2/4 = 0.8`,
- `exp(log(4) + log(0.6)) = 2.4` (medium effect) for group “B” which
  corresponds to a RR of `2.4/4 = 0.6`,
- `exp(log(4) + log(0.4)) = 1.6` (target effect) for group “C” which
  corresponds to a RR of `1.6/4 = 0.4`.

### Fixed designs

We now show how to use `batss.glm` to define the operating
characteristics of the fixed design when considering the same efficacy
and futiliy parameters as well as the same total maximum sample size of
$N = 260$ patients. The easiest way to fit a fixed design with **BATSS**
is to plan a single interim analysis defined so that no adaptations are
allowed. This is ensured by setting the `delta.eff` and `delta.fut`
values to `NA` for the interim (as well as by setting the `RAR` argument
to `NULL` if RAR was considered in the design). Note that the sample
size at the dummy interim should be a multiple of the number of arms in
the trial.

#### Scenario 1

Same as the adaptive Scenario 1 above, except `interim` is reduced to a
single dummy interim at 100 participants, and `delta.eff` and
`delta.fut` are set to `c(NA, 0)` and `c(NA, log(0.8))` respectively, so
that no decisions are made at the interim.

``` r
R = 25

scenario1_fixed = batss.glm(   
  model           = y~treatment,
  var             = list(y = rnbinom,          
                         treatment = treatalloc.fun),
  var.control     = list(y = list(size = 1/2)), 
  family          = "nbinomial",
  link            = "log",
  beta            = c(log(4),log(1),log(1),log(1)),
  which           = c(2:4),
  R               = R,
  alternative     = c("less"),
  RAR             = NULL,
  prob0           = c(control = .25, A = .25, B = .25, C = .25),
  N               = 260,
  interim         = list(recruited=c(100)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA,0), 
  eff.arm.control = list("b.eff"=0.009, "p.eff"=3),
  delta.fut       = c(NA,log(0.8)), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list("b.fut"=0.2025),
  computation     = "parallel",
  mc.cores        = 10,
  H0              = TRUE,
  extended        = 1) 
```

#### Scenario 2

Same as the adaptive Scenario 2 above, except `interim` is reduced to a
single dummy interim at 100 participants, and `delta.eff` and
`delta.fut` are set to `c(NA, 0)` and `c(NA, log(0.8))` respectively, so
that no decisions are made at the interim.

``` r
R = 25

scenario2_fixed = batss.glm(   
  model           = y~treatment,
  var             = list(y = rnbinom,          
                         treatment = treatalloc.fun),
  var.control     = list(y = list(size = 1/2)), 
  family          = "nbinomial",
  link            = "log",
  beta            = c(log(4),log(1),log(1),log(0.4)),
  which           = c(2:4),
  R               = R,
  alternative     = c("less"),
  RAR             = NULL,
  prob0           = c(control = .25, A = .25, B = .25, C = .25),
  N               = 260,
  interim         = list(recruited=c(100)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA,0), 
  eff.arm.control = list("b.eff"=0.009, "p.eff"=3),
  delta.fut       = c(NA,log(0.8)), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list("b.fut"=0.2025),
  computation     = "parallel",
  mc.cores        = 10,
  H0              = TRUE,
  extended        = 1) 
```

#### Scenario 3

Same as the adaptive Scenario 3 above, except `interim` is reduced to a
single dummy interim at 100 participants, and `delta.eff` and
`delta.fut` are set to `c(NA, 0)` and `c(NA, log(0.8))` respectively, so
that no decisions are made at the interim.

``` r
R = 25

scenario3_fixed = batss.glm(   
  model           = y~treatment,
  var             = list(y = rnbinom,          
                         treatment = treatalloc.fun),
  var.control     = list(y = list(size = 1/2)), 
  family          = "nbinomial",
  link            = "log",
  beta            = c(log(4),log(0.8),log(0.6),log(0.4)),
  which           = c(2:4),
  R               = R,
  alternative     = c("less"),
  RAR             = NULL,
  prob0           = c(control = .25, A = .25, B = .25, C = .25),
  N               = 260,
  interim         = list(recruited=c(100)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA,0), 
  eff.arm.control = list("b.eff"=0.009, "p.eff"=3),
  delta.fut       = c(NA,log(0.8)), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list("b.fut"=0.2025),
  computation     = "parallel",
  mc.cores        = 10,
  H0              = TRUE,
  extended        = 1) 
```

## References

Curtin, Francois, Herve Porchet, Robert Glanzman, Hans Martin Schneble,
Virginie Vidal, Marie-Laure Audoli-Inthavong, Estelle Lambert, and Hans
Peter Hartung. 2016. “A Placebo Randomized Controlled Study to Test the
Efficacy and Safety of GNbAC1, a Monoclonal Antibody for the Treatment
of Multiple Sclerosis – Rationale and Design.” *Multiple Sclerosis and
Related Disorders* 9: 95–100.
https://doi.org/<https://doi.org/10.1016/j.msard.2016.07.002>.

Gotmaker, Robert, Michael J Barrington, John Reynolds, Lorenzo Trippa,
and Stephane Heritier. 2019. “Bayesian Adaptive Design: The Future for
Regional Anesthesia Trials?” *Regional Anesthesia & Pain Medicine* 44
(6): 617–22. <https://doi.org/10.1136/rapm-2018-100248>.

Hartung, Hans-Peter, Tobias Derfuss, Bruce AC Cree, Maria Pia Sormani,
Krzysztof Selmaj, Jonathan Stutters, Ferran Prados, et al. 2022.
“Efficacy and Safety of Temelimab in Multiple Sclerosis: Results of a
Randomized Phase 2b and Extension Study.” *Multiple Sclerosis Journal*
28 (3): 429–40. <https://doi.org/10.1177/13524585211024997>.
