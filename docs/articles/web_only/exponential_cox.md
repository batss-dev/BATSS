# Exponential time-to-event endpoint

## Context

The case study is motivated by the *rEECur* trial, a seamless phase
II/III MAMS trial that investigated the optimal systemic anticancer
regimen for recurrent and refractory Ewing sarcoma (McCabe et al.
(2019)). As opposed to the original trial, this is designed as a
superiority trial where we are interested in demonstrating efficacy of
each intervention arm compared to a standard of care. No head-to-head
comparisons of interventions are made.

The primary endpoint is event-free survival time (EFS) - the time from
randomisation to the first failure event (progression, recurrence,
diagnosis of a second malignancy, or death) and is assumed to follow an
exponential distribution. In the original trial the 1-year EFS rate was
considered as an endpoint, here we look at a time-to-event endpoint.

## Design

The chosen design has the following characteristics:

- **treatment arms:** This case study has 4 treatment arms: standard of
  care (arm “control”), three investigational/novel drug regimen
  (referred to as arms “A”, “B”, and “C”).

- **alternative hypotheses:** This is a superiority trial where we are
  interested in demonstrating efficacy compared to the standard of care
  arm and intervention arms will not be compared to one another.

- **interim analyses and maximum sample size:** We will conduct an
  interim analyses after 1 year of recruitment and every 1 year (12
  months) thereafter until year 6. We will conduct a final analysis
  after the last participant has completed follow-up if no decision to
  stop the trial has been made earlier. Our maximum sample size is 800
  participants (across all tresatments).

- **participant accrual and censoring:** We assume exponentially
  distributed time intervals between participants’ enrolment into the
  study (i.e. a Poisson process) with a rate of 200 participants per
  year. This will result in an average enrollment period of 4 years. In
  a second scenario, we will look at a yearly-changing exponential
  recruitment rate with a slower start. Further, we will follow-up
  participants until the last participant experienced an event or 5
  years of follow-up after the last participant is recruited. No
  participant attrition will be considered for this example (i.e. no
  loss-to-follow-up).

- **endpoint conditional distribution:** We will use an exponential
  model for the data generating process and a Cox proportional hazards
  model for the analyses.

- **group allocation:** we will consider equal allocation probabilities
  per group (no response-adaptive randomisation).

- **efficacy stopping rule:** Early stopping of intervention arms for
  efficacy may occur if there is a high posterior probability of (any)
  benefit (i.e. ${HR} < 1$) at look $j$, i.e., when
  $${Prob}\left( \beta_{k} < \delta_{e}|{y,X} \right) > b_{e,i}$$ where
  $\beta_{k}$ denotes the linear predictor of the $k$th target
  parameter, $\delta_{e} = 0$ denotes the (efficacy-related) clinically
  meaningful treatment effect and $b_{e,i}$ the cut-off values used to
  declare futility (here, posterior probabilities). We will be using a
  stricter threshold $b_{e,1}$ during interim analyses and a more
  lenient value of $b_{e,2}$ during the final analysis.

- **futility stopping rule:** Early stopping of intervention arms for
  futility may occur if there is a low posterior probability of seeing
  any difference (i.e., ${HR} < 1$), i.e., when
  $${Prob}\left( \beta_{k} < \delta_{f}|{y,X} \right) < b_{f}$$ where
  $\delta_{f} = 0$ denotes the clinically meaningful treatment effect
  and $b_{f}$ is the cut-off value to declare futility.

- **trial stopping rule:** The trial will run until an efficacy decision
  has been reached for any intervention arm, a futility decision has
  been reached for each intervention arm or once the maximum sample size
  has been reached and follow-up period ended.

## Self-defined R functions

In the following, we define the group allocation, accrual, efficacy and
futility functions corresponding to the design described above.

### Group allocation

We will use the function `alloc.balanced` from the **BATSS** package,
which first allocates the largest possible number of units to the
different groups given their exact target probabilities and then assigns
randomly the remaining units to the different groups according to
multinomial draws.

### Participant accrual

In the standard scenario with constant accrual we can just use the base
R function `rexp` to get the exponentially distributed intervals,
i.e. `accr = rexp`, `accr.control = list(rate = 200)`.

For the scenario with a changing accrual rate, we need to generate a
function that allows us to change the rates of an exponential
distribution after certain time intervals. For simplicity, we will draw
from a exponential distribution with rate $\lambda_{1}$ and change the
rates for the first interval after the cumulative sum of previously
drawn random numbers exceeds the threshold value. The function
[`batss.surv()`](../../reference/batss.surv.md) will expect an input
value `n` indicating the number of generated numbers as the first
element of the function (as for any R base random number generating
function). The additional parameters can be passed on by using the
`accr.control` argument in `batss.surv`, in this example
`accr.control = list(rates = c(100, 180, 260), changes.at = c(1, 2))`.

``` r
# function
changing.accrual <- function(n,rates,changes.at){
    i <- 1
    res <- rexp(1,rate=rates[i])
    pos <- 2
    while (i <= length(changes.at)) {
      while(sum(res)<changes.at[i]) {
        res[pos] <- rexp(1,rate=rates[i])
        pos <- pos+1
      }
      i <- i+1
    }
    res <- c(res,rexp(n-(pos-1),rate=rates[i]))
    return(res)
} 

# test
changing.accrual(20, rates = c(5,20), changes.at = 1)
```

### Arm efficacy stopping rule

We need to generate a function that leads to a logical output and takes
as input

- the ingredients

  - `posterior` for the posterior probability of the target parameter
    being smaller than `delta.eff = 0`,
  - `curr.look` and `n.look`, respectively the number of the current
    interim analysis and the total number of interim analyses,

- the additional parameters (to be added to `eff.arm.control` in
  `batss.surv`)

  - $b_{e}$ that we will name `b`.

``` r
# function
efficacy.fun <- function(posterior, b, curr.look, n.look){
  if (curr.look != n.look) {
    posterior > b[1]
  } else {
    posterior > b[2]
  }
} 

# test
efficacy.fun(0.85, b = c(0.99, 0.95), curr.look = 2, n.look = 4)
```

### Trial efficacy stopping rule

We will use the function `eff.trial.any` from the **BATSS** package,
which will stop a trial if any target parameter reached efficacy
(indicated by the ingredient `eff.target`).

### Arm futility stopping rule

We can select the function `fut.arm.simple` from the package **BATSS**
to compare the ingredient `posterior`, the posterior probability of the
target parameter being smaller than `delta.fut = 0`, to the threshold
value `b = 0.05`.

### Trial futility stopping

We will stop the trial for futility if all target parameters were
declared ineffective. This corresponds to the (default) function
`fut.trial.all` available in the **BATSS** package.

## Monte Carlo Simulations

We will consider two scenarios for both the constant recruitment rate
and the changing recruitment rate:

- **Scenario 1 = ‘global null’:** the hazard ratio is equal to 1 for
  each intervention arm compared to the control.

- **Scenario 2 = ‘one treatment works’:** arms “B” and “C” have a hazard
  ratio of 1 compared to the control while the hazard ratio for
  treatment “A” equals 0.75. This corresponds with an increase of the
  1-year EFS rate from 0.2 to 0.3 (with our parametrisation of the data
  generating exponential distribution).

### Scenario 1 & 2

``` r
# number of trials
R <- 25

# simulation
rEECur.sim <- batss.surv(
  model           = inla.surv(time, status) ~ trt,
  family          = "coxph",
  surv            = simsurv,
  surv.control    = list(lambdas = -log(0.2), gammas = 1,
                         maxt = NULL),
  fup             = 5,
  var             = list(trt = alloc.balanced),
  accr            = rexp,
  accr.control    = list(rate = 200),
  accr.type       = "random",
  prob0           = c(control = 1, A = 1, B = 1, C = 1),
  hr              = c(0.75, 1, 1),
  which           = 1:3,
  alternative     = "less",
  RAR             = NULL,
  interim         = list(time = 1:6),
  eff.arm         = efficacy.fun,
  eff.arm.control = list(b = c(0.99, 0.95)),
  eff.trial       = eff.trial.any,
  delta.eff       = 0,
  fut.arm         = fut.arm.simple,
  fut.arm.control = list(b = 0.05),
  N               = 800,
  R               = R,
  H0              = TRUE,
  extended        = 1,
  computation     = "parallel",
  mc.cores        = parallel::detectCores() - 1,
  control.inla    = list(cmin = 0)
)
```

You can note

- The tuple (`time`, `status`) is generated via the `simsurv` function
  from the **simsurv** package. Arguments are passed via the
  `surv.control` argument. The default distribution in `simsurv` is
  Weibull. We set $\lambda = {log}(0.2)$ to get a 1-year EFS rate of
  0.2.

- the `family` argument is linked to the analytical model but shares the
  linear predictor with the data generating model,

- participant accrual is set in `accr`, `accr.control` and `accr.type`,

- the censoring of participants is done by a combination of
  `maxt = NULL` in `surv.control` (no individual administrative
  censoring) and `fup = 5` (the maximum follow-up time after recruitment
  of the last participant - ignored if `maxt` is set),

- `prob0` provides

  - the (equal) allocation probabilities at the start of the trial,
  - the names of the groups,

- the interim analyses schedule is set to time-based,
  i.e. `interim = list(time = 1:6)` for interim analyses after 1, 2,
  3,.. , 6 years, respectively, (See the vignette “Definition of the
  interim schedule in BATSS” to learn how to change this.)

- `eff.arm` is set to the function defined above (i.e., `efficacy.fun`)
  and the cut-off values for declaring an intervention efficacious are
  set to a probability of 0.99 during the trial and 0.95 at the final
  analysis in `eff.arm.control`.

- `fut.arm` is set to the standard option `fut.arm.simple` with
  additional parameters specified in `fut.arm.control`,

- `eff.trial` is set to `eff.trial.any` as discussed above,

- `fut.trial` is not specified and therefore equal to NULL (default)
  which leads to the desired behaviour in this case,

- `delta.fut` is set to 0 by default after setting `delta.eff = 0`,

- `H0 = TRUE` will additionally calculate the null hypotheses of all
  hazard ratios are equal to 1,

- `control.inla` passes information to the INLA routine via `...`,
  setting `cmin = 0` will stabilize the numerical optimization by
  setting the minimum value for the negative Hessian, but may bias
  estimates!

Attention: The number of repetitions is set to `R = 25` for
demonstration purposes only.

### Scenario 3 & 4 - changing accrual rates

``` r
# number of trials
R <- 25

# simulation
rEECur.sim.change <- batss.surv(
  model           = inla.surv(time, status) ~ trt,
  family          = "coxph",
  surv            = simsurv,
  surv.control    = list(lambdas = -log(0.2), gammas = 1,
                         maxt = NULL),
  fup             = 5,
  var             = list(trt = alloc.balanced),
  accr            = changing.accrual,
  accr.control    = list(rates = c(100, 180, 260), changes.at = c(1, 2)),
  accr.type       = "random",
  prob0           = c(control = 1, A = 1, B = 1, C = 1),
  hr              = c(0.75, 1, 1),
  which           = 1:3,
  alternative     = "less",
  RAR             = NULL,
  interim         = list(time = 1:6),
  eff.arm         = efficacy.fun,
  eff.arm.control = list(b = c(0.99, 0.95)),
  eff.trial       = eff.trial.any,
  delta.eff       = 0,
  fut.arm         = fut.arm.simple,
  fut.arm.control = list(b = 0.05),
  N               = 800,
  R               = R,
  H0              = TRUE,
  extended        = 1,
  computation     = "parallel",
  mc.cores        = parallel::detectCores() - 1,
  control.inla    = list(cmin = 0)
)
```

Same as above except for `accr` and `accr.control`.

## References

McCabe, Martin G., Veronica Moroz, Maria Khan, Uta Dirksen, Abigail
Evans, Nicola Fenwick, Nathalie Gaspar, et al. 2019. “Results of the
first interim assessment of rEECur, an international randomized
controlled trial of chemotherapy for the treatment of recurrent and
primary refractory Ewing sarcoma.” *Journal of Clinical Oncology* 37
(15_suppl): 11007–7.
[https://doi.org/10.1200/JCO.2019.37.15\\suppl.11007](https://doi.org/10.1200/JCO.2019.37.15%5C_suppl.11007).
