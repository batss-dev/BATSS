# Binomial endpoint

## Context

This example is motivated by the Personalised Immunotherapy Platform
(PIP) trial (Lo et al. 2022), which is a phase II randomised controlled
trial that is currently under development that will compare
immunotherapy combinations for advanced melanoma patients.

This study aims to assess the effectiveness of different novel drug
combinations compared to a control arm in treatment-naive patients that
are predicted to be unlikely to respond to standard first line
immunotherapy treatments. The primary outcome is the objective response
rate (ORR, i.e., the probability of observing a one/success) according
to the Response Evaluation Criteria in Solid Tumors (RECIST) (partial or
complete response vs stable disease or progressive disease) at 6 months
post randomisation.

## Design

The chosen design has the following characteristics:

- **treatment arms**: Our case study has 6 treatment arms: standard of
  care (arm “A”), five drug combinations (referred to as arms “B”, “C”,
  “D”, “E” and “F”).
- **alternative hypotheses**: This is a superiority trial where we are
  interested in demonstrating efficacy compared to the standard of care
  arm and intervention arms will not be compared to one another.
- **interim analyses and maximum sample size**: We will begin the
  interim analyses at 60 participants completing follow-up, and perform
  interim analyses every 12 participants completing follow-up
  thereafter. Our maximum sample size is 216 participants and so a total
  of 14 looks are planned.
- **endpoint conditional distribution**: the binary endpoint is assumed
  to follow a binomial distribution.
- **group allocation**: the trial will start with equal randomisation
  (1:1:1:1:1:1) and RAR will be used from interim analysis 1 onwards.
  The randomisation probabilities are calculated using a variation of
  Trippa et al. (2012) and Cellamare et al. (2017) where the allocation
  probabilities at look $j$ are proportional to:
  - control arm:
    $$\frac{\text{exp}\left( \text{max}\left( n_{kj} \right) - n_{0j} \right)^{\nu}}{K_{j}}$$
    where
    - $n_{kj}$ corresponds to the sample size of arm $k$ (where $k = 0$
      denotes the control arm) at stage $j$ and
      $\text{max}\left( n_{kj} \right)$ is evaluated across intervention
      arms only (i.e. for $k > 0$),
    - $K_{j}$ is the number of active arms (excluding control) at look
      $j$,
    - $\nu = 0.1$ is a tuning parameter.
  - intervention arms:
    $$\frac{\text{Prob}\left( \beta_{k} > \delta_{r}|y,X \right)^{h}}{\sum\limits_{k = 1}^{K}\text{Prob}\left( \beta_{k} > \delta_{r}|y,X \right)^{h}}$$
    where
    - $\beta_{k}$ denotes the $k$th target parameter (on the logit
      scale),
    - $\delta_{r} = 0$ denotes the (RAR-related) clinically meaningful
      treatment effect value,
    - and $h = \gamma\left( \sum_{k = 0}^{K}n_{kj}/N \right)^{\eta}$
      with
      - $n_{kj}$, the sample size of arm $k$ at stage $j$,
      - $N$, the maximum sample size,
      - $\gamma = 3$ and $\eta = 1.4$, two tuning parameters. Refer to
        Gotmaker et al. (2019) for details
- **efficacy stopping rule**: Stopping arms for efficacy is not
  permitted at the interim analyses as the main objective of the
  adaptations is to drop poorly performing arms.
- **end of trial efficacy rule**: An intervention may be declared
  effective at the end of the trial if the posterior probability of
  efficacy is above a certain (high) threshold:
  $$\text{Prob}\left( \beta_{k} > \delta_{\epsilon}|y,X \right) > 1 - b_{\epsilon}$$
  where
  - $\delta_{\epsilon} = 0$ denotes the (efficacy-related) clinically
    meaningful treatment effect value,
  - $b_{\epsilon} = 0.045$ is a tuning parameter.
- **futility stopping rule**: Early stopping of intervention arms for
  futility may occur if there is a low posterior probability of
  observing at least a 10% absolute increase (improvement) in the
  proportion of responders, i.e., when
  $$\text{Prob}\left( \beta_{k} > \delta_{f}|y,X \right) < b_{f}$$ where
  - $\delta_{f} = log(1.5)$ denotes the (futility-related) clinically
    meaningful treatment effect values assuming a control ORR of 40%,
  - $b_{f} = 0.1$ is the cut-off value used to declare futility.
- **trial stopping rule:** The trial will run until a futility decision
  has been reached for each intervention arm or once the maximum sample
  size has been reached.

Note that, in the above, the tuning parameter values have been optimised
to lead to suitable operating characteristics.

In the following, we define the group allocation, RAR, futility and
efficacy functions corresponding to the design described above

#### RAR

We need to generate a function that defines the allocation probabilities
of patients to the different active arms as a function of

- relevant *ingredients*, i.e.,
  - `posterior`, for the posterior probabilities of the target parameter
    being greater than `delta.RAR = 0`,
  - `n` and `N`, respectively the sample size per arm at the look of
    interest and the max sample size to define the information fraction,
  - `ref`, for a vector of logicals indicating which group is the
    reference one,
  - `active`, for a vector of logicals indicating which groups are
    active at the look of interest,
- relevant tuning parameters (to be added to `RAR.control` in
  `batss.glm`), namely $\gamma$, $\eta$, and $\nu$, that we will
  respectively call `gamma`, `eta` and `nu`

``` r
# function
prob.trippa = function(posterior,n,N,ref,active,gamma,eta,nu){
  g = sum(active)
  h = gamma*(sum(n)/N)^eta
  prob = rep(NA,g)
  # reference/control arm allocation
  prob[1] = (exp(max(n[!ref])-n[ref])^nu)/(g-1)
  # targets/interventions (that haven't been dropped)
  prob[2:g] = (posterior^h)/(sum(posterior^h))
  unlist(prob)   
}
# test after at the first look with two set of posterior probabilities
# 1/ c(.5,.5,.5,.5,.5)
# 2/ c(.5,.5,.5,.5,.6)
prob.trippa(c(.5,.5,.5,.5,.5),n=c(A=10,B=10,C=10,D=10,E=10,F=10),N=16,
            ref = c(TRUE,rep(FALSE,5)), active = rep(TRUE,6),
            gamma=3, eta=1.4, nu=0.1)
prob.trippa(c(.5,.5,.5,.5,.6),n=c(A=10,B=10,C=10,D=10,E=10,F=10),N=16,
            ref = c(TRUE,rep(FALSE,5)), active = rep(TRUE,6),
            gamma=3, eta=1.4, nu=0.1)
```

Note that this corresponds to the function `RAR.trippa` available in the
**BATSS** package.

#### Group allocation

We need to generate a function that randomises the `m` participants of
the next look according to the allocation ratios `prob`, where `m` and
`prob` are ingredients described in the [Ingredients
section](Ingredients.md) and where `prob` is the output of the
user-defined RAR function.

``` r
# function
treatalloc.fun = function(m,prob){
  prob = abs(prob)/sum(abs(prob)) 
  m0.g = floor(prob*m)
  m0   = sum(m0.g)
  factor(rep(names(prob),m0.g+rmultinom(1,m-m0,prob)),
         levels=names(prob))
}
# test on m = 60 patients and equal allocation per group
table(treatalloc.fun(m=60,prob=c(A=1,B=1,C=1,D=1,E=1,F=1)))
```

`treatalloc.fun` first allocates the largest possible number of units to
the different groups given their exact target probabilities and then
assigns randomly the remaining units to the different groups according
to multinomial draws.

Note that this corresponds to the function `alloc.balanced` available in
the **BATSS** package.

#### Efficacy rule

We need to generate a function that leads to a logical output when used
at the last look of the trial and takes as input

- the ingredient `posterior` for the posterior probability of the target
  parameter being greater than `delta.eff = 0`,
- the additional parameter $b_{\epsilon}$ that we will name `b.eff` (and
  that needs to be added to `eff.arm.control` in `batss.glm`).

``` r
# function
efficacy.arm.fun = function(posterior,b.eff){
  posterior > 1-b.eff 
}
# test for a parameter with a posterior = 0.999
efficacy.arm.fun(0.999, b.eff = 0.045)
# test for a parameter with a posterior = 0.95
efficacy.arm.fun(0.95, b.eff = 0.045)
```

#### Arm futility stopping rule

We need to generate a function that leads to a logical output and takes
as input

- the ingredient `posterior` for the posterior probability of the target
  parameter being greater than `delta.fut = log(1.5)`,
- the additional parameter $b_{f}$ that we will name `b.fut` (and that
  needs to be added to `fut.arm.control` in `batss.glm`)

``` r
# function
futility.arm.fun = function(posterior,b.fut){
  posterior < b.fut
}
# test 
futility.arm.fun(0.9, b.fut=0.1)
futility.arm.fun(0.075, b.fut=0.1)
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
futility.trial.fun(c(B=TRUE,C=TRUE,D=TRUE,E=TRUE,F=TRUE))
futility.trial.fun(c(B=TRUE,C=TRUE,D=TRUE,E=TRUE,F=FALSE))
```

Note that this corresponds to the function `fut.trial.all` available in
the **BATSS** package and that this behaviour is the default one of the
function `batss.glm` when `fut.trial = NULL`.

## Monte Carlo Simulations

We consider two scenarios:

- **Scenario 1**: each arm has an ORR of 0.4 (global null),
- **Scenario 2**:
  - arms “A”, “B” and “C” have an ORR of 0.4,
  - arm “D” has an ORR of 0.5,
  - arms “E” and “F” have an ORR of 0.7.

### Adaptive designs

We first show the code corresponding to the planned adaptive design for
each Scenario:

#### Scenario 1

``` r
# number of trials
R = 25

# logit function
logit = function(p){log(p/(1 - p))}

# simulation
scenario1 = batss.glm(   
  model           = y~group,
  var             = list(y = rbinom,
                         group = treatalloc.fun),
  var.control     = list(y = list(size = 1)), 
  family          = "binomial",
  link            = "logit",
  beta            = c(logit(0.4), rep(0,5)),
  which           = c(2:6),
  R               = R,
  alternative     = c("greater"),
  RAR             = prob.trippa,
  RAR.control     = list("gamma"=3, "eta"=1.4,"nu"=0.1),
  delta.RAR       = 0,
  prob0           = c(A=1,B=1,C=1,D=1,E=1,F=1),
  N               = 216, 
  interim         = list(recruited=list(m0=60, m=12)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(rep(NA,13), 0), 
  eff.arm.control = list(b.eff = 0.045),
  delta.fut       = log(1.5), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.1),
  computation     = "parallel",
  mc.cores        = 12,
  H0              = FALSE,
  extended        = 1)   
```

You can note that

- `y` is generated via the function `rbinom` and has an ORR of
  respectively
  - `expit(logit(0.4)) = 0.4` for the reference group, where
    `logit(0.4) = -0.405` corresponds to the log odds of ‘success’ for
    the control group,
  - `expit(logit(0.4+0)) = 0.4` for all other groups (as contrasts of
    type treatment are used for the factor `treatment` in the
    self-defined function `treatalloc.fun`), where `0` corresponds to
    the log odds ratio value related to groups “B” to “F” (i.e., no
    effect),
- the target parameters (corresponding to the shift in ORR of treatment
  arms “B” to “F” compared to the “control” on the logit scale) are in
  position 2 to 6 of the fitted coefficients obtained when using the
  formula `y~group`,
- `prob0` provides
  - the (equal) allocation probabilities at the start of the trial,
  - the names of the groups,
- `eff.arm` and `fut.arm` are set to the functions defined above (i.e.,
  `efficacy.arm.fun` and `futility.arm.fun`),
- `delta.eff` is set to `NA` for all 13 interim analyses and to `0` for
  the last look. This ensures that the trial can’t stop early for
  efficacy. This aim could be achieved using other strategies, like
  setting the output of the function indicated under `eff.trial` to
  `FALSE`,
- `fut.trial` is not specified and therefore equal to `NULL` (default)
  which leads to the behaviour wished in this case,
- the values of the additional parameters of `efficacy.arm.fun` and
  `futility.arm.fun` (i.e., `b.eff` and `b.fut`) are specified under
  `eff.arm.control` and `fut.arm.control`,
- `delta.eff` and `delta.fut` are respectively set to `0` and
  `log(1.5)`.

We chose here a low number of seeds/trials (`R=25`) to save time.

#### Scenario 2

``` r
# number of trials
R = 25

# logit function
logit = function(p){log(p/(1 - p))}

# simulation
scenario1 = batss.glm(   
  model           = y~group,
  var             = list(y = rbinom,
                         group = treatalloc.fun),
  var.control     = list(y = list(size = 1)), 
  family          = "binomial",
  link            = "logit",
  beta            = c(logit(0.4),                    # log odd of reference group
                      0,0,                           # log odds ratio of groups B and C       
                      logit(0.5)-logit(0.4),         # log odds ratio of group D
                      rep(logit(0.7)-logit(0.4),2)), # log odds ratio of groups E and F       
  which           = c(2:6),
  R               = R,
  alternative     = c("greater"),
  RAR             = prob.trippa,
  RAR.control     = list("gamma"=3, "eta"=1.4,"nu"=0.1),
  delta.RAR       = 0,
  prob0           = c(A=1,B=1,C=1,D=1,E=1,F=1),
  N               = 216, 
  interim         = list(recruited=list(m0=60, m=12)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(rep(NA,13), 0), 
  eff.arm.control = list(b.eff = 0.045),
  delta.fut       = log(1.5), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.1),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = FALSE,
  extended        = 1)  
```

Same as above except for `beta`.

### Fixed designs

We now show how to use `batss.glm` to define the operating
characteristics of the fixed design when considering the same efficacy
and futility parameters as well as the same total maximum sample size of
$N = 216$ patients. The easiest way to fit a fixed design with **BATSS**
is to plan a single interim analysis defined so that no adaptations are
allowed. This is ensured by setting the `delta.eff` and `delta.fut`
values to `NA` for the interim (as well as by setting the `RAR` argument
to `NULL` as treatment allocation after the interim analysis would be
based on posterior probabilities estimated on interim data otherwise).
Note that the sample size at the dummy interim should be a multiple of
the number of arms in the trial.

#### Scenario 1

Same as the adaptive Scenario 1 above, except `interim` is reduced to a
single dummy interim at 60 participants, and `delta.eff` and `delta.fut`
are set to `c(NA, 0)` and `c(NA, log(1.5))` respectively, and `RAR` is
set to `NULL` (and `RAR.control` and `delta.RAR` deleted), so that no
decisions are made at the interim and equal group allocation is used.

``` r
# number of trials
R = 25

# logit function
logit = function(p){log(p/(1 - p))}

# simulation
scenario1_fixed = batss.glm(   
  model           = y~group,
  var             = list(y = rbinom,
                         group = treatalloc.fun),
  var.control     = list(y = list(size = 1)), 
  family          = "binomial",
  link            = "logit",
  beta            = c(logit(0.4), rep(0,5)),
  which           = c(2:6),
  R               = R,
  alternative     = c("greater"),
  RAR             = NULL,
  prob0           = c(A=1,B=1,C=1,D=1,E=1,F=1),
  N               = 216, 
  interim         = list(recruited=list(60)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA, 0), 
  eff.arm.control = list(b.eff = 0.045),
  delta.fut       = c(NA,log(1.5)), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.1),
  computation     = "parallel",
  mc.cores        = 12,
  H0              = FALSE,
  extended        = 1)   
```

#### Scenario 2

Same as the adaptive Scenario 2 above, except `interim` is reduced to a
single dummy interim at 60 participants, and `delta.eff` and `delta.fut`
are set to `c(NA, 0)` and `c(NA, log(1.5))` respectively, and `RAR` is
set to `NULL` (and `RAR.control` and `delta.RAR` deleted), so that no
decisions are made at the interim and equal group allocation is used.

``` r
# number of trials
R = 25

# logit function
logit = function(p){log(p/(1 - p))}

# simulation
scenario2_fixed = batss.glm(   
  model           = y~group,
  var             = list(y = rbinom,
                         group = treatalloc.fun),
  var.control     = list(y = list(size = 1)), 
  family          = "binomial",
  link            = "logit",
  beta            = c(logit(0.4),                    # log odd of reference group
                      0,0,                           # log odds ratio of groups B and C       
                      logit(0.5)-logit(0.4),         # log odds ratio of group D
                      rep(logit(0.7)-logit(0.4),2)), # log odds ratio of groups E and F       
  which           = c(2:6),
  R               = R,
  alternative     = c("greater"),
  RAR             = NULL,
  prob0           = c(A=1,B=1,C=1,D=1,E=1,F=1),
  N               = 216, 
  interim         = list(recruited=list(60)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA, 0), 
  eff.arm.control = list(b.eff = 0.045),
  delta.fut       = c(NA,log(1.5)), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.1),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = FALSE,
  extended        = 1)  
```

## References

Cellamare, Matteo, Steffen Ventz, Elisabeth Baudin, Carole D. Mitnick,
and Lorenzo Trippa. 2017. “A Bayesian Response-Adaptive Trial in
Tuberculosis: The endTB Trial.” *Clinical Trials* 14 (1): 17–28.
<https://doi.org/10.1177/1740774516665090>.

Gotmaker, Robert, Michael J Barrington, John Reynolds, Lorenzo Trippa,
and Stephane Heritier. 2019. “Bayesian Adaptive Design: The Future for
Regional Anesthesia Trials?” *Regional Anesthesia & Pain Medicine* 44
(6): 617–22. <https://doi.org/10.1136/rapm-2018-100248>.

Lo, Serigne N., Tuba Nur Gide, Maria Gonzalez, Ines Silva, Alexander M.
Menzies, Matteo S. Carlino, Richard A. Scolyer, Stephane Heritier, James
S. Wilmott, and Georgina V. Long. 2022. “A Biomarker-Guided Bayesian
Response-Adaptive Phase II Trial for Metastatic Melanoma: The
Personalized Immunotherapy Platform (PIP) Trial Design.” *Journal of
Clinical Oncology* 40 (16): TPS9599–99.

Trippa, Lorenzo, Eudocia Q. Lee, Patrick Y. Wen, Tracy T. Batchelor,
Timothy Cloughesy, Giovanni Parmigiani, and Brian M. Alexander. 2012.
“Bayesian Adaptive Randomized Trial Design for Patients with Recurrent
Glioblastoma.” *Journal of Clinical Oncology* 30 (26): 3258–63.
<https://doi.org/10.1200/JCO.2011.39.8420>.
