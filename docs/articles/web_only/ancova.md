# ANCOVA

## Context

This example is motivated by the EMPACT trial - a phase IIb randomised
controlled clinical trial exploring the efficacy and safety of
psilocybin-assisted psychotherapy (PAP) in the treatment of people
diagnosed with treatment resistant major depressive disorder.

EMPACT is a Bayesian adaptive multi-arm multi-stage trial with
response-adaptive randomisation. It compares the average reversed
‘Hamilton Depression Rating Scale’ scores (R-HAM-D; where large R-HAM-D
correspond to better outcome) one month after the last PAP session for
patients allocated to an experimental group compared with the control
group. There are three experimental groups, each treated with different
Psilocybin doses; the control treatment consists of an active
comparator.

Here an ANCOVA model is considered to account for each patient’s
baseline R-HAM-D score, as this might lead to efficiency gains. Indeed,
the required sample size to detect any effect can be reduced by a factor
$(1 - \rho)$, where $\rho$ is the correlation between the baseline
R-HAM-D score and the R-HAM-D score at one month. For more details on
baseline covariate adjustment methods in RCTs, check, for example,
Shiyuan Zhang and Thabane (2014).

## Design

The chosen design has the following characteristics:

- **treatment arms**: Our study has 4 arms: an active control treatment
  (arm “Ctrl”) and 3 (experimental) groups treated with different
  Psilocybin doses (referred to as arms “D1”, “D2” and “D3”).
- **additional predictor:** to control for baseline R-HAM-D scores
  (continuous predictor). R-HAM-D scores at baseline are assumed to be
  associated with R-HAM-D score taken after one month, i.e., with the
  primary endpoint.
- **alternative hypotheses**: This is a superiority trial where we are
  interested in demonstrating efficacy for an experimental arm compared
  to the shared control group. Experimental arms will not be compared to
  one another.
- **interim analyses and maximum sample size**: We will begin the
  interim analyses at 50 participants completing follow-up, and perform
  interim analyses every 20 participants completing follow-up
  thereafter. Our maximum sample size is 130 participants with a total
  of 5 looks planned. We assume that recruitment is fairly slow in this
  trial.
- **endpoint conditional distribution**: \`R-HAM-D’ scores at one month
  are assumed to follow a Gaussian distribution.
- **group allocation**: the trial will start with equal randomisation
  (1:1:1:1) and RAR will be used from interim analysis 1 onwards. The
  randomisation probabilities are calculated using a variation of Trippa
  et al. (2012) and Cellamare et al. (2017) where the allocation
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
    - $\beta_{k}$ denotes the $k$th target parameter,
    - $\delta_{r} = 0$ denotes the (RAR-related) clinically meaningful
      treatment effect value,
    - and $h = \gamma\left( \sum_{k = 0}^{K}n_{kj}/N \right)^{\eta}$
      with
      - $n_{kj}$, the sample size of arm $k$ at stage $j$,
      - $N$, the maximum sample size,
      - $\gamma = 3$ and $\eta = 1.4$, two tuning parameters. Refer to
        Gotmaker et al. (2019) for details
- **arm efficacy stopping rule**: Early stopping of intervention arms
  for efficacy may occur if there is a high posterior probability of
  (any) benefit at look $j$, i.e., when
  $$\text{Prob}\left( \beta_{k} > \delta_{\epsilon}|y,X \right) > 1 - b_{\epsilon}\left( \sum\limits_{g = 0}^{G}n_{g}/N \right)^{p_{\epsilon}}$$
  where
  - $\delta_{\epsilon} = 0$ denotes the (efficacy-related) clinically
    meaningful treatment effect value,
  - $b_{\epsilon} = 0.0115$ and $p_{\epsilon} = 1.575$ are tuning
    parameters.
- **arm futility stopping rule**: Early stopping of intervention arms
  for futility may occur if there is a low posterior probability of
  observing a reduction of at least 3 in the mean R-HAM-D score (from
  baseline to 1 month), i.e., when
  $$\text{Prob}\left( \beta_{k} > \delta_{f}|y,X \right) < b_{f}$$ where
  - $\delta_{f} = 3$ denotes the (futility-related) clinically
    meaningful treatment effect values,
  - $b_{f} = 0.05$ is the cut-off value used to declare futility.
- **trial stopping rule:** The trial will run until an efficacy or
  futility decision has been reached for each intervention arm, or once
  the maximum sample size has been reached (separate stopping rules for
  each arm).

Note that, in the above, the tuning parameter values have been optimised
to lead to suitable operating characteristics, namely,

- a familywise error rate (FWER) of 0.05 under the global null
  hypothesis,
- a probability of efficacy per arm of 0.8 or higher when assuming a
  difference in mean of 5 compared to the control group,
- a probability of early futility stopping per arm of 0.7 under H0.

In the following, we define the group allocation, RAR, futility and
efficacy functions corresponding to the design described above.

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
trippa.fun = function(posterior,n,N,ref,active,gamma,eta,nu){
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
trippa.fun(c(.5,.5,.5,.5,.5),n=c(A=10,B=10,C=10,D=10,E=10,F=10),N=16,
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
# test on m = 50 patients and equal allocation per group 
table(treatalloc.fun(m=50,prob=c(Ctrl=1, D1=1, D2=1, D3=1)))
```

`treatalloc.fun` first allocates the largest possible number of units to
the different groups given their exact target probabilities and then
assigns randomly the remaining units to the different groups according
to multinomial draws.

Note that this corresponds to the function `alloc.balanced` available in
the **BATSS** package.

#### Arm efficacy stopping rule

We need to generate a function that leads to a logical output (`TRUE` if
the arm needs to be stopped for efficacy and `FALSE` otherwise) and
takes as input

- the ingredient `posterior` for the posterior probability of the target
  parameter being greater than `delta.eff = 0`,
- `n` and `N`, respectively the sample size per arm at the look of
  interest and the max sample size to define the information fraction,
- the additional parameters $b_{\epsilon}$ and $p_{\epsilon}$ that we
  will respectively name `b.eff` and `p.eff` (and that need to be added
  to `eff.arm.control` in `batss.glm`).

``` r
# function
efficacy.arm.fun = function(posterior,n,N,b.eff,p.eff){
  posterior > (1-(b.eff*(sum(n)/N)^p.eff))
}
# test for a parameter with a posterior = 0.999
efficacy.arm.fun(0.999, n=100, N=1000, b.eff = 0.045, p.eff=1.4)
# test for a parameter with a posterior = 0.95
efficacy.arm.fun(0.95, n=100, N=1000, b.eff = 0.045, p.eff=1.4)
```

Note that this corresponds to the function `eff.arm.infofract` available
in the **BATSS** package.

#### Arm futility stopping rule

We need to generate a function that leads to a logical output (`TRUE` if
the arm needs to be stopped for futility and `FALSE` otherwise) and
takes as input

- the ingredient `posterior` for the posterior probability of the target
  parameter being greater than `delta.fut = 3`,
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

#### Trial efficacy and futility stopping rule

We need to generate functions that, based on the ingredients
`eff.target` and `fut.target` (respectively indicating if efficacy and
futility were declared for each target parameter at the stage of
interest or before) leads to a logical output equal to `TRUE` if all
target parameters were respectively declared efficacious or futile, and
`FALSE` otherwise.

``` r
# function
efficacy.trial.fun = function(eff.target){
  all(eff.target)
}

# function
futility.trial.fun = function(fut.target){
  all(fut.target)
}

# test 
futility.trial.fun(c(D1=TRUE,D2=TRUE,D3=TRUE))
efficacy.trial.fun(c(D1=TRUE,D2=TRUE,D3=FALSE))
```

Note that these correspond to the functions `eff.trial.all` and
`fut.trial.all` available in the **BATSS** package which are used by
default and that this behaviour is the default one of the function
`batss.glm` when `fut.trial = NULL`.

## Monte Carlo Simulations

We consider three scenarios: a first one that does not control for
‘baseline R-HAM-D scores’ and two that do and assume different
correlation levels between R-HAM-D scores at baseline (predictor 4) and
after one month (primary endpoint). Thus,

- the primary predictor, $x_{1}$, is the 4-level treatment group factor
  with the control group as reference group,
- the secondary predictor, $x_{2}$, is the centered ‘R-HAM-D scores’ at
  baseline.

For the $i$th participant, our assumed model is:

$$y_{i} = \beta_{0} + \beta_{1}\iota\left( x_{i1} = \prime D1\prime \right) + \beta_{2}\iota\left( x_{i1} = \prime D2\prime \right) + \beta_{3}\iota\left( x_{i1} = \prime D3\prime \right) + \beta_{4}\left\lbrack x_{i2} - {\bar{x}}_{2} \right\rbrack + \epsilon_{i}$$,

where

- $x_{ip}$ denotes the value of the predictor $p$ for participant $i$,
- $\iota(.)$ denotes the indicator function that equals to 1 if the
  condition is satisfied and 0 otherwise,
- ${\bar{x}}_{2}$ denotes the average ‘R-HAM-D scores’ at baseline.

This parameterisation leads to the following parameter interpretation:

- $\beta_{0}$ corresponds to the average ‘R-HAM-D scores’ of patients at
  baseline, as well as after one month of treatment for patients of the
  “Ctrl” arm,
- $\beta_{1}$, $\beta_{2}$ and $\beta_{3}$ respectively correspond to
  the average shift in mean between the R-HAM-D scores at one month for
  each treatment arm (“D1”, “D2” and “D3”) compared to the control one
  (“Ctrl”),
- $\beta_{4}$ is the regression parameter related to participant
  ‘R-HAM-D’ scores at baseline.

We assume

- a mean of 5 in the control arm and then a mean change from control of
  5 in each intervention so that
  $\left\lbrack \beta_{0},\beta_{1},\beta_{2},\beta_{3} \right\rbrack^{T} = \lbrack 5,5,5,5\rbrack^{T}$  
- a standard deviation of the response *given* the 4-level ‘treatment’
  predictor, equal to 7, i.e., $\sigma_{y} = 7$.

When controlling for the (centered) ‘R-HAM-D scores’ at baseline, we
further need to define

- the value of $\beta_{4}$,
- the values of $\sigma_{x_{4}}$ and $\sigma_{\epsilon}$, the standard
  deviations of ‘baseline R-HAM-D scores’ and of the error term, so that
  $\rho$, the Pearson’s correlation level between the centered ‘baseline
  R-HAM-D scores’ and ‘R-HAM-D scores at one month’ *given* the 4-level
  ‘treatment’ predictor, equals the target correlation level. This is
  achieved, in the ANCOVA setting, by choosing values respecting the
  following equalities:
- $\rho = \beta_{4}\frac{\sigma_{x_{4}}}{\sigma_{y}}$,
- $\sigma_{y}^{2} = \beta_{4}^{2}\sigma_{x_{4}}^{2} + \sigma_{\epsilon}^{2}$,
  where the former comes from the relationship between (simple) linear
  regression and Pearson’s correlation coefficient and the later comes
  from the law of total variance.

### Adaptive designs

We first show the code corresponding to the planned adaptive design for
each Scenario:

#### Scenario 0

In Scenario 0, we don’t correct for baseline R-HAM-D scores.

``` r
library(BATSS)
# number of trials
R = 100


# simulation
scenario0 = batss.glm(   
  model           = y~group,
  var             = list(y = rnorm,
                         group = treatalloc.fun),
  var.control     = list(y = list(sd = 7)), 
  family          = "gaussian",
  link            = "identity",
  beta            = c(5,5,5,5),
  which           = c(2:4),
  R               = R,
  alternative     = c("greater"),
  RAR             = trippa.fun,
  RAR.control     = list(gamma=3,eta=1.4,nu=0.1),
  delta.RAR       = 0,
  prob0           = c(Ctrl=1,D1=1,D2=1,D3=1),
  N               = 130, 
  interim         = list(recruited=list(m0=50, m=20)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = 0, 
  eff.arm.control = list(b.eff = 0.0115, p.eff=1.575),
  delta.fut       = 3, 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.05),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = TRUE,
  extended        = 1)   
```

You can note that

- `beta`, the regression parameter vector equals
  $\lbrack 5,5,5,5\rbrack^{T}$,

- `y` is generated via the function `rnorm` as well and has an average
  of respectively

  - `5` for the reference group,
  - `5+5 = 10` for all other groups (as contrasts of type treatment are
    used for the factor `treatment` in the self-defined function
    `treatalloc.fun`),

- the target parameters (average shift in mean between the R-HAM-D
  scores at one month for each treatment arm (“D1”, “D2” and “D3”)
  compared to the control arm) are in position 2 to 4 of the fitted
  coefficients obtained when using the formula `y~group`,

- `prob0` provides

  - the (equal) allocation probabilities at the start of the trial,
  - the names of the groups,

- `eff.arm` and `fut.arm` are set to the functions defined above (i.e.,
  `efficacy.arm.fun` and `futility.arm.fun`),

- `fut.trial` is not specified and therefore equal to `NULL` (default)
  which leads to the behaviour wished in this case,

- the values of the additional parameters of `efficacy.arm.fun` and
  `futility.arm.fun` (i.e., `b.eff` and `b.fut`) are specified under
  `eff.arm.control` and `fut.arm.control`,

- `delta.eff` and `delta.fut` are respectively set to `0` and `3`.

We chose here a low number of seeds/trials (`R=100`) to save time. Note
that, as `H0=TRUE`, the results under the global null are also defined.

#### Scenario 1

In Scenario 1, we assume a correlation of 0 between R-HAM-D scores at
baseline and after one month so that

- $\beta_{4} = 0$,
- $\sigma_{\epsilon} = 7$. This corresponds to a very unlikely scenario
  allowing us to see the drop in power induced by adding to the model a
  predictor unrelated to the response.

``` r
library(BATSS)
# number of trials
R = 100


# simulation
scenario1 = batss.glm(   
  model           = y~group+baseline,
  var             = list(y = rnorm,
                         group = treatalloc.fun,
                         baseline = rnorm),
  var.control     = list(y = list(sd = 7)), 
  family          = "gaussian",
  link            = "identity",
  beta            = c(5,5,5,5,0),
  which           = c(2:4),
  R               = R,
  alternative     = c("greater"),
  RAR             = trippa.fun,
  RAR.control     = list(gamma=3,eta=1.4,nu=0.1),
  delta.RAR       = 0,
  prob0           = c(Ctrl=1,D1=1,D2=1,D3=1),
  N               = 130, 
  interim         = list(recruited=list(m0=50, m=20)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = 0, 
  eff.arm.control = list(b.eff = 0.0115, p.eff=1.575),
  delta.fut       = 3, 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.05),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = TRUE,
  extended        = 1)   
```

Compared to Scenario 0, you can note that

- `y`, the outcome vector, now depends on the R-HAM-D scores at
  baseline,
- `beta`, the regression parameter vector equals
  $\lbrack 5,5,5,5,0\rbrack^{T}$,
- `baseline` is generated with the function `rnorm` (without further
  arguments in `var.control` so that a standard normal distribution is
  assumed, which is irrelevant as the corresponding beta coefficient is
  equal to 0).

#### Scenario 2

In Scenario 2, we assume that $\rho = 0.6$ (a value mentioned in the
EMPACT protocol), $\sigma_{x_{4}} = 3.5$ while still requiring that
$\sigma_{\epsilon} = 7$, so that

$$\beta_{4} = \rho\frac{\sigma_{y}}{\sigma_{x_{4}}} = 0.6 \times \frac{7}{3.5} = 1.2,$$$$\sigma_{\epsilon}^{2} = \sigma_{y}^{2} - \beta_{4}^{2}\sigma_{x_{4}}^{2} = 7^{2} - 1.2^{2} \times 3.5^{2} = 5.6.$$

The following code first defines the values of $\beta_{4}$ and
$\sigma_{\epsilon}$, then checks these values by simulation and finally
defines the operating characteristics of the design.

``` r
##
## calculate the SIGMA.e and BETA values outside of the batss.glm function
## 

SIGMA.x = 3.5
SIGMA.y = 7
RHO     = 0.6
BETA    = c(5,5,5,5,RHO*SIGMA.y/SIGMA.x)
SIGMA.e = sqrt(SIGMA.y^2 - BETA[5]^2*SIGMA.x^2)
n       = 1e+5
X       = cbind(kronecker(matrix(1,nrow=n/4,ncol=1),cbind(1,diag(4)[,-1])),NA)

##
## check values by simulation
##

sim.fun = function(seed, Xmat, BETA, SIGMA.x,SIGMA.e){
    # seed=2
    set.seed(seed)
    Xmat[,5] = rnorm(nrow(Xmat),0,SIGMA.x)
    e = rnorm(nrow(Xmat),0,SIGMA.e)
    y = Xmat%*%BETA+e
    #
    fitw0 = lm.fit(Xmat[,1:4],y)
    fitw1 = lm.fit(Xmat,y)
    c(coef(fitw1),rho=cor(resid(fitw0),Xmat[,5]),sigma_ytot=sqrt(var(y)),
      sigma_y = mean(sqrt(tapply(y,rep(1:4,nrow(Xmat)/4),var))),
      sigma_e = sqrt(var(resid(fitw1))))
}

library(future.apply)
plan(multisession)
toto = t(future_sapply(1:1000, sim.fun, Xmat=X, BETA=BETA, 
         SIGMA.x=SIGMA.x,SIGMA.e=SIGMA.e, future.seed=TRUE))
boxplot(toto-matrix(rep(c(BETA,RHO,NA,SIGMA.y,SIGMA.e),nrow(toto)),nrow=nrow(toto),byrow=TRUE))
abline(h=0,col=c("blue"))


##
## define operating characteristics 
##

# number of trials
R = 100

# simulation
scenario2 = batss.glm(   
  model           = y~group+baseline,
  var             = list(y = rnorm,
                         group = treatalloc.fun,
                         baseline = rnorm),
  var.control     = list(y = list(sd = SIGMA.e), baseline = list(sd = SIGMA.x)), 
  family          = "gaussian",
  link            = "identity",
  beta            = BETA,
  which           = c(2:4),
  R               = R,
  alternative     = c("greater"),
  RAR             = trippa.fun,
  RAR.control     = list(gamma=3,eta=1.4,nu=0.1),
  delta.RAR       = 0,
  prob0           = c(Ctrl=1,D1=1,D2=1,D3=1),
  N               = 130, 
  interim         = list(recruited=list(m0=50, m=20)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = 0, 
  eff.arm.control = list(b.eff = 0.0115, p.eff=1.575),
  delta.fut       = 3, 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.05),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = TRUE,
  extended        = 1)  
```

Compared to ‘Scenario 1’, we can note the following differences:

- `beta`, now taking `1.2` as last value instead of `0`,
- `var.control`, now indicating
  - an amended value for the standard deviation of the error term,
  - the standard deviation of the variable `baseline`.

### Fixed designs

We now show how to use `batss.glm` to obtain the operating
characteristics of the corresponding fixed designs, considering the same
efficacy and futility parameters as well as the same maximum sample size
of $N = 130$ patients. The easiest way to fit a fixed design with
**BATSS** is to plan a single interim analysis defined so that no
adaptations are allowed. This is ensured by setting the `delta.eff` and
`delta.fut` values to `NA` for the interim, as well as by setting the
`RAR` argument to `NULL` (as treatment allocation after the interim
analysis would otherwise be based on posterior probabilities estimated
on interim data). Note that the sample size at the dummy interim should
ideally be a multiple of the number of arms in the trial.

#### Scenario 0 (fixed)

Same as the adaptive Scenario 0 above, except `interim` is reduced to a
single dummy interim at 52 participants, `delta.eff` and `delta.fut` are
set to `c(NA, 0)` and `c(NA, 3)` respectively so that no decisions are
made at the interim, and `RAR` is set to `NULL` (with `RAR.control` and
`delta.RAR` removed) so that equal group allocation is used throughout.

``` r
library(BATSS)
# number of trials
R = 100

# simulation
scenario0_fixed = batss.glm(   
  model           = y~group,
  var             = list(y = rnorm,
                         group = treatalloc.fun),
  var.control     = list(y = list(sd = 7)), 
  family          = "gaussian",
  link            = "identity",
  beta            = c(5,5,5,5),
  which           = c(2:4),
  R               = R,
  alternative     = c("greater"),
  RAR             = NULL,
  prob0           = c(Ctrl=1,D1=1,D2=1,D3=1),
  N               = 130, 
  interim         = list(recruited=c(52)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA, 0), 
  eff.arm.control = list(b.eff = 0.0115, p.eff=1.575),
  delta.fut       = c(NA, 3), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.05),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = TRUE,
  extended        = 1)   
```

#### Scenario 1 (fixed)

Same as the fixed Scenario 0 above, but now including `baseline` as in
the adaptive Scenario 1 (with $\beta_{4} = 0$).

``` r
# number of trials
R = 100

# simulation
scenario1_fixed = batss.glm(   
  model           = y~group+baseline,
  var             = list(y = rnorm,
                         group = treatalloc.fun,
                         baseline = rnorm),
  var.control     = list(y = list(sd = 7)), 
  family          = "gaussian",
  link            = "identity",
  beta            = c(5,5,5,5,0),
  which           = c(2:4),
  R               = R,
  alternative     = c("greater"),
  RAR             = NULL,
  prob0           = c(Ctrl=1,D1=1,D2=1,D3=1),
  N               = 130, 
  interim         = list(recruited=c(52)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA, 0), 
  eff.arm.control = list(b.eff = 0.0115, p.eff=1.575),
  delta.fut       = c(NA, 3), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.05),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = TRUE,
  extended        = 1)   
```

#### Scenario 2 (fixed)

Same as the fixed Scenario 1 above, but with `beta`, `var.control` and
`SIGMA.e`/`SIGMA.x` values as defined in the adaptive Scenario 2
($\rho = 0.6$).

``` r
# number of trials
R = 100

# simulation
scenario2_fixed = batss.glm(   
  model           = y~group+baseline,
  var             = list(y = rnorm,
                         group = treatalloc.fun,
                         baseline = rnorm),
  var.control     = list(y = list(sd = SIGMA.e), baseline = list(sd = SIGMA.x)), 
  family          = "gaussian",
  link            = "identity",
  beta            = BETA,
  which           = c(2:4),
  R               = R,
  alternative     = c("greater"),
  RAR             = NULL,
  prob0           = c(Ctrl=1,D1=1,D2=1,D3=1),
  N               = 130, 
  interim         = list(recruited=c(52)),
  eff.arm         = efficacy.arm.fun,
  delta.eff       = c(NA, 0), 
  eff.arm.control = list(b.eff = 0.0115, p.eff=1.575),
  delta.fut       = c(NA, 3), 
  fut.arm         = futility.arm.fun,
  fut.arm.control = list(b.fut = 0.05),
  computation     = "parallel",
  mc.cores        = 9,
  H0              = TRUE,
  extended        = 1)   
```

#### Scenario comparison

Based on a simulation considering 10,000 Monte Carlo trials,

- For *Scenario 0*, the (estimated) FWER under the global null
  hypothesis and power per arm under the alternative hypothesis
  respectively equal 0.0498 and 0.8011,
- For *Scenario 1*, now controlling for a further predictor not
  associated with the primary endpoint, the FWER slightly increases to
  0.0527 and the power per arm marginally decreases (0.7975),
- For *Scenario 2*, now controlling for a further predictor strongly
  associated with the primary endpoint, the FWER slightly increases to
  0.0550 and the power per arm jumps to 0.9424.

This suggests that adding a predictor has little downside as:  
\* it only marginally reduces the power when the endpoint and new
predictor are not associated, \* it strongly increases the power when
the endpoint and new predictor are associated, whilst maintaining very
reasonable FWER control.

## References

Cellamare, Matteo, Steffen Ventz, Elisabeth Baudin, Carole D. Mitnick,
and Lorenzo Trippa. 2017. “A Bayesian Response-Adaptive Trial in
Tuberculosis: The endTB Trial.” *Clinical Trials* 14 (1): 17–28.
<https://doi.org/10.1177/1740774516665090>.

Gotmaker, Robert, Michael J Barrington, John Reynolds, Lorenzo Trippa,
and Stephane Heritier. 2019. “Bayesian Adaptive Design: The Future for
Regional Anesthesia Trials?” *Regional Anesthesia & Pain Medicine* 44
(6): 617–22. <https://doi.org/10.1136/rapm-2018-100248>.

Shiyuan Zhang, Manyat Nantha-Aree, James Paul, and Lehana Thabane. 2014.
“Empirical Comparison of Four Baseline Covariate Adjustment Methods in
Analysis of Continuous Outcomes in Randomized Controlled Trials.” *Clin
Epidemiol.* 6: 227–35. <https://doi.org/10.1177/0193841X0102500101>.

Trippa, Lorenzo, Eudocia Q. Lee, Patrick Y. Wen, Tracy T. Batchelor,
Timothy Cloughesy, Giovanni Parmigiani, and Brian M. Alexander. 2012.
“Bayesian Adaptive Randomized Trial Design for Patients with Recurrent
Glioblastoma.” *Journal of Clinical Oncology* 30 (26): 3258–63.
<https://doi.org/10.1200/JCO.2011.39.8420>.
