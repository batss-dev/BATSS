# Prior distributions

``` r
library(BATSS)
```

## Introduction

Laplace approximation provides an analytical method to estimate the
posterior distribution of parameters in complex models by assuming that
the posterior distributions can be approximated by Gaussian densities
centered at their mode (or maximum). The calculation of integrals
related to the posteriors are therefore greatly simplified, typically
leading to computational gain, especially for complex models.

R has several implementations of Laplace approximations for Bayesian
inference, such as in the
[TMB](https://cran.r-project.org/web/packages/TMB/) (Template Model
Builder),
[LaplacesDemon](https://cran.r-project.org/web/packages/LaplacesDemon/)
and [mgcv](https://cran.r-project.org/web/packages/mgcv/) packages, for
example.

`BATSS` builds on integrated nested Laplace approximation (INLA) and
R-INLA, its implementation in R. INLA methodology focuses on models that
can be expressed as latent Gaussian Markov random fields (GMRF).
Consider the model

$$\eta_{i} = g\left( \mu_{i} \right) = \alpha + \sum\limits_{j = 1}^{p}\beta_{j}z_{ij} + \sum\limits_{k = 1}^{q}f_{k}\left( u_{ik} \right) + \epsilon_{i}$$

where

- $\eta$ is the linear predictor,
- $g(.)$ is the link function,
- $\beta$ are the $p$ coefficients for the linear effects of covariates
  $z$,
- $f_{k}(.)$ are the $q$ unknown functions of the covariates $u$  
- $\epsilon$ are unstructured terms.

Latent Gaussian models are the subset of Bayesian additive models that
assign a Gaussian prior to the latent field
$\mathcal{X} = \{\alpha,\beta,f(.)\}$

$$\mathcal{X}|\theta \sim N\left( 0,\mathcal{Q}^{- 1}(\theta) \right)$$

where $\mathcal{Q}$ denotes the precision matrix and $\theta$ a vector
of hyperparameters. Then the joint density of the latent field,
hyperparameters and the data is

$$\pi\left( \mathcal{X},\theta|y \right) \propto \pi(\theta)\pi\left( \mathcal{X}|\theta \right)\prod\limits_{i = 1}^{n}\pi\left( y_{i}|(A\mathcal{X})_{i},\theta \right)$$

where $\eta = A\mathcal{X}$. From here the posterior mariginals can be
approximated. Refer to Rue, Martino, and Chopin (2009) and Van Niekerk
et al. (2023) for details.

## Priors for fixed effects

Hence, priors for fixed effects in INLA are normally distributed (with a
default mean $\mu = 0$ and precision $\tau = 0.001$, and since
$\sigma^{2} = 1/\tau$, $N(0,\, 1000)$).

The `control.fixed` argument of the `inla()` function allows users to
change $\mu$ or $\tau$. The
[`batss.glm()`](../../reference/batss.glm.md) function allows control
elements of the `inla()` function to be passed on to the embedded INLA
routine via the `...` argument. `inla()` differentiates between the
intercept, `mean.intercept` and `prec.intercept`, and all other fixed
effects, `mean` and `prec`, e.g. `control.fixed = list(mean=1)` would
assign a prior mean of $1$ to all fixed effects but the intercept.

Individual prior means or precision can be set with named lists in the
`control.fixed` argument,
e.g. `control.fixed = list(prec=list(A=1,B=2,default=0.1))` would set
the prior precision to $1$ for fixed effect A, to $2$ for fixed effect B
and to $0.1$ for all other fixed effects.

Take note that in `BATSS` nomenclature of fixed effects is a combination
of names in the right hand side of the `model` argument and factor
levels from the allocation function in the `var` argument (which have to
match the names of the vector `prob0`) as is standard for the
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) function -
the first level of a factor being ignored for naming purposes as
reference level, e.g. 

``` r
y <- rnorm(30)
treatment <- factor(rep(c("A","B","C"),10))
levels(treatment)
```

    ## [1] "A" "B" "C"

``` r
model <- y ~ treatment
mm <- model.matrix(model)
colnames(mm)
```

    ## [1] "(Intercept)" "treatmentB"  "treatmentC"

The following example shows a
[`batss.glm()`](../../reference/batss.glm.md) call without and with
changes to the default prior for the fixed effects (Note: For all
examples, we are only running `R = 10` simulated trials each for
demonstration purposes.):

``` r
efficacy.fun = function(posterior,efficacy.bound){
  posterior > efficacy.bound
  }
group.fun = function(m,prob){
  prob = abs(prob)/sum(abs(prob))
  m0.g = floor(prob*m)
  m0   = sum(m0.g)
  factor(rep(names(prob),m0.g+rmultinom(1,m-m0,prob)),
         levels=names(prob))
  }
 
#batss.glm() using INLA's default normal priors, N(0,1000)
sim1 = batss.glm(  
          model             = y ~ trt,   
          var               = list(y   = rnorm,
                                   trt = group.fun),
          var.control       = list(y = list(sd = 2)),
          beta              = c(1, 1, 2),
          which             = c(2:3),
          R                 = 10,
          N                 = 200,
          interim           = list(recruited=100),
          prob0             = c(C=1/3, T1=1/3, T2=1/3),
          eff.arm           = efficacy.fun,
          eff.arm.control   = list(efficacy.bound = 0.975),
          delta.eff         = c(NA,0),
          fut.arm           = NULL,
          computation       = "parallel",
          H0                = FALSE,
          mc.cores          = parallel::detectCores()-1)
 
t(apply(sim1$H1$estimate[,3,],1,summary))

#batss.glm() using normal priors with N(1,4) for the fixed effect of intervention T1
#and the intercept and N(0,4) for intervention T2 and all other fixed effects.
sim2 = batss.glm(  
          model             = y ~ trt,   
          var               = list(y   = rnorm,
                                   trt = group.fun),
          var.control       = list(y = list(sd = 2)),
          beta              = c(1, 1, 2),
          which             = c(2:3),
          R                 = 10,
          N                 = 200,
          interim           = list(recruited=100),
          prob0             = c(C=1/3, T1=1/3, T2=1/3),
          eff.arm           = efficacy.fun,
          eff.arm.control   = list(efficacy.bound = 0.975),
          delta.eff         = c(NA,0),
          fut.arm           = NULL,
          computation       = "parallel",
          H0                = FALSE,
          mc.cores          = parallel::detectCores()-1,
          control.fixed     = list(mean.intercept=1,
                                   prec.intercept=0.25,
                                   mean = list(trtT1 = 1, default = 0), 
                                   prec = 0.25)
          ) 
 
t(apply(sim2$H1$estimate[,3,],1,summary))
```

And the same for a [`batss.surv()`](../../reference/batss.surv.md) call:

``` r
#batss.surv() using INLA's default normal priors, N(0,1000)
sim1.surv = batss.surv(
  model           = inla.surv(time, status) ~ trt,
  family          = "coxph",
  surv            = simsurv,
  surv.control    = list(lambdas = 2, gammas = 2,
                         maxt = 3),
  var             = list(trt = alloc.balanced),
  accr            = rexp,
  accr.control    = list(rate = 300),
  accr.type       = "random",
  prob0           = c(control = 1, A = 1, B = 1, C = 1),
  hr              = c(0.7, 1, 1),
  which           = 1:3,
  alternative     = "less",
  interim         = list(time = 1:3),
  eff.arm         = eff.arm.simple,
  eff.arm.control = list(b = c(0.98)),
  eff.trial       = eff.trial.any,
  delta.eff       = 0,
  fut.arm         = fut.arm.simple,
  fut.arm.control = list(b = 0.05),
  N               = 600,
  R               = 10,
  H0              = FALSE,
  extended        = 1,
  computation     = "parallel",
  mc.cores        = parallel::detectCores() - 1
)

exp(t(apply(sim1.surv$H1$estimate[,3,],1,summary)))

#batss.surv() using normal priors with N(1,4) for the fixed effect of intervention A
#and N(0,4) for the intercept and all other fixed effects.
sim2.surv = batss.surv(
  model           = inla.surv(time, status) ~ trt,
  family          = "coxph",
  surv            = simsurv,
  surv.control    = list(lambdas = 2, gammas = 2,
                         maxt = 3),
  var             = list(trt = alloc.balanced),
  accr            = rexp,
  accr.control    = list(rate = 300),
  accr.type       = "random",
  prob0           = c(control = 1, A = 1, B = 1, C = 1),
  hr              = c(0.7, 1, 1),
  which           = 1:3,
  alternative     = "less",
  interim         = list(time = 1:3),
  eff.arm         = eff.arm.simple,
  eff.arm.control = list(b = c(0.98)),
  eff.trial       = eff.trial.any,
  delta.eff       = 0,
  fut.arm         = fut.arm.simple,
  fut.arm.control = list(b = 0.05),
  N               = 600,
  R               = 10,
  H0              = FALSE,
  extended        = 1,
  computation     = "parallel",
  mc.cores        = parallel::detectCores() - 1,
  control.fixed   = list(mean.intercept=0,
                         prec.intercept=0.25,
                         mean = list(treatmentA = 1, default = 0), 
                         prec = 0.25)
)

exp(t(apply(sim2.surv$H1$estimate[,3,],1,summary)))
```

## Priors for hyperparameters

| **Distribution**  | **Hyperparameter** | **Hyperparameter default prior** |
|-------------------|--------------------|----------------------------------|
| normal            | log precision      | log-gamma                        |
| binomial          | none               | \-                               |
| negative binomial | log size           | penalised complexity gamma       |
| poisson           | none               | \-                               |

INLA hyperparameters and their default priors currently available in
`batss.glm`

| **Distribution** | **Hyperparameter** | **Hyperparameter default prior** |
|------------------|--------------------|----------------------------------|
| exponentialsurv  | none               | \-                               |
| weibullsurv      | log alpha          | penalised complexity alpha       |
| gompertzsurv     | shape              | normal                           |
| loglogisticsurv  | log alpha          | log-gamma                        |
| lognormalsurv    | log precision      | log-gamma                        |
| coxph            | log precision for  | log-gamma                        |
|                  | the piecewise      |                                  |
|                  | constant hazard    |                                  |

INLA hyperparameters and their default priors currently available in
`batss.surv`

Based on INLA’s defaults, family, hyperparameter and the associated
default prior used in `BATSS` are listed in Table 1 and 2 (see Simpson
et al. (2017) for an in depth description of complexity priors).
Detailed information can be viewed by typing
`inla.models()$likelihood$[FAMILYNAME]` including default parameters,
e.g. for the normal distribution:

``` r
INLA::inla.models()$likelihood$gaussian
```

    ## $doc
    ## [1] "The Gaussian likelihoood"
    ## 
    ## $hyper
    ## $hyper$theta1
    ## $hyper$theta1$hyperid
    ## [1] 65001
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$name
    ## [1] "log precision"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$short.name
    ## [1] "prec"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$output.name
    ## [1] "Precision for the Gaussian observations"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$output.name.intern
    ## [1] "Log precision for the Gaussian observations"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$initial
    ## [1] 4
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$fixed
    ## [1] FALSE
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$prior
    ## [1] "loggamma"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$param
    ## [1] 1e+00 5e-05
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta1$to.theta
    ## function (x) 
    ## log(x)
    ## <bytecode: 0x10e120550>
    ## <environment: 0x10e11a9b0>
    ## attr(,"inla.read.only")
    ## [1] TRUE
    ## 
    ## $hyper$theta1$from.theta
    ## function (x) 
    ## exp(x)
    ## <bytecode: 0x10e120668>
    ## <environment: 0x10e11a9b0>
    ## attr(,"inla.read.only")
    ## [1] TRUE
    ## 
    ## 
    ## $hyper$theta2
    ## $hyper$theta2$hyperid
    ## [1] 65002
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$name
    ## [1] "log precision offset"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$short.name
    ## [1] "precoffset"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$output.name
    ## [1] "NOT IN USE"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$output.name.intern
    ## [1] "NOT IN USE"
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$initial
    ## [1] 72.08731
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$fixed
    ## [1] TRUE
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$prior
    ## [1] "none"
    ## attr(,"inla.read.only")
    ## [1] TRUE
    ## 
    ## $hyper$theta2$param
    ## numeric(0)
    ## attr(,"inla.read.only")
    ## [1] FALSE
    ## 
    ## $hyper$theta2$to.theta
    ## function (x) 
    ## log(x)
    ## <bytecode: 0x10e120550>
    ## <environment: 0x10e11a9b0>
    ## attr(,"inla.read.only")
    ## [1] TRUE
    ## 
    ## $hyper$theta2$from.theta
    ## function (x) 
    ## exp(x)
    ## <bytecode: 0x10e120668>
    ## <environment: 0x10e11a9b0>
    ## attr(,"inla.read.only")
    ## [1] TRUE
    ## 
    ## 
    ## 
    ## $survival
    ## [1] FALSE
    ## 
    ## $discrete
    ## [1] FALSE
    ## 
    ## $link
    ## [1] "default"   "identity"  "logit"     "loga"      "cauchit"   "log"      
    ## [7] "logoffset"
    ## 
    ## $pdf
    ## [1] "gaussian"

In case of the normal distribution the default prior for the
hyperparameter $\theta = log(\tau)$ is a log-gamma distribution with
parameters $a = 1$ and $b = 0.00005$. To change the prior distribution
or its parameters use the `control.family` argument of the `inla()`
function.

``` r
#batss.glm() using normal priors with N(1,4) for the fixed effect of intervention T1  
#and the intercept and N(0,4) for intervention T2 and all other fixed effects and a 
#log-gamma(0.5,0.01) prior for the log precision hyperparameter.
sim3 = batss.glm(  
          model             = y ~ trt,   
          var               = list(y   = rnorm,
                                   trt = group.fun),
          var.control       = list(y = list(sd = 2)),
          beta              = c(1, 1, 2),
          which             = c(2:3),
          R                 = 10,
          N                 = 200,
          interim           = list(recruited=100),
          prob0             = c(C=1/3, T1=1/3, T2=1/3),
          eff.arm           = efficacy.fun,
          eff.arm.control   = list(efficacy.bound = 0.975),
          delta.eff         = c(NA,0),
          fut.arm           = NULL,
          computation       = "parallel",
          H0                = FALSE,
          mc.cores          = parallel::detectCores()-1,
          control.fixed     = list(mean.intercept=1,
                                   prec.intercept=0.25,
                                   mean = list(trtT1 = 1, default = 0), 
                                   prec = 0.25),
          control.family    = list(hyper = list(prec = list(param = c(0.5,0.01))))
          ) 
 
t(apply(sim3$H1$estimate[,3,],1,summary))
```

And for a [`batss.surv()`](../../reference/batss.surv.md) call:

``` r
#batss.surv() using normal priors with N(1,4) for the fixed effect of intervention A
#and N(0,4) for the intercept and all other fixed effects. The hyperparameter for the
#Cox PH model refers to the baseline hazard, in this case a random walk of order 1 with 
#20 intervals and a log-gamma(0.001,0.001) prior for the log precision therein. 
sim3.surv = batss.surv(
  model           = inla.surv(time, status) ~ trt,
  family          = "coxph",
  surv            = simsurv,
  surv.control    = list(lambdas = 2, gammas = 2,
                         maxt = 3),
  var             = list(trt = alloc.balanced),
  accr            = rexp,
  accr.control    = list(rate = 300),
  accr.type       = "random",
  prob0           = c(control = 1, A = 1, B = 1, C = 1),
  hr              = c(0.7, 1, 1),
  which           = 1:3,
  alternative     = "less",
  interim         = list(time = 1:3),
  eff.arm         = eff.arm.simple,
  eff.arm.control = list(b = c(0.98)),
  eff.trial       = eff.trial.any,
  delta.eff       = 0,
  fut.arm         = fut.arm.simple,
  fut.arm.control = list(b = 0.05),
  N               = 600,
  R               = 10,
  H0              = FALSE,
  extended        = 1,
  computation     = "parallel",
  mc.cores        = parallel::detectCores() - 1,
  control.fixed   = list(mean.intercept=0,
                         prec.intercept=0.25,
                         mean = list(treatmentA = 1, default = 0), 
                         prec = 0.25),
  control.hazard  = list(model="rw1", 
                         n.intervals=20,
                         hyper = list(prec = list(param = c(0.001, 0.001))))
)

exp(t(apply(sim3.surv$H1$estimate[,3,],1,summary)))
```

Priors for hyperparameters can either be a built-in prior or a
user-defined function. Available predefined priors in INLA:

``` r
INLA::inla.list.models("prior")
```

    ## Section [prior]
    ##      betacorrelation               Beta prior for the correlation          
    ##      dirichlet                     Dirichlet prior                         
    ##      expression:                   A generic prior defined using expressions
    ##      flat                          A constant prior                        
    ##      gamma                         Gamma prior                             
    ##      gaussian                      Gaussian prior                          
    ##      invalid                       Void prior                              
    ##      jeffreystdf                   Jeffreys prior for the doc              
    ##      laplace                       Laplace prior                           
    ##      linksnintercept               Skew-normal-link intercept-prior        
    ##      logflat                       A constant prior for log(theta)         
    ##      loggamma                      Log-Gamma prior                         
    ##      logiflat                      A constant prior for log(1/theta)       
    ##      logitbeta                     Logit prior for a probability           
    ##      logtgaussian                  Truncated Gaussian prior                
    ##      logtnormal                    Truncated Normal prior                  
    ##      minuslogsqrtruncnormal        (obsolete)                              
    ##      mvnorm                        A multivariate Normal prior             
    ##      none                          No prior                                
    ##      normal                        Normal prior                            
    ##      pc                            Generic PC prior                        
    ##      pc.alphaw                     PC prior for alpha in Weibull           
    ##      pc.ar                         PC prior for the AR(p) model            
    ##      pc.cor0                       PC prior correlation, basemodel cor=0   
    ##      pc.cor1                       PC prior correlation, basemodel cor=1   
    ##      pc.dof                        PC prior for log(dof-2)                 
    ##      pc.egptail                    PC prior for the tail in the EGP likelihood
    ##      pc.fgnh                       PC prior for the Hurst parameter in FGN 
    ##      pc.gamma                      PC prior for a Gamma parameter          
    ##      pc.gammacount                 PC prior for the GammaCount likelihood  
    ##      pc.gevtail                    PC prior for the tail in the GEV likelihood
    ##      pc.matern                     PC prior for the Matern SPDE            
    ##      pc.mgamma                     PC prior for a Gamma parameter          
    ##      pc.prec                       PC prior for log(precision)             
    ##      pc.prw2.range                 PCprior for the range in PRW2           
    ##      pc.range                      PC prior for the range in the Matern SPDE
    ##      pc.sn                         PC prior for the skew-normal            
    ##      pc.spde.GA                    (experimental)                          
    ##      pom                           #classes-dependent prior for the POM model
    ##      ref.ar                        Reference prior for the AR(p) model, p<=3
    ##      rprior:                       A R-function defining the prior         
    ##      table:                        A generic tabulated prior               
    ##      wishart1d                     Wishart prior dim=1                     
    ##      wishart2d                     Wishart prior dim=2                     
    ##      wishart3d                     Wishart prior dim=3                     
    ##      wishart4d                     Wishart prior dim=4                     
    ##      wishart5d                     Wishart prior dim=5                     
    ##      wishartkd                     Wishart prior

To change a prior to one of the predefined INLA priors simply use the
`control.family` argument of the `inla()` function. For example,
`control.family = list(hyper = list(prec = list(prior = "logtnormal")))`
will change the prior distribution of the precision to a truncated
normal prior. (Note that the name will depend on the hyperparameter,
here `prec`.)

``` r
#batss.glm() with default normal priors for the fixed effects and a truncated normal 
#prior for the precision hyperparameter
sim4 = batss.glm(  
          model             = y ~ trt,   
          var               = list(y   = rnorm,
                                   trt = group.fun),
          var.control       = list(y = list(sd = 2)),
          beta              = c(1, 1, 2),
          which             = c(2:3),
          R                 = 10,
          N                 = 200,
          interim           = list(recruited=100),
          prob0             = c(C=1/3, T1=1/3, T2=1/3),
          eff.arm           = efficacy.fun,
          eff.arm.control   = list(efficacy.bound = 0.975),
          delta.eff         = c(NA,0),
          fut.arm           = NULL,
          computation       = "parallel",
          H0                = FALSE,
          mc.cores          = parallel::detectCores()-1,
          control.family    = list(hyper = list(prec = list(prior = "logtnormal")))
          ) 
 
t(apply(sim4$H1$estimate[,3,],1,summary))
```

There are two options to define priors that are not predefined in
INLA: 1. Providing a tabulated distribution via `table:`

``` r
#manually define the log-gamma prior
# either as a 'table'
loggamma.function <-  function(log_precision) {
  a <- 0.05;
  b <- 0.01;
  precision <- exp(log_precision);
  logdens <- log(b^a) - lgamma(a) + (a-1)*log_precision - b*precision;
  log_jacobian <- log_precision;
  return(logdens + log_jacobian)
}
lprec <- seq(-10, 10, len=10000)
loggamma.tab <- paste(c("table:", cbind(lprec, loggamma.function(lprec))), 
                      sep = "", collapse = " ")
loggamma.table <- list(prec = list(prior = loggamma.tab))

#batss.glm() using normal priors with N(1,4) for the fixed effect of intervention T1 
#and the intercept and N(0,4) for intervention T2 and all other fixed effects and a 
#tabulated log-gamma(0.5,0.01) prior for the log precision hyperparameter.
sim5 = batss.glm(  
          model             = y ~ trt,   
          var               = list(y   = rnorm,
                                   trt = group.fun),
          var.control       = list(y = list(sd = 2)),
          beta              = c(1, 1, 2),
          which             = c(2:3),
          R                 = 10,
          N                 = 200,
          interim           = list(recruited=100),
          prob0             = c(C=1/3, T1=1/3, T2=1/3),
          eff.arm           = efficacy.fun,
          eff.arm.control   = list(efficacy.bound = 0.975),
          delta.eff         = c(NA,0),
          fut.arm           = NULL,
          computation       = "parallel",
          H0                = FALSE,
          mc.cores          = parallel::detectCores()-1,
          control.fixed     = list(mean.intercept=1,
                                   prec.intercept=0.25,
                                   mean = list(trtT1 = 1, default = 0), 
                                   prec = 0.25),
          control.family    = list(hyper = loggamma.table)
          ) 
 
t(apply(sim5$H1$estimate[,3,],1,summary))
```

or 2. providing the function via the `expression:` statement.

``` r
#or 'expression'
loggamma = "expression:
              a = 0.5;
              b = 0.01;
              precision = exp(log_precision);
              logdens = log(b^a) - lgamma(a) + (a-1)*log_precision - b*precision;
              log_jacobian = log_precision;
              return(logdens + log_jacobian);"
loggamma.expression = list(prec = list(prior = loggamma))

#batss.glm() using normal priors with N(1,4) for the fixed effect of intervention T1 
#and the intercept and N(0,4) for intervention T2 and all other fixed effects and a 
#user-defined log-gamma(0.5,0.01) prior for the log precision hyperparameter.
sim6 = batss.glm(  
          model             = y ~ trt,   
          var               = list(y   = rnorm,
                                   trt = group.fun),
          var.control       = list(y = list(sd = 2)),
          beta              = c(1, 1, 2),
          which             = c(2:3),
          R                 = 10,
          N                 = 200,
          interim           = list(recruited=100),
          prob0             = c(C=1/3, T1=1/3, T2=1/3),
          eff.arm           = efficacy.fun,
          eff.arm.control   = list(efficacy.bound = 0.975),
          delta.eff         = c(NA,0),
          fut.arm           = NULL,
          computation       = "parallel",
          H0                = FALSE,
          mc.cores          = parallel::detectCores()-1,
          control.fixed     = list(mean.intercept=1,
                                   prec.intercept=0.25,
                                   mean = list(trtT1 = 1, default = 0), 
                                   prec = 0.25),
          control.family    = list(hyper = loggamma.expression)
          ) 
 
t(apply(sim6$H1$estimate[,3,],1,summary))

#check if it matches the predefined loggamma prior
t(apply(sim3$H1$estimate[,3,],1,summary))
```

## Further reading

Further details (including parametrisation of distributions) and
examples for use of priors in `INLA` (and by extension in `BATSS`) can
be found in the INLA documentation `inla.doc("[TOPIC]")` or online
(<https://www.r-inla.org/>).

## References

Rue, Håvard, Sara Martino, and Nicolas Chopin. 2009. “Approximate
Bayesian Inference for Latent Gaussian Models by Using Integrated Nested
Laplace Approximations.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 71 (2): 319–92.
<https://doi.org/10.1111/j.1467-9868.2008.00700.x>.

Simpson, Daniel, Håvard Rue, Andrea Riebler, Thiago G. Martins, and
Sigrunn H. Sørbye. 2017. “Penalising Model Component Complexity: A
Principled, Practical Approach to Constructing Priors.” *Statistical
Science* 32 (1): 1–28. <https://doi.org/10.1214/16-STS576>.

Van Niekerk, Janet, Elias Krainski, Denis Rustand, and Håvard Rue. 2023.
“A New Avenue for Bayesian Inference with INLA.” *Computational
Statistics & Data Analysis* 181: 107692.
<https://doi.org/10.1016/j.csda.2023.107692>.
