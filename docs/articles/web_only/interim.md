# Definition of the interim schedule

## Introduction

Setting a viable and effective interim schedule in multi-arm,
multi-stage designs is paramount to an efficient trial design. With
`BATSS` version \>=2.0.0 and the introduction of a time component to
simulated trials, the scheduling of interim analyses has gained a lot of
flexibility. (Note: this feature is currently only available in the
[`batss.surv()`](../../reference/batss.surv.md) function but will be
extended to [`batss.glm()`](../../reference/batss.glm.md) in the near
future.)

## Available options for the scheduling of interim analyses in `BATSS`

At the present, `BATSS` allows the scheduling of interim analyses to
depend on the number of participants recruited, the number of events
observed and calendar time. In the case of event-based interim analyses,
the “number of events” to trigger an interim analysis can be the count
of events across all treatment arms (control arm included), the number
of events in the control group only, the minimum number of events across
all arms (control group included), or a combination of the number of
events in the control and treatment arm $i$ (number of events in a pair
of arms, where one of the pair is always the control). This allows for
the following interim analysis (timing) strategies based on:

- number of participants recruited (across the entire trial; i.e. across
  all treatment arms)
- time units
- total number of events observed
- events observed in the control group
- the last arm reaching a minimum number of events
- events observed in the control and treatment arm $i$ combined (number
  of events in a pair of arms)
- a specific number of events to be reached for the first interim
  (reference group options are limited to all, control and minimum),
  then after specified time units

In general, any interim strategy will be passed to `batss.surv` through
the `interim` argument, a named list item that can contain one or more
of the following named items:

- `recruited` (numeric vector)
- `time` (numeric vector)
- `event` (numeric vector)
- `event.type` (character string)

## Basic options

The length of numeric vectors provided in `interim` will set the total
number of interim analyses. An additional final analysis will be
performed after the last participant has finished follow-up (if no trial
ending criteria have been met beforehand). For example, setting
`interim = list(recruited = c(100, 200, 300))` will trigger an interim
analysis after the 100th, 200th and 300th participant has been recruited
into the trial (across all arms).

In the following, we will focus on trials that will only evaluate for
efficacy at the final analysis, so that the scheduling of interim
analyses can be seen without the trial being stopped early
(i.e. interims are performed but no decisions are made / no statistical
triggers are implemented).

### Number of participants recruited

Let’s look at this in an example and call the function `batss.surv` with
argument `interim = list(recruited = c(100, 200, 300)`. The maximum
sample size $N$ is set at 400. (Note: only 2 trials are run for
demonstration purposes of the interim specification).

``` r
sim.surv.recruited = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(recruited = c(100, 200, 300)),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

We can then look at the summary for the first trial to check the timing
of the interims:

``` r
sim.surv.recruited$H1$trial[[1]]$look
```

We can see that interim analyses have in fact been performed after
`n = c(100, 200, 300)` and a final analysis after the last patient
completed follow-up (at `n = 400`)

``` r
sim.surv.recruited$H1$trial[[1]]$look$n
```

at times `t`.

``` r
sim.surv.recruited$H1$trial[[1]]$look$t
```

### Time units

Changing the `interim` argument in the call to
[`batss.surv()`](../../reference/batss.surv.md) to
`interim = list(time = 1:3)` will result in the interim analysis being
instead performed after 1, 2 and 3 time units (here years), respectively
(and a final analysis again after the last participant completed
follow-up.)

``` r
sim.surv.time      = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(time = 1:3),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

Again looking at

``` r
sim.surv.time$H1$trial[[1]]$look
```

will show us that interim and final analysis have been performed after
time `t`

``` r
sim.surv.time$H1$trial[[1]]$look$t
```

and with `n` participants recruited, respectively.

``` r
sim.surv.time$H1$trial[[1]]$look$n
```

### Events

For interim analyses based on events, the `interim` argument has to be
changed to `interim = list(event = c(50, 100, 150))`. However, the
additional list item `event.type` will specify which (arms’) events will
be counted.

#### Overall events

To count all events occuring, we simply set `event.type = "all"`,
e.g. `interim = list(event = c(80, 160, 240), event.type = "all")` will
result in an interim analysis being conducted after a total of 80, 160
and 240 events in all arms combined.

``` r
sim.surv.event.all = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(event = c(80, 160, 240), event.type = "all"),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

The single trial summary

``` r
sim.surv.event.all$H1$trial[[1]]$look
```

will show us that interim and final analysis have been performed after
overall events `ev(n)`

``` r
sim.surv.event.all$H1$trial[[1]]$look$"ev(n)"
```

We can see that these have occurred at times `t`

``` r
sim.surv.event.all$H1$trial[[1]]$look$t
```

and with `n` participants recruited, respectively.

``` r
sim.surv.event.all$H1$trial[[1]]$look$n
```

#### Events in the control group

Setting `event.type = "control"`,
e.g. `interim = list(event = c(20, 40, 60), event.type = "control")`
will achieve interim analyses being performed at the occurence of 20,
40, and 60 events in the control group.

``` r
sim.surv.event.ctr = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(event = c(20, 40, 60), event.type = "control"),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

Again, the single trial summary

``` r
sim.surv.event.ctr$H1$trial[[1]]$look
```

will show us that interim and final analysis have been performed after
events in the control group, `ev(control)`. Please note that this call
will depend on the naming of the vector `prob0`.

``` r
sim.surv.event.ctr$H1$trial[[1]]$look$"ev(control)"
```

These have been reached times `t`

``` r
sim.surv.event.ctr$H1$trial[[1]]$look$t
```

and with `n` participants recruited, respectively.

``` r
sim.surv.event.ctr$H1$trial[[1]]$look$n
```

#### The last arm reaches event threshold (minimum number of events)

Setting `event.type = "min"`,
e.g. `interim = list(event = c(20, 40, 60), event.type = "min")`, will
achieve interim analyses being performed when the last arm reaches 20,
40, and 60 events, respectively. This will ensure that interim analyses
are only performed once each treatment arm has a certain number of
events. All arms will be analysed once the last arm has reached the
defined number of events

``` r
sim.surv.event.min = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(event = c(20, 40, 60), event.type = "min"),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

Again, the single trial summary

``` r
sim.surv.event.min$H1$trial[[1]]$look
```

will show us that interim and final analysis have been performed after
the last arm reaches the respective threshold, we can check by looking
at the minimum number of events (across arms) at each interim

``` r
apply(sim.surv.event.min$H1$trial[[1]]$look[,16:19], 1, min)
```

Again, we can check the times `t`

``` r
sim.surv.event.min$H1$trial[[1]]$look$t
```

and participants recruited `n` at each step.

``` r
sim.surv.event.min$H1$trial[[1]]$look$n
```

#### Events in control and treatment $i$

The last basic schedule is to perform an interim analysis every time an
intervention arm and control reach a pre-specified threshold value
(number of events in each pair of arms). This will result in multiple
interim analyses per threshold value (one for each active/remaining(!)
arm). This is a different approach than thos outlined above because only
treatment $i$ will be compared with the control and the exact number of
interim analysis will not be known at the beginning of the trial (only
the maximum number). To perform this scheduling the option
`event.type = "cplusone"` has to be chosen,
e.g. `interim = list(event = c(40, 80, 120), event.type = "cplusone")`
will perform a maximum of 9 interim analysis (one for each intervention
arm (plus control) at each threshold value). Again, please note that we
only enforce decision rules (for efficacy) at the final analysis and
therefore will reach all 9 interim analyses in the data.

``` r
sim.surv.event.cp1 = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(event = c(40, 80, 120), event.type = "cplusone"),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

The single trial summary in this case

``` r
sim.surv.event.cp1$H1$trial[[1]]$look
```

will show us that interim and final analyses have been performed after
the last arm reaches the respective threshold. We can check this by
adding the number of events in the control group to the matrix of events
in the intervention groups

``` r
sim.surv.event.cp1$H1$trial[[1]]$look$`ev(control)` + 
  sim.surv.event.cp1$H1$trial[[1]]$look[,17:19]
```

Again, we can check the times `t`

``` r
sim.surv.event.cp1$H1$trial[[1]]$look$t
```

and participants recruited `n` at each step.

``` r
sim.surv.event.cp1$H1$trial[[1]]$look$n
```

## Combinations of events and time

As a further feature it is possible to combine event-based and
time-based interim schedules. Using both `event` and `time` in the
`interim` argument will allow us to set a time-based schedule after a
certain threshold/number of events has been observed. Setting
`interim = list(event = 20, time = c(NA, 1, 2), event.type = "min")`
will schedule an interim analysis 1 and 2 years after the last arm has
reached 20 events. The leading `NA` in the vector `time` indicates that
an analysis based on events is to be performed prior to the time
dependent analyses. (`event.type` can be set to either `all`, `control`,
or `min`.)

``` r
sim.surv.ev.t      = batss.surv(
  model            = inla.surv(time, status) ~ trt,
  surv             = simsurv,
  surv.control     = list(lambdas = 2.5, gammas = 1, maxt = 2),
  var              = list(trt = alloc.balanced),
  accr             = rexp,
  accr.control     = list(rate = 100),
  accr.type        = "random",
  family           = "weibullsurv",
  control.family   = list(variant = 0),
  hr               = c(0.6667, 1, 1.5),
  which            = 1:3,
  R                = 2,
  alternative      = "less",
  RAR              = NULL,
  N                = 400,
  prob0            = c(control = 1, A = 1, B = 1, C = 1),
  interim          = list(event = 20, time = c(NA, 1, 2), event.type = "min"),
  eff.arm          = eff.arm.simple,
  eff.arm.control  = list(b = 0.012),
  delta.eff        = c(NA, NA, NA, 0),
  fut.arm          = NULL,
  computation      = "parallel",
  mc.cores         = parallel::detectCores() - 1,
  H0               = FALSE,
  extended         = 1,
  control.inla     = list(cmin = 0))
```

Looking at the summary of a single trial again

``` r
sim.surv.ev.t$H1$trial[[1]]$look
```

we can see that an interim is being performed after the last arm reaches
20 events

``` r
sim.surv.ev.t$H1$trial[[1]]$look[1,16:19]
```

Checking the times `t` we see that the following interim analyses are
conducted exactly 1 and 2 years thereafter

``` r
sim.surv.ev.t$H1$trial[[1]]$look$t
```

And finally, we can check the participants recruited `n` at each step.

``` r
sim.surv.ev.t$H1$trial[[1]]$look$n
```
