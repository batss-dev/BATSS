# trial futility stop

allows stopping the trial for efficacy if *all* active treatment reached
futility at the look of interest or before.

## Usage

``` r
fut.trial.all(fut.target)
```

## Arguments

- fut.target:

  the 'BATSS' ingredient '`fut.target`' corresponding to a
  [logical](https://rdrr.io/r/base/logical.html) vector of the same
  length as argument `which` (i.e., the number of target parameters)
  indicating if futility was declared for each target parameter at that
  stage or at a previous stage.

## Value

fut.trial.all returns a logical constant.
