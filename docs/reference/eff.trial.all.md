# trial efficacy stop

allows stopping the trial for efficacy if *all* target parameters
reached efficacy at the look of interest or before.

## Usage

``` r
eff.trial.all(eff.target)
```

## Arguments

- eff.target:

  the 'BATSS' ingredient '`eff.target`' corresponding to a
  [logical](https://rdrr.io/r/base/logical.html) vector of the same
  length as argument `which` (i.e., the number of target parameters)
  indicating if efficacy was reached for each target parameter at that
  stage or at a previous stage.

## Value

eff.trial.all returns a logical constant.
