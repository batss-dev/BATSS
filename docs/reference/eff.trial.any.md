# trial efficacy stop

allows stopping the trial for efficacy if *at least one* target
parameter reached efficacy at the look of interest.

## Usage

``` r
eff.trial.any(eff.target)
```

## Arguments

- eff.target:

  the 'BATSS' ingredient '`eff.target`' corresponding to a
  [logical](https://rdrr.io/r/base/logical.html) vector of the same
  length as argument `which` (i.e., the number of target parameters)
  indicating if efficacy was reached for each target parameter at that
  stage or at a previous stage.

## Value

eff.trial.any returns a logical constant.
