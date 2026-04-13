# 'Optimal' control allocation

technically not response adaptive but keeps allocation ratio to control
at the square root of active intervention arms

## Usage

``` r
RAR.optimal(active)
```

## Arguments

- active:

  the 'BATSS' ingredient '`active`' corresponding to a
  [logical](https://rdrr.io/r/base/logical.html) vector of the same
  length and order as `'prob0'` (i.e., number of arms initially included
  in the study including the reference group)) and indicating if each
  arm is active at the look of interest.

## Value

RAR.optimal returns a vector of probabilities with length of active.
