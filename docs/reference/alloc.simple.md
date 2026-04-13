# Simple allocation function

alloc.simple independently randomises each unit to a group (i.e., flips
a coin for each unit) so that the observed allocation probabilities may
be far from the target ones. This strategy is often considered to be a
poor choice.

## Usage

``` r
alloc.simple(m, prob)
```

## Arguments

- m:

  the 'BATSS' ingredient '`m`', a scalar corresponding to the number of
  participants to be allocated.

- prob:

  the 'BATSS' ingredient '`prob`', a named vector of allocation ratios
  or probabilities.

## Value

alloc.simple returns an object of class
[factor](https://rdrr.io/r/base/factor.html) of length '`m`' with levels
matching the names of the vector '`prob`'.

## See also

[`alloc.balanced()`](alloc.balanced.md), another group allocation
function.

## Examples

``` r
alloc.simple(100, prob = c(A=.4,B=.6))
#>   [1] B B B B A A A B B B B B B B B A B A A B A A B B A A B B B A B B A A B B A
#>  [38] B A A B A A A A B B B B B B A B A B B A B B B B B A B B B B B B B B B B B
#>  [75] A B B A B B B A B B A B A B B A A B B A B B A A B B
#> Levels: A B
table(alloc.simple(100, prob = c(A=.4,B=.6)))
#> 
#>  A  B 
#> 49 51 
table(alloc.simple(100, prob = c(A=.4,B=.6)))
#> 
#>  A  B 
#> 47 53 
```
