# Balanced allocation function

alloc.balanced first allocates the largest possible number of units to
the different groups given their exact target probabilities and then
assigns randomly the remaining units to the different groups according
to multinomial draws. This method leads to observed allocation
probabilities matching the target ones when m\*prob is an integer for
each group and to observed allocation probabilities (on average) closer
to the target ones compared to [alloc.simple](alloc.simple.md).

## Usage

``` r
alloc.balanced(m, prob)
```

## Arguments

- m:

  the 'BATSS' ingredient '`m`', a scalar corresponding to the number of
  participants to be allocated.

- prob:

  the 'BATSS' ingredient '`prob`', a named vector of allocation ratios
  or probabilities.

## Value

alloc.balanced returns an object of class
[factor](https://rdrr.io/r/base/factor.html) of length '`m`' with levels
matching the names of the vector '`prob`'.

## See also

[`alloc.simple()`](alloc.simple.md), another group allocation function.

## Examples

``` r
alloc.balanced(100, prob = c(A=.4,B=.6))
#>   [1] A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A
#>  [38] A A A B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B
#>  [75] B B B B B B B B B B B B B B B B B B B B B B B B B B
#> Levels: A B
table(alloc.balanced(100, prob = c(A=.4,B=.6)))
#> 
#>  A  B 
#> 40 60 
table(alloc.balanced(100, prob = c(A=.4,B=.6)))
#> 
#>  A  B 
#> 40 60 
```
