# Simple arm efficacy stop

allows stopping an arm for efficacy at a given look when the probability
of the corresponding target parameter being greater or smaller
(depending on the argument `'alternative'` of [batss.glm](batss.glm.md))
than `delta.eff` is greater than a fixed value `b`.

## Usage

``` r
eff.arm.simple(posterior, b)
```

## Arguments

- posterior:

  the 'BATSS' ingredient '`posterior`' corresponding, in this context,
  to the (posterior) probability of the target parameter being greater
  or smaller (depending on the argument `'alternative'` of
  [batss.glm](batss.glm.md)) than '`delta.eff`'.

- b:

  the cut-off value used to declare efficacy (to be defined in
  `eff.arm.control`).

## Value

eff.arm.simple returns a logical constant.
