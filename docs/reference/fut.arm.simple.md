# arm futility stop

allows stopping an arm for futility when the probability of the
corresponding target parameter being greater or smaller (depending on
the argument `'alternative'` of [batss.glm](batss.glm.md)) than
'`delta.fut`' is smaller than a fixed value '`b`'

## Usage

``` r
fut.arm.simple(posterior, b)
```

## Arguments

- posterior:

  the 'BATSS' ingredient '`posterior`' corresponding, in this context,
  to the (posterior) probability of the target parameter being greater
  or smaller (depending on the argument `'alternative'` of
  [batss.glm](batss.glm.md)) than '`delta.fut`'.

- b:

  the cut-off value used to declare futility (to be defined in
  `fut.arm.control`).

## Value

fut.arm.simple returns a logical constant.
