# Summary function for 'BATSS' outputs

Summary method function for objects of class 'batss'.

## Usage

``` r
# S3 method for class 'batss'
summary(object, extended = NULL, ...)
```

## Arguments

- object:

  An object of class 'batss' (i.e., output of the function
  [batss.glm](batss.glm.md)).

- extended:

  A logical indicating if a standard (extended = FALSE, default) or
  extended output (extended = TRUE) should be returned. Default to
  `NULL` in which case the input of the argument `extended` chosen when
  generating `object` with [`batss.glm()`](batss.glm.md) or
  [`batss.surv()`](batss.surv.md) is used.

- ...:

  For future use

## Value

Object of class 'summary.batss'.

The function summary.batss returns an S3 list of class 'summary.batss'
with available print functions. The list elements are

- beta - A data frame providing information related to the beta
  parameter vector, such as parameter names and values, for example.

- look - A data frame providing information related to looks, like
  sample size of a given interim (m) and cumulative sample size at a
  given interim (n), for example.

- par - A list providing different information, like the used seeds
  (seed) and the groups (group), for example.

- H1 - A list providing trial aggregated results under the alternative,
  like the probability of efficacy, futility, or both, per arm or
  globally (`object$H1$target`), the probability of stopping early for
  efficacy (`object$H1$efficacy`) and futility (`object$H1$futility`),
  the sample size expectation, standard deviation, and quantiles 0.1,
  0.5 and 0.9, per group and overall (`object$H1$summary.sample.sizes`),
  the probabilities associated to each combination of efficacy and
  futility per group (scenario).

- H0 - A list providing trial aggregated results under the global null
  hypothesis (same structure as H1).

- call - The matched call.

- type - The type of 'BATSS' analysis (currently either 'glm' or
  'surv').

## See also

[`batss.glm()`](batss.glm.md), [`batss.surv()`](batss.surv.md), the
functions generating S3 objects of class 'batss'.
