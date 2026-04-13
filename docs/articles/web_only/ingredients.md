# Ingredients

At each look, the function `batss.glm` defines a number of *ingredients*
or derived quantities which can be used as inputs to the user-defined
functions specified in arguments `eff.arm`, `fut.arm`, `eff.trial`,
`fut.trial`, and `RAR` to perform adaptations. This document describes
each *ingredient*.

## `n`

- **object type**: a vector of integers of the same length and order as
  argument `prob0`.

- **description**: `n` provides the number of participants that have
  completed follow up in each treatment arm at the look of interest.
  This includes control and intervention arms, both active and dropped
  (inactive). `sum(n)` therefore corresponds to the total number of
  participants that have completed follow-up to-date in the trial.

## `m`

- **object type**: an integer scalar.

- **description**: `m` corresponds to the number of participants to be
  randomised to the different groups (including control) at a given
  look. It is therefore needed in `RAR` and group allocation functions .

## `N`

- **object type**: an integer scalar.

- **description**: `N` corresponds to the maximum sample size of the
  trial as specified in argument `N` of `batss.glm`.

## `posterior`

Depending on the context - efficacy, futility or RAR - `posterior`
corresponds to (possibly) different probabilities and takes different
lengths.

- Efficacy-related *ingredient*
  - **object type**: a numerical scalar.
  - **description**: in this context, `posterior` corresponds to the
    (posterior) probability of the target parameter being greater or
    smaller (depending on the argument `'alternative'` of \[batss.glm\])
    than ‘`delta.eff`’.
- Futility-related *ingredient*
  - **object type**: a numerical scalar.
  - **description**: in this context, `posterior` corresponds to the
    (posterior) probability of the target parameter being greater or
    smaller (depending on the argument `'alternative'` of \[batss.glm\])
    than ‘`delta.fut`’.
- RAR-related *ingredient*
  - **object type**: a numerical vector of length ‘number of active
    target parameters’ at a given look.
  - **description**: in this context, `posterior` corresponds to the
    vector of (posterior) probabilities of the active target parameters
    being greater or smaller (depending on the argument `'alternative'`
    of \[batss.glm\]) than ‘`delta.RAR`’.

## `active`

- **object type**: vector of the same length and order as `'prob0'`.

- **description**: `active` indicates if the different arms (including
  the control group) are active (`TRUE`) or not (`FALSE`) at the look of
  interest.

## `ref`

- **object type**: a vector of the same length and order as `'prob0'`.

- **description**: `ref` indicates which group is the reference one
  (`TRUE`). It is typically the first one as the reference group is
  typically not dropped and has position 1 in `'prob0'`.

## `prob`

- **object type**: a named numerical vector of length ‘number of active
  arms (including the reference group)’.

- **description**: `prob` is the output of the under-defined function
  provided in ‘`RAR`’ and and is used in defining the vector of
  allocation ratios or probabilities for all active arms (including the
  reference group). Note that, if `RAR = NULL`, `prob` equals `prob0`
  (with appropriate subsetting for arms that have been dropped).

## `eff.target`

- **object type**: a vector of the same length as argument `which`
  (i.e., the number of target parameters)

- **description**: `eff.target` indicates if efficacy was reached for
  each target parameter at the stage of interest or before.

## `fut.target`

- **object type**: a vector of the same length as argument `which`
  (i.e., the number of target parameters)

- **description**: `fut.target` indicates if futility was declared for
  each target parameter at the stage of interest or before.
