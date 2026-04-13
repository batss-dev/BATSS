# RAR of Trippa et al. (2012)

define the group allocation probabilities based on the response adaptive
randomisation rule of Trippa et al. (2012)

## Usage

``` r
RAR.trippa(posterior, n, N, ref, active, gamma, eta, nu)
```

## Arguments

- posterior:

  the 'BATSS' ingredient '`posterior`' corresponding, in this context,
  to the (posterior) probability of the active target parameters being
  greater or smaller (depending on the argument `'alternative'` of
  [batss.glm](batss.glm.md)) than '`delta.RAR`'.

- n:

  the 'BATSS' ingredient '`n`' corresponding to the vector of number of
  recruited participants per arm including the control group at the look
  of interest.

- N:

  the 'BATSS' ingredient '\`N' corresponding to the maximum (planned)
  sample size.

- ref:

  the 'BATSS' ingredient '`ref`' corresponding to a
  [logical](https://rdrr.io/r/base/logical.html) vector of the same
  length and order as `'prob0'` (i.e., number of arms initially included
  in the study including the reference group)) and indicating which
  group is the reference one.

- active:

  the 'BATSS' ingredient '`active`' corresponding to a
  [logical](https://rdrr.io/r/base/logical.html) vector of the same
  length and order as `'prob0'` (i.e., number of arms initially included
  in the study including the reference group)) and indicating if each
  arm is active at the look of interest.

- gamma:

  a scaling factor (to be defined in `RAR.arm.control`).

- eta:

  a scaling factor (to be defined in `RAR.arm.control`).

- nu:

  a scaling factor (to be defined in `RAR.arm.control`).

## Value

RAR.trippa returns a vector of probabilities with length of active.
