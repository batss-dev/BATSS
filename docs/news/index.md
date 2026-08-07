# Changelog

## Version 2.0.6

- Fixes a bug in `batss.surv` (reported by Stina Zetterstrom) that
  required intervention groups to be listed in alphabetical order in
  `prob0` for the corresponding allocation probabilities to be correctly
  assigned when using RAR.

## Version 2.0.5

- Resolution of bugs affecting censoring in some scenarios.

## Version 2.0.4

- Minor improvement to `summary.batss` to handles case with only one
  treatment arm (thanks to Trinh, N.-D. for the feedback)

## Version 2.0.3

- Resolution of bugs affecting censoring in some scenarios.

## Version 2.0.2

- Add LICENCE.txt and README.md files

## Version 2.0.1

- Minor improvements to internal/generic functions

## Version 2.0.0

- New `batss.surv` function for time-to-event endpoints: `batss.surv`
  creates `batss` class objects of type `surv` with dediacted `print`,
  `summary` and `plots` methods.
- Update of the website:
  - New ‘Survival endpoint’ example
  - New ‘Interim schedule’ page
  - Improved ‘Priors’ page now including information specific to
    survival endpoints

## Version 1.2.1

- Fixes a bug in `batss.glm` (reported by Stina Zetterstrom) that
  required intervention groups to be listed in alphabetical order in
  `prob0` for the corresponding allocation probabilities to be correctly
  assigned when using RAR.

## Version 1.2.0

CRAN release: 2026-05-28

- Improvements to `batss.glm`:
  - `interim = NA` now specifies a fixed (non-adaptive) design with a
    single look at the maximum sample size `N`
  - better warnings
- Fixes to minor bugs detected by Claude (Opus 4.6) in version 1.1.1 of
  `batss.glm`, `batss.combine`, `summary.batss`, `plot.batss` and
  internal functions, including a fix to `eff.trial.control` and
  `fut.trial.control` being ignored in trial-level stopping decisions

## Version 1.1.0

CRAN release: 2025-09-21

- Improvements to `summary.batss`, that now
  - works with additional predictors to the treatment effect
  - allows to save outputs (as objects of class `summary.batss` with
    corresponding `print` function)
- Improvements to `plot.batss`, that now includes
  - a violin plot of sample sizes per group,
  - a barplot of the probability of stopping the trial at each look
- Update of the website: new ‘ANCOVA’ and ‘Parallelisation’ pages

## Version 1.0.1

CRAN release: 2025-08-27

- New website address (<https://batss-stable.github.io/BATSS/>)

## Version 1.0.0

CRAN release: 2024-10-02

- Improvements to checks
- Improvements to internal/generic functions

## Version 0.7.15

- Resolution of a bug affecting designs with a single target parameter

## Version 0.7.14

CRAN release: 2024-07-14

- Improvements to functions `plot.batss` and `summary.batss` related to
  cases in which efficacy and futility criteria are simulatneously met
  for an arm
- Minor improvements (help pages, printed messages)

## Version 0.7.13

- Minor improvements to function `plot.batss`
- Minor changes to the vignette dedicated to a binomial endpoint

## Version 0.7.12

- Initial release of BATSS with pkgdown
