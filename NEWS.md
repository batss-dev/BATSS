## Version 2.0.0

-   New `batss.surv` function for time-to-event endpoints: 
    `batss.surv` creates `batss` class objects of type `surv` with 
    dediacted `print`, `summary` and `plots` methods. 
-   Update of the website: 
    -   New 'Survival endpoint' example
    -   Improved 'Priors' page now including information specific
        to survival endpoints

## Version 1.1.0

-   Improvements to `summary.batss`, that now 
    -   works with additional predictors to the treatment effect 
    -   allows to save outputs (as objects of class `summary.batss`
        with corresponding `print` function)
-   Improvements to `plot.batss`, that now includes
    -   a violin plot of sample sizes per group,
    -   a barplot of the probability of stopping the trial at each look
-   Update of the website: new 'ANCOVA' and 'Parallelisation' pages 

## Version 1.0.1

-   New website address (https://batss-stable.github.io/BATSS/)

## Version 1.0.0

-   Improvements to checks
-   Improvements to internal/generic functions

## Version 0.7.15

-   Resolution of a bug affecting designs with a 
    single target parameter

## Version 0.7.14

-   Improvements to functions `plot.batss` and `summary.batss` 
    related to cases in which efficacy and futility criteria 
    are simulatneously met for an arm
-   Minor improvements (help pages, printed messages)

## Version 0.7.13

-   Minor improvements to function `plot.batss`
-   Minor changes to the vignette dedicated to a binomial endpoint

## Version 0.7.12

-   Initial release of BATSS with pkgdown
