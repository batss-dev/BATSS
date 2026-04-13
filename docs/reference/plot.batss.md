# Plot function for 'BATSS' outputs

Plot for objects of class 'batss'

## Usage

``` r
# S3 method for class 'batss'
plot(
  x,
  which = if (x$type != "surv") 1:4 else 1:5,
  ask = TRUE,
  hypothesis = ifelse(is.null(x$H1), "H0", "H1"),
  title = TRUE,
  legend = TRUE,
  ess = TRUE,
  percentage = TRUE,
  beta = TRUE,
  col = c("#008B0040", "#8B3A3A40", "#8B897040", "#FF990075"),
  smooth = 1,
  ...
)
```

## Arguments

- x:

  An object of class 'batss' (i.e., output of the function
  [batss.glm](batss.glm.md)).

- which:

  An integer vector indicating the list of desired plots. If a subset of
  the plots is required, specify a subset of the numbers `1:4`. By
  default, all plots are provided, i.e., `which=1:4`. Plot `1` displays
  the boxplot of the total sample size as well as the boxplot of the
  sample sizes per group, Plot `2` displays a barplot of the probability
  of stopping at each look, plot `3` displays the violin plot of the
  sample size per group, and plot `4` displays the Monte Carlo trial
  target estimates as a function of the sample size.

- ask:

  A logical. If `TRUE`, the user is prompted to hit the Enter key before
  each plot. See [`par`](https://rdrr.io/r/graphics/par.html)`(ask=.)`.
  The default is `ask=TRUE`.

- hypothesis:

  A character string indicating which alternative hypothesis to use for
  analyses considering both "H0" and "H1", with options "H1" (default)
  and "H0".

- title:

  Either a [logical](https://rdrr.io/r/base/logical.html) indicating if
  a title should be added or a string (of class
  [character](https://rdrr.io/r/base/character.html)) indicating the
  title to be added. If `title` equals `TRUE` (default), the title
  'Under 'H1' or 'Under 'H0' (depending on the argument `hypothesis`) is
  added to the outer margin of the plot. No outer margin space is added
  if `title = FALSE`.

- legend:

  A [logical](https://rdrr.io/r/base/logical.html) (with default set to
  `TRUE`) if a legend should be added at the bottom of plots 3 and 4, or
  a list with names `height`, `cex` and `pt` respectively indicating i/
  the fraction of the plot to be used for the legend as a numeric
  (default is `.15`), ii/ the character expansion factor relative to
  current `par("cex")` as a numeric (default to `1.25`), and iii/ the
  expansion factor(s) for the points as a numeric (default to `2`). The
  input `legend = TRUE` is equivalent to
  `legend = c(height=.15, cex=1, pt=2)`.

- ess:

  A [logical](https://rdrr.io/r/base/logical.html) (with default set to
  `TRUE`) indicating if the expected sample size should be displayed in
  plots 2, 3 and 4, or a list with names `col`, `cex` and `bg`
  respectively indicating i/ the colour of the label as a character
  (plots 2, 3, and 4), ii/ the text expension level as a numerical value
  (plots 2, 3 and 4) and iii/ the background colour as a character (plot
  3). The input `ess = TRUE` is equivalent to
  `ess = list(col="blue", cex=1, bg="#FFD70070")`.

- percentage:

  A [logical](https://rdrr.io/r/base/logical.html) (with default set to
  `TRUE`) indicating if the probability of stopping at each look should
  be displayed in plots 2 (as a percentage), or a list with names `col`,
  and `cex` indicating i/ the colour of the label as a character, ii/
  the text expension level as a numerical value. The input
  `percentage = TRUE` is equivalent to
  `percentage = list(col="violet", cex=1)`.

- beta:

  A [logical](https://rdrr.io/r/base/logical.html) (with default set to
  `TRUE`) indicating if the assumed target parameter value(s) should be
  displayed in plot 4, or a list with names `col`, `lwd`, and `lty`
  respectively indicating i/ the colour of the target parameter value
  horizontal line(s) as a character, ii/ the thickness of the line(s) as
  a numerical value and iii/ the type of line as an integer (see
  [par](https://rdrr.io/r/graphics/par.html) for details). The input
  `beta = TRUE` is equivalent to
  `beta = list(col="blue", lwd=1.5, lty=1)`.

- col:

  A vector of length 4 specifying the colours to be used. Default to
  `c("#8B897040","#008B0040","#8B3A3A40","#FF990075")` where the last
  two digits of the hexadecimal strings specify the level of
  transparency. Refer to the Section 'colour specification' in
  [par](https://rdrr.io/r/graphics/par.html) for details. If the length
  of `col` equals 1, the same colour is used for all cases. For plots
  `1` and `2`, the 3rd colour of the vector `col` is used to display the
  barplot and boxplots.

- smooth:

  A numerical (\>0) indictating the level of smoothing of the violin
  plots (for plot `3`). Default to `1`. When `smooth=NULL`, the
  smoothing value is optimised in the
  [sm.density](https://rdrr.io/pkg/sm/man/sm.density.html) function.

- ...:

  Additional arguments affecting the plot produced, like ylim and ylab.

## Value

Generates graphical displays of results for objects of class 'batss'.

## See also

[`batss.glm()`](batss.glm.md), [`batss.surv()`](batss.surv.md), the
function generating S3 objects of class 'batss'.
