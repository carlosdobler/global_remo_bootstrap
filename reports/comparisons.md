Comparisons
================
Carlos Dobler
2022-10-27

Here I compare different ensembles of REMO/RegCM model output
corresponding to days above 32 C in a 2C warming level over Europe.

The following three set of maps compare “raw” output between REMO (3
members), RegCM (3 members), and REMO+RegCM (6 members) ensembles. I
show differences in the mean count of days above threshold, as well as
percentiles 5th, 50th (median) and 95th.

![](comparisons_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

![](comparisons_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

![](comparisons_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

In the next set of maps I compare some of the previous maps with two
bootstrapped experiments conducted with REMO. The first experiment
represents a regular bootstrap. For each cell, I resampled the original
sample of 63 observations (21 years x 3 REMO models) with replacement
1000 times. This means: imagine I have a bag of 63 balls, each one
representing one observation (1 observation = the count of days above
32C in 1 year for 1 model). I randomly pick one ball, register its
value, and return it to the bag. I do that 63 times. Then I calculate
the mean and percentiles 5th, 50th and 95th out of those 63
observations. Then I re-do that 1000 times. Lastly, I calculate the mean
of those 1000 means, 5th, 50th, and 95th percentiles. (Note: those 1000
means (or any of the percentiles) are somewhat normally distributed.
This means that I could also calculate the 2.5th and 97.5th percentiles
of those 1000 means to obtain my 95% confidence interval. In other
words, I could say something like, “based on REMO, the mean count of
days above 32C is between X (the 2.5th perc.) and Y (the 97.5th perc.)
with 95% confidence”).

![](comparisons_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

![](comparisons_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

As a second experiment I conducted a parametric bootstrapping. For each
cell, I fitted a binomial distribution to its 63 values. I’m using
binomial because our underlying variable is binomial: days either exceed
the threshold (=1), or they don’t (=0). Once fitted, I generate 1000
samples of 63 observations based on that distribution. And then, same as
before, for each of the 1000 iterations, I calculate the mean and
percentiles of the 63 observations, to then calculate the mean of those
1000 statistics.

Note: when fitting a distribution we assume our data (the 63
observations) is “incomplete”; it’s just a sample out of a larger
“population”. A fitted distribution essentially represents how that
population looks like. The following figure shows nine examples
comparing the empirical (“observed”) distribution of the 63 observations
found in a cell (solid line) and their fitted distribution (dashed).

![](comparisons_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

As mentioned above, from the fitted distribution I drew 1000 random
samples of 63 observations and calculated the means of the 1000 means
and percentiles. See below how that ensemble compares with the raw REMO
ensemble.

![](comparisons_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

![](comparisons_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
