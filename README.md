

# rBiasCorrection

<!-- badges: start -->

[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://www.r-pkg.org/badges/version/rBiasCorrection)](https://cran.r-project.org/package=rBiasCorrection)
[![CRAN
checks](https://badges.cranchecks.info/worst/rBiasCorrection.svg)](https://cran.r-project.org/web/checks/check_results_rBiasCorrection.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/rBiasCorrection?color=blue)](https://cran.r-project.org/package=rBiasCorrection)
[![](http://cranlogs.r-pkg.org/badges/last-month/rBiasCorrection?color=blue)](https://cran.r-project.org/package=rBiasCorrection)
[![Dependencies](https://tinyverse.netlify.app/badge/rBiasCorrection)](https://cran.r-project.org/package=rBiasCorrection)
[![R build
status](https://github.com/kapsner/rBiasCorrection/workflows/R%20CMD%20Check%20via%20%7Btic%7D/badge.svg)](https://github.com/kapsner/rBiasCorrection/actions)
[![R build
status](https://github.com/kapsner/rBiasCorrection/workflows/lint/badge.svg)](https://github.com/kapsner/rBiasCorrection/actions)
[![R build
status](https://github.com/kapsner/rBiasCorrection/workflows/test-coverage/badge.svg)](https://github.com/kapsner/rBiasCorrection/actions)
[![](https://codecov.io/gh/https://github.com/kapsner/rBiasCorrection/branch/master/graph/badge.svg)](https://codecov.io/gh/https://github.com/kapsner/rBiasCorrection)
[![](https://img.shields.io/badge/doi-10.1002/ijc.33681-yellow.svg)](https://doi.org/10.1002/ijc.33681)

<!-- badges: end -->

`rBiasCorrection` is published in *‘BiasCorrector: fast and accurate
correction of all types of experimental biases in quantitative DNA
methylation data derived by different technologies’ (2021)* in the
*International Journal of Cancer* (DOI:
[https://onlinelibrary.wiley.com/doi/10.1002/ijc.33681](https://doi.org/10.1002/ijc.33681)).

`rBiasCorrection` is the R implementation with minor modifications of
the algorithms described by Moskalev et al. in their research article
*‘Correction of PCR-bias in quantitative DNA methylation studies by
means of cubic polynomial regression’*, published 2011 in *Nucleic acids
research, Oxford University Press* (DOI:
<https://doi.org/10.1093/nar/gkr213>).

# Installation

## CRAN version

You can install `rBiasCorrection` simply with via R’s `install.packages`
interface:

``` r
install.packages("rBiasCorrection")
```

## Development version

If you want to use the latest development version, you can install the
github version of `rBiasCorrection` with:

``` r
install.packages("remotes")
remotes::install_github("kapsner/rBiasCorrection")
```

## Example

This is a basic example which shows you how to correct PCR-bias in
quantitative DNA methylation data:

``` r
library(rBiasCorrection)

# define input file paths
experimental <- file.path(tempdir(), "/experimental_data.csv")
calibration <- file.path(tempdir(), "/calibration_data.csv")

# create example files from provided example dataset
data.table::fwrite(
  rBiasCorrection::example.data_experimental$dat,
  experimental
)
data.table::fwrite(
  rBiasCorrection::example.data_calibration$dat,
  calibration
)

# run bias correction algorithm
biascorrection(
  experimental = experimental,
  calibration = calibration,
  samplelocusname = "BRAF"
)
```

More detailed information on how to use the package `rBiasCorrection`
can be found in the
[vignette](https://cran.r-project.org/web/packages/rBiasCorrection/vignettes/rBiasCorrection_howto.html)
and the
[FAQs](https://github.com/kapsner/rBiasCorrection/blob/master/FAQ.md).

## Available Fitting Options (TODO)

There are three fitting options available for fitting the non-linear
least squares (nls) algorithm with `rBiasCorrection`. The default method
(used in the publication) is to fit nls with the Gauss-Newton algorithm
and define for each parameter that should be optimized a random grid
between -1000 and 1000 for initializing the starting estimates
(`options(rBiasCorrection.nls_implementation = "GN.paper")`.  
For making a better guess on the starting estimates when fitting nls
with the Gauss-Newton algorithm
(`options(rBiasCorrection.nls_implementation = "GN.guess")`), the
estimates of a linear model (for both hyperbolic corrections) and of a
cubic model (for the cubic correction with defined minimum- and maximum
values (`minmax = TRUE`)) are computed for initializing the nls (see
details below).  
The third option is to fit nls with the Levenberg-Marquardt algorithm
(using the implementation from the `minpack.lm` R package). In this
case, the start estimates of the nls model are also guessed using either
a linear or a cubic model (as previously described).

### `GN.paper`

Algorithm: Gauss-Newton

Parameterizing `nls2::nls2()` with starting values:

- hyperbolic equation: a = b = d = c(-1000, 1000)
- hyperbolic equation (minmax): b = c(-1000, 1000)
- cubic equation (minmax): a, b = c(-1000, 1000)

``` r
options(rBiasCorrection.nls_implementation = "GN.paper")
```

### `GN.guess`

Algorithm: Gauss-Newton

Parameterizing `nls2::nls2()` with starting values:

- hyperbolic equation: fitting a linear regression and taking the
  intercept and the beta as starting values and defaulting `d` to `1000`
- hyperbolic equation (minmax): fitting a linear regression and taking
  the beta as starting value
- cubic equation (minmax): fitting a cubic regression and taking the
  betas for the cubic and the quadratic term as starting values

``` r
options(rBiasCorrection.nls_implementation = "GN.guess")
```

### `LM`

Algorithm: Levenberg-Marquardt

Parameterizing `minpack.lm::nlsLM()` with starting values: same as
guessing starting values for option `GN.guess`

``` r
options(rBiasCorrection.nls_implementation = "LM")
```

## BiasCorrector

The GUI `BiasCorrector` provides the functionality implemented in
`rBiasCorrection` in a web application. For further information please
visit <https://github.com/kapsner/BiasCorrector>.

## FAQ

For further information, please refer to the [frequently asked
questions](https://github.com/kapsner/rBiasCorrection/blob/master/FAQ.md).

## Citation

L.A. Kapsner, M.G. Zavgorodnij, S.P. Majorova, A. Hotz‐Wagenblatt, O.V.
Kolychev, I.N. Lebedev, J.D. Hoheisel, A. Hartmann, A. Bauer, S. Mate,
H. Prokosch, F. Haller, and E.A. Moskalev, BiasCorrector: fast and
accurate correction of all types of experimental biases in quantitative
DNA methylation data derived by different technologies, Int. J. Cancer.
(2021) ijc.33681.
doi:[10.1002/ijc.33681](https://onlinelibrary.wiley.com/doi/10.1002/ijc.33681).

``` bibtex
@article{kapsner2021,
  title = {{{BiasCorrector}}: Fast and Accurate Correction of All Types of Experimental Biases in Quantitative {{DNA}} Methylation Data Derived by Different Technologies},
  author = {Kapsner, Lorenz A. and Zavgorodnij, Mikhail G. and Majorova, Svetlana P. and Hotz-Wagenblatt, Agnes and Kolychev, Oleg V. and Lebedev, Igor N. and Hoheisel, J{\"o}rg D. and Hartmann, Arndt and Bauer, Andrea and Mate, Sebastian and Prokosch, Hans-Ulrich and Haller, Florian and Moskalev, Evgeny A.},
  year = {2021},
  month = may,
  pages = {ijc.33681},
  issn = {0020-7136, 1097-0215},
  doi = {10.1002/ijc.33681},
  journal = {International Journal of Cancer},
  language = {en}
}
```

## More Infos

- Original work by Moskalev et al.: https://doi.org/10.1093/nar/gkr213
