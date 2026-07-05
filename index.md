# mize

[![codecov](https://codecov.io/github/jlmelville/mize/graph/badge.svg?token=v14HaDPD2g)](https://codecov.io/github/jlmelville/mize)
[![R-CMD-check](https://github.com/jlmelville/mize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jlmelville/mize/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/mize)](https://cran.r-project.org/package=mize)
[![CRAN monthly
downloads](https://cranlogs.r-pkg.org/badges/mize)](https://cran.r-project.org/package=mize)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/mize)](https://cran.r-project.org/package=mize)
[![Last
commit](https://img.shields.io/github/last-commit/jlmelville/mize)](https://github.com/jlmelville/mize)

Unconstrained numerical optimization algorithms implemented in R.

`mize` can be used as a standalone optimizer like
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), or it can be
integrated into other packages by creating a stateful optimizer and
handling iterations, convergence, and logging externally.

`mize` includes Broyden-Fletcher-Goldfarb-Shanno (BFGS), limited-memory
BFGS (L-BFGS), conjugate gradient (CG), Nesterov accelerated gradient
(NAG), momentum-based methods, and related line searches.

## Installation

``` r

install.packages("mize")
```

Install the development version from GitHub with:

``` r

# install.packages("devtools")
devtools::install_github("jlmelville/mize")
```

## Quick Start

``` r

library(mize)

rosenbrock_fg <- list(
  fn = function(x) {
    100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
  },
  gr = function(x) {
    c(
      -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
      200 * (x[2] - x[1] * x[1])
    )
  }
)

res <- mize(c(-1.2, 1), rosenbrock_fg, method = "L-BFGS")
res$par
res$f
```

For stateful optimization, create an optimizer and call
[`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md)
manually:

``` r

opt <- make_mize(method = "L-BFGS", par = c(-1.2, 1), fg = rosenbrock_fg)
par <- c(-1.2, 1)

for (i in seq_len(30)) {
  step <- mize_step(opt, par, rosenbrock_fg)
  par <- step$par
  opt <- step$opt
}
```

## Documentation

The package website is <https://jlmelville.github.io/mize/>.

Start with:

- [Getting
  started](https://jlmelville.github.io/mize/articles/mize.html)
- [Convergence](https://jlmelville.github.io/mize/articles/convergence.html)
- [Stateful
  optimization](https://jlmelville.github.io/mize/articles/stateful.html)
- [Metric MDS
  example](https://jlmelville.github.io/mize/articles/mmds.html)
- [Function
  reference](https://jlmelville.github.io/mize/reference/index.html)

The CRAN package also includes the main vignettes.

## See Also

The Wolfe line searches use conversions of Mark Schmidt’s [minFunc
routines](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html), Carl
Edward Rasmussen’s [Matlab
code](https://gaussianprocess.org/gpml/code/matlab/doc/), and Dianne
O’Leary’s Matlab translation of the [More-Thuente line
search](https://www.cs.umd.edu/users/oleary/software/) algorithm from
[MINPACK](https://www.netlib.org/minpack/).

I also maintain the
[funconstrain](https://github.com/jlmelville/funconstrain) package,
which contains test problems for numerical optimization.

## License

[BSD 2-Clause](https://opensource.org/licenses/BSD-2-Clause).

## Acknowledgments

I am grateful to Hans Werner Borchers, who provided assistance and
encouragement in getting mize onto CRAN.
