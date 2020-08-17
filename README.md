# mize

[![Travis-CI Build Status](https://travis-ci.org/jlmelville/mize.svg?branch=master)](https://travis-ci.org/jlmelville/mize)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/jlmelville/mize?branch=master&svg=true)](https://ci.appveyor.com/project/jlmelville/mize)
[![Coverage status](https://codecov.io/gh/jlmelville/mize/branch/master/graph/badge.svg)](https://codecov.io/github/jlmelville/mize?branch=master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mize)](https://cran.r-project.org/package=mize)
[![Dependencies](https://tinyverse.netlify.com/badge/mize)](https://cran.r-project.org/package=mize)
[![CRAN Monthly Downloads](https://cranlogs.r-pkg.org/badges/mize)](https://cran.r-project.org/package=mize)
![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mize)
[![Last Commit](https://img.shields.io/github/last-commit/jlmelville/mize)](https://github.com/jlmelville/mize)

Unconstrained Numerical Optimization Algorithms.

`mize` can be used as a standalone function like the `stats::optim` function,
or can be integrated into other packages by creating a stateful optimizer and
handling the iterations, convergence, logging and so on externally.

`mize` knows how to do Broyden-Fletcher-Goldfarb-Shanno (BFGS),
the limited-memory BFGS (L-BFGS), various flavors of Conjugate Gradient (CG),
Nesterov Accelerated Gradient (NAG) and momentum-based methods, among others.

## Installing

```R
# Install from CRAN:
install.packages("mize")

# Or install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("jlmelville/mize")
```

## Documentation

```R
?mize
```

There are also some vignettes:

* `mize.Rmd`, which goes through many of the options available.
* `mmds.Rmd`, which does a simple, but non-trivial, application of `mize` to
carry out metric Multi-Dimensional Scaling on the `eurodist` data set.
* `stateful.Rmd`, which demonstrates how to use `mize` statefully, so you can
manually and externally invoke each iteration step.

## Examples

```R
# Make a list containing the function and gradient:
rosenbrock_fg <- list(
   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
                          200 *        (x[2] - x[1] * x[1])) })
# Starting point:
rb0 <- c(-1.2, 1)

# Minimize using L-BFGS
res <- mize(rb0, rosenbrock_fg, method = "L-BFGS")
# Optimized parameters are in res$par

# Or create an optimizer and then loop manually
opt <- make_mize(method = "L-BFGS")
opt <- mize_init(opt, rb0, rosenbrock_fg)

par <- rb0
done <- FALSE
iter <- 0
while (!done) {
  iter <- iter + 1
  res <- mize_step(opt, par, rosenbrock_fg)
  par <- res$par
  opt <- res$opt
  # Look at res$f for current function value
  # you get to (i.e. have to) decide when to stop
  done <- iter > 30
}
```

## See also

The Wolfe line searches use conversion of Mark Schmidt's
[minFunc routines](http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html),
Carl Edward Rasmussen's
[Matlab code](http://learning.eng.cam.ac.uk/carl/code/minimize/) and Dianne
O'Leary's Matlab translation of the
[Moré-Thuente line search](http://www.cs.umd.edu/users/oleary/software/)
algorithm from [MINPACK](http://www.netlib.org/minpack/).

I also maintain the [funconstrain](https://github.com/jlmelville/funconstrain) package, which contains a large number of test
problems for numerical optimization. See this [gist](https://gist.github.com/jlmelville/2cb8905edd0dbc23806d3122a7a05c5d) for functions
to use mize with funconstrain.

## License

[BSD 2-Clause](https://opensource.org/licenses/BSD-2-Clause).

## Acknowledgments

I am grateful to Hans Werner Borchers, who provided assistance and
encouragement in getting mize onto CRAN.
