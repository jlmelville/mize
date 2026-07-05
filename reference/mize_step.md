# One Step of Optimization

Performs one iteration of optimization using a specified optimizer.

## Usage

``` r
mize_step(opt, par, fg)
```

## Arguments

- opt:

  Optimizer, created by
  [`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md).

- par:

  Vector of initial values for the function to be optimized over.

- fg:

  Function and gradient list. See the documentation of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

## Value

Result of the current optimization step, a list with components:

- `opt`: Updated version of the optimizer passed to the `opt` argument.
  Should be passed as the `opt` argument in the next iteration.

- `par`: Updated version of the parameters passed to the `par` argument.
  Should be passed as the `par` argument in the next iteration.

- `nf`: Running total number of function evaluations carried out since
  iteration 1.

- `ng`: Running total number of gradient evaluations carried out since
  iteration 1.

- `f`: Optional. The new value of the function, evaluated at the
  returned value of `par`. Only present if calculated as part of the
  optimization step (e.g. during a line search calculation).

- `g`: Optional. The gradient vector, evaluated at the returned value of
  `par`. Only present if the gradient was calculated as part of the
  optimization step (e.g. during a line search calculation.)

## Details

This function returns both the (hopefully) optimized vector of
parameters, and an updated version of the optimizer itself. This is
intended to be used when you want more control over the optimization
process compared to the more black box approach of the
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md)
function. In return for having to manually call this function every time
you want the next iteration of optimization, you gain the ability to do
your own checks for convergence, logging and so on, as well as take
other action between iterations, e.g. visualization.

Normally calling this function should return a more optimized vector of
parameters than the input, or at least leave the parameters unchanged if
no improvement was found, although this is determined by how the
optimizer was configured by
[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md).
It is very possible to create an optimizer that can cause a solution to
diverge. It is the responsibility of the caller to check that the result
of the optimization step has actually reduced the value returned from
the function being optimized.

Details of the `fg` list can be found in the 'Details' section of
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

## See also

[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md)
to create a value to pass to `opt`,
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
to initialize `opt` before passing it to this function for the first
time. [`mize()`](https://jlmelville.github.io/mize/reference/mize.md)
creates an optimizer and carries out a full optimization with it.

## Examples

``` r
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
rb0 <- c(-1.2, 1)

opt <- make_mize(
  method = "SD", line_search = "const", step0 = 0.0001,
  par = rb0, fg = rosenbrock_fg
)
par <- rb0
for (iter in 1:3) {
  res <- mize_step(opt, par, rosenbrock_fg)
  par <- res$par
  opt <- res$opt
}
```
