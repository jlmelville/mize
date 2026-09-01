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

- `nf`: Total number of function evaluations over the optimizer's
  lifetime. This count persists across
  [`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
  calls.

- `ng`: Total number of gradient evaluations over the optimizer's
  lifetime. This count persists across
  [`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
  calls.

- `nh`: Total number of accepted Hessian callback evaluations over the
  optimizer's lifetime. This count persists across
  [`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
  calls.

- `nhi`: Total number of accepted inverse-Hessian callback evaluations
  over the optimizer's lifetime. This count persists across
  [`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
  calls.

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

A stateful optimization loop divides stopping work across three calls.
`mize_step()` can stop immediately when a hard callback budget is
exhausted or an optimization method fails. Ordinary numerical tolerances
and `max_iter` are applied by
[`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md).
After each active step, call
[`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md),
retain its returned `opt`, and then call
[`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md)
if the summary did not terminate the optimizer. The examples below show
this sequence. A candidate parameter vector containing `Inf`, `NaN`, or
`NA` terminates the optimizer with `opt$terminate$what = "par_inf"`; the
returned `par` is rolled back to its value at the start of the step.

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
  par = rb0, fg = rosenbrock_fg, max_iter = 3
)
par <- rb0
while (!opt$is_terminated) {
  par_old <- par
  step_result <- mize_step(opt, par, rosenbrock_fg)
  opt <- step_result$opt
  par <- step_result$par
  if (opt$is_terminated) {
    break
  }

  step_info <- mize_step_summary(opt, par, rosenbrock_fg, par_old)
  opt <- step_info$opt
  if (opt$is_terminated) {
    break
  }
  opt <- check_mize_convergence(step_info)
}
```
