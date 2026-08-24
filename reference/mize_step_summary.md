# Mize Step Summary

Produces a result summary for an optimization iteration. Information
such as function value, gradient norm and step size may be returned.

## Usage

``` r
mize_step_summary(opt, par, fg, par_old = NULL, calc_fn = NULL, calc_gr = NULL)
```

## Arguments

- opt:

  Optimizer to generate summary for, from return value of
  [`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md).

- par:

  Vector of parameters at the end of the iteration, from return value of
  [`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md).

- fg:

  Function and gradient list. See the documentation of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- par_old:

  (Optional). Vector of parameters at the end of the previous iteration.
  Used to calculate step size.

- calc_fn:

  (Optional). If `TRUE`, force calculation of function if not already
  cached in `opt`, even if it would not be needed for convergence
  checking.

- calc_gr:

  (Optional). If `TRUE`, force calculation of gradient if not already
  cached in `opt`, even if it would not be needed for convergence
  checking.

## Value

A list containing the available items below. Items whose owning method
or calculation did not supply a value are omitted.

- `opt`: Optimizer with updated state (e.g. function and gradient
  counts).

- `iter`: Iteration number.

- `f`: Function value at `par`.

- `g2n`: 2-norm of the gradient at `par`.

- `ginfn`: Infinity-norm of the gradient at `par`.

- `nf`: Number of function evaluations so far.

- `ng`: Number of gradient evaluations so far.

- `step`: Size of the step between `par_old` and `par`, if `par_old` is
  provided.

- `alpha`: Step length of the gradient descent part of the step.

- `mu`: Momentum coefficient for this iteration.

- `alpha_init`: Initial line-search step length after safeguards.

- `slope_init`: Directional derivative at the start of the line search.

- `ls_reason`: Reason the line search stopped.

- `ls_outcome`: Kind of point selected by the line search.

- `ls_nf`, `ls_ng`: Function and gradient callbacks owned by the line
  search.

- `direction_reason`: Exact-Newton direction provenance.

See the 'Progress' section of
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
diagnostic value meanings.

## Details

By default, convergence tolerance parameters will be used to determine
what function and gradient data is returned. The function value will be
returned if it was already calculated and cached in the optimization
iteration. Otherwise, it will be calculated only if a non-null absolute
or relative tolerance value was asked for. A gradient norm will be
returned only if a non-null gradient tolerance was specified, even if
the gradient is available.

If a function value is required by a tolerance but is not cached for the
relevant value of `par`, it will be calculated here. The calculation
contributes to the total function count and is cached for potential use
in the next iteration. The same rule applies when a gradient tolerance
requires a gradient calculation. Function and gradient calculation can
also be forced here by setting the `calc_fn` and `calc_gr` parameters,
respectively, to `TRUE`.

Requested objective and gradient calculations honor the optimizer's hard
callback budgets. If a request would exceed a budget, the corresponding
field is omitted and the returned `opt` records the termination. Always
retain this updated optimizer before continuing. This summary does not
apply ordinary numerical tolerances or `max_iter`; pass it to
[`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md)
for those checks. Non-finite observations are likewise interpreted
there.

## Examples

``` r
rb_fg <- list(
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

opt <- make_mize(method = "BFGS", par = rb0, fg = rb_fg, max_iter = 30)
mize_res <- mize_step(opt = opt, par = rb0, fg = rb_fg)
# Get info about first step, use rb0 to compare new par with initial value
step_info <- mize_step_summary(mize_res$opt, mize_res$par, rb_fg, rb0)
opt <- step_info$opt
if (!opt$is_terminated) {
  opt <- check_mize_convergence(step_info)
}
```
