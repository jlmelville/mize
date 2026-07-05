# Check Optimization Convergence

Updates the optimizer with information about convergence or termination,
signaling if the optimization process should stop.

## Usage

``` r
check_mize_convergence(mize_step_info)
```

## Arguments

- mize_step_info:

  Step info for this iteration, created by
  [`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)

## Value

`opt` updated with convergence and termination data. See 'Details'.

## Details

On returning from this function, the updated value of `opt` will
contain:

- A boolean value `is_terminated` which is `TRUE` if termination has
  been indicated, and `FALSE` otherwise.

- A list `terminate` if `is_terminated` is `TRUE`. This contains two
  items: `what`, a short string describing what caused the termination,
  and `val`, the value of the termination criterion that caused
  termination. This list will not be present if `is_terminated` is
  `FALSE`.

Convergence criteria are only checked here. To set these criteria, use
[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md)
or
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md).

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
step_info <- mize_step_summary(mize_res$opt, mize_res$par, rb_fg, rb0)
# check convergence by looking at opt$is_terminated
opt <- check_mize_convergence(step_info)
```
