# Check Gradient Consistency

Compares an analytic gradient with a finite-difference approximation of
the objective function at a parameter vector.

## Usage

``` r
check_mize_gradient(
  fg,
  par,
  method = c("central", "forward"),
  rel_eps = NULL,
  abs_eps = 0,
  ...
)
```

## Arguments

- fg:

  Function and gradient list. See the documentation of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- par:

  Parameter vector where the gradient should be checked.

- method:

  Finite-difference method: either `"central"` or `"forward"`.

- rel_eps:

  Relative finite-difference step. If `NULL`, a method-specific default
  is used.

- abs_eps:

  Non-negative scalar added to each finite-difference step.

- ...:

  Additional arguments passed to `fg$fn`, `fg$gr`, and `fg$fg` when
  present.

## Value

A list with components:

- `par`: Parameter vector checked.

- `method`: Finite-difference method used.

- `step`: Per-coordinate finite-difference steps.

- `fn`: Function value at `par`.

- `gr`: Gradient returned by `fg$gr(par)`.

- `gr_fd`: Finite-difference gradient approximation.

- `error`: Difference `gr - gr_fd`.

- `abs_error`: Absolute coordinate-wise errors.

- `rel_error`: Relative coordinate-wise errors, scaled by the larger
  coordinate magnitude of `gr` and `gr_fd`.

- `max_abs_error`: Maximum absolute error.

- `max_rel_error`: Maximum relative error.

- `worst_abs`: One-row data frame for the coordinate with the largest
  absolute error.

- `worst_rel`: One-row data frame for the coordinate with the largest
  relative error.

- `by_coord`: Data frame with per-coordinate gradient comparisons.

- `fg_consistency`: List describing whether `fg$fg()` was checked and,
  when available, its function and gradient consistency with `fg$fn()`
  and `fg$gr()`.

## Details

The `fg` argument uses the same list shape as
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md): it must
contain `fn` and `gr` functions. If it also contains an `fg` function,
the combined function-and-gradient result is checked for consistency
with the separate `fn` and `gr` functions at `par`.

The finite-difference step for each coordinate is
`abs_eps + rel_eps * pmax(abs(par), 1)`. Central differences are usually
more accurate and use twice as many function evaluations as forward
differences.

## See also

[`mize()`](https://jlmelville.github.io/mize/reference/mize.md)

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

check <- check_mize_gradient(rb_fg, c(-1.2, 1))
check$max_abs_error
#> [1] 2.208253e-08
```
