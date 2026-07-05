# Initialize the Optimizer.

Prepares the optimizer for use with a specific function and starting
point.

## Usage

``` r
mize_init(
  opt,
  par,
  fg,
  max_iter = Inf,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  abs_tol = NULL,
  rel_tol = abs_tol,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
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

- max_iter:

  (Optional). Maximum number of iterations. See the 'Convergence'
  section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- max_fn:

  (Optional). Maximum number of function evaluations. See the
  'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- max_gr:

  (Optional). Maximum number of gradient evaluations. See the
  'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- max_fg:

  (Optional). Maximum number of function or gradient evaluations. See
  the 'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- abs_tol:

  (Optional). Absolute tolerance for comparing two function evaluations.
  See the 'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- rel_tol:

  (Optional). Relative tolerance for comparing two function evaluations.
  See the 'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- grad_tol:

  (Optional). Absolute tolerance for the length (l2-norm) of the
  gradient vector. See the 'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- ginf_tol:

  (Optional). Absolute tolerance for the infinity norm (maximum absolute
  component) of the gradient vector. See the 'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- step_tol:

  (Optional). Absolute tolerance for the size of the parameter update.
  See the 'Convergence' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

## Value

Initialized optimizer.

## Details

Should be called after creating an optimizer with
[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md)
and before beginning any optimization with
[`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md).
Alternatively, if `fg` and `par` are available when calling
[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md),
they can be passed to that function and the returned optimizer will
already be initialized.
[`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md)
requires an initialized optimizer and does not carry out initialization
itself.

Optional convergence parameters may also be passed here, for use with
[`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md).
They are optional if you do your own convergence checking.

Details of the `fg` list can be found in the 'Details' section of
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

## Examples

``` r

# Create an optimizer
opt <- make_mize(method = "L-BFGS")

# Function to optimize and starting point defined after creating optimizer
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

# Initialize with function and starting point before commencing optimization
opt <- mize_init(opt, rb0, rosenbrock_fg)

# Finally, can commence the optimization loop
par <- rb0
for (iter in 1:3) {
  res <- mize_step(opt, par, rosenbrock_fg)
  par <- res$par
  opt <- res$opt
}
```
