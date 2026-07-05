# Create an Optimizer

Factory function for creating a (possibly uninitialized) optimizer.

## Usage

``` r
make_mize(
  method = "L-BFGS",
  norm_direction = FALSE,
  scale_hess = TRUE,
  memory = 5,
  cg_update = "PR+",
  preconditioner = "",
  tn_init = 0,
  tn_exit = "curvature",
  nest_q = 0,
  nest_convex_approx = FALSE,
  nest_burn_in = 0,
  step_up = 1.1,
  step_up_fun = c("*", "+"),
  step_down = NULL,
  dbd_weight = 0.1,
  line_search = "More-Thuente",
  c1 = 1e-04,
  c2 = NULL,
  step0 = NULL,
  step_next_init = NULL,
  try_newton_step = NULL,
  ls_max_fn = 20,
  ls_max_gr = Inf,
  ls_max_fg = Inf,
  ls_max_alpha_mult = Inf,
  ls_max_alpha = Inf,
  ls_safe_cubic = FALSE,
  strong_curvature = NULL,
  approx_armijo = NULL,
  mom_type = NULL,
  mom_schedule = NULL,
  mom_init = NULL,
  mom_final = NULL,
  mom_switch_iter = NULL,
  mom_linear_weight = FALSE,
  use_init_mom = FALSE,
  restart = NULL,
  restart_wait = 10,
  par = NULL,
  fg = NULL,
  max_iter = 100,
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

- method:

  Optimization method. See 'Details' of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- norm_direction:

  If `TRUE`, then the steepest descent direction is normalized to unit
  length. Useful for adaptive step size methods where the previous step
  size is used to initialize the next iteration.

- scale_hess:

  if `TRUE`, the approximation to the inverse Hessian is scaled
  according to the method described by Nocedal and Wright (approximating
  an eigenvalue). Applies only to the methods `BFGS` (where the scaling
  is applied only during the first step) and `L-BFGS` (where the scaling
  is applied during every iteration). Ignored otherwise.

- memory:

  The number of updates to store if using the `L-BFGS` method. Ignored
  otherwise. Must be a positive integer.

- cg_update:

  Type of update to use for the `"CG"` method. For details see the "CG"
  subsection of the "Optimization Methods" section. Ignored if `method`
  is not `"CG"`.

- preconditioner:

  Type of preconditioner to use in Truncated Newton. Leave blank or set
  to `"L-BFGS"` to use a limited memory BFGS preconditioner. Use the
  `"memory"` parameter to control the number of updates to store.
  Applies only if `method = "TN"`, or `"CG"`, ignored otherwise.

- tn_init:

  Type of initialization to use in inner loop of Truncated Newton. Use
  `0` to use the zero vector (the usual TN initialization), or
  `"previous"` to use the final result from the previous iteration, as
  suggested by Martens (2010). Applies only if `method = "TN"`, ignored
  otherwise.

- tn_exit:

  Type of exit criterion to use when terminating the inner CG loop of
  Truncated Newton method. Either `"curvature"` to use the standard
  negative curvature test, or `"strong"` to use the modified "strong"
  curvature test in TNPACK (Xie and Schlick, 1999). Applies only if
  `method = "TN"`, ignored otherwise.

- nest_q:

  Strong convexity parameter for the `"NAG"` method's momentum term.
  Must take a value between 0 (strongly convex) and 1 (results in
  steepest descent).Ignored unless the `method` is `"NAG"` and
  `nest_convex_approx` is `FALSE`.

- nest_convex_approx:

  If `TRUE`, then use an approximation due to Sutskever for calculating
  the momentum parameter in the NAG method. Only applies if `method` is
  `"NAG"`.

- nest_burn_in:

  Number of iterations to wait before using a non-zero momentum. Only
  applies if using the `"NAG"` method or setting the `momentum_type` to
  "Nesterov".

- step_up:

  Value by which to increase the step size for the `"bold"` step size
  method or the `"DBD"` method.

- step_up_fun:

  Operator to use when combining the current step size with `step_up`.
  Can be one of `"*"` (to multiply the current step size with `step_up`)
  or `"+"` (to add).

- step_down:

  Multiplier to reduce the step size by if using the `"DBD"` method or
  the `"bold"`. Can also be used with the `"back"` line search method,
  but is optional. Should be a positive value less than 1.

- dbd_weight:

  Weighting parameter used by the `"DBD"` method only, and only if no
  momentum scheme is provided. Must be an integer between 0 and 1.

- line_search:

  Type of line search to use. See 'Details' of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- c1:

  Sufficient decrease parameter for Wolfe-type line searches. Should be
  a value between 0 and 1.

- c2:

  Sufficient curvature parameter for line search for Wolfe-type line
  searches. Should be a value between `c1` and 1.

- step0:

  Initial value for the line search on the first step. See 'Details' of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- step_next_init:

  For Wolfe-type line searches only, how to initialize the line search
  on iterations after the first. See 'Details' of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- try_newton_step:

  For Wolfe-type line searches only, try the line step value of 1 as the
  initial step size whenever `step_next_init` suggests a step size \> 1.
  Defaults to `TRUE` for quasi-Newton methods such as BFGS and L-BFGS,
  `FALSE` otherwise.

- ls_max_fn:

  Maximum number of function evaluations allowed during a line search.

- ls_max_gr:

  Maximum number of gradient evaluations allowed during a line search.

- ls_max_fg:

  Maximum number of function or gradient evaluations allowed during a
  line search.

- ls_max_alpha_mult:

  The maximum value that can be attained by the ratio of the initial
  guess for alpha for the current line search, to the final value of
  alpha of the previous line search. Used to stop line searches
  diverging due to very large initial guesses. Only applies for
  Wolfe-type line searches.

- ls_max_alpha:

  Maximum value of alpha allowed during line search. Only applies for
  `line_search = "more-thuente"`.

- ls_safe_cubic:

  (Optional). If `TRUE`, check that cubic interpolation in the Wolfe
  line search does not produce too small a value. Only applies for
  `line_search = "more-thuente"`.

- strong_curvature:

  (Optional). If `TRUE` use the strong curvature condition in Wolfe line
  search. See the 'Line Search' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- approx_armijo:

  (Optional). If `TRUE` use the approximate Armijo condition in Wolfe
  line search. See the 'Line Search' section of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for
  details.

- mom_type:

  Momentum type, either `"classical"` or `"nesterov"`.

- mom_schedule:

  Momentum schedule. See 'Details' of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- mom_init:

  Initial momentum value.

- mom_final:

  Final momentum value.

- mom_switch_iter:

  For `mom_schedule` `"switch"` only, the iteration when `mom_init` is
  changed to `mom_final`.

- mom_linear_weight:

  If `TRUE`, the gradient contribution to the update is weighted using
  momentum contribution.

- use_init_mom:

  If `TRUE`, then the momentum coefficient on the first iteration is
  non-zero. Otherwise, it's zero. Only applies if using a momentum
  schedule.

- restart:

  Momentum restart type. Can be one of "fn" or "gr". See 'Details' of
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md).

- restart_wait:

  Number of iterations to wait between restarts. Ignored if `restart` is
  `NULL`.

- par:

  (Optional) Initial values for the function to be optimized over.

- fg:

  (Optional). Function and gradient list. See 'Details' of
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

## Details

If the function to be optimized and starting point are not present at
creation time, then the optimizer should be initialized using
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
before being used with
[`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md).

See the documentation to
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) for an
explanation of all the parameters.

Details of the `fg` list containing the function to be optimized and its
gradient can be found in the 'Details' section of
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md). It is
optional for this function, but if it is passed to this function, along
with the vector of initial values, `par`, the optimizer will be returned
already initialized for this function. Otherwise,
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
must be called before optimization begins.

Additionally, optional convergence parameters may also be passed here,
for use with
[`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md).
They are optional here if you plan to call
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
later, or if you want to do your own convergence checking.

## Examples

``` r
# Function to optimize and starting point
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

# Create an optimizer and initialize it for use with the Rosenbrock function
opt <- make_mize(method = "L-BFGS", par = rb0, fg = rosenbrock_fg)

# Create optimizer without initialization
opt <- make_mize(method = "L-BFGS")

# Need to call mize_init separately:
opt <- mize_init(opt, rb0, rosenbrock_fg)
```
