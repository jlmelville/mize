# Numerical Optimization

Numerical optimization including conjugate gradient,
Broyden-Fletcher-Goldfarb-Shanno (BFGS), and the limited memory BFGS.

## Usage

``` r
mize(
  par,
  fg,
  method = "L-BFGS",
  norm_direction = FALSE,
  memory = 5,
  scale_hess = TRUE,
  cg_update = "PR+",
  preconditioner = "",
  tn_init = 0,
  tn_exit = "curvature",
  nest_q = 0,
  nest_convex_approx = FALSE,
  nest_burn_in = 0,
  step_up = 1.1,
  step_up_fun = "*",
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
  max_iter = 100,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  abs_tol = sqrt(.Machine$double.eps),
  rel_tol = abs_tol,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = sqrt(.Machine$double.eps),
  check_conv_every = 1,
  log_every = check_conv_every,
  verbose = FALSE,
  store_progress = FALSE
)
```

## Arguments

- par:

  Initial values for the function to be optimized over.

- fg:

  Function and gradient list. See 'Details'.

- method:

  Optimization method. See 'Details'.

- norm_direction:

  If `TRUE`, then the steepest descent direction is normalized to unit
  length. Useful for adaptive step size methods where the previous step
  size is used to initialize the next iteration.

- memory:

  The number of updates to store if using the `L-BFGS` method. Ignored
  otherwise. Must be a positive integer.

- scale_hess:

  if `TRUE`, the approximation to the inverse Hessian is scaled
  according to the method described by Nocedal and Wright (approximating
  an eigenvalue). Applies only to the methods `BFGS` (where the scaling
  is applied only during the first step) and `L-BFGS` (where the scaling
  is applied during every iteration). Ignored otherwise.

- cg_update:

  Type of update to use for the `"CG"` method. For details see the "CG"
  subsection of the "Optimization Methods" section. Ignored if `method`
  is not `"CG"`.

- preconditioner:

  Type of preconditioner to use in Truncated Newton. Leave blank or set
  to `"L-BFGS"` to use a limited memory BFGS preconditioner. Use the
  `"memory"` parameter to control the number of updates to store.
  Applies only if `method = "TN"` or `"CG"`, ignored otherwise.

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

  Strong convexity parameter for the NAG momentum term. Must take a
  value between 0 (strongly convex) and 1 (zero momentum). Only applies
  using the NAG method or a momentum method with Nesterov momentum
  schedule. Also does nothing if `nest_convex_approx` is `TRUE`.

- nest_convex_approx:

  If `TRUE`, then use an approximation due to Sutskever for calculating
  the momentum parameter in the NAG method. Only applies using the NAG
  method or a momentum method with Nesterov momentum schedule.

- nest_burn_in:

  Number of iterations to wait before using a non-zero momentum. Only
  applies using the NAG method or a momentum method with Nesterov
  momentum schedule.

- step_up:

  Value by which to increase the step size for the `"bold"` step size
  method or the `"DBD"` method.

- step_up_fun:

  Operator to use when combining the current step size with `step_up`.
  Can be one of `"*"` (to multiply the current step size with `step_up`)
  or `"+"` (to add).

- step_down:

  Multiplier to reduce the step size by if using the `"DBD"` method or
  the `"bold"` line search method. Should be a positive value less
  than 1. Also optional for use with the `"back"` line search method.

- dbd_weight:

  Weighting parameter used by the `"DBD"` method only, and only if no
  momentum scheme is provided. Must be an integer between 0 and 1.

- line_search:

  Type of line search to use. See 'Details'.

- c1:

  Sufficient decrease parameter for Wolfe-type line searches. Should be
  a value between 0 and 1.

- c2:

  Sufficient curvature parameter for line search for Wolfe-type line
  searches. Should be a value between `c1` and 1.

- step0:

  Initial value for the line search on the first step. See 'Details'.

- step_next_init:

  For Wolfe-type line searches only, how to initialize the line search
  on iterations after the first. See 'Details'.

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
  line search does not produce too small a value, using method of Xie
  and Schlick (2002). Only applies for `line_search = "more-thuente"`.

- strong_curvature:

  (Optional). If `TRUE` use the strong curvature condition in Wolfe line
  search. See the 'Line Search' section for details.

- approx_armijo:

  (Optional). If `TRUE` use the approximate Armijo condition in Wolfe
  line search. See the 'Line Search' section for details.

- mom_type:

  Momentum type, either `"classical"` or `"nesterov"`. See 'Details'.

- mom_schedule:

  Momentum schedule. See 'Details'.

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

  Momentum restart type. Can be one of "fn", "gr" or "speed". See
  'Details'. Ignored if no momentum scheme is being used.

- restart_wait:

  Number of iterations to wait between restarts. Ignored if `restart` is
  `NULL`.

- max_iter:

  Maximum number of iterations to optimize for. Defaults to 100. See the
  'Convergence' section for details.

- max_fn:

  Maximum number of function evaluations. See the 'Convergence' section
  for details.

- max_gr:

  Maximum number of gradient evaluations. See the 'Convergence' section
  for details.

- max_fg:

  Maximum number of function or gradient evaluations. See the
  'Convergence' section for details.

- abs_tol:

  Absolute tolerance for comparing two function evaluations. See the
  'Convergence' section for details.

- rel_tol:

  Relative tolerance for comparing two function evaluations. See the
  'Convergence' section for details.

- grad_tol:

  Absolute tolerance for the length (l2-norm) of the gradient vector.
  See the 'Convergence' section for details.

- ginf_tol:

  Absolute tolerance for the infinity norm (maximum absolute component)
  of the gradient vector. See the 'Convergence' section for details.

- step_tol:

  Absolute tolerance for the size of the parameter update. See the
  'Convergence' section for details.

- check_conv_every:

  Positive integer indicating how often to check convergence. Default is
  1, i.e. every iteration. See the 'Convergence' section for details.

- log_every:

  Positive integer indicating how often to log convergence results to
  the console. Ignored if `verbose` is `FALSE`. If not an integer
  multiple of `check_conv_every`, it will be set to `check_conv_every`.

- verbose:

  If `TRUE`, log information about the progress of the optimization to
  the console.

- store_progress:

  If `TRUE` store information about the progress of the optimization in
  a data frame, and include it as part of the return value.

## Value

A list with components:

- `par`: Optimized parameters. Normally, this is the best set of
  parameters seen during optimization, i.e. the set that produced the
  minimum function value. This requires that convergence checking with
  is carried out, including function evaluation where necessary. See the
  'Convergence' section for details.

- `nf`: Total number of function evaluations carried out. This includes
  any extra evaluations required for convergence calculations. Also, a
  function evaluation may be required to calculate the value of `f`
  returned in this list (see below). Additionally, if the `verbose`
  parameter is `TRUE`, then function and gradient information for the
  initial value of `par` will be logged to the console. These values are
  cached for subsequent use by the optimizer.

- `ng`: Total number of gradient evaluations carried out. This includes
  any extra evaluations required for convergence calculations using
  `grad_tol`. As with `nf`, additional gradient calculations beyond what
  you're expecting may have been needed for logging, convergence and
  calculating the value of `g2n` or `ginfn` (see below).

- `f`: Value of the function, evaluated at the returned value of `par`.

- `best_par`: The best parameters returned by `mize()`. This is
  currently the same value as `par`, and is provided so callers can use
  an explicit best-vs-last naming convention.

- `best_f`: Value of the function at `best_par`. This is currently the
  same value as `f`.

- `last_par`: Parameters from the last optimizer state before any final
  best-result restoration. This is the same as `par` unless `mize()`
  returns an earlier best point.

- `last_f`: Value of the function at `last_par` when it is known without
  requiring an extra function evaluation. If `last_par` differs from
  `par` and the function value was not already available, this is
  `NA_real_`.

- `g2n`: Optional: the length (Euclidean or l2-norm) of the gradient
  vector, evaluated at the returned value of `par`. Calculated only if
  `grad_tol` is non-null.

- `ginfn`: Optional: the infinity norm (maximum absolute component) of
  the gradient vector, evaluated at the returned value of `par`.
  Calculated only if `ginf_tol` is non-null.

- `iter`: The number of iterations the optimization was carried out for.

- `terminate`: List containing items: `what`, indicating what
  convergence criterion was met, and `val` specifying the value at
  convergence. See the 'Convergence' section for more details.

- `converged`: Logical value indicating whether `terminate$what` is one
  of the tolerance-based convergence criteria.

- `status`: Short string classifying the termination reason. One of
  `"converged"`, `"budget_exhausted"`, `"failed"`, or `"terminated"`.

- `message`: Human-readable summary of the termination reason.

- `progress`: Optional data frame containing information on the value of
  the function, gradient, momentum, and step sizes evaluated at each
  iteration where convergence is checked. Only present if
  `store_progress` is set to `TRUE`. Could get quite large if the
  optimization is long and the convergence is checked regularly.

## Details

The function to be optimized should be passed as a list to the `fg`
parameter. This should consist of:

- `fn`. The function to be optimized. Takes a vector of parameters and
  returns a scalar.

- `gr`. The gradient of the function. Takes a vector of parameters and
  returns a vector with the same length as the input parameter vector.

- `fg`. (Optional) function which calculates the function and gradient
  in the same routine. Takes a vector of parameters and returns a list
  containing the function result as `fn` and the gradient result as
  `gr`.

- `hs`. (Optional) Hessian of the function. Takes a vector of parameters
  and returns a square matrix with dimensions the same as the length of
  the input vector, containing the second derivatives. For `"NEWTON"`,
  it may also return a vector, which is treated as the diagonal of the
  Hessian. Required by `"PHESS"` and by `"NEWTON"` when no `hi` function
  is provided.

- `hi`. (Optional) inverse Hessian of the function. Takes a vector of
  parameters and returns either a square matrix with dimensions the same
  as the length of the input vector or a vector, which is treated as the
  diagonal of the inverse Hessian. Used by `"NEWTON"` when `hs` is not
  provided, by `"BFGS"` and `"SR1"` to initialize their inverse-Hessian
  approximation, and by `"L-BFGS"` as the initial inverse-Hessian
  approximation in its two-loop recursion.

The `fg` function is optional, but for some methods (e.g. line search
methods based on the Wolfe criteria), both the function and gradient
values are needed for the same parameter value. Calculating them in the
same function can save time if there is a lot of shared work.

## Optimization Methods

The `method` specifies the optimization method:

- `"SD"` is plain steepest descent. Not very effective on its own, but
  can be combined with various momentum approaches.

- `"BFGS"` is the Broyden-Fletcher-Goldfarb-Shanno quasi-Newton method.
  This stores an approximation to the inverse of the Hessian of the
  function being minimized, which requires storage proportional to the
  square of the length of `par`, so is unsuitable for large problems.

- `"SR1"` is the Symmetric Rank 1 quasi-Newton method, incorporating the
  safeguard given by Nocedal and Wright. Even with the safeguard, the
  SR1 method is not guaranteed to produce a descent direction. If this
  happens, the BFGS update is used for that iteration instead. Note that
  I have not done any research into the theoretical justification for
  combining BFGS with SR1 like this, but empirically (comparing to BFGS
  results with the datasets in the funconstrain package
  <https://github.com/jlmelville/funconstrain>), it works competitively
  with BFGS, particularly with a loose line search.

- `"L-BFGS"` is the Limited memory Broyden-Fletcher-Goldfarb-Shanno
  quasi-Newton method. This does not store the inverse Hessian
  approximation directly and so can scale to larger-sized problems than
  `"BFGS"`. The amount of memory used can be controlled with the
  `memory` parameter.

- `"CG"` is the conjugate gradient method. The `cg_update` parameter
  allows for different methods for choosing the next direction:

  - `"FR"` The method of Fletcher and Reeves.

  - `"PR"` The method of Polak and Ribiere.

  - `"PR+"` The method of Polak and Ribiere with a restart to steepest
    descent if conjugacy is lost. The default.

  - `"HS"` The method of Hestenes and Stiefel.

  - `"DY"` The method of Dai and Yuan.

  - `"HZ"` The method of Hager and Zhang.

  - `"HZ+"` The method of Hager and Zhang with restart, as used in
    CG_DESCENT.

  - `"PRFR"` The modified PR-FR method of Gilbert and Nocedal (1992).
    The `"PR+"` and `"HZ+"` are likely to be most robust in practice.
    Other updates are available more for curiosity purposes.

- `"TN"` is the Truncated Newton method, which approximately solves the
  Newton step without explicitly calculating the Hessian (at the expense
  of extra gradient calculations).

- `"NEWTON"` is Newton's method using a Hessian or inverse Hessian
  supplied in `fg`. It is supported, but is less commonly used than the
  quasi-Newton and truncated Newton methods.

- `"PHESS"` is a partial-Hessian Newton variant that reuses a Hessian
  factorization. It is accepted by `make_mize`, but should be treated as
  experimental/internal and may change without the same compatibility
  guarantees as the other methods.

- `"NAG"` is the Nesterov Accelerated Gradient method. The exact form of
  the momentum update in this method can be controlled with the
  following parameters:

  - `nest_q`: Strong convexity parameter. Must take a value between 0
    (strongly convex) and 1 (zero momentum). Ignored if
    `nest_convex_approx` is `TRUE`.

  - `nest_convex_approx`: If `TRUE`, then use an approximation due to
    Sutskever for calculating the momentum parameter.

  - `nest_burn_in`: Number of iterations to wait before using a non-zero
    momentum.

- `"DBD"` is the Delta-Bar-Delta method of Jacobs.

- `"Momentum"` is steepest descent with momentum. See below for momentum
  options.

For more details on gradient-based optimization in general, and the
BFGS, L-BFGS and CG methods, see Nocedal and Wright.

## Line Search

The parameter `line_search` determines the line search to be carried
out:

- `"Rasmussen"` carries out a line search using the strong Wolfe
  conditions as implemented by Carl Edward Rasmussen's minimize.m
  routines.

- `"More-Thuente"` carries out a line search using the strong Wolfe
  conditions and the method of More-Thuente. Can be abbreviated to
  `"MT"`.

- `"Schmidt"` carries out a line search using the strong Wolfe
  conditions as implemented in Mark Schmidt's minFunc routines.

- `"Backtracking"` carries out a back tracking line search using the
  sufficient decrease (Armijo) condition. By default, cubic
  interpolation using function and gradient values is used to find an
  acceptable step size. A constant step size reduction can be used by
  specifying a value for `step_down` between 0 and 1 (e.g. step size
  will be halved if `step_down` is set to `0.5`). If a constant step
  size reduction is used then only function evaluations are carried out
  and no extra gradient calculations are made.

- `"Bold Driver"` carries out a back tracking line search until a
  reduction in the function value is found.

- `"Constant"` uses a constant line search, the value of which should be
  provided with `step0`. Note that this value will be multiplied by the
  magnitude of the direction vector used in the gradient descent method.
  For method `"SD"` only, setting the `norm_direction` parameter to
  `TRUE` will scale the direction vector so it has unit length.

If using one of the methods: `"BFGS"`, `"L-BFGS"`, `"CG"` or `"NAG"`,
one of the Wolfe line searches: `"Rasmussen"` or `"More-Thuente"`,
`"Schmidt"` or `"Hager-Zhang"` should be used, otherwise very poor
performance is likely to be encountered. The following parameters can be
used to control the line search:

- `c1`: The sufficient decrease condition. Normally left at its default
  value of 1e-4.

- `c2`: The sufficient curvature condition. Defaults to 0.9 if using the
  methods `"BFGS"` and `"L-BFGS"`, and to 0.1 for every other method,
  more or less in line with the recommendations given by Nocedal and
  Wright. The smaller the value of `c2`, the stricter the line search,
  but it should not be set to smaller than `c1`.

- `step0`: Initial value for the line search on the first step. If a
  positive numeric value is passed as an argument, that value is used
  as-is. Otherwise, by passing a string as an argument, a guess is made
  based on values of the gradient, function or parameters, at the
  starting point:

  - `"rasmussen"`: As used by Rasmussen in `minimize.m`:
    \$\$\frac{1}{1+\left\|g\right\|^2}\$\$

  - `"scipy"`: As used in `optimize.py` in the python library Scipy.
    \$\$\frac{1}{\left\|g\right\|}\$\$

  - `"schmidt"`: As used by Schmidt in `minFunc.m` (the reciprocal of
    the l1 norm of g) \$\$\frac{1}{\left\|g\right\|\_1}\$\$

  - `"hz"`: The method suggested by Hager and Zhang (2006) for the
    CG_DESCENT software.

  These arguments can be abbreviated.

- `step_next_init`: How to initialize alpha value of subsequent line
  searches after the first, using results from the previous line search:

  - `"slope ratio"`: Slope ratio method.

  - `"quadratic"`: Quadratic interpolation method.

  - `"hz"`: The QuadStep method of Hager and Zhang (2006) for the
    CG_DESCENT software.

  - scalar numeric: Set to a numeric value (e.g. `step_next_init = 1`)
    to explicitly set alpha to this value initially.

  These arguments can be abbreviated. Details on the first two methods
  are provided by Nocedal and Wright.

- `try_newton_step`: For quasi-Newton methods (e.g. `"TN"`, `"BFGS"` and
  `"L-BFGS"`), setting this to `TRUE` will try the "natural" step size
  of 1, whenever the `step_next_init` method suggests an initial step
  size larger than that.

- `strong_curvature`: If `TRUE`, then the strong curvature condition
  will be used to check termination in Wolfe line search methods. If
  `FALSE`, then the standard curvature condition will be used. The
  default is `NULL` which lets the different Wolfe line searches choose
  whichever is their default behavior. This option is ignored if not
  using a Wolfe line search method.

- `approx_armijo`: If `TRUE`, then the approximate Armijo sufficient
  decrease condition (Hager and Zhang, 2005) will be used to check
  termination in Wolfe line search methods. If `FALSE`, then the exact
  curvature condition will be used. The default is `NULL` which lets the
  different Wolfe line searches choose whichever is their default
  behavior. This option is ignored if not using a Wolfe line search
  method.

For the Wolfe line searches, the methods of `"Rasmussen"`, `"Schmidt"`
and `"More-Thuente"` default to using the strong curvature condition and
the exact Armijo condition to terminate the line search (i.e. Strong
Wolfe conditions). The default step size initialization methods use the
Rasmussen method for the first iteration and quadratic interpolation for
subsequent iterations.

The `"Hager-Zhang"` Wolfe line search method defaults to the standard
curvature condition and the approximate Armijo condition (i.e.
approximate Wolfe conditions). The default step size initialization
methods are those used by Hager and Zhang (2006) in the description of
CG_DESCENT.

If the `"DBD"` is used for the optimization `"method"`, then the
`line_search` parameter is ignored, because this method controls both
the direction of the search and the step size simultaneously. The
following parameters can be used to control the step size:

- `step_up`: The amount by which to increase the step size in a
  direction where the current step size is deemed to be too short. This
  should be a positive scalar.

- `step_down`: The amount by which to decrease the step size in a
  direction where the currents step size is deemed to be too long. This
  should be a positive scalar smaller than 1. Default is 0.5.

- `step_up_fun`: How to increase the step size: either the method of
  Jacobs (addition of `step_up`) or Janet and co-workers (multiplication
  by `step_up`). Note that the step size decrease `step_down` is always
  a multiplication.

The `"bold driver"` line search also uses the `step_up` and `step_down`
parameters with similar meanings to their use with the `"DBD"` method:
the backtracking portion reduces the step size by a factor of
`step_down`. Once a satisfactory step size has been found, the line
search for the next iteration is initialized by multiplying the
previously found step size by `step_up`.

## Momentum

For `method` `"Momentum"`, momentum schemes can be accessed through the
momentum arguments:

- `mom_type`: Momentum type, either `"classical"` or `"nesterov"` (case
  insensitive, can be abbreviated). Using `"nesterov"` applies the
  momentum step before the gradient descent as suggested by Sutskever,
  emulating the behavior of the Nesterov Accelerated Gradient method.

- `mom_schedule`: How the momentum changes over the course of the
  optimization:

  - If a numerical scalar is provided, a constant momentum will be
    applied throughout.

  - `"nsconvex"`: Use the momentum schedule from the Nesterov
    Accelerated Gradient method suggested for non-strongly convex
    functions. Parameters which control the NAG momentum can also be
    used in combination with this option.

  - `"switch"`: Switch from one momentum value (specified via
    `mom_init`) to another (`mom_final`) at a a specified iteration
    (`mom_switch_iter`).

  - `"ramp"`: Linearly increase from one momentum value (`mom_init`) to
    another (`mom_final`).

  - If a function is provided, this will be invoked to provide a
    momentum value. It must take one argument (the current iteration
    number) and return a scalar.

  String arguments are case insensitive and can be abbreviated.

The `restart` parameter provides a way to restart the momentum if the
optimization appears to be not be making progress, inspired by the
method of O'Donoghue and Candes (2013) and Su and co-workers (2014).
There are three strategies:

- `"fn"`: A restart is applied if the function does not decrease on
  consecutive iterations.

- `"gr"`: A restart is applied if the direction of the optimization is
  not a descent direction.

- `"speed"`: A restart is applied if the update vector is not longer (as
  measured by Euclidean 2-norm) in consecutive iterations.

The effect of the restart is to "forget" any previous momentum update
vector, and, for those momentum schemes that change with iteration
number, to effectively reset the iteration number back to zero. If the
`mom_type` is `"nesterov"`, the gradient-based restart is not available.
The `restart_wait` parameter controls how many iterations to wait after
a restart, before allowing another restart. Must be a positive integer.
Default is 10, as used by Su and co-workers (2014). Setting this too low
could cause premature convergence. These methods were developed
specifically for the NAG method, but can be employed with any momentum
type and schedule.

If `method` type `"momentum"` is specified with no other values, the
momentum scheme will default to a constant value of `0.9`.

## Convergence

There are several ways for the optimization to terminate. The type of
termination is communicated by a two-item list `terminate` in the return
value, consisting of `what`, a short string describing what caused the
termination, and `val`, the value of the termination criterion that
caused termination.

The `converged`, `status`, and `message` return values summarize
`terminate` without replacing it. `converged` is `TRUE` only for
tolerance based termination (`"abs_tol"`, `"rel_tol"`, `"grad_tol"`,
`"ginf_tol"`, or `"step_tol"`). The `status` value is `"converged"` for
those tolerance exits, `"budget_exhausted"` for `"max_iter"`,
`"max_fn"`, `"max_gr"`, or `"max_fg"`, `"failed"` for `"fn_inf"` or
`"gr_inf"`, and `"terminated"` for any other termination reason.

The following parameters control various stopping criteria:

- `max_iter`: Maximum number of iterations to calculate. Reaching this
  limit is indicated by `terminate$what` being `"max_iter"`.

- `max_fn`: Maximum number of function evaluations allowed. Indicated by
  `terminate$what` being `"max_fn"`.

- `max_gr`: Maximum number of gradient evaluations allowed. Indicated by
  `terminate$what` being `"max_gr"`.

- `max_fg`: Maximum number of gradient evaluations allowed. Indicated by
  `terminate$what` being `"max_fg"`.

- `abs_tol`: Absolute tolerance of the function value. If the absolute
  value of the function falls below this threshold, `terminate$what`
  will be `"abs_tol"`. Will only be triggered if the objective function
  has a minimum value of zero.

- `rel_tol`: Relative tolerance of the function value, comparing
  consecutive function evaluation results. Indicated by `terminate$what`
  being `"rel_tol"`.

- `grad_tol`: Absolute tolerance of the l2 (Euclidean) norm of the
  gradient. Indicated by `terminate$what` being `"grad_tol"`. Note that
  the gradient norm is not a very reliable stopping criterion (see
  Nocedal and co-workers 2002), but is quite commonly used, so this
  might be useful for comparison with results from other optimization
  software.

- `ginf_tol`: Absolute tolerance of the infinity norm (maximum absolute
  component) of the gradient. Indicated by `terminate$what` being
  `"ginf_tol"`.

- `step_tol`: Absolute tolerance of the step size, i.e. the Euclidean
  distance between values of `par` fell below the specified value.
  Indicated by `terminate$what` being `"step_tol"`. For those
  optimization methods which allow for abandoning the result of an
  iteration and restarting using the previous iteration's value of `par`
  an iteration, `step_tol` will not be triggered.

Convergence is checked between specific iterations. How often is
determined by the `check_conv_every` parameter, which specifies the
number of iterations between each check. By default, this is set for
every iteration.

Be aware that if `abs_tol` or `rel_tol` are non-`NULL`, this requires
the function to have been evaluated at the current position at the end
of each iteration. If the function at that position has not been
calculated, it will be calculated and will contribute to the total
reported in the `counts` list in the return value. The calculated
function value is cached for use by the optimizer in the next iteration,
so if the optimizer would have needed to calculate the function anyway
(e.g. use of the strong Wolfe line search methods), there is no
significant cost accrued by calculating it earlier for convergence
calculations. However, for methods that don't use the function value at
that location, this could represent a lot of extra function evaluations.
On the other hand, not checking convergence could result in a lot of
extra unnecessary iterations. Similarly, if `grad_tol` or `ginf_tol` is
non-`NULL`, then the gradient will be calculated if needed.

If extra function or gradient evaluations is an issue, set
`check_conv_every` to a higher value, but be aware that this can cause
convergence limits to be exceeded by a greater amount.

Note also that if the `verbose` parameter is `TRUE`, then a summary of
the results so far will be logged to the console whenever a convergence
check is carried out. If the `store_progress` parameter is `TRUE`, then
the same information will be returned as a data frame in the return
value. For a long optimization this could be a lot of data, so by
default it is not stored.

Other ways for the optimization to terminate is if an iteration
generates a non-finite (i.e. `Inf` or `NaN`) gradient or function value.
Some, but not all, line-searches will try to recover from the latter, by
reducing the step size, but a non-finite gradient calculation during the
gradient descent portion of optimization is considered catastrophic by
mize, and it will give up. Termination under non-finite gradient or
function conditions will result in `terminate$what` being `"gr_inf"` or
`"fn_inf"` respectively. Unlike the convergence criteria, the
optimization will detect these error conditions and terminate even if a
convergence check would not be carried out for this iteration.

The value of `par` in the return value should be the parameters which
correspond to the lowest value of the function that has been calculated
during the optimization. As discussed above however, determining which
set of parameters requires a function evaluation at the end of each
iteration, which only happens if either the optimization method
calculates it as part of its own operation or if a convergence check is
being carried out during this iteration. Therefore, if your method does
not carry out function evaluations and `check_conv_every` is set to be
so large that no convergence calculation is carried out before
`max_iter` is reached, then the returned value of `par` is the last
value encountered.

## References

Gilbert, J. C., & Nocedal, J. (1992). Global convergence properties of
conjugate gradient methods for optimization. *SIAM Journal on
optimization*, *2*(1), 21-42.

Hager, W. W., & Zhang, H. (2005). A new conjugate gradient method with
guaranteed descent and an efficient line search. *SIAM Journal on
Optimization*, *16*(1), 170-192.

Hager, W. W., & Zhang, H. (2006). Algorithm 851: CG_DESCENT, a conjugate
gradient method with guaranteed descent. *ACM Transactions on
Mathematical Software (TOMS)*, *32*(1), 113-137.

Jacobs, R. A. (1988). Increased rates of convergence through learning
rate adaptation. *Neural networks*, *1*(4), 295-307.

Janet, J. A., Scoggins, S. M., Schultz, S. M., Snyder, W. E., White, M.
W., & Sutton, J. C. (1998, May). Shocking: An approach to stabilize
backprop training with greedy adaptive learning rates. In *1998 IEEE
International Joint Conference on Neural Networks Proceedings.* (Vol. 3,
pp. 2218-2223). IEEE.

Martens, J. (2010, June). Deep learning via Hessian-free optimization.
In *Proceedings of the International Conference on Machine Learning.*
(Vol. 27, pp. 735-742).

More', J. J., & Thuente, D. J. (1994). Line search algorithms with
guaranteed sufficient decrease. *ACM Transactions on Mathematical
Software (TOMS)*, *20*(3), 286-307.

Nocedal, J., Sartenaer, A., & Zhu, C. (2002). On the behavior of the
gradient norm in the steepest descent method. *Computational
Optimization and Applications*, *22*(1), 5-35.

Nocedal, J., & Wright, S. (2006). Numerical optimization. Springer
Science & Business Media.

O'Donoghue, B., & Candes, E. (2013). Adaptive restart for accelerated
gradient schemes. *Foundations of computational mathematics*, *15*(3),
715-732.

Schmidt, M. (2005). minFunc: unconstrained differentiable multivariate
optimization in Matlab.
<https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html>

Su, W., Boyd, S., & Candes, E. (2014). A differential equation for
modeling Nesterov's accelerated gradient method: theory and insights. In
*Advances in Neural Information Processing Systems* (pp. 2510-2518).

Sutskever, I. (2013). *Training recurrent neural networks* (Doctoral
dissertation, University of Toronto).

Sutskever, I., Martens, J., Dahl, G., & Hinton, G. (2013). On the
importance of initialization and momentum in deep learning. In
*Proceedings of the 30th international conference on machine learning
(ICML-13)* (pp. 1139-1147).

Xie, D., & Schlick, T. (1999). Remark on Algorithm 702 - The updated
truncated Newton minimization package. *ACM Transactions on Mathematical
Software (TOMS)*, *25*(1), 108-122.

Xie, D., & Schlick, T. (2002). A more lenient stopping rule for line
search algorithms. *Optimization Methods and Software*, *17*(4),
683-700.

## Examples

``` r
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

# Minimize using L-BFGS
res <- mize(rb0, rosenbrock_fg, method = "L-BFGS")

# Conjugate gradient with Fletcher-Reeves update, tight Wolfe line search
res <- mize(rb0, rosenbrock_fg, method = "CG", cg_update = "FR", c2 = 0.1)

# Steepest decent with constant momentum = 0.9
res <- mize(rb0, rosenbrock_fg, method = "MOM", mom_schedule = 0.9)

# Steepest descent with constant momentum in the Nesterov style as described
# in papers by Sutskever and Bengio
res <- mize(rb0, rosenbrock_fg,
  method = "MOM", mom_type = "nesterov",
  mom_schedule = 0.9
)

# Nesterov momentum with adaptive restart comparing function values
res <- mize(rb0, rosenbrock_fg,
  method = "MOM", mom_type = "nesterov",
  mom_schedule = 0.9, restart = "fn"
)
```
