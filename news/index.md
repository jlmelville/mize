# Changelog

## mize 0.3.0

This release improves optimizer robustness (mainly bug fixes when
pathological conditions are encountered) and observability.

### New features

- New function:
  [`check_mize_gradient()`](https://jlmelville.github.io/mize/reference/check_mize_gradient.md),
  which compares an analytic gradient with a finite-difference
  approximation.
- [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) now
  returns status fields: `converged`, `status`, and `message`, plus
  explicit best/last result fields (`best_par`, `best_f`, `last_par`,
  and `last_f`).
- [`mize()`](https://jlmelville.github.io/mize/reference/mize.md),
  [`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md),
  [`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md),
  and stored progress now report accepted Hessian and inverse-Hessian
  callback counts as `nh` and `nhi`.
- With `store_progress = TRUE`,
  [`mize()`](https://jlmelville.github.io/mize/reference/mize.md) now
  exposes optional line-search reason, selected-point provenance, local
  callback counts, initial scale, and exact-Newton direction provenance.

### Bug fixes and minor improvements

- DBD now returns structured `fn_inf` or `gr_inf` failures when a later
  objective or gradient is non-finite, validates its initial step
  controls, and safeguards named initial-step estimates using the
  current parameter scale.
- Momentum schedules now reject malformed or non-finite configured
  values and function results. The `"ramp"` and `"switch"` schedules
  require numeric `mom_init` and `mom_final` values, and adaptive
  restart safely rejects non-finite comparisons instead of raising a
  control-flow error.
- Stateful optimization now handles lifecycle dependencies consistently,
  and global function and gradient evaluation limits are enforced across
  optimizer and line-search callbacks, including truncated Newton inner
  iterations.
- Optimizer inputs and objective, gradient, Hessian, and inverse-Hessian
  callback results now receive consistent early validation and clearer
  errors.
- Quasi-Newton updates and exact-Newton directions now use safer
  fallbacks when curvature information or Hessian factorization is
  unsuitable.
- NEWTON and L-BFGS now preserve plain vector parameters across multiple
  iterations when a supplied inverse-Hessian callback returns a full
  matrix.
- Wolfe line searches now require an explicitly numeric `step0` to be a
  positive finite scalar; string initializers are unchanged.
- Line searches using weak Wolfe curvature now accept equality at the
  curvature boundary. The optional Hager-Zhang initializer probe now
  counts toward the line search’s local function and combined evaluation
  limits. Hager-Zhang initializer arithmetic also safely handles
  non-finite values and uses the specified Euclidean gradient norm.
- A line search that selects no usable step and produces no complete
  optimizer transition now reports `line_search_failed` instead of
  tolerance convergence.

## mize 0.2.5

CRAN release: 2026-02-09

Bug fix release.

### Bug fixes

- Removed `LazyData` from the `DESCRIPTION` file.
- Some flaky unit tests were removed.

## mize 0.2.4

CRAN release: 2020-08-30

Bug fix release.

### Bug fixes.

- If using `line_search = "backtracking"` with a specified `step_down`
  parameter, an incorrectly large number of gradient calculations was
  being reported.
- The documentation now specifies that if you *don’t* provide a
  `step_down` argument with `line_search = "backtracking"`,
  interpolation using function and gradient evaluations is carried out.
  To get a typical Armijo-style backtracking line search, specify a
  value for `step_down` (e.g. `step_down = 0.5` to halve the step size),
  and only function evaluations are used.

## mize 0.2.3

CRAN release: 2019-12-05

A patch release to fix an incompatibility with R-devel.

### Bug fixes

- Fixed a bug where `class` was being checked directly and a scalar
  value was assumed. The correct behavior is to use
  [`methods::is`](https://rdrr.io/r/methods/is.html).

## mize 0.2.2

CRAN release: 2019-07-11

A patch release for a bug fix.

### Bug fixes

- Fixed a bug where if the maximum number of function evaluations for
  the Schmidt line search is exceeded (controlled by the `ls_max_fn`
  parameter), a `'bracket_step' not found` error could result. Thank you
  to reporter Charles Driver.
- Fixed a couple of vignette links that were missing the “<http://>” at
  the front.

## mize 0.2.1

CRAN release: 2019-04-14

A patch release to fix an incompatibility with R-devel.

### Bug fixes

- Fixed an error with bold driver and back-tracking line search where
  the stage was being incorrectly checked.

## mize 0.2.0

CRAN release: 2018-09-14

### New features

- New method: Truncated Newton (`method = "TN"`). Can be controlled
  using the `tn_init` and `tn_exit` options.
- New method: SR1 (`method = "SR1"`), falling back to the BFGS direction
  if a descent direction is not found.
- New option `preconditioner`, which applies to the conjugate gradient
  and truncated newton methods. The only value currently available is
  `preconditioner = "L-BFGS"` which uses L-BFGS to estimate the inverse
  Hessian for preconditioning. The number of updates to store for this
  preconditioner is controlled by the `memory` parameter, just as if you
  were using `method = "L-BFGS"`.
- BFGS, SR1, L-BFGS methods will now make use of a user-supplied inverse
  Hessian function if provided. In the input `fg` list, supply a
  function `hi`, that takes the `par` vector as input. The function can
  return a matrix (obviously not a great idea for memory use), or a
  vector, the latter of which is assumed to be the diagonal of the
  matrix.
- `ls_max_alpha` (for `line_search = "More-Thuente"` only): sets maximum
  value of alpha that can be attained during line search.
- `ls_max_alpha_mult` (for Wolfe-type line search only): sets maximum
  value that can be attained by the ratio of the initial guess for alpha
  for the current line search, to the final value of alpha of the
  previous line search. Used to stop line searches diverging due to very
  large initial guesses.
- `ls_safe_cubic` (for `line_search = "More-Thuente"` only): if `TRUE`,
  use the safe-guarded cubic modification suggested by Xie and Schlick.
- `cg_update = "prfr"`, the “PR-FR” (Polak-Ribiere/Fletcher-Reeves)
  conjugate gradient update suggested by Gilbert and Nocedal.

### Bug fixes

- An error occurred when checking if a step size was finite during line
  search.
- DBD method didn’t use momentum when asked to.
- Fix incorrectly specified conjugate gradient descent methods:
  Hestenes-Stiefel (`cg_update = "hs"`), Conjugate Descent
  (`cg_update = "cd"`), Dai-Yuan (`cg_update = "dy"`) and Liu-Storey
  (`cg_update = "ls"`).

## mize 0.1.1

CRAN release: 2017-07-14

Initial release.
