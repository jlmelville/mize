# mize 0.2.3

A patch release to fix an incompatibility with R-devel.

## Bug fixes

* Fixed a bug where `class` was being checked directly and a scalar value was
assumed. The correct behavior is to use `methods::is`.

# mize 0.2.2

A patch release for a bug fix.

## Bug fixes

* Fixed a bug where if the maximum number of function evaluations for the
Schmidt line search is exceeded (controlled by the `ls_max_fn` parameter), a 
`'bracket_step' not found` error could result. Thank you to reporter Charles 
Driver.
* Fixed a couple of vignette links that were missing the "http://" at the front.

# mize 0.2.1

A patch release to fix an incompatibility with R-devel.

## Bug fixes

* Fixed an error with bold driver and back-tracking line search where the
stage was being incorrectly checked.

# mize 0.2.0

## New features

* New method: Truncated Newton (`method = "TN"`). Can be controlled using the
`tn_init` and `tn_exit` options.
* New method: SR1 (`method = "SR1"`), falling back to the BFGS direction if a 
descent direction is not found.
* New option `preconditioner`, which applies to the conjugate gradient and
truncated newton methods. The only value currently available is `preconditioner
= "L-BFGS"` which uses L-BFGS to estimate the inverse Hessian for
preconditioning. The number of updates to store for this preconditioner is
controlled by the `memory` parameter, just as if you were using `method =
"L-BFGS"`.
* BFGS, SR1, L-BFGS methods will now make use of a user-supplied inverse Hessian
function if provided. In the input `fg` list, supply a function `hi`, that takes
the `par` vector as input. The function can return a matrix (obviously not a
great idea for memory use), or a vector, the latter of which is assumed to be
the diagonal of the matrix.
* `ls_max_alpha` (for `line_search = "More-Thuente"` only): sets maximum value
of alpha that can be attained during line search.
* `ls_max_alpha_mult` (for Wolfe-type line search only): sets maximum value that
can be attained by the ratio of the initial guess for alpha for the current line
search, to the final value of alpha of the previous line search. Used to stop
line searches diverging due to very large initial guesses.
* `ls_safe_cubic` (for `line_search = "More-Thuente"` only): if `TRUE`,
use the safe-guarded cubic modification suggested by Xie and Schlick.
* `cg_update = "prfr"`, the "PR-FR" (Polak-Ribiere/Fletcher-Reeves) conjugate 
gradient update suggested by Gilbert and Nocedal.

## Bug fixes

* Error occurred when checking if a step size was finite during line search.
* DBD method didn't use momentum when asked to.
* Fix incorrectly specified conjugate gradient descent methods: 
Hestenes-Steifel (`cg_udpate = "hs"`), Conjugate Descent (`cg_udpate = "cd"`), 
Dai-Yuan (`cg_udpate = "dy"`) and Liu-Storey (`cg_udpate = "ls"`).

# mize 0.1.1

Initial release.
