# Convergence

An optimization can stop for three importantly different reasons:

| Outcome | `status` | `converged` | Typical termination reasons |
|:---|:---|:---|:---|
| A tolerance was satisfied | `"converged"` | `TRUE` | `abs_tol`, `rel_tol`, `grad_tol`, `ginf_tol`, `step_tol` |
| A hard limit was exhausted | `"budget_exhausted"` | `FALSE` | `max_iter`, `max_fn`, `max_gr`, `max_fg` |
| The optimization failed | `"failed"` | `FALSE` | a non-finite callback or a failed line search |

The result contract also reserves `status = "terminated"` for other
terminal reasons.

A limit often marks the intended end of a bounded run. It leaves
convergence unestablished. Inspect `status` and `converged` for the
broad outcome and `terminate$what` and `terminate$val` for the detailed
reason. This article shows how to choose among the controls and
interpret tolerance, budget, and failure outcomes.

## Choose stopping controls

Each control answers a different question:

| Control | What it detects | Main limitation | Quantity required |
|:---|:---|:---|:---|
| `abs_tol` / `rel_tol` | Little change in the objective | Stagnation need not mean stationarity | Objective |
| `grad_tol` / `ginf_tol` | A small gradient | Stationary points include saddles and maxima; absolute thresholds depend on scaling | Gradient |
| `step_tol` | Little movement in the parameters | May reflect poor scaling or a rejected step | Parameter update |
| `max_iter` | An outer iteration limit | Iterations can contain different amounts of work | None |
| `max_fn` / `max_gr` / `max_fg` | A callback cost limit | Says nothing about convergence | None |

Out of the box,
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) checks
for small absolute and relative changes in the objective and for a small
parameter update; all three tolerances use `sqrt(.Machine$double.eps)`.
Gradient-based stopping is disabled. A run is capped at 100 outer
iterations, while `max_fn`, `max_gr`, and `max_fg` do not impose finite
run-wide limits unless you set them. Those run-wide callback limits are
separate from the controls applied within an individual line search.

For a general smooth problem, enable a gradient criterion when its scale
is meaningful, retain callback or iteration limits as safety bounds, and
treat objective-change or step tolerances as evidence of stagnation
rather than proof of stationarity.

We’ll use the two-dimensional Rosenbrock function, which has a minimum
at `c(1, 1)` where the function equals zero.

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
```

## Iteration limit

`max_iter` bounds the number of optimizer iterations. Disable the
tolerances when the purpose of an example is to demonstrate that limit:

``` r

res <- mize(
  rb0,
  rb_fg,
  max_iter = 10,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>             status terminate iter nf ng objective
#> 1 budget_exhausted  max_iter   10 18 18    2.5526
```

Iteration counts are method-specific: one method may perform much more
work inside an iteration than another. For comparisons, report callback
counts and elapsed time as appropriate. Equal iteration limits do not
imply equal computational effort.

## Function-change tolerances

`abs_tol` compares consecutive function values and is reached when
`abs(fn_old - fn_new) < abs_tol`.

``` r

res <- mize(
  rb0,
  rb_fg,
  abs_tol = 1e-8,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>      status terminate iter nf ng   objective
#> 1 converged   abs_tol   38 50 50 4.60648e-11
```

At a [`mize()`](https://jlmelville.github.io/mize/reference/mize.md)
call, `rel_tol` defaults to the value supplied for `abs_tol`. Set
`rel_tol = NULL` explicitly when you want the absolute criterion alone;
similarly, set `abs_tol = NULL` when isolating relative change.

Relative change is calculated as
`abs(fn_old - fn_new) / min(abs(fn_old), abs(fn_new))`. If the
denominator is zero, the relative change is zero when the values are
equal and `Inf` otherwise.

Consecutive objective values can become similar because a run has
approached a solution or because its steps have become too small. The
following one-dimensional quadratic isolates that distinction. Its
minimizer is 10, and a deliberately small constant step produces little
relative change while the parameter remains far away:

``` r

slow_minimum <- 10
slow_fg <- list(
  fn = function(x) (x - slow_minimum)^2,
  gr = function(x) 2 * (x - slow_minimum)
)

relative_result <- mize(
  0,
  slow_fg,
  method = "SD",
  line_search = "constant",
  step0 = 1e-4,
  max_iter = 5,
  abs_tol = NULL,
  rel_tol = 1e-3,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
```

| Reason  | Relative change | Parameter | Distance from 10 |
|:--------|----------------:|----------:|-----------------:|
| rel_tol |           4e-04 |     0.004 |            9.996 |

The relative criterion has answered its specific question: recent
objective values changed little. The distance column shows why that
observation alone does not establish that the parameter is near the
minimizer. Pair it with a stationarity criterion when that is the
conclusion the run must support.

## Gradient tolerances

`grad_tol` tests the Euclidean (2-)norm of the current gradient.

``` r

res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = 1e-3,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>      status terminate iter nf ng   objective
#> 1 converged  grad_tol   37 49 49 4.52167e-10
```

A small gradient indicates approximate stationarity. It does not
distinguish a minimum from a saddle point or maximum, and an absolute
threshold depends on the scaling of both the parameters and the
objective. Gradient norms are still commonly used when comparing
optimizers; see, for example, [Nocedal and
co-workers](https://doi.org/10.1023/A:1014897230089).

`ginf_tol` instead tests the infinity norm: the largest absolute
component of the gradient. This is another common choice, particularly
for larger problems; see the [conjugate-gradient paper by Hager and
Zhang](https://doi.org/10.1137/030601880).

``` r

res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = 1e-3,
  step_tol = NULL
)
result_summary(res)
#>      status terminate iter nf ng   objective
#> 1 converged  ginf_tol   37 49 49 4.52167e-10
```

Checking either gradient tolerance requires a gradient at the current
position. A cached current gradient avoids an additional callback;
otherwise the check requests one. The method and line search determine
whether the value is already available.

## Step tolerance

`step_tol` detects when the Euclidean distance between consecutive
parameter vectors is small.

``` r

res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = 1e-4
)
result_summary(res)
#>      status terminate iter nf ng   objective
#> 1 converged  step_tol   38 50 50 4.60648e-11
```

Some methods can reject an iteration and restart from the same
parameters with different optimizer state. `step_tol` ignores that
restart. Because the criterion uses an absolute Euclidean distance,
rescaling one parameter can also change when it triggers even when the
underlying path describes the same problem.

## Evaluation budgets

When function or gradient callbacks dominate runtime, their counts can
be a useful proxy for work. The relative costs still depend on the
problem, and methods and line searches differ in which callbacks they
need and can reuse. Accepted Hessian and inverse-Hessian callbacks are
reported separately as `nh` and `nhi`; they are observational fields and
are not governed by `max_fn`, `max_gr`, or `max_fg`.

`max_fn` caps function evaluations, `max_gr` caps gradient evaluations,
and `max_fg` caps their combined total. The following runs use the same
problem and disable every ordinary tolerance so that each budget
supplies the terminal condition:

``` r

budget_names <- c("max_fn", "max_gr", "max_fg")

run_with_budget <- function(budget) {
  arguments <- list(
    par = rb0,
    fg = rb_fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  arguments[[budget]] <- 10
  do.call(mize, arguments)
}

budget_runs <- setNames(
  lapply(budget_names, run_with_budget),
  budget_names
)
```

| Control | Reason | Function calls | Gradient calls | Total calls |
|:--------|:-------|---------------:|---------------:|------------:|
| max_fn  | max_fn |             10 |             10 |          20 |
| max_gr  | max_gr |             10 |             10 |          20 |
| max_fg  | max_fg |              5 |              5 |          10 |

The separate limits allow ten evaluations of each callback in this
example; the combined limit stops after ten callbacks in total. The
table reports the actual work because an iteration can consume different
numbers of callbacks.

Each actual callback consumes the applicable budget; reuse of a cached
value does not. Once a callback budget is exhausted, `mize` makes no
later callback covered by that cap. This is why a final objective or
requested gradient diagnostic can be absent when its value was not
already available.

## When several criteria are active

All non-`NULL` criteria are active. `terminate$what` identifies the
terminal condition reported by the optimizer; it does not imply that
every other measure was far from its threshold. Hard callback limits are
checked before a covered callback is made, so the configured cap is
respected. When the final permitted callback itself establishes
convergence or produces a non-finite result, that observed tolerance or
failure is reported. The [Convergence section of the `mize()`
reference](https://jlmelville.github.io/mize/reference/mize.html#convergence)
gives the complete contract.

## Failure result

A non-finite callback produces a failure outcome. Here the objective is
finite, but the gradient callback returns `NaN`:

``` r

bad_fg <- list(
  fn = function(x) sum(x^2),
  gr = function(x) rep(NaN, length(x))
)

failure_result <- mize(c(1, 1), bad_fg)
failure_result[c("status", "message", "terminate")]
#> $status
#> [1] "failed"
#> 
#> $message
#> [1] "Failed: gr_inf encountered."
#> 
#> $terminate
#> $terminate$what
#> [1] "gr_inf"
#> 
#> $terminate$val
#> [1] Inf
```

The `gr_inf` reason identifies the non-finite gradient as the immediate
problem. Begin failure triage with `message` and `terminate`, then
inspect available progress or line-search diagnostics for the callback
and iteration that produced it.

## Checking tolerances less often

By default, `mize` checks tolerances every iteration. Function-change
tolerances need a function value at the current parameters, and gradient
tolerances need a current gradient. When the method or line search has
not already produced the needed value, checking can add a callback; when
it has, the cached value is reused.

`check_conv_every` changes the cadence of tolerance detection. For
example, this run tests the gradient criterion only every five
iterations:

``` r

res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = 1e-3,
  ginf_tol = NULL,
  step_tol = NULL,
  check_conv_every = 5
)
result_summary(res)
#>      status terminate iter nf ng objective
#> 1 converged  grad_tol   40 52 52 4.147e-18
```

A larger interval can delay recognition of a satisfied tolerance and
therefore overshoot it. Callback budgets remain hard: `max_fn`,
`max_gr`, and `max_fg` are enforced when a covered callback would be
made. Non-finite callback results and line-search no-step failures are
also handled independently of the ordinary tolerance cadence.

If `verbose = TRUE` or `store_progress = TRUE`, `log_every` controls how
often the corresponding summary is emitted or stored. It defaults to
`check_conv_every` and must be an integer multiple of it; otherwise
`mize` uses `check_conv_every` for both.

## Best and last results

`mize` returns explicit best-versus-last fields: `best_par` and `best_f`
name the best observed result, while `last_par` and `last_f` describe
the final optimizer state. The traditional `par` and `f` fields
currently match the best result.

Cached function or gradient information can make the best and last
results differ even with every tolerance disabled. A hard callback
budget can also prevent a final observation, leaving an objective field
absent or `last_f` as `NA`. The named fields identify the returned
points; stopping settings alone leave them ambiguous.

## See also

- [Getting started](https://jlmelville.github.io/mize/articles/mize.md)
  introduces status and result fields at the first optimization call.
- [Choosing methods and
  tuning](https://jlmelville.github.io/mize/articles/methods.md)
  explains method-specific work and line-search controls.
- [Stateful
  optimization](https://jlmelville.github.io/mize/articles/stateful.md)
  applies these termination controls in a caller-owned loop.
