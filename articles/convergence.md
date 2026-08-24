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
reason. This article concentrates on tolerances and limits. For a failed
run, start with `message` and the callback or line-search diagnostics.

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
equal and `Inf` otherwise. A large relative tolerance can declare
convergence after an iteration that made little progress:

``` r

res <- mize(
  rb0,
  rb_fg,
  abs_tol = NULL,
  rel_tol = 1e-2,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>      status terminate iter nf ng objective
#> 1 converged   rel_tol    5  7  7   4.13905
```

Here the result is still far from the known minimum even though the
relative criterion was satisfied. Efficient methods such as L-BFGS can
make little progress on an individual iteration, so choose `rel_tol`
with care.

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

The gradient is zero at a differentiable local minimum, even when the
function value is nonzero, and gradient norms are commonly used when
comparing optimizers. A small gradient can still mislead; see, for
example, [Nocedal and
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
restart.

## Evaluation budgets

When function or gradient callbacks dominate runtime, their counts can
be a useful proxy for work. The relative costs still depend on the
problem, and methods and line searches differ in which callbacks they
need and can reuse.

`max_fn` is a hard function-evaluation budget:

``` r

res <- mize(
  rb0,
  rb_fg,
  max_fn = 10,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>             status terminate iter nf ng objective
#> 1 budget_exhausted    max_fn    8 10 10   4.09761
```

`max_gr` applies the same rule to gradient evaluations:

``` r

res <- mize(
  rb0,
  rb_fg,
  max_gr = 10,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>             status terminate iter nf ng objective
#> 1 budget_exhausted    max_gr    8 10 10   4.09761
```

`max_fg` bounds the combined total:

``` r

res <- mize(
  rb0,
  rb_fg,
  max_fg = 10,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
result_summary(res)
#>             status terminate iter nf ng objective
#> 1 budget_exhausted    max_fg    3  5  5   4.28063
```

Each actual callback consumes the applicable budget; reuse of a cached
value does not. Once a callback budget is exhausted, `mize` makes no
later callback covered by that cap. This is why a final objective or
requested gradient diagnostic can be absent when its value was not
already available.

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
