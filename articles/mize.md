# Getting started

`mize` performs unconstrained numerical optimization when you can
provide an objective function and its gradient. This guide follows one
complete first-use path: define the problem, check the gradient, run the
optimizer, and interpret the result. It ends with pointers to the
specialist guides for convergence, stateful use, and method tuning.

## Define the objective and gradient

[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) receives
the objective and gradient together in a list. The `fn` function must
return one numeric value; `gr` must return a numeric vector with the
same length as its input. Here is the two-dimensional Rosenbrock
function:

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
```

For a more involved problem, a factory function can return this list so
that `fn` and `gr` share data through a closure. The list may also
contain an `fg` function that computes both values together when they
share expensive work. The [Metric
MDS](https://jlmelville.github.io/mize/articles/mmds.md) article
develops that pattern in a realistic example.

## Choose a starting point and check the gradient

The parameters are a numeric vector. Before optimizing, check an
analytic gradient at a representative, non-special starting point:

``` r

rb0 <- c(-1.2, 1)
gradient_check <- check_mize_gradient(rb_fg, rb0)

gradient_summary <- data.frame(
  "Maximum absolute error" = formatC(
    gradient_check$max_abs_error,
    format = "e",
    digits = 2
  ),
  "Maximum relative error" = formatC(
    gradient_check$max_rel_error,
    format = "e",
    digits = 2
  ),
  check.names = FALSE
)
knitr::kable(gradient_summary, align = c("r", "r"))
```

| Maximum absolute error | Maximum relative error |
|-----------------------:|-----------------------:|
|               2.21e-08 |               1.02e-10 |

[`check_mize_gradient()`](https://jlmelville.github.io/mize/reference/check_mize_gradient.md)
compares `gr` with a finite-difference approximation of `fn`. Large
errors are a strong reason to fix the gradient before trusting an
optimization. At this point, the analytic and finite-difference
gradients agree to approximately 1.0e-10 in relative terms. The full
result also identifies the worst coordinate and, when an optional
combined `fg` function is present, checks it against the separate
callbacks.

## Run the optimization

The default method is limited-memory BFGS (`"L-BFGS"`), a useful general
starting point for a smooth problem. Set `store_progress = TRUE` here so
that we can inspect the trajectory later:

``` r

res <- mize(rb0, rb_fg, store_progress = TRUE)
```

Start with the fitted parameters, objective, and high-level outcome:

``` r

res[c("par", "f", "status", "message")]
#> $par
#> [1] 0.9999972 0.9999923
#> 
#> $f
#> [1] 4.521673e-10
#> 
#> $status
#> [1] "converged"
#> 
#> $message
#> [1] "Converged: abs_tol reached."
```

The Rosenbrock function has a known minimum at `c(1, 1)`, where its
objective is zero. The returned parameters are close to that point, the
objective is close to zero, and `status = "converged"` says that a
configured numerical tolerance was reached.

### Read the termination record

`status` gives the quickest classification of any run:

- `"converged"` means a configured tolerance was reached;
- `"budget_exhausted"` means an iteration or callback limit was reached;
- `"failed"` identifies failures such as a non-finite callback or
  line-search failure; and
- `"terminated"` covers another recognized stopping reason.

`converged` is `TRUE` only for the tolerance-based class, while
`message` is a human-readable summary. Use the detailed record when code
needs the precise reason:

``` r

res[c("converged", "terminate")]
#> $converged
#> [1] TRUE
#> 
#> $terminate
#> $terminate$what
#> [1] "abs_tol"
#> 
#> $terminate$val
#> [1] 1.06609e-08
```

## Inspect progress

With `store_progress = TRUE`, `res$progress` contains the observations
and diagnostics available at stored iterations. Plotting the objective
against the cumulative number of function evaluations gives a quick view
of the run:

![Line plot of the Rosenbrock objective decreasing from about 24 to near
zero as cumulative function evaluations increase. The objective axis is
logarithmic.](mize_files/figure-html/progress-plot-1.png)

Rosenbrock objective values over cumulative function evaluations. The
vertical axis uses a logarithmic scale.

The stored objective falls by several orders of magnitude during this
run. For exact values, select the fields useful for the question and
inspect a few rows from the beginning and end:

``` r

progress_view <- res$progress[, c("f", "step", "nf", "ng"), drop = FALSE]
```

| Iteration | Objective |      Step | Function calls | Gradient calls |
|----------:|----------:|----------:|---------------:|---------------:|
|         0 |      24.2 |         0 |              1 |              0 |
|         1 |      19.5 |   0.02169 |              3 |              3 |
|         2 |     11.57 |   0.04729 |              4 |              4 |
|        35 | 2.304e-06 |   0.01136 |             47 |             47 |
|        36 | 1.111e-08 |  0.003086 |             48 |             48 |
|        37 | 4.522e-10 | 0.0002131 |             49 |             49 |

The initial row represents the starting point. Here `f` is the current
objective, `step` is the size of the parameter update, and `nf` and `ng`
are cumulative function and gradient callback counts. The available
columns depend on the method and calculations requested; see the
[Progress section of the `mize()`
reference](https://jlmelville.github.io/mize/reference/mize.html#progress)
for the complete dynamic schema.

For live console output, use `verbose = TRUE`; `log_every` controls how
often progress is printed or stored. Keeping every iteration can produce
a large result on a long run.

> **Returned points and hard budgets.**
> [`mize()`](https://jlmelville.github.io/mize/reference/mize.md)
> normally returns the best evaluated point. For nonmonotone runs,
> `best_*` and `last_*` distinguish that point from the final iterate.
> Under an exact callback budget, some objective fields may be
> unavailable. The
> [Convergence](https://jlmelville.github.io/mize/articles/convergence.md)
> article explains these cases, the status classes, and the stopping
> controls in detail.

## Where to go next

The default L-BFGS method is a useful baseline for a smooth problem.
Problem size and available derivative information may suggest another
method:

| Situation | Method to consider | Main trade-off |
|----|----|----|
| General smooth problem | `L-BFGS` | Limited-memory quasi-Newton; the default |
| Smaller smooth problem | `BFGS` | Dense inverse-Hessian approximation uses quadratic storage |
| Very tight memory budget | `CG` | Low storage; sensitive to scaling and line-search choices |
| Extra gradient work is acceptable | `TN` | Approximates a Newton step without forming a Hessian |
| Hessian or inverse Hessian is available | `NEWTON` | Requires an `hs` or `hi` callback in the function list |

Treat these as starting points and benchmark promising methods on the
problem at hand. Continue with the guide that matches your next task:

- [Convergence](https://jlmelville.github.io/mize/articles/convergence.md)
  explains tolerances, limits, failures, and result interpretation.
- [Stateful
  optimization](https://jlmelville.github.io/mize/articles/stateful.md)
  shows how to step, monitor, checkpoint, and resume an optimizer.
- [Metric MDS](https://jlmelville.github.io/mize/articles/mmds.md) is a
  complete case study using closures, a checked gradient, and combined
  callbacks.
- [Choosing methods and
  tuning](https://jlmelville.github.io/mize/articles/methods.md) covers
  algorithm selection and advanced controls, including line searches,
  step sizes, and momentum.
