# Stateful optimization

A stateful optimizer separates configuration from iteration. You keep
the current parameters and optimizer state, decide when to advance them,
and can log, visualize, intervene, or checkpoint between steps.

This is useful when a single blocking call to
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) does not
leave enough room for application-specific work. The price of that
control is that your loop owns the lifecycle: it must retain every
updated optimizer, check terminal state, and treat function and gradient
observations as optional.

We will use the two-dimensional Rosenbrock function throughout:

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

## Why retaining state matters

A tempting workaround with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) is to run a small
batch, take its parameters, and start another optimizer:

``` r

first_batch <- stats::optim(
  par = rb0,
  fn = rb_fg$fn,
  gr = rb_fg$gr,
  method = "BFGS",
  control = list(maxit = 10)
)
second_batch <- stats::optim(
  par = first_batch$par,
  fn = rb_fg$fn,
  gr = rb_fg$gr,
  method = "BFGS",
  control = list(maxit = 10)
)
c(first = first_batch$value, second = second_batch$value)
#>     first    second 
#> 1.3673831 0.2010733
```

The second call starts a new BFGS run. It receives the parameters and
discards the first run’s inverse-Hessian approximation and other
transient method state. Shorter batches provide more interruptions while
throwing away more of the information that makes a stateful method
effective.

## The three lifecycle operations

The stateful API has three main operations:

| Operation | Role |
|:---|:---|
| [`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md) | Configure an optimizer; optionally initialize it when `par` and `fg` are available. |
| [`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md) | Bind an optimizer to a starting point and callbacks, and initialize a new run. |
| [`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md) | Advance one iteration and return updated `opt` and `par`. |

Configuration and initialization can be separate:

``` r

opt <- make_mize(method = "BFGS")
opt <- mize_init(
  opt,
  par = rb0,
  fg = rb_fg,
  max_iter = 3,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
```

Or pass `par` and `fg` to
[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md)
to initialize immediately. In either form, always use the updated `opt`
and `par` returned by each step.

## Optional observations

[`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md)
always returns `opt`, `par`, `nf`, and `ng`. It returns `f` or `g` only
when that value was calculated at the returned parameters. For example,
constant-step steepest descent needs a gradient at the old parameters
but does not observe either value at the new parameters:

``` r

observation_opt <- make_mize(
  method = "SD",
  line_search = "constant",
  step0 = 0.0001,
  par = rb0,
  fg = rb_fg,
  max_iter = 3,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
observation_step <- mize_step(observation_opt, rb0, rb_fg)
names(observation_step)
#> [1] "opt" "par" "nf"  "ng"
```

Code such as `if (step$f < best_f)` is therefore unsafe. An absent
observation is normal and leaves the step status unchanged.

## A minimal terminal-aware loop

A stateful run does not check its ordinary convergence controls until
you summarize a step and call
[`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md).
A minimal safe loop retains state and checks for termination after every
call that can update the optimizer:

``` r

par <- rb0

while (!opt$is_terminated) {
  par_old <- par
  step <- mize_step(opt, par, rb_fg)
  opt <- step$opt
  par <- step$par
  if (opt$is_terminated) {
    break
  }

  step_info <- mize_step_summary(
    opt,
    par,
    rb_fg,
    par_old = par_old,
    calc_fn = FALSE,
    calc_gr = FALSE
  )
  opt <- step_info$opt
  if (opt$is_terminated) {
    break
  }

  opt <- check_mize_convergence(step_info)
}

opt[c("status", "converged", "message", "terminate")]
#> $status
#> [1] "budget_exhausted"
#> 
#> $converged
#> [1] FALSE
#> 
#> $message
#> [1] "Budget exhausted: max_iter reached."
#> 
#> $terminate
#> $terminate$what
#> [1] "max_iter"
#> 
#> $terminate$val
#> [1] 3
```

This run deliberately ends at `max_iter = 3`: its status is
`budget_exhausted`, and `converged` is `FALSE`.

## Deliberate observations and diagnostics

[`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)
is where a loop can request observations. Set `calc_fn = TRUE` when
logging or best-value tracking needs the objective, and `calc_gr = TRUE`
when gradient norms are needed. The returned optimizer must be retained
because a requested callback updates the lifetime counts and can itself
terminate at a hard budget.

The summary schema is dynamic. Common fields include `iter`, `f`, `nf`,
`ng`, `step`, and `alpha`; line-search methods can also expose
`alpha_init`, `slope_init`, `ls_reason`, `ls_outcome`, `ls_nf`, and
`ls_ng`. Some directions add `direction_reason`, and momentum methods
can add `mu`. See
[`?mize_step_summary`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)
for the complete current schema.

Here is one canonical production loop. Its `state` list includes the
required `opt` and `par`, plus caller-owned best and progress records.
The best value is the best post-step objective this loop has actually
observed; it is distinct from the optimizer’s current state.

``` r

value_or_na <- function(x) {
  if (is.null(x)) NA else x
}

advance_stateful <- function(state, fg, max_steps = Inf) {
  steps <- 0L

  while (!state$opt$is_terminated && steps < max_steps) {
    par_old <- state$par
    step <- mize_step(state$opt, state$par, fg)
    state$opt <- step$opt
    state$par <- step$par
    if (state$opt$is_terminated) {
      break
    }

    step_info <- mize_step_summary(
      state$opt,
      state$par,
      fg,
      par_old = par_old,
      calc_fn = TRUE
    )
    state$opt <- step_info$opt
    if (state$opt$is_terminated) {
      break
    }

    if (!is.null(step_info$f) && is.finite(step_info$f)) {
      if (step_info$f < state$best_f) {
        state$best_f <- step_info$f
        state$best_par <- state$par
      }

      state$progress <- rbind(
        state$progress,
        data.frame(
          iter = step_info$iter,
          f = step_info$f,
          nf = step_info$nf,
          ng = step_info$ng,
          step = step_info$step,
          alpha = value_or_na(step_info$alpha),
          ls_outcome = value_or_na(step_info$ls_outcome)
        )
      )
    }

    state$opt <- check_mize_convergence(step_info)
    steps <- steps + 1L
  }

  state
}
```

We will use a bounded BFGS run so the terminal classification is
unambiguous:

``` r

controls <- list(
  max_iter = 30,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)

new_bfgs_state <- function() {
  opt <- do.call(
    make_mize,
    c(list(method = "BFGS", par = rb0, fg = rb_fg), controls)
  )
  list(
    opt = opt,
    par = rb0,
    best_par = NULL,
    best_f = Inf,
    progress = data.frame()
  )
}

uninterrupted <- advance_stateful(new_bfgs_state(), rb_fg)
tail(uninterrupted$progress, 3)
#>    iter           f nf ng       step alpha ls_outcome
#> 28   28 0.026689274 38 38 0.12722348     1      wolfe
#> 29   29 0.015185570 39 39 0.13148806     1      wolfe
#> 30   30 0.005002828 40 40 0.04709997     1      wolfe
uninterrupted$opt[c("status", "converged", "message", "terminate")]
#> $status
#> [1] "budget_exhausted"
#> 
#> $converged
#> [1] FALSE
#> 
#> $message
#> [1] "Budget exhausted: max_iter reached."
#> 
#> $terminate
#> $terminate$what
#> [1] "max_iter"
#> 
#> $terminate$val
#> [1] 30
```

The selected progress columns are useful for this run. Test
method-specific columns for presence before application code relies on
them.

### Hard budgets can deny a requested observation

`max_fn` remains hard when `calc_fn = TRUE`. With no function
evaluations available, the summary returns no `f` and its updated
optimizer reports the exhausted budget:

``` r

budget_opt <- make_mize(
  method = "SD",
  line_search = "constant",
  step0 = 0.0001,
  par = rb0,
  fg = rb_fg,
  max_fn = 0,
  abs_tol = NULL,
  rel_tol = NULL,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
)
budget_info <- mize_step_summary(
  budget_opt,
  rb0,
  rb_fg,
  calc_fn = TRUE
)
list(
  f_available = !is.null(budget_info$f),
  nf = budget_info$nf,
  status = budget_info$opt$status,
  terminate = budget_info$opt$terminate
)
#> $f_available
#> [1] FALSE
#> 
#> $nf
#> [1] 0
#> 
#> $status
#> [1] "budget_exhausted"
#> 
#> $terminate
#> $terminate$what
#> [1] "max_fn"
#> 
#> $terminate$val
#> [1] 0
```

This is why the production loop both retains `step_info$opt` and guards
`step_info$f`. Calling `fg$fn()` directly bypasses `mize`’s counters and
limits. Keep direct callback calls outside hard-budget workflows.

## Checkpoint and resume

Checkpoint the complete caller-owned state needed to continue. Here that
means the optimizer and parameters, plus the best observation and
progress log. The callbacks are ordinary program code and are supplied
again after restoration.

``` r

partial <- advance_stateful(
  new_bfgs_state(),
  rb_fg,
  max_steps = 10
)

checkpoint_file <- tempfile(fileext = ".rds")
saveRDS(partial, checkpoint_file)
rm(partial)

resumed <- readRDS(checkpoint_file)
unlink(checkpoint_file)
resumed <- advance_stateful(resumed, rb_fg)
```

To show what preservation buys us, compare the resumed run with both an
uninterrupted stateful run and a one-shot
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) call
using the same method and termination controls:

``` r

one_shot <- do.call(
  mize,
  c(list(par = rb0, fg = rb_fg, method = "BFGS"), controls)
)

comparison <- data.frame(
  run = c("one-shot", "stateful", "checkpoint/resume"),
  objective = c(one_shot$f, uninterrupted$best_f, resumed$best_f),
  nf = c(one_shot$nf, uninterrupted$opt$counts$fn, resumed$opt$counts$fn),
  ng = c(one_shot$ng, uninterrupted$opt$counts$gr, resumed$opt$counts$gr),
  status = c(
    one_shot$status,
    uninterrupted$opt$status,
    resumed$opt$status
  ),
  terminate = c(
    one_shot$terminate$what,
    uninterrupted$opt$terminate$what,
    resumed$opt$terminate$what
  )
)
comparison
#>                 run   objective nf ng           status terminate
#> 1          one-shot 0.005002828 40 40 budget_exhausted  max_iter
#> 2          stateful 0.005002828 40 40 budget_exhausted  max_iter
#> 3 checkpoint/resume 0.005002828 40 40 budget_exhausted  max_iter

comparison_tolerance <- 1e-10
equivalent <- c(
  one_shot_best = isTRUE(all.equal(
    one_shot$par,
    uninterrupted$best_par,
    tolerance = comparison_tolerance
  )),
  one_shot_objective = isTRUE(all.equal(
    one_shot$f,
    uninterrupted$best_f,
    tolerance = comparison_tolerance
  )),
  one_shot_last = isTRUE(all.equal(
    one_shot$last_par,
    uninterrupted$par,
    tolerance = comparison_tolerance
  )),
  resumed_current = isTRUE(all.equal(
    resumed$par,
    uninterrupted$par,
    tolerance = comparison_tolerance
  )),
  resumed_best_par = isTRUE(all.equal(
    resumed$best_par,
    uninterrupted$best_par,
    tolerance = comparison_tolerance
  )),
  resumed_best = isTRUE(all.equal(
    resumed$best_f,
    uninterrupted$best_f,
    tolerance = comparison_tolerance
  ))
)
equivalent
#>      one_shot_best one_shot_objective      one_shot_last    resumed_current 
#>               TRUE               TRUE               TRUE               TRUE 
#>   resumed_best_par       resumed_best 
#>               TRUE               TRUE
```

For this deterministic problem and implementation, all comparisons pass
at `1e-10`. The comparison shows that serializing and restoring the
optimizer preserves method state. Bit-for-bit identity across platforms
remains outside the contract.

## Reinitialization starts a new run

Calling
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md)
on an existing optimizer deliberately starts a new run. It resets the
iteration counter, terminal fields, caches, and transient algorithm
state. Function and gradient counts remain optimizer-lifetime totals:

``` r

counts_before <- uninterrupted$opt$counts
reinitialized <- mize_init(
  uninterrupted$opt,
  rb0,
  rb_fg
)

data.frame(
  iter_before = uninterrupted$opt$iter,
  iter_after = reinitialized$iter,
  counts_preserved = identical(counts_before, reinitialized$counts),
  terminated_after = reinitialized$is_terminated
)
#>   iter_before iter_after counts_preserved terminated_after
#> 1          30          0             TRUE            FALSE
```

Use reinitialization when a new run is intended; use serialization when
the same run must continue. See
[`?mize_init`](https://jlmelville.github.io/mize/reference/mize_init.md),
[`?mize_step`](https://jlmelville.github.io/mize/reference/mize_step.md),
[`?mize_step_summary`](https://jlmelville.github.io/mize/reference/mize_step_summary.md),
and
[`?check_mize_convergence`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md)
for the complete API.

## See also

- [Getting started](https://jlmelville.github.io/mize/articles/mize.md)
  introduces the one-shot workflow and result contract.
- [Convergence](https://jlmelville.github.io/mize/articles/convergence.md)
  explains how to choose and interpret termination controls.
- [Metric MDS](https://jlmelville.github.io/mize/articles/mmds.md)
  develops a callback closure suitable for either one-shot or stateful
  use.
