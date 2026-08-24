# Stateful optimization

A stateful optimizer separates configuration from iteration. You keep
the current parameters and optimizer state, decide when to advance them,
and can log, visualize, apply a caller-owned stopping policy, or
checkpoint between steps.

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

> **Important.**
> [`mize_step()`](https://jlmelville.github.io/mize/reference/mize_step.md)
> advances the optimizer and may stop on a hard callback budget or an
> immediate method failure. A requested observation from
> [`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)
> can also stop at a hard callback budget. Ordinary numerical
> tolerances, `max_iter`, and non-finite summary observations are
> handled by
> [`check_mize_convergence()`](https://jlmelville.github.io/mize/reference/check_mize_convergence.md).
> A caller-owned loop must retain the optimizer returned by each
> operation, stop if either of the first two operations terminates it,
> and otherwise call the convergence check after every step.

The run below has disabled every ordinary numerical tolerance and
demonstrates only `max_iter`. Its summary therefore sets
`calc_fn = FALSE` and `calc_gr = FALSE` to avoid requesting observations
that the convergence check does not need. A loop that enables function-
or gradient-based tolerances should let
[`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)
request the corresponding observations.

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

## Work safely between steps

The optimizer’s stored state describes the current problem and
trajectory. Reading that state and deciding whether to continue leaves
it coherent. Changing the problem, the current parameters, or controls
that shape stored algorithm history begins a new run and should go
through
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md).

| Action | Examples | Treatment |
|:---|:---|:---|
| Observe or decide | Log or plot state, inspect diagnostics, save a checkpoint, update caller-owned records, apply an external stopping rule | Continue with the returned `opt` and `par`. |
| Change the problem | Replace callbacks or captured data, change dimension | Reinitialize. |
| Change the parameter state | Replace `par` independently, or change parameter meaning or scale | Reinitialize. |
| Change trajectory-owning controls | Change method; momentum or restart settings; L-BFGS memory; CG or TN preconditioning; line-search strategy or history-dependent step initialization | Reinitialize. |

Curvature pairs, cached values, momentum, and line-search history all
refer to the trajectory that produced them. Replacing part of that
trajectory while retaining the rest can make the state internally
inconsistent. Supply stopping controls through
[`make_mize()`](https://jlmelville.github.io/mize/reference/make_mize.md)
or
[`mize_init()`](https://jlmelville.github.io/mize/reference/mize_init.md);
this guide does not rely on direct mutation of `opt$convergence`.

## Deliberate observations and diagnostics

[`mize_step_summary()`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)
is where a loop can request observations. Set `calc_fn = TRUE` when
logging or best-value tracking needs the objective, and `calc_gr = TRUE`
when gradient norms are needed. The returned optimizer must be retained
because a requested callback updates the lifetime counts and can itself
terminate at a hard budget.

The summary schema is dynamic. A field is absent when its calculation
was not requested or its method did not produce it:

| Field group | Examples | Availability |
|:---|:---|:---|
| Core | `iter`, `nf`, `ng`, `step` | Generally available |
| Requested observations | `f`, `g2n`, `ginfn` | When calculated, cached, or requested |
| Line search | `alpha`, `ls_reason`, `ls_outcome`, `ls_nf`, `ls_ng` | Searches that report these diagnostics |
| Direction diagnostics | `direction_reason` | Methods that report direction provenance |
| Momentum | `mu` | Methods with a momentum stage |

Line searches may also report initialization diagnostics such as
`alpha_init` and `slope_init`. See
[`?mize_step_summary`](https://jlmelville.github.io/mize/reference/mize_step_summary.md)
for the complete field definitions, and test optional names before
application code relies on them.

The monitored loop follows the same lifecycle on every pass:

1.  Advance one step and retain its optimizer and parameters.
2.  Stop immediately if the step terminated.
3.  Request the observations needed by caller-owned records, and retain
    that optimizer too.
4.  Update the best value and append one typed progress record.
5.  Apply the ordinary convergence checks.

The starting objective is requested through the same budget-aware
summary API before the first step. The complete helper is available
below for copying.

Complete executable state helpers

``` r

numeric_or_na <- function(x) {
  if (is.null(x)) NA_real_ else as.numeric(x)
}

character_or_na <- function(x) {
  if (is.null(x)) NA_character_ else as.character(x)
}

progress_record <- function(step_info) {
  data.frame(
    iter = as.integer(step_info$iter),
    f = as.numeric(step_info$f),
    nf = as.integer(step_info$nf),
    ng = as.integer(step_info$ng),
    step = as.numeric(step_info$step),
    alpha = numeric_or_na(step_info$alpha),
    ls_outcome = character_or_na(step_info$ls_outcome)
  )
}

bind_progress <- function(records) {
  if (length(records) == 0L) {
    return(data.frame(
      iter = integer(),
      f = double(),
      nf = integer(),
      ng = integer(),
      step = double(),
      alpha = double(),
      ls_outcome = character()
    ))
  }
  do.call(rbind, records)
}

initialize_stateful <- function(opt, par, fg) {
  initial_info <- mize_step_summary(
    opt,
    par,
    fg,
    calc_fn = TRUE,
    calc_gr = FALSE
  )
  state <- list(
    opt = initial_info$opt,
    par = par,
    best_par = NULL,
    best_f = Inf,
    progress = list()
  )

  if (state$opt$is_terminated || is.null(initial_info$f)) {
    return(state)
  }
  if (!is.finite(initial_info$f)) {
    state$opt <- check_mize_convergence(initial_info)
    return(state)
  }

  state$best_par <- par
  state$best_f <- initial_info$f
  state$progress[[1L]] <- progress_record(initial_info)
  state
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

      next_record <- length(state$progress) + 1L
      state$progress[[next_record]] <- progress_record(step_info)
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
  initialize_stateful(opt, rb0, rb_fg)
}

initial_state <- new_bfgs_state()
initial_observation <- data.frame(
  `Starting objective` = initial_state$best_f,
  `Function calls` = initial_state$opt$counts$fn,
  `Progress records` = length(initial_state$progress),
  check.names = FALSE
)
knitr::kable(initial_observation)
```

| Starting objective | Function calls | Progress records |
|-------------------:|---------------:|-----------------:|
|               24.2 |              1 |                1 |

``` r


uninterrupted <- advance_stateful(initial_state, rb_fg)
uninterrupted_progress <- bind_progress(uninterrupted$progress)
knitr::kable(
  tail(uninterrupted_progress, 3),
  digits = 6,
  row.names = FALSE
)
```

| iter |        f |  nf |  ng |     step | alpha | ls_outcome |
|-----:|---------:|----:|----:|---------:|------:|:-----------|
|   28 | 0.026689 |  38 |  38 | 0.127223 |     1 | wolfe      |
|   29 | 0.015186 |  39 |  39 | 0.131488 |     1 | wolfe      |
|   30 | 0.005003 |  40 |  40 | 0.047100 |     1 | wolfe      |

``` r

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

For this BFGS configuration, initialization has not evaluated the
objective. The first summary therefore consumes one function call and
records the starting point as the caller-owned best. A repeated request
at the same point reuses the cached value. Progress remains a list of
typed one-row data frames until `bind_progress()` prepares it for
display; method-specific columns should still be tested for presence
before application code relies on them.

### Hard budgets can deny a requested observation

`max_fn` remains hard when `calc_fn = TRUE`. With no function
evaluations available, initialization records neither a best value nor a
progress row, and its updated optimizer reports the exhausted budget:

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
budget_state <- initialize_stateful(
  budget_opt,
  rb0,
  rb_fg
)
budget_result <- data.frame(
  `Best recorded` = !is.null(budget_state$best_par),
  `Function calls` = budget_state$opt$counts$fn,
  Status = budget_state$opt$status,
  Reason = budget_state$opt$terminate$what,
  check.names = FALSE
)
knitr::kable(budget_result)
```

| Best recorded | Function calls | Status           | Reason |
|:--------------|---------------:|:-----------------|:-------|
| FALSE         |              0 | budget_exhausted | max_fn |

This is why both initialization and the production loop retain the
optimizer returned by a summary and guard the requested value. Calling
`fg$fn()` directly bypasses `mize`’s counters and limits. Keep direct
callback calls outside hard-budget workflows.

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
checkpoint_progress_is_list <- is.list(resumed$progress)
checkpoint_records_are_typed <- all(vapply(
  resumed$progress,
  is.data.frame,
  logical(1)
))
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
knitr::kable(
  comparison,
  digits = 8,
  col.names = c(
    "Run",
    "Objective",
    "Function calls",
    "Gradient calls",
    "Status",
    "Reason"
  )
)
```

| Run | Objective | Function calls | Gradient calls | Status | Reason |
|:---|---:|---:|---:|:---|:---|
| one-shot | 0.00500283 | 40 | 40 | budget_exhausted | max_iter |
| stateful | 0.00500283 | 40 | 40 | budget_exhausted | max_iter |
| checkpoint/resume | 0.00500283 | 40 | 40 | budget_exhausted | max_iter |

``` r


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
```

The resumed run reaches the same current and best parameters as the
uninterrupted run to numerical tolerance, and both agree with the
one-shot run. This is the practical guarantee supplied by checkpointing
the complete optimizer state; exact byte-for-byte identity across
different R or platform builds is not required.

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
