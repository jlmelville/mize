state_hook_quadratic <- function() {
  list(
    fn = function(x) {
      sum(x^2)
    },
    gr = function(x) {
      2 * x
    }
  )
}

register_test_hooks <- function(opt, hooks) {
  for (hook in hooks) {
    opt <- register_hook(opt, hook)
  }
  opt
}

make_trace_opt_hook <- function(name, event) {
  hook <- function(opt, par, fg, iter, ...) {
    opt$trace <- c(opt$trace, name)
    opt
  }
  attr(hook, "name") <- name
  attr(hook, "event") <- event
  hook
}

make_trace_stage_hook <- function(name, event) {
  hook <- function(opt, stage, par, fg, iter, ...) {
    opt$trace <- c(opt$trace, name)
    list(opt = opt)
  }
  attr(hook, "name") <- name
  attr(hook, "event") <- event
  hook
}

test_that("stateful steps run lifecycle hooks in optimizer order", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )
  opt <- register_test_hooks(
    opt,
    list(
      make_trace_opt_hook("before_step", "before step"),
      make_trace_stage_hook(
        "before_gradient_descent",
        "before gradient_descent"
      ),
      make_trace_stage_hook("after_stage", "after stage"),
      make_trace_opt_hook("before_validation", "before validation"),
      make_trace_opt_hook("during_validation", "during validation"),
      make_trace_opt_hook("after_step", "after step")
    )
  )

  res <- mize_step(opt, c(1), fg)

  expect_equal(
    res$opt$trace,
    c(
      "before_step",
      "before_gradient_descent",
      "after_stage",
      "before_validation",
      "during_validation",
      "after_step"
    )
  )
  expect_equal(res$par, 0.8)
})

test_that("stage and sub-stage hooks persist returned state mutations", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  mark_stage <- function(opt, stage, par, fg, iter, ...) {
    stage$marked_by_hook <- iter
    list(stage = stage)
  }
  attr(mark_stage, "name") <- "mark_stage"
  attr(mark_stage, "event") <- "before gradient_descent"

  double_direction <- function(opt, stage, sub_stage, par, fg, iter, ...) {
    sub_stage$value <- 2 * sub_stage$value
    sub_stage$marked_by_hook <- iter
    list(sub_stage = sub_stage)
  }
  attr(double_direction, "name") <- "double_direction"
  attr(double_direction, "event") <- "after gradient_descent direction"

  opt <- register_test_hooks(opt, list(mark_stage, double_direction))

  res <- mize_step(opt, c(1), fg)
  stage <- res$opt$stages$gradient_descent

  expect_equal(res$par, 0.6)
  expect_equal(stage$marked_by_hook, 1)
  expect_equal(stage$direction$marked_by_hook, 1)
  expect_equal(stage$direction$value, -4)
  expect_equal(stage$result, -0.4)
})

test_that("termination from a stage hook short-circuits later lifecycle hooks", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  terminate_after_stage <- function(opt, stage, par, fg, iter, ...) {
    opt$trace <- c(opt$trace, "after_stage")
    opt <- set_mize_termination(
      opt,
      list(what = "hook_stop", val = iter)
    )
    list(opt = opt)
  }
  attr(terminate_after_stage, "name") <- "terminate_after_stage"
  attr(terminate_after_stage, "event") <- "after stage"

  opt <- register_test_hooks(
    opt,
    list(
      terminate_after_stage,
      make_trace_opt_hook("before_validation", "before validation"),
      make_trace_opt_hook("after_step", "after step")
    )
  )

  res <- mize_step(opt, c(1), fg)

  expect_equal(res$par, 1)
  expect_equal(res$opt$trace, "after_stage")
  expect_equal(res$opt$terminate$what, "hook_stop")
  expect_equal(res$opt$status, "terminated")
})

test_that("validation veto rolls back parameters before after-step hooks", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  veto_validation <- function(opt, par, fg, iter, par0, update) {
    opt$trace <- c(opt$trace, "veto")
    opt$ok <- FALSE
    opt
  }
  attr(veto_validation, "name") <- "veto_validation"
  attr(veto_validation, "event") <- "during validation"

  after_step_probe <- function(opt, par, fg, iter, par0, update) {
    opt$trace <- c(
      opt$trace,
      "after_step",
      paste0("after_par=", par[1]),
      paste0("update=", update[1]),
      paste0("ok=", opt$ok)
    )
    opt
  }
  attr(after_step_probe, "name") <- "after_step_probe"
  attr(after_step_probe, "event") <- "after step"

  opt <- register_test_hooks(opt, list(veto_validation, after_step_probe))

  res <- mize_step(opt, c(1), fg)

  expect_equal(res$par, 1)
  expect_equal(
    res$opt$trace,
    c("veto", "after_step", "after_par=1", "update=-0.2", "ok=FALSE")
  )
  expect_false(res$opt$ok)
})
