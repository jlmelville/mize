#' Mize Step Summary
#'
#' Produces a result summary for an optimization iteration. Information such as
#' function value, gradient norm and step size may be returned.
#'
#' By default, convergence tolerance parameters will be used to determine what
#' function and gradient data is returned. The function value will be returned if
#' it was already calculated and cached in the optimization iteration. Otherwise,
#' it will be calculated only if a non-null absolute or relative tolerance value
#' was asked for. A gradient norm will be returned only if a non-null gradient
#' tolerance was specified, even if the gradient is available.
#'
#' Note that if a function tolerance was specified, but was not calculated for
#' the relevant value of `par`, they will be calculated here and the
#' calculation does contribute to the total function count (and will be cached
#' for potential use in the next iteration). The same applies for gradient
#' tolerances and gradient calculation. Function and gradient calculation can
#' also be forced here by setting the `calc_fn` and `calc_gr`
#' (respectively) parameters to `TRUE`.
#'
#' @param opt Optimizer to generate summary for, from return value of
#'  [mize_step()].
#' @param par Vector of parameters at the end of the iteration, from return value
#'  of [mize_step()].
#' @param fg Function and gradient list. See the documentation of
#'  [mize()].
#' @param par_old (Optional). Vector of parameters at the end of the previous
#'  iteration. Used to calculate step size.
#' @param calc_fn (Optional). If `TRUE`, force calculation of function if
#'  not already cached in `opt`, even if it would not be needed for
#'  convergence checking.
#' @param calc_gr (Optional). If `TRUE`, force calculation of gradient if
#'  not already cached in `opt`, even if it would not be needed for
#'  convergence checking.
#' @return A list with the following items:
#'
#' * `opt`: Optimizer with updated state (e.g. function and gradient
#'  counts).
#' * `iter`: Iteration number.
#' * `f`: Function value at `par`.
#' * `g2n`: 2-norm of the gradient at `par`.
#' * `ginfn`: Infinity-norm of the gradient at `par`.
#' * `nf`: Number of function evaluations so far.
#' * `ng`: Number of gradient evaluations so far.
#' * `step`: Size of the step between `par_old` and `par`,
#'  if `par_old` is provided.
#' * `alpha`: Step length of the gradient descent part of the step.
#' * `mu`: Momentum coefficient for this iteration.
#' @export
#' @examples
#' rb_fg <- list(
#'   fn = function(x) {
#'     100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
#'   },
#'   gr = function(x) {
#'     c(
#'       -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'       200 * (x[2] - x[1] * x[1])
#'     )
#'   }
#' )
#' rb0 <- c(-1.2, 1)
#'
#' opt <- make_mize(method = "BFGS", par = rb0, fg = rb_fg, max_iter = 30)
#' mize_res <- mize_step(opt = opt, par = rb0, fg = rb_fg)
#' # Get info about first step, use rb0 to compare new par with initial value
#' step_info <- mize_step_summary(mize_res$opt, mize_res$par, rb_fg, rb0)
mize_step_summary <- function(
  opt,
  par,
  fg,
  par_old = NULL,
  calc_fn = NULL,
  calc_gr = NULL
) {
  iter <- opt$iter
  # An internal flag useful for unit tests: if FALSE, doesn't count any
  # fn/gr calculations towards their counts. Can still get back fn/gr values
  # without confusing the issue of the expected number of fn/gr evaluations
  if (!is.null(opt$count_res_fg)) {
    count_fg <- opt$count_res_fg
  } else {
    count_fg <- TRUE
  }

  # Whether and what convergence info to calculate if fn/gr calculation not
  # explicitly asked for
  if (is.null(calc_fn)) {
    calc_fn <- is_finite_numeric(opt$convergence$abs_tol) ||
      is_finite_numeric(opt$convergence$rel_tol)
  }

  gr_norms <- c()
  if (is_finite_numeric(opt$convergence$grad_tol)) {
    gr_norms <- c(gr_norms, 2)
  }
  if (is_finite_numeric(opt$convergence$ginf_tol)) {
    gr_norms <- c(gr_norms, Inf)
  }
  if (isTRUE(calc_gr)) {
    gr_norms <- unique(c(gr_norms, 2, Inf))
  }
  if (is.null(calc_gr)) {
    calc_gr <- length(gr_norms) > 0
  } else {
    calc_gr <- isTRUE(calc_gr)
  }

  f <- NULL
  if (calc_fn || has_fn_curr(opt, iter + 1)) {
    if (!has_fn_curr(opt, iter + 1)) {
      f <- fg$fn(par)
      if (count_fg) {
        opt <- set_fn_curr(opt, f, iter + 1)
        opt$counts$fn <- opt$counts$fn + 1
      }
    } else {
      f <- opt$cache$fn_curr
    }
  }

  g2n <- NULL
  ginfn <- NULL
  if (calc_gr || has_gr_curr(opt, iter + 1)) {
    if (!has_gr_curr(opt, iter + 1)) {
      g <- fg$gr(par)
      if (grad_is_first_stage(opt) && count_fg) {
        opt <- set_gr_curr(opt, g, iter + 1)
        opt$counts$gr <- opt$counts$gr + 1
      }
    } else {
      g <- opt$cache$gr_curr
    }
    if (2 %in% gr_norms) {
      g2n <- norm2(g)
    }
    if (Inf %in% gr_norms) {
      ginfn <- norm_inf(g)
    }
  }

  if (!is.null(par_old)) {
    step_size <- norm2(par - par_old)
  } else {
    step_size <- 0
  }

  alpha <- 0
  if (!is.null(opt$stages[["gradient_descent"]]$step_size$value)) {
    alpha <- norm2(opt$stages[["gradient_descent"]]$step_size$value)
    if (is.null(alpha)) {
      alpha <- 0
    }
  }

  res <- list(
    opt = opt,
    f = f,
    g2n = g2n,
    ginfn = ginfn,
    nf = opt$counts$fn,
    ng = opt$counts$gr,
    step = step_size,
    alpha = alpha,
    iter = iter
  )

  if ("momentum" %in% names(opt$stages)) {
    res$mu <- opt$stages[["momentum"]]$step_size$value
    if (is.null(res$mu)) {
      res$mu <- 0
    }
  }

  Filter(Negate(is.null), res)
}


#' Check Optimization Convergence
#'
#' Updates the optimizer with information about convergence or termination,
#' signaling if the optimization process should stop.
#'
#' On returning from this function, the updated value of `opt` will
#' contain:
#'
#' * A boolean value `is_terminated` which is `TRUE` if
#' termination has been indicated, and `FALSE` otherwise.
#' * A list `terminate` if `is_terminated` is `TRUE`. This
#' contains two items: `what`, a short string describing what caused the
#' termination, and `val`, the value of the termination criterion that
#' caused termination. This list will not be present if `is_terminated` is
#' `FALSE`.
#'
#' Convergence criteria are only checked here. To set these criteria, use
#' [make_mize()] or [mize_init()].
#'
#' @param mize_step_info Step info for this iteration, created by
#'   [mize_step_summary()]
#' @return `opt` updated with convergence and termination data. See
#'   'Details'.
#' @export
#' @examples
#' rb_fg <- list(
#'   fn = function(x) {
#'     100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
#'   },
#'   gr = function(x) {
#'     c(
#'       -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'       200 * (x[2] - x[1] * x[1])
#'     )
#'   }
#' )
#' rb0 <- c(-1.2, 1)
#'
#' opt <- make_mize(method = "BFGS", par = rb0, fg = rb_fg, max_iter = 30)
#' mize_res <- mize_step(opt = opt, par = rb0, fg = rb_fg)
#' step_info <- mize_step_summary(mize_res$opt, mize_res$par, rb_fg, rb0)
#' # check convergence by looking at opt$is_terminated
#' opt <- check_mize_convergence(step_info)
check_mize_convergence <- function(mize_step_info) {
  opt <- mize_step_info$opt

  convergence <- opt$convergence

  terminate <- check_counts(
    opt,
    convergence$max_fn,
    convergence$max_gr,
    convergence$max_fg
  )
  if (!is.null(terminate)) {
    opt$terminate <- terminate
    opt$is_terminated <- TRUE
    return(opt)
  }

  terminate <- check_step_conv(
    opt,
    opt$iter,
    mize_step_info$step,
    convergence$step_tol
  )
  if (!is.null(terminate)) {
    opt$terminate <- terminate
    opt$is_terminated <- TRUE
    return(opt)
  }

  terminate <- check_gr_conv(opt, convergence$grad_tol, convergence$ginf_tol)
  if (!is.null(terminate)) {
    opt$terminate <- terminate
    opt$is_terminated <- TRUE
    return(opt)
  }

  if (!is.null(opt$cache$fn_curr)) {
    fn_new <- opt$cache$fn_curr
    fn_old <- convergence$fn_new
    convergence$fn_new <- fn_new
    opt$convergence <- convergence

    terminate <- check_fn_conv(
      opt,
      opt$iter,
      fn_old,
      fn_new,
      convergence$abs_tol,
      convergence$rel_tol
    )
    if (!is.null(terminate)) {
      opt$is_terminated <- TRUE
      opt$terminate <- terminate
      return(opt)
    }
  }

  if (opt$iter == convergence$max_iter) {
    opt$is_terminated <- TRUE
    opt$terminate <- list(what = "max_iter", val = convergence$max_iter)
  }

  opt
}
