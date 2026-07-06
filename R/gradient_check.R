#' Check Gradient Consistency
#'
#' Compares an analytic gradient with a finite-difference approximation of the
#' objective function at a parameter vector.
#'
#' The `fg` argument uses the same list shape as [mize()]: it must contain
#' `fn` and `gr` functions. If it also contains an `fg` function, the
#' combined function-and-gradient result is checked for consistency with the
#' separate `fn` and `gr` functions at `par`.
#'
#' The finite-difference step for each coordinate is
#' `abs_eps + rel_eps * pmax(abs(par), 1)`. Central differences are usually
#' more accurate and use twice as many function evaluations as forward
#' differences.
#'
#' @param fg Function and gradient list. See the documentation of [mize()].
#' @param par Parameter vector where the gradient should be checked.
#' @param method Finite-difference method: either `"central"` or `"forward"`.
#' @param rel_eps Relative finite-difference step. If `NULL`, a method-specific
#'   default is used.
#' @param abs_eps Non-negative scalar added to each finite-difference step.
#' @param ... Additional arguments passed to `fg$fn`, `fg$gr`, and `fg$fg` when
#'   present.
#' @return A list with components:
#'
#' * `par`: Parameter vector checked.
#' * `method`: Finite-difference method used.
#' * `step`: Per-coordinate finite-difference steps.
#' * `fn`: Function value at `par`.
#' * `gr`: Gradient returned by `fg$gr(par)`.
#' * `gr_fd`: Finite-difference gradient approximation.
#' * `error`: Difference `gr - gr_fd`.
#' * `abs_error`: Absolute coordinate-wise errors.
#' * `rel_error`: Relative coordinate-wise errors, scaled by the larger
#'   coordinate magnitude of `gr` and `gr_fd`.
#' * `max_abs_error`: Maximum absolute error.
#' * `max_rel_error`: Maximum relative error.
#' * `worst_abs`: One-row data frame for the coordinate with the largest
#'   absolute error.
#' * `worst_rel`: One-row data frame for the coordinate with the largest
#'   relative error.
#' * `by_coord`: Data frame with per-coordinate gradient comparisons.
#' * `fg_consistency`: List describing whether `fg$fg()` was checked and, when
#'   available, its function and gradient consistency with `fg$fn()` and
#'   `fg$gr()`.
#' @seealso [mize()]
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
#'
#' check <- check_mize_gradient(rb_fg, c(-1.2, 1))
#' check$max_abs_error
check_mize_gradient <- function(
  fg,
  par,
  method = c("central", "forward"),
  rel_eps = NULL,
  abs_eps = 0,
  ...
) {
  method <- match.arg(method)
  mize_check_gradient_inputs(fg, par)

  par_names <- names(par)
  par <- as.numeric(par)
  names(par) <- par_names

  rel_eps <- mize_gradient_rel_eps(rel_eps, method)
  abs_eps <- mize_gradient_abs_eps(abs_eps)
  step <- abs_eps + rel_eps * pmax(abs(par), 1)
  names(step) <- par_names

  fn <- fg[["fn"]]
  gr <- fg[["gr"]]

  f0 <- mize_scalar_result(fn(par, ...), "fg$fn(par)")
  analytic <- mize_gradient_result(gr(par, ...), length(par), "fg$gr(par)")
  names(analytic) <- par_names

  fd <- mize_fd_gradient(par, fn, f0, step, method, ...)
  comparison <- mize_error_by_coord(
    par = par,
    lhs = analytic,
    rhs = fd,
    lhs_name = "gr",
    rhs_name = "gr_fd"
  )

  err <- comparison$error
  abs_error <- comparison$abs_error
  rel_error <- comparison$rel_error
  names(err) <- par_names
  names(abs_error) <- par_names
  names(rel_error) <- par_names

  list(
    par = par,
    method = method,
    step = step,
    fn = f0,
    gr = analytic,
    gr_fd = fd,
    error = err,
    abs_error = abs_error,
    rel_error = rel_error,
    max_abs_error = max(abs_error),
    max_rel_error = max(rel_error),
    worst_abs = mize_worst_coordinate(comparison, "abs_error"),
    worst_rel = mize_worst_coordinate(comparison, "rel_error"),
    by_coord = comparison,
    fg_consistency = mize_fg_consistency(fg, par, f0, analytic, ...)
  )
}

mize_check_gradient_inputs <- function(fg, par) {
  if (!is.list(fg)) {
    stop("fg must be a list", call. = FALSE)
  }
  if (!is.function(fg[["fn"]])) {
    stop("fg$fn must be a function", call. = FALSE)
  }
  if (!is.function(fg[["gr"]])) {
    stop("fg$gr must be a function", call. = FALSE)
  }
  if (!is.numeric(par) || !is.null(dim(par)) || length(par) == 0) {
    stop("par must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(par))) {
    stop("par must contain only finite values", call. = FALSE)
  }
}

mize_gradient_rel_eps <- function(rel_eps, method) {
  if (is.null(rel_eps)) {
    if (method == "central") {
      return(.Machine$double.eps^(1 / 3))
    }
    return(sqrt(.Machine$double.eps))
  }
  if (
    !is.numeric(rel_eps) ||
      length(rel_eps) != 1 ||
      !is.finite(rel_eps) ||
      rel_eps <= 0
  ) {
    stop("rel_eps must be a positive finite numeric scalar", call. = FALSE)
  }
  rel_eps
}

mize_gradient_abs_eps <- function(abs_eps) {
  if (
    !is.numeric(abs_eps) ||
      length(abs_eps) != 1 ||
      !is.finite(abs_eps) ||
      abs_eps < 0
  ) {
    stop("abs_eps must be a non-negative finite numeric scalar", call. = FALSE)
  }
  abs_eps
}

mize_scalar_result <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value)) {
    stop(label, " must return a finite numeric scalar", call. = FALSE)
  }
  as.numeric(value)
}

mize_gradient_result <- function(value, n, label) {
  if (!is.numeric(value) || length(value) != n || any(!is.finite(value))) {
    stop(
      label,
      " must return a finite numeric vector with length matching par",
      call. = FALSE
    )
  }
  as.numeric(value)
}

mize_fd_gradient <- function(par, fn, f0, step, method, ...) {
  fd <- numeric(length(par))

  for (i in seq_along(par)) {
    par_plus <- par
    par_plus[i] <- par_plus[i] + step[i]
    fplus <- mize_scalar_result(
      fn(par_plus, ...),
      paste0("fg$fn(par + step)[", i, "]")
    )

    if (method == "central") {
      par_minus <- par
      par_minus[i] <- par_minus[i] - step[i]
      fminus <- mize_scalar_result(
        fn(par_minus, ...),
        paste0("fg$fn(par - step)[", i, "]")
      )
      fd[i] <- (fplus - fminus) / (2 * step[i])
    } else {
      fd[i] <- (fplus - f0) / step[i]
    }
  }

  names(fd) <- names(par)
  fd
}

mize_relative_error <- function(error, lhs, rhs) {
  abs(error) / pmax(abs(lhs), abs(rhs), .Machine$double.eps)
}

mize_error_by_coord <- function(par, lhs, rhs, lhs_name, rhs_name) {
  error <- lhs - rhs
  res <- data.frame(
    index = seq_along(par),
    name = mize_par_names(par),
    par = as.numeric(par),
    lhs = as.numeric(lhs),
    rhs = as.numeric(rhs),
    error = as.numeric(error),
    abs_error = abs(error),
    rel_error = mize_relative_error(error, lhs, rhs),
    row.names = NULL
  )
  names(res)[names(res) == "lhs"] <- lhs_name
  names(res)[names(res) == "rhs"] <- rhs_name
  res
}

mize_par_names <- function(par) {
  par_names <- names(par)
  if (is.null(par_names)) {
    return(rep(NA_character_, length(par)))
  }
  ifelse(nzchar(par_names), par_names, NA_character_)
}

mize_worst_coordinate <- function(comparison, error_name) {
  res <- comparison[which.max(comparison[[error_name]]), , drop = FALSE]
  row.names(res) <- NULL
  res
}

mize_fg_consistency <- function(fg, par, fn, gr, ...) {
  if (is.null(fg[["fg"]])) {
    return(list(checked = FALSE))
  }
  if (!is.function(fg[["fg"]])) {
    stop("fg$fg must be a function", call. = FALSE)
  }

  combined <- fg[["fg"]](par, ...)
  if (!is.list(combined)) {
    stop("fg$fg(par) must return a list", call. = FALSE)
  }

  fg_fn <- mize_scalar_result(combined[["fn"]], "fg$fg(par)$fn")
  fg_gr <- mize_gradient_result(
    combined[["gr"]],
    length(par),
    "fg$fg(par)$gr"
  )
  names(fg_gr) <- names(par)

  fn_error <- fg_fn - fn
  gr_comparison <- mize_error_by_coord(
    par = par,
    lhs = fg_gr,
    rhs = gr,
    lhs_name = "fg_gr",
    rhs_name = "gr"
  )
  gr_error <- gr_comparison$error
  gr_abs_error <- gr_comparison$abs_error
  gr_rel_error <- gr_comparison$rel_error
  names(gr_error) <- names(par)
  names(gr_abs_error) <- names(par)
  names(gr_rel_error) <- names(par)

  list(
    checked = TRUE,
    fn = list(
      fn = fn,
      fg_fn = fg_fn,
      error = fn_error,
      abs_error = abs(fn_error),
      rel_error = mize_relative_error(fn_error, fg_fn, fn)
    ),
    gr = list(
      gr = gr,
      fg_gr = fg_gr,
      error = gr_error,
      abs_error = gr_abs_error,
      rel_error = gr_rel_error,
      max_abs_error = max(gr_abs_error),
      max_rel_error = max(gr_rel_error),
      worst_abs = mize_worst_coordinate(gr_comparison, "abs_error"),
      worst_rel = mize_worst_coordinate(gr_comparison, "rel_error"),
      by_coord = gr_comparison
    )
  )
}
