# Returns a termination list if step falls below step_tol
# A zero step is allowed if this is a restart step
check_step_conv <- function(opt, iter, step = NULL, step_tol = NULL) {
  if (is.null(step) || is.null(step_tol) || is_restart_iter(opt, iter) ||
      step >= step_tol) {
    return()
  }
  list(what = "step_tol", val = step)
}

# Return a termination list if maximum number of function and/or gradient
# calls has been exceeded
check_counts <- function(opt, max_fn, max_gr, max_fg) {
  terminate <- NULL
  if (opt$counts$fn >= max_fn) {
    terminate <- list(
      what = "max_fn",
      val = opt$counts$fn
    )
  }
  else if (opt$counts$gr >= max_gr) {
    terminate <- list(
      what = "max_gr",
      val = opt$counts$gr
    )
  }
  else if (opt$counts$fn + opt$counts$gr >= max_fg) {
    terminate <- list(
      what = "max_fg",
      val = opt$counts$fn + opt$counts$gr
    )
  }
  terminate
}

# Return a termination list if the gradient 2 norm tolerance (grad_tol) or
# infinity norm tolerance is reached. Termination is also indicated if
# any element of the gradient vector is not finite. Requires the gradient
# have already been calculated - this routine does NOT calculate it if it's
# not present
check_gr_conv <- function(opt, grad_tol, ginf_tol) {
  if (is.null(opt$cache$gr_curr)) {
    return()
  }

  if (any(!is.finite(opt$cache$gr_curr))) {
    return(list(what = "gr_inf", val = Inf))
  }

  if (!is.null(grad_tol)) {
    gtol <- norm2(opt$cache$gr_curr)
    if (gtol <= grad_tol) {
      return(list(what = "grad_tol", val = gtol))
    }
  }

  if (!is.null(ginf_tol)) {
    gitol <- norm_inf(opt$cache$gr_curr)
    if (gitol <= ginf_tol) {
      return(list(what = "ginf_tol", val = gitol))
    }
  }
}

# Return a termination list if the absolute or relative tolerance is reached
# for the difference between fn_old and fn_new. Termination is also indicated
# if fn_new is non-finite. Tolerance is not checked if this is a restart
# iteration.
check_fn_conv <- function(opt, iter, fn_old, fn_new, abs_tol, rel_tol) {
  if (!is.finite(fn_new)) {
    return(list(what = "fn_inf", val = fn_new))
  }

  if (is.null(fn_old)) {
    return()
  }

  if (is_restart_iter(opt, iter)) {
    return()
  }

  if (!is.null(abs_tol)) {
    atol <- abs(fn_old - fn_new)
    if (atol < abs_tol) {
      return(list(what = "abs_tol", val = atol))
    }
  }

  if (!is.null(rel_tol)) {
    rtol <- abs(fn_old - fn_new) / min(abs(fn_new), abs(fn_old))
    if (rtol < rel_tol) {
      return(list(what = "rel_tol", val = rtol))
    }
  }
}

# True if this iteration was marked as a restart
# Zero step size and function difference is allowed under these circumstances.
is_restart_iter <- function(opt, iter) {
  !is.null(opt$restart_at) && opt$restart_at == iter
}

