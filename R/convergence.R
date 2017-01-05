# if any of the termination criteria are fulfilled, a list is returned
# saying which one and what the value was that triggered the termination.
# Otherwise, an empty list is returned
# Gradient and Function-based termination (abs_tol, rel_tol and grad_tol)
# are checked only if the function and gradient values were calculated
# in the optimization step.
check_termination <- function(terminate, opt, iter, step = NULL,
                              max_fn, max_gr, max_fg,
                              abs_tol, rel_tol, grad_tol, ginf_tol, step_tol) {
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
  else if (step_converged(opt, iter, step, step_tol)) {
    terminate <- list(
      what = "step_tol",
      val = step
    )
  }
  if (!is.null(grad_tol) && !is.null(opt$cache$gr_curr)) {
    if (any(!is.finite(opt$cache$gr_curr))) {
      terminate$what <- "gr_inf"
      terminate$val <- Inf
      return(terminate)
    }
    gtol <- norm2(opt$cache$gr_curr)
    if (gtol <= grad_tol) {
      terminate <- list(
        what = "grad_tol",
        val = gtol
      )
    }
  }
  if (!is.null(ginf_tol) && !is.null(opt$cache$gr_curr)) {
    if (any(!is.finite(opt$cache$gr_curr))) {
      terminate$what <- "gr_inf"
      terminate$val <- Inf
      return(terminate)
    }
    gitol <- norm_inf(opt$cache$gr_curr)
    if (gitol <= ginf_tol) {
      terminate <- list(
        what = "ginf_tol",
        val = gitol
      )
    }
  }
  if (!is.null(rel_tol) || !is.null(abs_tol)) {
    if (!is.null(opt$cache$fn_curr)) {
      fn_new <- opt$cache$fn_curr
      if (!is.finite(fn_new)) {
        terminate$what <- "fn_inf"
        terminate$val <- fn_new
        return(terminate)
      }

      if (!is.null(terminate$fn_new)) {
        fn_old <- terminate$fn_new
        atol <- abs(fn_old - fn_new)
        if (is.null(opt$restart_at) || (opt$restart_at != iter)) {
          if (!is.null(abs_tol) && atol < abs_tol) {
            terminate$what <- "abs_tol"
            terminate$val <- atol
          }
          else {
            if (!is.null(rel_tol)) {
              rtol <- abs(fn_old - fn_new) / min(abs(fn_new), abs(fn_old))
              if (rtol < rel_tol) {
                terminate$what <- "rel_tol"
                terminate$val <- rtol
              }
            }
          }
        }
      }

      terminate$fn_new <- fn_new
    }
  }
  terminate
}

step_converged <- function(opt, iter, step = NULL, step_tol = NULL) {
  !is.null(step) &&
  !is.null(step_tol) &&
  step < step_tol &&
  (is.null(opt$restart_at) || opt$restart_at != iter)
}
