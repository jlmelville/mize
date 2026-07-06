# Gradient Dependencies ------------------------------------------------------------

# Calculate the gradient at par.
require_gradient <- function(opt, stage, par, fg, iter) {
  if (!has_gr_curr(opt, iter)) {
    opt <- calc_gr_curr(opt, par, fg$gr, iter)

    if (any(!is.finite(opt$cache$gr_curr))) {
      opt <- set_mize_termination(
        opt,
        list(
          what = "gr_inf",
          val = Inf
        )
      )
    }
  }

  list(opt = opt)
}
attr(require_gradient, "event") <- "before gradient_descent"
attr(require_gradient, "name") <- "gradient"

# Caches the gradient at the current step.
require_gradient_old <- function(opt, par, fg, iter, par0, update) {
  opt$cache$gr_old <- opt$cache$gr_curr
  opt
}
attr(require_gradient_old, "event") <- "after step"
attr(require_gradient_old, "name") <- "gradient_old"
