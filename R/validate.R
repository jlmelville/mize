# Validate ----------------------------------------------------------------

# Checks that the function value has decreased over the step
require_validate_fn <- function(opt, par, fg, iter, par0, update) {
  if (can_restart(opt, iter)) {
    opt$ok <- opt$cache$fn_new < opt$cache$fn_curr
  }
  opt
}
attr(require_validate_fn, 'name') <- 'validate_fn'
attr(require_validate_fn, 'event') <- 'during validation'
attr(require_validate_fn, 'depends') <- 'fn_new fn_curr save_cache_on_failure'

# Checks that the gradient is a descent direction
# This relies on the gradient being calculated in the "classical" location
# i.e. not using the implementation of Nesterov Acceleration
require_validate_gr <- function(opt, par, fg, iter, par0, update) {
  if (can_restart(opt, iter)) {
    opt$ok <- dot(opt$cache$gr_curr, update) < 0
  }
  opt
}
attr(require_validate_gr, 'name') <- 'validate_gr'
attr(require_validate_gr, 'event') <- 'during validation'
attr(require_validate_gr, 'depends') <- 'gradient save_cache_on_failure'

# Checks that the update vector is getting larger
require_validate_speed <- function(opt, par, fg, iter, par0, update) {
  if (can_restart(opt, iter)) {
    opt$ok <- sqnorm2(update) > sqnorm2(opt$cache$update_old)
  }
  opt
}
attr(require_validate_speed, 'name') <- 'validate_speed'
attr(require_validate_speed, 'event') <- 'during validation'
attr(require_validate_speed, 'depends') <- 'save_cache_on_failure'

# Validate Dependencies ------------------------------------------------------------

# Caches the current fn value
require_fn_curr <- function(opt, par, fg, iter, par0, update) {
  if (!has_fn_curr(opt, iter)) {
    opt <- calc_fn_curr(opt, par, fg$fn, iter)
  }
  opt
}
attr(require_fn_curr, 'name') <- 'fn_curr'
attr(require_fn_curr, 'event') <- 'before step'
attr(require_fn_curr, 'depends') <- 'update_fn_cache'

# Caches the new fn value
require_fn_new <- function(opt, par, fg, iter, par0, update) {
  if (!has_fn_new(opt, iter)) {
    opt <- calc_fn_new(opt, par, fg$fn, iter)
  }
  opt
}
attr(require_fn_new, 'name') <- 'fn_new'
attr(require_fn_new, 'event') <- 'before validation'

# Caches the new fn value as the current value for the next iteration
require_update_fn_cache <- function(opt, par, fg, iter, par0, update) {
  if (opt$ok && has_fn_new(opt, iter)) {
    opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
  }
  opt
}
attr(require_update_fn_cache, 'name') <- 'update_fn_cache'
attr(require_update_fn_cache, 'event') <- 'after step'

# Keep the old cached values around in the event of failure
require_save_cache_on_failure <- function(opt, par, fg, iter, par0, update) {
  # not safe to re-use gr_curr and fn_curr unless gradient calc is the first
  # stage: Nesterov results in moving par via momentum before grad calc.
  # Different result will occur after restart
  if (!opt$ok && opt$stages[[1]]$type == "gradient_descent") {
    cache <- opt$cache
    for (name in names(cache)) {
      if (endsWith(name, "_curr")) {
        iter_name <- paste0(name, "_iter")
        cache_iter <- cache[[iter_name]]
        if (!is.null(cache_iter) && cache_iter == iter) {
          cache[[iter_name]] <- cache_iter + 1
        }
      }
    }
    opt$cache <- cache
  }
  opt
}
attr(require_save_cache_on_failure, 'name') <- 'save_cache_on_failure'
attr(require_save_cache_on_failure, 'event') <- 'after step'
