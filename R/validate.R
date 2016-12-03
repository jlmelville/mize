# Validate ----------------------------------------------------------------

require_validate_fn <- function(opt, par, fn, gr, iter, par0, update) {
  #message("validating fn")
  #message("fn_new = ", formatC(opt$cache$fn_new), " fn_curr = ", formatC(opt$cache$fn_curr))
  opt$ok <- opt$cache$fn_new < opt$cache$fn_curr
  opt
}
attr(require_validate_fn, 'name') <- 'validate_fn'
attr(require_validate_fn, 'event') <- 'during validation'
attr(require_validate_fn, 'depends') <- 'fn_new fn_curr save_cache_on_failure'


# This relies on the gradient being calculated in the "classical" location
# i.e. not using the implementation of Nesterov Acceleration
# Wants: gradient
require_validate_gr <- function(opt, par, fn, gr, iter, par0, update) {
  opt$ok <- dot(opt$cache$gr_curr, update) < 0
  # message("validating by gradient: ", "g = ", vec_formatC(opt$cache$gr_curr),
  #   " v = ", vec_formatC(update), " g.v = ", formatC(dot(opt$cache$gr_curr, update)),
  #   " ok = ", opt$ok)

  opt
}
attr(require_validate_gr, 'name') <- 'validate_gr'
attr(require_validate_gr, 'event') <- 'during validation'
attr(require_validate_gr, 'depends') <- 'gradient save_cache_on_failure'

# Validate Dependencies ------------------------------------------------------------

require_fn_curr <- function(opt, par, fn, gr, iter, par0, update) {
  #message("requiring fn_curr")
  if (!has_fn_curr(opt, iter)) {
    #message("require fn curr: calculating fn_curr ", iter)
    opt <- calc_fn_curr(opt, par, fn, iter)
  }
  #else {
  #message("require fn curr: already have fn_curr")
  #}

  opt
}
attr(require_fn_curr, 'name') <- 'fn_curr'
attr(require_fn_curr, 'event') <- 'before step'
attr(require_fn_curr, 'depends') <- 'update_fn_cache'

require_fn_new <- function(opt, par, fn, gr, iter, par0, update) {
  #message("require fn_new")
  if (!has_fn_new(opt, iter)) {
    #message("require fn new: calculating fn_new ", iter)
    opt <- calc_fn_new(opt, par, fn, iter)
  }
  else {
    #message("require fn new: already have fn_new")
  }
  opt
}
attr(require_fn_new, 'name') <- 'fn_new'
attr(require_fn_new, 'event') <- 'before validation'

require_update_fn_cache <- function(opt, par, fn, gr, iter, par0, update) {
  if (opt$ok && has_fn_new(opt, iter)) {
    #message("prestoring f_curr for iter ", iter + 1)
    opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
  }
  opt
}
attr(require_update_fn_cache, 'name') <- 'update_fn_cache'
attr(require_update_fn_cache, 'event') <- 'after step'

require_save_cache_on_failure <- function(opt, par, fn, gr, iter, par0, update) {
  # not safe to re-use gr_curr and fn_curr unless gradient calc is the first
  # stage: Nesterov results in moving par via momentum before grad calc.
  # Different result will occur after restart
  if (!opt$ok && opt$stages[[1]]$type == "gradient_descent") {
    #    message("Saving 'curr' values in cache due to validation failure")
    cache <- opt$cache
    for (name in names(cache)) {
      if (endsWith(name, "_curr")) {
        iter_name <- paste0(name, "_iter")
        cache_iter <- cache[[iter_name]]
        if (!is.null(cache_iter) && cache_iter == iter) {
          #message("Saving '", name, "'")
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
