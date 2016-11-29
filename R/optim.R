
# Test Data ---------------------------------------------------------------

out0 <- c(-1.2, 1)
# taken from the optim man page
rosenbrock_fg <- list(
  fn = function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
  },
  gr = function(x) {
    x1 <- x[1]
    x2 <- x[2]
    c(
     -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
      200 *      (x2 - x1 * x1))
  },
  fg <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    a <- (x2 - x1 * x1)
    b <- 1 - x1
    list(
       f = 100 * a * a + b * b,
       g = c(
         -400 * x1 * a - 2 * b,
          200 * a
        )
    )
  },
  n = 2
)

wrap_fg <- function(fg) {
  nc <- fg$n
  nr <- 1

  list(
    nc = nc,
    nr = nr,
    method = list(
      cost = list(
        fn = function(inp, out, method) {
          fg$fn(out$ym)
        }
      ),
      eps = .Machine$double.eps
    ),
    grad_func = function(inp, out, method, mat_name) {
      list(gm = matrix(fg$gr(out[[mat_name]]), nrow = nr, ncol = nc))
    }
  )
}
rosenbrock <- wrap_fg(rosenbrock_fg)

# Function / Gradient ----------------------------------------------------------------

# Uncached function evaluation for arbitrary values of par
calc_fn <- function(opt, par, fn) {
  opt$fn <- fn(par)
  opt$counts$fn <- opt$counts$fn + 1
  opt
}

# Cached function evaluation for par value after finding a step size
# (possibly re-usable)
calc_fn_new <- function(opt, par, fn, iter) {
  if (is.null(opt$cache$fn_new_iter) || opt$cache$fn_new_iter != iter) {
    opt <- set_fn_new(opt, fn(par), iter)
    opt$counts$fn <- opt$counts$fn + 1
  }
  opt
}

set_fn_new <- function(opt, val, iter) {
  opt$cache$fn_new <- val
  opt$cache$fn_new_iter <- iter
  opt
}

# Cached function evaluation for par at starting point
# (possibly re-usable)
calc_fn_curr <- function(opt, par, fn, iter) {
  if (is.null(opt$cache$fn_curr_iter) || opt$cache$fn_curr_iter != iter) {
    opt <- set_fn_curr(opt, fn(par), iter)
    opt$counts$fn <- opt$counts$fn + 1
    #message("Calculated and cached fn = ", formatC(opt$cache$fn_curr))
  }
  opt
}

set_fn_curr <- function(opt, val, iter) {
  opt$cache$fn_curr <- val
  opt$cache$fn_curr_iter <- iter
  opt
}

# Cached gradient evaluation for par value at start of iteration
# (possibly re-usable)
calc_gr_curr <- function(opt, par, gr, iter) {
  if (is.null(opt$cache$gr_curr_iter) || opt$cache$gr_curr_iter != iter) {
    opt <- set_gr_curr(opt, gr(par), iter)
    opt$counts$gr <- opt$counts$gr + 1
  }
  opt
}

set_gr_curr <- function(opt, val, iter) {
  opt$cache$gr_curr <- val
  opt$cache$gr_curr_iter <- iter
  opt
}

# Uncached gradient evaluation for arbitrary values of par
calc_gr <- function(opt, par, gr) {
  opt$gr <- gr(par)
  opt$counts$gr <- opt$counts$gr + 1
  opt
}

# Optimizer ---------------------------------------------------------------

optloop <- function(opt, par, fn, gr, max_iter = 10, verbose = FALSE,
                    store_progress = FALSE, invalidate_cache = FALSE) {

  opt <- opt_init(opt, par, fn, gr, 0)

  progress <- data.frame()

  if (verbose || store_progress) {
    res <- opt_results(opt, par, fn, gr, 0)
    if (store_progress) {
      progress <- update_progress(opt_res = res, progress = progress)
    }
    if (verbose) {
      opt_report(res, print_time = TRUE, print_par = TRUE)
    }
  }

  if (max_iter < 1) {
    if (!verbose) {
      res <- opt_results(opt, par, fn, gr, 0)
      if (store_progress) {
        res$progress <- progress
      }
    }
    return(res)
  }

  for (iter in 1:max_iter) {

    if (invalidate_cache) {
      for (name in names(opt$cache)) {
        iter_name <- paste0(name, "_iter")
        if (!is.null(opt$cache[[iter_name]])) {
          opt$cache[[iter_name]] <- "invalid"
        }
      }
    }

    par0 <- par

    step_res <- optimize_step(opt, par, fn, gr, iter)
    opt <- step_res$opt
    par <- step_res$par


    if (opt$verbose || store_progress) {
      res <- opt_results(opt, par, fn, gr, iter, par0)
      if (store_progress) {
        progress <- update_progress(opt_res = res, progress = progress)
      }
      if (verbose) {
        opt_report(res, print_time = TRUE, print_par = TRUE)
      }
    }
  }

  if (verbose) {
    res <- opt_results(opt, par, fn, gr, iter)
  }
  if (store_progress) {
    res$progress <- progress
  }
  res
}

# One Step of Optimization
#
optimize_step <- function(opt, par, fn, gr, iter) {
  opt <- life_cycle_hook("step", "before", opt, par, fn, gr, iter)

  par0 <- par
  step_result <- NULL
  counts <- make_counts()
  for (i in 1:length(opt$stages)) {

    stage <- opt$stages[[i]]
    opt <- life_cycle_hook(stage$type, "before", opt, par, fn, gr, iter)
    opt <- life_cycle_hook(stage$type, "during", opt, par, fn, gr, iter)
    opt <- life_cycle_hook(stage$type, "after", opt, par, fn, gr, iter)

    stage <- opt$stages[[i]]
    if (is.null(step_result)) {
      step_result <- stage$result
    }
    else {
      step_result <- step_result + stage$result
    }

    if (opt$eager_update) {
      par <- par + stage$result
    }
    message(iter, " ", substr(stage$type, 1, 2)
            ," par = ", vec_formatC(par)
            ," p = ", vec_formatC(stage$direction$value)
            , " a = ", formatC(stage$step_size$value)
            , " ap = ", vec_formatC(stage$result)
            , " f = ", formatC(fn(par)))
  }

  if (!opt$eager_update) {
    par <- par + step_result
  }

  # intercept whether we want to accept the new solution
  opt <- life_cycle_hook("validation", "before", opt, par, fn, gr, iter,
                         par0, step_result)
  opt$ok <- TRUE
  opt <- life_cycle_hook("validation", "during", opt, par, fn, gr, iter,
                         par0, step_result)

  # If the this solution was vetoed, roll back to the previous one.
  if (!opt$ok) {
    par <- par0
  }

  opt <- life_cycle_hook("step", "after", opt, par, fn, gr, iter, par0,
                         step_result)
  opt$counts <- update_counts(opt$counts, counts)

  list(opt = opt, par = par)
}


opt_init <- function(opt, par, fn, gr, iter) {
  opt <- register_hooks(opt)
  opt <- life_cycle_hook("opt", "init", opt, par, fn, gr, iter)
  list_hooks(opt)

  opt
}

opt_results <- function(opt, par, fn, gr, iter, par0 = NULL) {

  f <- fn(par)
  g <- gr(par)
  g2n <- norm2(g)
  step_size <- norm2(par - par0)

  list(
    f = f,
    g2n = g2n,
    nf = opt$counts$fn,
    ng = opt$counts$gr,
    par = par,
    step = step_size,
    iter = iter
  )
}

opt_report <- function(opt_result, print_time = FALSE, print_par = FALSE) {
  msg <- paste0("iter ", opt_result$iter
                , " f = ", formatC(opt_result$f)
                , " |g| = ", formatC(opt_result$g2n)
                , " nf = ", opt_result$nf
                , " ng = ", opt_result$ng
  )

  if (print_time) {
    msg <- paste(format(Sys.time(), "%H:%M:%S"), msg, collapse = " ")
  }

  if (print_par) {
    msg <- paste0(msg, " par = ", vec_formatC(opt_result$par))
  }

  message(msg)
}

update_progress <- function(opt_res, progress) {
  progress <- rbind(progress, c(opt_res$f, opt_res$g2n, opt_res$nf, opt_res$ng,
                                opt_res$step))
  # Probably not a major performance issue to regenerate column names each time
  colnames(progress) <- c("f", "g2n", "nf", "ng", "step")
  rownames(progress)[nrow(progress)] <- opt_res$iter
  progress
}

# Constructor -------------------------------------------------------------

make_opt <- function(stages,
                     verbose = FALSE) {
  opt <- list(
    init = function(opt, par, fn, gr, iter) {
      for (i in 1:length(opt$stages)) {
        opt <- life_cycle_hook(opt$stages[[i]]$type, "init", opt, par, fn, gr, iter)
      }
      opt
    },
    cache = list(),
    stages = stages,
    counts = make_counts(),
    eager_update = FALSE,
    verbose = verbose
  )

  if (!is.null(opt$init)) {
    attr(opt$init, 'event') <- 'init opt'
    attr(opt$init, 'name') <- 'init opt'
  }
  opt
}

make_stage <- function(type, direction, step_size, depends = NULL) {

  stage <- list(
    type = type,
    direction = direction,
    step_size = step_size,
    init = function(opt, stage, par, fn, gr, iter) {
      for (sub_stage_name in c("direction", "step_size")) {
        phase <- paste0(stage$type, " ", sub_stage_name)
        opt <- life_cycle_hook(phase, "init", opt, par, fn, gr, iter)
      }

      list(opt = opt)
    },
    calculate = function(opt, stage, par, fn, gr, iter) {
      for (sub_stage_name in c("direction", "step_size")) {
        phase <- paste0(stage$type, " ", sub_stage_name)
        opt <- life_cycle_hook(phase, "during", opt, par, fn, gr, iter)
      }

      list(opt = opt)
    },
    after_stage = function(opt, stage, par, fn, gr, iter) {
      #message("After stage: Calculating stage result")
      stage$result <- stage$direction$value * stage$step_size$value
      list(stage = stage)
    },
    counts = make_counts()
  )

  if (!is.null(depends)) {
    stage$depends <- depends
  }

  if (!is.null(stage$init)) {
    attr(stage$init, 'event') <- paste0('init ', type)
    attr(stage$init, 'name') <- paste0(type,' init')
  }
  if (!is.null(stage$calculate)) {
    attr(stage$calculate, 'event') <- paste0('during ', type)
    attr(stage$calculate, 'name') <- paste0(type, ' calculate')
  }
  if (!is.null(stage$after_stage)) {
    attr(stage$after_stage, 'event') <- paste0('after ', type)
    attr(stage$after_stage, 'name') <- paste0(type, ' after stage')
  }
  if (!is.null(stage$after_step)) {
    attr(stage$after_step, 'event') <- 'after step'
    attr(stage$after_step, 'name') <- paste0(type, ' after step')
  }

  res <- list()
  res[[type]] <- stage
  res
}

gradient_stage <- function(direction, step_size) {
  make_stage(type = "gradient_descent", direction, step_size,
             depends = c('gradient'))
}

momentum_stage <- function(direction = momentum_direction(normalize = FALSE),
                           step_size) {
  make_stage(type = "momentum", direction, step_size)
}

momentum_correction_stage <- function(
  direction = momentum_correction_direction(),
  step_size = momentum_correction_step()) {
  make_stage(type = "momentum_correction", direction, step_size)
}

make_stages <- function(...) {
  stages <- list()
  varargs <- list(...)
  for (arg in varargs) {
    for (i in names(arg)) {
      stages[[i]] <- arg[[i]]
    }
  }
  stages
}


make_counts <- function() {
  list(
    fn = 0,
    gr = 0
  )
}

update_counts <- function(counts, new_counts) {
  counts$fn <- counts$fn + new_counts$fn
  counts$gr <- counts$gr + new_counts$gr
  counts
}

# Wrappers ----------------------------------------------------------------

prepend_stage <- function(opt, stage) {
  opt$stages <- append(stage, opt$stages)
  opt
}

append_stage <- function(opt, stage) {
  opt$stages <- append(opt$stages, stage)
  opt
}




# Momentum Dependencies ------------------------------------------------------------

# After step
require_update_old <- function(opt, par, fn, gr, iter, par0, update) {
  #message("update_old: saving old update")
  opt$cache$update_old <- update
  opt
}
attr(require_update_old, 'event') <- 'after step'
attr(require_update_old, 'name') <- 'update_old'
attr(require_update_old, 'depends') <- 'update_old_init'

require_update_old_init <- function(opt, stage, sub_stage, par, fn, gr, iter) {
  message("Initializing update_old")
  opt$cache$update_old <- rep(0, length(par))
  list(opt = opt)
}
attr(require_update_old_init, 'event') <- 'init momentum direction'
attr(require_update_old_init, 'name') <- 'update_old_init'

require_adaptive_restart <- function(opt, par, fn, gr, iter, par0, update) {
  if (!opt$ok) {
    opt <- life_cycle_hook("momentum", "init", opt, par, fn, gr, iter)
    message("adaptive restart: clearing old update")
  }
  else {
    opt$cache$update_old <- update
    message("adaptive restart: saving old update")
  }
  opt
}
attr(require_adaptive_restart, 'event') <- 'after step'
# Should have the same name as normal update old: we want to replace that hook
attr(require_adaptive_restart, 'name') <- 'update_old'
attr(require_adaptive_restart, 'depends') <- 'update_old_init'



# Gradient Dependencies ------------------------------------------------------------

require_gradient <- function(opt, par, fn, gr, iter) {
  if (!has_gr_curr(opt, iter)) {
    #message("require gradient: calculating gr_curr ", iter)
    opt <- calc_gr_curr(opt, par, gr, iter)

    if (any(is.nan(opt$cache$gr_curr))) {
      stop("NaN in grad. descent at iter ", iter)
    }
  }
  else {
    #message("require gradient: already have gr_curr")
  }
  opt
}
attr(require_gradient, 'event') <- 'before gradient_descent'
attr(require_gradient, 'name') <- 'gradient'

require_gradient_old <- function(opt, par, fn, gr, iter, par0, update) {
#  message("saving old gradient")
  opt$cache$gr_old <- opt$cache$gr_curr
  opt
}
attr(require_gradient_old, 'event') <- 'after step'
attr(require_gradient_old, 'name') <- 'gradient_old'

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
  message("validating by gradient: ", "g = ", vec_formatC(opt$cache$gr_curr),
    " v = ", vec_formatC(update), " g.v = ", formatC(dot(opt$cache$gr_curr, update)),
    " ok = ", opt$ok)

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
    message("Saving 'curr' values in cache due to validation failure")
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

# Predicates --------------------------------------------------------------

is_single_stage <- function(opt) {
  length(opt$stages) == 1
}

is_first_stage <- function(opt, stage) {
  stage$type == opt$stages[[1]]$type
}

is_last_stage <- function(opt, stage) {
  stage$type == opt$stages[[length(opt$stages)]]$type
}

has_fn_new <- function(opt, iter) {
  (!is.null(opt$cache$fn_new)
   && !is.null(opt$cache$fn_new_iter)
   && opt$cache$fn_new_iter == iter)
}

# Predicate for whether fn is already cached
has_fn_curr <- function(opt, iter) {
  (!is.null(opt$cache$fn_curr)
   && !is.null(opt$cache$fn_curr_iter)
   && opt$cache$fn_curr_iter == iter)
}

has_gr_curr <- function(opt, iter) {
  (!is.null(opt$cache$gr_curr)
   && !is.null(opt$cache$gr_curr_iter)
   && opt$cache$gr_curr_iter == iter)
}
