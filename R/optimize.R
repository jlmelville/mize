# Optimizer ---------------------------------------------------------------

optloop <- function(opt, par, fn, gr, max_iter = 10, verbose = FALSE,
                    store_progress = FALSE, invalidate_cache = FALSE,
                    ret_opt = FALSE) {

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


    if (verbose || store_progress) {
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
  if (ret_opt) {
    res$opt <- opt
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
    opt$stage_i <- i
    stage <- opt$stages[[i]]
    opt <- life_cycle_hook(stage$type, "before", opt, par, fn, gr, iter)
    opt <- life_cycle_hook(stage$type, "during", opt, par, fn, gr, iter)
    opt <- life_cycle_hook(stage$type, "after", opt, par, fn, gr, iter)

    # should run "after stage"
    stage <- opt$stages[[i]]

    if (is.null(step_result)) {
      step_result <- stage$result
    }
    else {
      step_result <- step_result + stage$result
    }

    # should run "after stage"
    if (opt$eager_update) {
      par <- par + stage$result
    }

    opt <- life_cycle_hook("stage", "after", opt, par, fn, gr, iter)
  }

  # should run "after stages"
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
  #  list_hooks(opt)

  opt <- life_cycle_hook("opt", "init", opt, par, fn, gr, iter)

  opt
}

opt_results <- function(opt, par, fn, gr, iter, par0 = NULL) {
  f <- fn(par)
  g <- gr(par)
  g2n <- norm2(g)
  step_size <- norm2(par - par0)

  res <- list(
    f = f,
    g2n = g2n,
    nf = opt$counts$fn,
    ng = opt$counts$gr,
    par = par,
    step = step_size,
    iter = iter
  )

  if ("momentum" %in% names(opt$stages)) {
    res$mu <- opt$stages[["momentum"]]$step_size$value
    if (is.null(res$mu)) {
      res$mu <- 0
    }
  }

  res
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
  res_names <- c("f", "g2n", "nf", "ng", "step")
  if (!is.null(opt_res$mu)) {
    res_names <- c(res_names, "mu")
  }
  progress <- rbind(progress, opt_res[res_names])

  # Probably not a major performance issue to regenerate column names each time
  colnames(progress) <- res_names
  rownames(progress)[nrow(progress)] <- opt_res$iter
  progress
}

# Constructor -------------------------------------------------------------

make_opt <- function(stages,
                     verbose = FALSE) {
  opt <- list(
    init = function(opt, par, fn, gr, iter) {
      opt <- default_handler("opt", "init", opt, par, fn, gr, iter)
      for (i in 1:length(opt$stages)) {
        opt$stage_i <- i
        opt <- life_cycle_hook(opt$stages[[i]]$type, "init", opt, par, fn, gr, iter)
      }
      opt
    },
    cache = list(),
    stages = stages,
    counts = make_counts(),
    hooks = list(),
    handlers = list(),
    eager_update = FALSE,
    verbose = verbose
  )

  if (!is.null(opt$init)) {
    attr(opt$init, 'event') <- 'init opt'
    attr(opt$init, 'name') <- 'handler'
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
        #      message("emitting during ", phase)

        opt <- life_cycle_hook(phase, "during", opt, par, fn, gr, iter)
      }

      list(opt = opt)
    },
    after_stage = function(opt, stage, par, fn, gr, iter) {
      for (sub_stage_name in c("direction", "step_size")) {
        phase <- paste0(stage$type, " ", sub_stage_name)
        #      message("emitting after ", phase)
        opt <- life_cycle_hook(phase, "after", opt, par, fn, gr, iter)
      }
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
    attr(stage$init, 'name') <- 'handler'
  }
  if (!is.null(stage$calculate)) {
    attr(stage$calculate, 'event') <- paste0('during ', type)
    attr(stage$calculate, 'name') <- 'handler'
  }
  if (!is.null(stage$after_stage)) {
    attr(stage$after_stage, 'event') <- paste0('after ', type)
    attr(stage$after_stage, 'name') <- 'handler'
  }
  if (!is.null(stage$after_step)) {
    attr(stage$after_step, 'event') <- 'after step'
    attr(stage$after_step, 'name') <- paste0(type, ' after step')
  }

  res <- list()
  res[[type]] <- stage
  res
}

make_sub_stage <- function(sub_stage, type) {
  sub_stage$type <- type
  if (!is.null(sub_stage$init)) {
    attr(sub_stage$init, 'event') <- paste0('init ', sub_stage$type)
    attr(sub_stage$init, 'name') <- 'handler'
  }
  if (!is.null(sub_stage$calculate)) {
    attr(sub_stage$calculate, 'event') <- paste0('during ', sub_stage$type)
    attr(sub_stage$calculate, 'name') <- 'handler'
  }
  if (!is.null(sub_stage$after_step)) {
    attr(sub_stage$after_step, 'event') <- 'after step'
    attr(sub_stage$after_step, 'name') <-  paste0(sub_stage$type, ' after step')
  }
  sub_stage
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

add_stage <- function(opt, stage) {
  opt$stages <- c(opt$stages, stage)
  opt
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