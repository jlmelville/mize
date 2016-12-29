# Optimizer ---------------------------------------------------------------

# Repeatedly minimizes par using opt until one of the termination conditions
# is met
opt_loop <- function(opt, par, fg, max_iter = 10, verbose = FALSE,
                    store_progress = FALSE, invalidate_cache = FALSE,
                    max_fn = Inf, max_gr = Inf, max_fg = Inf,
                    abs_tol = sqrt(.Machine$double.eps),
                    rel_tol = abs_tol, grad_tol = NULL, ginf_tol = NULL,
                    step_tol = .Machine$double.eps,
                    check_conv_every = 1, log_every = check_conv_every,
                    ret_opt = FALSE, count_res_fg = TRUE) {

  # log_every must be an integer multiple of check_conv_every
  if (!is.null(check_conv_every) && log_every %% check_conv_every != 0) {
    log_every <- check_conv_every
  }

  opt <- mize_init(opt, par, fg)
  opt$counts$max_fn <- max_fn
  opt$counts$max_gr <- max_gr
  opt$counts$max_fg <- max_fg

  progress <- data.frame()
  terminate <- list()

  # Whether and what function convergence info to calculate
  calc_fn <- (is.numeric(abs_tol) && is.finite(abs_tol)) ||
             (is.numeric(rel_tol) && is.finite(rel_tol))
  # Whether and what gradient convergence info to calculate
  calc_gr <- (is.numeric(grad_tol) && is.finite(grad_tol)) ||
             (is.numeric(ginf_tol) && is.finite(ginf_tol))
  gr_norms <- c()
  if (is.numeric(grad_tol) && is.finite(grad_tol)) {
    gr_norms <- c(gr_norms, 2)
  }
  if (is.numeric(ginf_tol) && is.finite(ginf_tol)) {
    gr_norms <- c(gr_norms, Inf)
  }

  res <- NULL

  if (verbose || store_progress) {
    res <- opt_results(opt, par, fg, 0, count_fg = count_res_fg,
                       calc_fn = calc_fn,
                       calc_gr = calc_gr, gr_norms = gr_norms)
    opt <- res$opt
    if (store_progress) {
      progress <- update_progress(opt_res = res, progress = progress)
    }
    if (verbose) {
      opt_report(res, print_time = TRUE, print_par = FALSE)
    }
  }

  best_crit <- NULL
  best_fn <- Inf
  best_grn <- Inf
  best_par <- NULL
  if (!is.null(opt$cache$fn_curr)) {
    best_crit <- "fn"
    best_fn <- opt$cache$fn_curr
    best_par <- par
  }
  else if (!is.null(opt$cache$gr_curr)) {
    best_crit <- "gr"
    best_grn <- norm_inf(opt$cache$gr_curr)
    best_par <- par
  }

  iter <- 0
  par0 <- par
  if (max_iter > 0) {
    for (iter in 1:max_iter) {

      if (invalidate_cache) {
        opt <- opt_clear_cache(opt)
      }

      par0 <- par

      # We're going to use this below to guess whether our optimization
      # requires function evaluations (this is only useful if max_fn or max_fg
      # is specified, but not really time consuming)
      if (iter == 1) {
        fn_count_before <- opt$counts$fn
      }
      step_res <- mize_step(opt, par, fg, iter)
      opt <- step_res$opt
      par <- step_res$par
      if (!is.null(opt$error)) {
        terminate$what <- opt$error
        terminate$val <- "Error"
        break
      }

      # After the first iteration, if we don't have the function available for
      # the current value of par, we probably won't have it at future iterations
      # So if we are limiting the number of function evaluations, we need to keep
      # one spare to evaluate fn after the loop finishes for when we return par
      if (iter == 1) {
        # message("has fn curr? ", has_fn_curr(opt, iter + 1))
        # message("counts before = ", fn_count_before, " counts after = ", opt$counts$fn)
        if (!has_fn_curr(opt, iter + 1)) {
          if (fn_count_before != opt$counts$fn) {
            opt$counts$max_fn <- opt$counts$max_fn - 1
          }
          opt$counts$max_fg <- opt$counts$max_fg - 1
        }
      }

      # Check termination conditions
      if (!is.null(check_conv_every) && iter %% check_conv_every == 0) {
        res <- opt_results(opt, par, fg, iter, par0, count_fg = count_res_fg,
                           calc_fn = calc_fn,
                           calc_gr = calc_gr, gr_norms = gr_norms)
        opt <- res$opt

        if (store_progress && iter %% log_every == 0) {
          progress <- update_progress(opt_res = res, progress = progress)
        }
        if (verbose && iter %% log_every == 0) {
          opt_report(res, print_time = TRUE, print_par = FALSE)
        }

        terminate <- check_termination(terminate, opt, iter = iter,
                                       step = res$step,
                                       max_fn = opt$counts$max_fn,
                                       max_gr = opt$counts$max_gr,
                                       max_fg = opt$counts$max_fg,
                                       abs_tol = abs_tol, rel_tol = rel_tol,
                                       grad_tol = grad_tol, ginf_tol = ginf_tol,
                                       step_tol = step_tol)
      }

      # might not have worked out which criterion to use on iteration 0
      if (has_fn_curr(opt, iter + 1)) {
        if (is.null(best_crit)) {
          best_crit <- "fn"
        }
        if (best_crit == "fn" && opt$cache$fn_curr < best_fn) {
          best_fn <- opt$cache$fn_curr
          best_par <- par
        }
      }
      else if (has_gr_curr(opt, iter + 1)) {
        if (is.null(best_crit)) {
          best_crit <- "gr"
        }
        if (best_crit == "gr" && norm_inf(opt$cache$gr_curr) < best_grn) {
          best_grn <- norm_inf(opt$cache$gr_curr)
          best_par <- par
        }
      }

      if (!is.null(terminate$what)) {
        break
      }
    }
  }

  # If we were keeping track of the best result and that's not currently par:
  if (!is.null(best_par)
      && ((best_crit == "fn" && best_fn != opt$cache$fn_curr) ||
          (best_crit == "gr" && best_grn != norm_inf(opt$cache$gr_curr)))) {
    par <- best_par
    opt <- opt_clear_cache(opt)
    opt <- set_fn_curr(opt, best_fn, iter + 1)
    # recalculate result for this iteration
    res <- opt_results(opt, par, fg, iter, par0, count_fg = count_res_fg,
                       calc_fn = calc_fn,
                       calc_gr = calc_gr, gr_norms = gr_norms)
    if (verbose) {
      message("Returning best result found")
    }
  }

  if (is.null(res) || res$iter != iter) {
    res <- opt_results(opt, par, fg, iter, par0, count_fg = count_res_fg,
                       calc_fn = calc_fn,
                       calc_gr = calc_gr, gr_norms = gr_norms)
    opt <- res$opt
  }
  if (verbose && iter %% log_every != 0) {
    opt_report(res, print_time = TRUE, print_par = FALSE)
  }
  if (store_progress && iter %% log_every != 0) {
    progress <- update_progress(opt_res = res, progress = progress)
  }

  if (store_progress) {
    res$progress <- progress
  }
  if (!ret_opt) {
    res["opt"] <- NULL
  }
  if (is.null(terminate$what)) {
    terminate <- list(what = "max_iter", val = max_iter)
  }
  res$terminate <- terminate[c("what", "val")]
  Filter(Negate(is.null), res)
}


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
  else if (!is.null(step) && !is.null(step_tol) && step < step_tol
            && (is.null(opt$restart_at) || opt$restart_at != iter)) {
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
            rtol <- abs(fn_old - fn_new) / min(abs(fn_new), abs(fn_old))
            if (rtol < rel_tol) {
              terminate$what <- "rel_tol"
              terminate$val <- rtol
            }
          }
        }
      }

      terminate$fn_new <- fn_new
    }
  }
  terminate
}

# Clears the cache. Results should be identical whether a cache is used or not.
opt_clear_cache <- function(opt) {
  for (name in names(opt$cache)) {
    iter_name <- paste0(name, "_iter")
    if (!is.null(opt$cache[[iter_name]])) {
      opt$cache[[iter_name]] <- "invalid"
    }
  }
  opt
}

# Creates a result object.
# If the function and gradient were not calculated as part of the optimization
# step, they WILL be calculated here, and do contribute to the total
# fn or gr count reported.
# if calc_gr is TRUE then the gradient will be calculated if it isn't
# available.
# gr_norms is a vector containing zero or more of the norms to be calculated:
#   2 for the l2 (Euclidean) norm
#   Inf for the infinity norm (max absolute component)
# Other reported results: alpha is the step size portion of the gradient
# descent stage (i.e. the result of the line search). Step is the total step
# size taken during the optimization step, including momentum.
# If a momentum stage is present, the value of the momentum is stored as 'mu'.
opt_results <- function(opt, par, fg, iter, par0 = NULL, count_fg = TRUE,
                        calc_fn = FALSE, calc_gr = FALSE, gr_norms = c()) {

  f <- NULL
  if (calc_fn || has_fn_curr(opt, iter + 1)) {
    if (!has_fn_curr(opt, iter + 1)) {
      f <- fg$fn(par)
      if (count_fg) {
        opt <- set_fn_curr(opt, f, iter + 1)
        opt$counts$fn <- opt$counts$fn + 1
      }
    }
    else {
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
    }
    else {
      g <- opt$cache$gr_curr
    }
    if (2 %in% gr_norms) {
     g2n <- norm2(g)
    }
    if (Inf %in% gr_norms) {
      ginfn <- norm_inf(g)
    }
  }

  if (!is.null(par0)) {
    step_size <- norm2(par - par0)
  }
  else {
    step_size <- 0
  }

  alpha <- NULL
  if (length(opt$stages) == 1 && opt$stages[[1]]$type == "gradient_descent") {
    alpha <- norm2(opt$stages[[1]]$step_size$value)
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
    par = par,
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

  res
}

# Prints information about the current optimization result
opt_report <- function(opt_result, print_time = FALSE, print_par = FALSE) {

  fmsg <- ""
  if (!is.null(opt_result$f)) {
    fmsg <- paste0(fmsg, " f = ", formatC(opt_result$f))
  }
  if (!is.null(opt_result$g2n)) {
    fmsg <- paste0(fmsg, " g2 = ", formatC(opt_result$g2n))
  }
  if (!is.null(opt_result$ginfn)) {
    fmsg <- paste0(fmsg, " ginf = ", formatC(opt_result$ginfn))
  }

  msg <- paste0("iter ", opt_result$iter
                , fmsg
                , " nf = ", opt_result$nf
                , " ng = ", opt_result$ng
                , " step = ", formatC(opt_result$step)
  )

  if (print_time) {
    msg <- paste(format(Sys.time(), "%H:%M:%S"), msg, collapse = " ")
  }

  if (print_par) {
    msg <- paste0(msg, " par = ", vec_formatC(opt_result$par))
  }

  message(msg)
}

# Transfers data from the result object to the progress data frame
update_progress <- function(opt_res, progress) {
  res_names <- c("f", "g2n", "ginf", "nf", "ng", "step", "alpha", "mu")
  res_names <- Filter(function(x) { !is.null(opt_res[[x]]) }, res_names)

  progress <- rbind(progress, opt_res[res_names])

  # Probably not a major performance issue to regenerate column names each time
  colnames(progress) <- res_names
  rownames(progress)[nrow(progress)] <- opt_res$iter
  progress
}

# Constructor -------------------------------------------------------------

# Creates an optimizer
make_opt <- function(stages,
                     verbose = FALSE) {
  opt <- list(
    init = function(opt, par, fg, iter) {
      opt <- default_handler("opt", "init", opt, par, fg, iter)
      for (i in 1:length(opt$stages)) {
        opt$stage_i <- i
        opt <- life_cycle_hook(opt$stages[[i]]$type, "init", opt, par, fg, iter)
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

# Creates a stage of the optimizer: a gradient_descent or momentum stage
# normally
make_stage <- function(type, direction, step_size, depends = NULL) {

  stage <- list(
    type = type,
    direction = direction,
    step_size = step_size,
    init = function(opt, stage, par, fg, iter) {
      for (sub_stage_name in c("direction", "step_size")) {
        phase <- paste0(stage$type, " ", sub_stage_name)
        opt <- life_cycle_hook(phase, "init", opt, par, fg, iter)
      }

      list(opt = opt)
    },
    calculate = function(opt, stage, par, fg, iter) {
      for (sub_stage_name in c("direction", "step_size")) {
        phase <- paste0(stage$type, " ", sub_stage_name)
        opt <- life_cycle_hook(phase, "during", opt, par, fg, iter)
      }

      list(opt = opt)
    },
    after_stage = function(opt, stage, par, fg, iter) {
      for (sub_stage_name in c("direction", "step_size")) {
        phase <- paste0(stage$type, " ", sub_stage_name)
        opt <- life_cycle_hook(phase, "after", opt, par, fg, iter)
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

# Creates a sub stage: a direction or a step size
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

# Creates a gradient_descent stage
gradient_stage <- function(direction, step_size) {
  make_stage(type = "gradient_descent", direction, step_size,
             depends = c('gradient'))
}

# Creates a momentum stage
momentum_stage <- function(direction = momentum_direction(normalize = FALSE),
                           step_size) {
  make_stage(type = "momentum", direction, step_size)
}

# Creates a momentum "correction" stage. If linear weighting is asked for, then
# mu * the gradient direction is substracted from the result.
momentum_correction_stage <- function(
  direction = momentum_correction_direction(),
  step_size = momentum_correction_step()) {
  make_stage(type = "momentum_correction", direction, step_size)
}

# Creates stages from a passed list. Should contain lists created from calling
# a specific stage function like momentum_stage or gradient_stage
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

# Add a stage to the end of an optimizer stage list
append_stage <- function(opt, stage) {
  opt$stages <- c(opt$stages, stage)
  opt
}

# Add a stage to the beginning of an optimizer stage list
prepend_stage <- function(opt, stage) {
  opt$stages <- c(stage, opt$stages)
  opt
}

# Initialize a list to store the number of times the function and gradient
# is called.
make_counts <- function() {
  list(
    fn = 0,
    gr = 0
  )
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

# Store val as fn_new for the specified iteration
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
  }
  opt
}

# Store val as fn_curr for the specified iteration
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

# Store val as gr_curr for the specified iteration
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

# Does the optimizer only have one stage (e.g. a gradient-only approach like
# BFGS)
is_single_stage <- function(opt) {
  length(opt$stages) == 1
}

# Is stage the first stage in the optimizers list of stages
is_first_stage <- function(opt, stage) {
  stage$type == opt$stages[[1]]$type
}

# Is stage the last stage in the optimizers list of stages
is_last_stage <- function(opt, stage) {
  stage$type == opt$stages[[length(opt$stages)]]$type
}

# Is the first stage of optimization gradient descent
# (i.e. not nesterov momentum)
grad_is_first_stage <- function(opt) {
  is_first_stage(opt, opt$stages[["gradient_descent"]])
}

# Has fn_new already been calculated for the specified iteration
has_fn_new <- function(opt, iter) {
  (!is.null(opt$cache$fn_new)
   && !is.null(opt$cache$fn_new_iter)
   && opt$cache$fn_new_iter == iter)
}

# Has fn_curr already been calculated for the specified iteration
has_fn_curr <- function(opt, iter) {
  (!is.null(opt$cache$fn_curr)
   && !is.null(opt$cache$fn_curr_iter)
   && opt$cache$fn_curr_iter == iter)
}

# Has gr_curr already been calculated for the specified iteration
has_gr_curr <- function(opt, iter) {
  (!is.null(opt$cache$gr_curr)
   && !is.null(opt$cache$gr_curr_iter)
   && opt$cache$gr_curr_iter == iter)
}
