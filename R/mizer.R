mizer <- function(par, fn, gr,
                  method = "SD",
                  memory = 10,
                  # Nesterov
                  nest_q = 0, # 1 - SD,
                  nest_convex_approx = FALSE,
                  # Line Search configuration
                  line_search = "MT",
                  c1 = 1e-4,
                  c2 = 0.1,
                  ls_initializer = "q",
                  # Adaptive Restart
                  restart = NULL, # one of "fn" or "gr"
                  # Termination criterion
                  max_iter = 100,
                  verbose = FALSE,
                  store_progress = FALSE) {

  opt <- make_mizer(method = method,
                    memory = memory,
                    nest_q = nest_q, nest_convex_approx = nest_convex_approx,
                    line_search = line_search, c1 = c1, c2 = c2,
                    ls_initializer = ls_initializer,
                    restart = restart,
                    verbose = verbose)

  optloop(opt, par, fn, gr,
          max_iter = max_iter,
          store_progress = store_progress,
          verbose = verbose)
}


make_mizer <- function(method = "L-BFGS",
                       memory = 10,
                       nest_q = 0,
                       nest_convex_approx = FALSE,
                       line_search = "MT",
                       c1 = 1e-4, c2 = 0.1,
                       ls_initializer = "q",
                       restart = NULL,
                       verbose = FALSE) {
  dir_type <- NULL
  method <- toupper(method)
  if (method == "SD") {
    dir_type <- sd_direction()
  }
  else if (method == "CG") {
    dir_type <- cg_direction()
  }
  else if (method == "BFGS") {
    dir_type <- bfgs_direction()
  }
  else if (method == "L-BFGS") {
    dir_type <- lbfgs_direction(memory = memory)
  }
  else if (method == "NAG") {
    dir_type <- sd_direction()
  }
  else {
    stop("Unknown method: '", method, "'")
  }

  step_type <- NULL
  if (line_search == "MT") {
    step_type <- more_thuente_ls(c1 = c1, c2 = c2, initializer = tolower(ls_initializer))
  }
  else if (line_search == "RAS") {
    step_type <- rasmussen_ls(c1 = c1, c2 = c2, initializer = tolower(ls_initializer))
  }
  else {
    stop("Unknown line search method: '", line_search, "'")
  }

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = dir_type,
        step_size = step_type),
      verbose = FALSE))

  if (method == "NAG") {
    if (nest_convex_approx) {
      nest_step <- nesterov_convex_approx_step()
    }
    else {
      nest_step <- nesterov_convex_step(q = nest_q)
    }
    opt <- add_stage(
      opt,
      momentum_stage(
        direction = nesterov_momentum_direction(),
        step_size = nest_step
    ))
  }

  if (!is.null(restart)) {
    restart <- tolower(restart)
    if (restart %in% c("fn", "gr")) {
      opt <- adaptive_restart(opt, restart)
    }
    else {
      stop("Unknown adaptive restart type: '", restart, "'")
    }
  }

  opt
}
