mizer <- function(par, fn, gr,
                  method = "SD",
                  norm_direction = FALSE,
                  # L-BFGS
                  memory = 10,
                  scale_hess = TRUE,
                  # CG
                  cg_update = "PR+",
                  # Nesterov
                  nest_q = 0, # 1 - SD,
                  nest_convex_approx = FALSE,
                  # Line Search configuration
                  line_search = "MT",
                  c1 = 1e-4,
                  c2 = 0.1,
                  step0 = 1,
                  ls_initializer = "q",
                  # Momentum
                  mom_type = "classical",
                  mom_schedule = NULL,
                  mom_init = NULL,
                  mom_final = NULL,
                  mom_switch_iter = NULL,
                  mom_linear_weight = FALSE,
                  # Adaptive Restart
                  restart = NULL, # one of "fn" or "gr"
                  # Termination criterion
                  max_iter = 100,
                  max_fn = Inf,
                  max_gr = Inf,
                  max_fg = Inf,
                  rel_tol = sqrt(.Machine$double.eps),
                  grad_tol = 1e-5,
                  verbose = FALSE,
                  store_progress = FALSE) {

  opt <- make_mizer(method = method,
                    norm_direction = norm_direction,
                    scale_hess = scale_hess,
                    memory = memory,
                    cg_update = cg_update,
                    nest_q = nest_q, nest_convex_approx = nest_convex_approx,
                    line_search = line_search, step0 = step0, c1 = c1, c2 = c2,
                    ls_initializer = ls_initializer,
                    mom_type = mom_type,
                    mom_schedule = mom_schedule,
                    mom_init = mom_init,
                    mom_final = mom_final,
                    mom_switch_iter = mom_switch_iter,
                    mom_linear_weight = mom_linear_weight,
                    max_iter = max_iter,
                    restart = restart,
                    verbose = verbose)

  optloop(opt, par, fn, gr,
          max_iter = max_iter,
          max_fn = max_fn, max_gr = max_gr, max_fg = max_fg,
          rel_tol = rel_tol, grad_tol = grad_tol,
          store_progress = store_progress,
          verbose = verbose)
}


make_mizer <- function(method = "L-BFGS",
                       norm_direction = FALSE,
                       scale_hess = TRUE,
                       memory = 10,
                       cg_update = "PR+",
                       nest_q = 0,
                       nest_convex_approx = FALSE,
                       line_search = "MT",
                       c1 = 1e-4, c2 = 0.1,
                       step0 = 1,
                       ls_initializer = "q",
                       mom_type = "classical",
                       mom_schedule = NULL,
                       mom_init = NULL,
                       mom_final = NULL,
                       mom_switch_iter = NULL,
                       mom_linear_weight = FALSE,
                       max_iter = NULL,
                       restart = NULL,
                       verbose = FALSE) {
  dir_type <- NULL
  method <- toupper(method)
  if (method == "SD") {
    dir_type <- sd_direction(normalize = norm_direction)
  }
  else if (method == "CG") {
    cg_update_fn <- NULL
    cg_update <- toupper(cg_update)
    if (cg_update == "PR+") {
      cg_update_fn <- pr_plus_update
    }
    else if (cg_update == "PR") {
      cg_update_fn <- pr_update
    }
    else if (cg_update == "FR") {
      cg_update_fn <- fr_update
    }
    else if (cg_update == "DY") {
      cg_update_fn <- dy_update
    }
    else {
      stop("Unknown CG update method '", cg_update, "'")
    }
    dir_type <- cg_direction(cg_update = cg_update_fn)
  }
  else if (method == "BFGS") {
    dir_type <- bfgs_direction(scale_inverse = scale_hess)
  }
  else if (method == "L-BFGS") {
    dir_type <- lbfgs_direction(memory = memory, scale_inverse = scale_hess)
  }
  else if (method == "NAG") {
    dir_type <- sd_direction()
  }
  else if (method == "MOM") {
    dir_type <- sd_direction(normalize = norm_direction)
  }
  else if (method == "DBD") {
    dir_type <- sd_direction(normalize = norm_direction)
  }
  else {
    stop("Unknown method: '", method, "'")
  }

  step_type <- NULL
  line_search <- toupper(line_search)
  if (method == "DBD") {
    eps_init <- 1
    if (is.numeric(step0)) {
      eps_init <- step0
    }
    step_type <- delta_bar_delta(epsilon = eps_init,
                                 use_momentum = is.null(mom_schedule))
  }
  else {
    if (line_search == "MT") {
      step_type <- more_thuente_ls(c1 = c1, c2 = c2, initializer = tolower(ls_initializer))
    }
    else if (line_search == "RAS") {
      step_type <- rasmussen_ls(c1 = c1, c2 = c2, initializer = tolower(ls_initializer))
    }
    else if (line_search == "BOLD") {
      step_type <- bold_driver()
    }
    else {
      stop("Unknown line search method: '", line_search, "'")
    }
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
    opt <- append_stage(
      opt,
      momentum_stage(
        direction = nesterov_momentum_direction(),
        step_size = nest_step
      ))
  }
  if (!is.null(mom_schedule)) {
    if (is.numeric(mom_schedule)) {
      mom_step <- constant_step_size(value = mom_schedule)
    }
    else {
      mom_schedule <- tolower(mom_schedule)
      if (mom_schedule == "ramp") {
        mom_step <- make_momentum_step(
          make_ramp(max_iter = max_iter,
                    momentum_init = mom_init,
                    momentum_final = mom_final))
      }
      else if (mom_schedule == "switch") {
        mom_step <- make_momentum_step(
          make_switch(
            init_value = mom_init,
            final_value = mom_final,
            switch_iter = mom_switch_iter))
      }
      else if (mom_schedule == "nesterov") {
        if (nest_convex_approx) {
          mom_step <- nesterov_convex_approx_step()
        }
        else {
          mom_step <- nesterov_convex_step(q = nest_q)
        }
      }
    }

    if (mom_type == "classical") {
      opt <- append_stage(
        opt,
        momentum_stage(
          direction = momentum_direction(),
          step_size = mom_step
        ))
    }
    else {
      # Nesterov Momentum
      opt <- prepend_stage(
        opt,
        momentum_stage(
          direction = momentum_direction(),
          step_size = mom_step
        ))
      opt$eager_update <- TRUE
    }
    if (mom_linear_weight) {
      opt <- append_stage(opt, momentum_correction_stage())
    }

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
