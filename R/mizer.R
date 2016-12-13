mizer <- function(par, fg,
                  method = "SD",
                  norm_direction = FALSE,
                  # L-BFGS
                  memory = 10,
                  scale_hess = TRUE,
                  # CG
                  cg_update = "PR+",
                  # NAG
                  nest_q = 0, # 1 - SD,
                  nest_convex_approx = FALSE,
                  nest_burn_in = 0, use_nest_mu_zero = FALSE,
                  # DBD
                  kappa = 1.1,
                  kappa_fun = "*",
                  phi = 0.5,
                  theta = 0.1,
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
                  abs_tol = sqrt(.Machine$double.eps),
                  rel_tol = abs_tol,
                  grad_tol = 1e-5,
                  verbose = FALSE,
                  store_progress = FALSE) {

  opt <- make_mizer(method = method,
                    norm_direction = norm_direction,
                    scale_hess = scale_hess,
                    memory = memory,
                    cg_update = cg_update,
                    nest_q = nest_q, nest_convex_approx = nest_convex_approx,
                    nest_burn_in = nest_burn_in,
                    use_nest_mu_zero = use_nest_mu_zero,
                    kappa = kappa,
                    kappa_fun = kappa_fun,
                    phi = phi,
                    theta = theta,
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

  res <- opt_loop(opt, par, fg,
          max_iter = max_iter,
          max_fn = max_fn, max_gr = max_gr, max_fg = max_fg,
          abs_tol = abs_tol, rel_tol = rel_tol, grad_tol = grad_tol,
          store_progress = store_progress,
          verbose = verbose)

  res[c("f", "g2n", "nf", "ng", "par", "iter", "terminate")]
}

make_mizer <- function(method = "L-BFGS",
                       norm_direction = FALSE,
                       # BFGS
                       scale_hess = TRUE,
                       memory = 10,
                       # CG
                       cg_update = "PR+",
                       # NAG
                       nest_q = 0,
                       nest_convex_approx = FALSE,
                       nest_burn_in = 0, use_nest_mu_zero = FALSE,
                       # DBD
                       kappa = 1.1,
                       kappa_fun = "*",
                       phi = 0.5,
                       theta = 0.1,
                       # Line Search
                       line_search = "MT",
                       c1 = 1e-4, c2 = 0.1,
                       step0 = 1,
                       ls_initializer = "q",
                       # Momentum
                       mom_type = "classical",
                       mom_schedule = NULL,
                       mom_init = NULL,
                       mom_final = NULL,
                       mom_switch_iter = NULL,
                       mom_linear_weight = FALSE,
                       max_iter = NULL,
                       restart = NULL,
                       verbose = FALSE,
                       par = NULL,
                       fg = NULL) {
  dir_type <- NULL
  method <- toupper(method)
  if (method == "SD") {
    dir_type <- sd_direction(normalize = norm_direction)
  }
  else if (method == "NEWTON") {
    dir_type <- newton_direction()
  }
  else if (method == "PHESS") {
    dir_type <- partial_hessian_direction()
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
    if (kappa_fun == "*") {
      kappa_fun <- `*`
    }
    else if (kappa_fun == "+") {
      kappa_fun <- `+`
    }
    else {
      stop("Unknown delta-bar-delta kappa function '", kappa_fun, "'")
    }
    step_type <- delta_bar_delta(epsilon = eps_init,
                                 kappa = kappa, kappa_fun = kappa_fun,
                                 phi = phi, theta = theta,
                                 use_momentum = is.null(mom_schedule))
  }
  else {
    if (line_search == "MT") {
      step_type <- more_thuente_ls(c1 = c1, c2 = c2,
                                   initializer = tolower(ls_initializer),
                                   initial_step_length = step0)
    }
    else if (line_search == "RAS") {
      step_type <- rasmussen_ls(c1 = c1, c2 = c2,
                                initializer = tolower(ls_initializer),
                                initial_step_length = step0)
    }
    else if (line_search == "BOLD") {
      step_type <- bold_driver(inc_mult = kappa, dec_mult = phi)
    }
    else if (line_search == "BACK") {
      step_type <- backtracking(rho = phi, c1 = c1)
    }
    else if (line_search == "CONST") {
      step_type <- constant_step_size(value = step0)
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
      nest_step <- nesterov_convex_approx_step(burn_in = nest_burn_in,
                                               use_mu_zero = use_nest_mu_zero)
    }
    else {
      nest_step <- nesterov_convex_step(q = nest_q, burn_in = nest_burn_in)
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
                    init_value = mom_init,
                    final_value = mom_final))
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
          mom_step <- nesterov_convex_approx_step(burn_in = nest_burn_in,
                                                  use_mu_zero = use_nest_mu_zero)
        }
        else {
          mom_step <- nesterov_convex_step(q = nest_q, burn_in = nest_burn_in)
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

  # Initialize for specific dataset if par and fg are provided
  if (!is.null(par) && !is.null(fg)) {
    opt <- mizer_init(opt, par, fg)
  }

  opt
}


#' One Step of Optimization
#'
#' Performs one iteration of optimization using a specified optimizer.
#'
#' This function returns both the (hopefully) optimized vector of parameters,
#' and an updated version of the optimizer itself. This is intended to be used
#' when you want more control over the optimization process compared to the more
#' black box approach of the \code{\link{mizer}} function. In return for having
#' to manually call this function every time you want the next iteration of
#' optimization, you gain the ability to do your own checks for convergence,
#' logging and so on, as well as take other action between iterations, e.g.
#' visualization.
#'
#' Normally callng this function should return a more optimized vector of
#' parameters than the input, or at  least leave the parameters unchanged if no
#' improvement was found, although this is determined by how the optimizer was
#' configured by \code{\link{make_mizer}}. It is very possible to create an
#' optimizer that can cause a solution to diverge. It is the responsibility of
#' the caller to check that the result of the optimization step has actually
#' reduced the value returned from function being optimized.
#'
#' The function to be optimized should be passed as a list to the \code{fg}
#' parameter. This should consist of:
#' \itemize{
#' \item{\code{fn}}. The function to be optimized. Takes a vector of parameters
#'   and returns a scalar.
#' \item{\code{gr}}. The gradient of the function. Takes a vector of parameters and
#'   and returns a vector with the same length as the input parameter vector.
#' \item{\code{fg}}. Optional function which calculates the function and gradient in the
#' same routine. Takes a vector of parameters and returns a list containing
#' the function result as \code{fn} and the gradient result as \code{gr}.
#' }
#'
#' The \code{fg} function is optional, but for some methods (e.g. line search
#' methods based on the Wolfe criteria), both the function and gradient values
#' are needed for the same parameter value. Calculating them in the same
#' function can save time if there is a lot of shared work.
#'
#' @param opt Optimizer, created by \code{\link{make_mizer}}.
#' @param par Initial values for the function to be optimized over.
#' @param fg Function and gradient list. See 'Details'.
#' @param iter Current iteration number. Should increase by one each time this
#'   function is invoked.
#' @return Result of the current optimization step, a list with components:
#'\itemize{
#'  \item{\code{opt}}. Updated version of the optimizer passed to the \code{opt}
#'    argument Should be passed as the \code{opt} argument in the next
#'    iteration.
#'  \item{\code{par}}. Updated version of the parameters passed to the \code{par}
#'    argument. Should be passed as the \code{par} argument in the next
#'    iteration.
#'  \item{\code{nf}}. Running total number of function evaluations carried out since
#'    iteration 1.
#'  \item{\code{ng}}. Running total number of gradient evaluations carried out since
#'    iteration 1.
#'  \item{\code{f}}. Optional. The new value of the function, evaluated at the returned
#'    value of \code{par}. Only present if calculated as part of the
#'    optimization step (e.g. during a line search calculation).
#'  \item{\code{g2n}}. Optional. The length (2-norm) of the gradient vector, evaluated
#'    at the returned value of \code{par}. Only present if the gradient was
#'    calculated as part of the optimization step (e.g. during a line search
#'    calculation.)
#'}
#' @seealso \code{\link{make_mizer}} to create a value to pass to \code{opt},
#' \code{\link{mizer_init}} to initialize \code{opt} before passing it to this
#' function for the first time. \code{\link{mizer}} creates an optimizer and
#' carries out a full optimization with it.
#' @export
mizer_step <- function(opt, par, fg, iter) {
  opt <- life_cycle_hook("step", "before", opt, par, fg, iter)

  par0 <- par
  step_result <- NULL
  for (i in 1:length(opt$stages)) {
    opt$stage_i <- i
    stage <- opt$stages[[i]]
    opt <- life_cycle_hook(stage$type, "before", opt, par, fg, iter)
    opt <- life_cycle_hook(stage$type, "during", opt, par, fg, iter)
    opt <- life_cycle_hook(stage$type, "after", opt, par, fg, iter)

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

    opt <- life_cycle_hook("stage", "after", opt, par, fg, iter)
  }

  if (!opt$eager_update) {
    par <- par + step_result
  }

  # intercept whether we want to accept the new solution
  opt <- life_cycle_hook("validation", "before", opt, par, fg, iter,
                         par0, step_result)
  opt$ok <- TRUE
  opt <- life_cycle_hook("validation", "during", opt, par, fg, iter,
                         par0, step_result)

  # If the this solution was vetoed, roll back to the previous one.
  if (!opt$ok) {
    par <- par0
  }

  opt <- life_cycle_hook("step", "after", opt, par, fg, iter, par0,
                         step_result)

  res <- list(opt = opt, par = par, nf = opt$counts$fn, ng = opt$counts$gr)
  if (has_fn_curr(opt, iter + 1)) {
    res$f <- opt$cache$fn_curr
  }
  if (has_gr_curr(opt, iter + 1)) {
    res$g2n <- norm2(opt$cache$gr_curr)
  }

  res
}

#' Initialize the Optimizer.
#'
#' Prepares the optimizer for use with a specific function and starting point.
#'
#' The function to be optimized should be passed as a list to the \code{fg}
#' parameter. This should consist of:
#' \itemize{
#' \item{\code{fn}}. The function to be optimized. Takes a vector of parameters
#'   and returns a scalar.
#' \item{\code{gr}}. The gradient of the function. Takes a vector of parameters and
#'   and returns a vector with the same length as the input parameter vector.
#' \item{\code{fg}}. Optional function which calculates the function and gradient in the
#' same routine. Takes a vector of parameters and returns a list containing
#' the function result as \code{fn} and the gradient result as \code{gr}.
#' }
#'
#' The \code{fg} function is optional, but for some methods (e.g. line search
#' methods based on the Wolfe criteria), both the function and gradient values
#' are needed for the same parameter value. Calculating them in the same
#' function can save time if there is a lot of shared work.
#'
#' @param opt Optimizer, created by \code{\link{make_mizer}}.
#' @param par Initial values for the function to be optimized over.
#' @param fg Function and gradient list. See 'Details'.
#' @return Initialized optimizer.
#' @export
mizer_init <- function(opt, par, fg) {
  opt <- register_hooks(opt)
  opt <- life_cycle_hook("opt", "init", opt, par, fg, 0)
  opt
}
