#' One Step of Optimization
#'
#' Performs one iteration of optimization using a specified optimizer.
#'
#' This function returns both the (hopefully) optimized vector of parameters, and
#' an updated version of the optimizer itself. This is intended to be used when
#' you want more control over the optimization process compared to the more black
#' box approach of the [mize()] function. In return for having to
#' manually call this function every time you want the next iteration of
#' optimization, you gain the ability to do your own checks for convergence,
#' logging and so on, as well as take other action between iterations, e.g.
#' visualization.
#'
#' Normally calling this function should return a more optimized vector of
#' parameters than the input, or at least leave the parameters unchanged if no
#' improvement was found, although this is determined by how the optimizer was
#' configured by [make_mize()]. It is very possible to create an
#' optimizer that can cause a solution to diverge. It is the responsibility of
#' the caller to check that the result of the optimization step has actually
#' reduced the value returned from the function being optimized.
#'
#' Details of the `fg` list can be found in the 'Details' section of
#' [mize()].
#'
#' @param opt Optimizer, created by [make_mize()].
#' @param par Vector of initial values for the function to be optimized over.
#' @param fg Function and gradient list. See the documentation of
#'  [mize()].
#' @return Result of the current optimization step, a list with components:
#'
#' * `opt`: Updated version of the optimizer passed to the `opt`
#'  argument. Should be passed as the `opt` argument in the next iteration.
#' * `par`: Updated version of the parameters passed to the
#'  `par` argument. Should be passed as the `par` argument in the next
#'  iteration.
#' * `nf`: Running total number of function evaluations carried out
#'  since iteration 1.
#' * `ng`: Running total number of gradient evaluations carried out
#'  since iteration 1.
#' * `f`: Optional. The new value of the function, evaluated at the
#'  returned value of `par`. Only present if calculated as part of the
#'  optimization step (e.g. during a line search calculation).
#' * `g`: Optional. The gradient vector, evaluated at the returned
#'  value of `par`. Only present if the gradient was calculated as part of
#'  the optimization step (e.g. during a line search calculation.)
#'
#' @seealso [make_mize()] to create a value to pass to `opt`,
#'  [mize_init()] to initialize `opt` before passing it to this
#'  function for the first time. [mize()] creates an optimizer and
#'  carries out a full optimization with it.
#' @examples
#' rosenbrock_fg <- list(
#'   fn = function(x) {
#'     100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
#'   },
#'   gr = function(x) {
#'     c(
#'       -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'       200 * (x[2] - x[1] * x[1])
#'     )
#'   }
#' )
#' rb0 <- c(-1.2, 1)
#'
#' opt <- make_mize(
#'   method = "SD", line_search = "const", step0 = 0.0001,
#'   par = rb0, fg = rosenbrock_fg
#' )
#' par <- rb0
#' for (iter in 1:3) {
#'   res <- mize_step(opt, par, rosenbrock_fg)
#'   par <- res$par
#'   opt <- res$opt
#' }
#' @export
mize_step <- function(opt, par, fg) {
  opt$iter <- opt$iter + 1
  iter <- opt$iter
  opt <- life_cycle_hook("step", "before", opt, par, fg, iter)

  par0 <- par
  step_result <- NULL

  for (i in 1:length(opt$stages)) {
    opt$stage_i <- i
    stage <- opt$stages[[i]]
    opt <- life_cycle_hook(stage$type, "before", opt, par, fg, iter)
    if (!is.null(opt$terminate)) {
      break
    }
    opt <- life_cycle_hook(stage$type, "during", opt, par, fg, iter)
    if (!is.null(opt$terminate)) {
      break
    }
    opt <- life_cycle_hook(stage$type, "after", opt, par, fg, iter)
    if (!is.null(opt$terminate)) {
      break
    }

    stage <- opt$stages[[i]]

    if (is.null(step_result)) {
      step_result <- stage$result
    } else {
      step_result <- step_result + stage$result
    }

    if (opt$eager_update) {
      par <- par + stage$result
    }

    opt <- life_cycle_hook("stage", "after", opt, par, fg, iter)
    if (!is.null(opt$terminate)) {
      break
    }
  }

  if (is.null(opt$terminate)) {
    opt$ok <- TRUE
    if (!opt$eager_update) {
      par <- par + step_result
    }

    # intercept whether we want to accept the new solution
    opt <- life_cycle_hook(
      "validation",
      "before",
      opt,
      par,
      fg,
      iter,
      par0,
      step_result
    )
    opt <- life_cycle_hook(
      "validation",
      "during",
      opt,
      par,
      fg,
      iter,
      par0,
      step_result
    )
  }
  # If the this solution was vetoed or the catastrophe happened,
  # roll back to the previous one.
  if (!is.null(opt$terminate) || !opt$ok) {
    par <- par0
  }

  if (is.null(opt$terminate)) {
    opt <- life_cycle_hook(
      "step",
      "after",
      opt,
      par,
      fg,
      iter,
      par0,
      step_result
    )
  }

  res <- list(opt = opt, par = par, nf = opt$counts$fn, ng = opt$counts$gr)
  if (has_fn_curr(opt, iter + 1)) {
    res$f <- opt$cache$fn_curr
  }
  if (has_gr_curr(opt, iter + 1)) {
    res$g <- opt$cache$gr_curr
  }
  res
}

#' Initialize the Optimizer.
#'
#' Prepares the optimizer for use with a specific function and starting point.
#'
#' Should be called after creating an optimizer with [make_mize()] and
#' before beginning any optimization with [mize_step()]. Alternatively,
#' if `fg` and `par` are available when calling
#' [make_mize()], they can be passed to that function and the returned
#' optimizer will already be initialized. [mize_step()] requires an
#' initialized optimizer and does not carry out initialization itself.
#'
#' Optional convergence parameters may also be passed here, for use with
#' [check_mize_convergence()]. They are optional if you do your own
#' convergence checking.
#'
#' Details of the `fg` list can be found in the 'Details' section of
#' [mize()].
#'
#' @param opt Optimizer, created by [make_mize()].
#' @param par Vector of initial values for the function to be optimized over.
#' @param fg Function and gradient list. See the documentation of
#'   [mize()].
#' @param max_iter (Optional). Maximum number of iterations. See the
#'   'Convergence' section of [mize()] for details.
#' @param max_fn (Optional). Maximum number of function evaluations. See the
#'   'Convergence' section of [mize()] for details.
#' @param max_gr (Optional). Maximum number of gradient evaluations. See the
#'   'Convergence' section of [mize()] for details.
#' @param max_fg (Optional). Maximum number of function or gradient evaluations.
#'   See the 'Convergence' section of [mize()] for details.
#' @param abs_tol (Optional). Absolute tolerance for comparing two function
#'   evaluations. See the 'Convergence' section of [mize()] for
#'   details.
#' @param rel_tol (Optional). Relative tolerance for comparing two function
#'   evaluations. See the 'Convergence' section of [mize()] for
#'   details.
#' @param grad_tol (Optional). Absolute tolerance for the length (l2-norm) of
#'   the gradient vector. See the 'Convergence' section of [mize()]
#'   for details.
#' @param ginf_tol (Optional). Absolute tolerance for the infinity norm (maximum
#'   absolute component) of the gradient vector. See the 'Convergence' section
#'   of [mize()] for details.
#' @param step_tol (Optional). Absolute tolerance for the size of the parameter
#'   update. See the 'Convergence' section of [mize()] for details.
#' @return Initialized optimizer.
#' @export
#' @examples
#'
#' # Create an optimizer
#' opt <- make_mize(method = "L-BFGS")
#'
#' # Function to optimize and starting point defined after creating optimizer
#' rosenbrock_fg <- list(
#'   fn = function(x) {
#'     100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
#'   },
#'   gr = function(x) {
#'     c(
#'       -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'       200 * (x[2] - x[1] * x[1])
#'     )
#'   }
#' )
#' rb0 <- c(-1.2, 1)
#'
#' # Initialize with function and starting point before commencing optimization
#' opt <- mize_init(opt, rb0, rosenbrock_fg)
#'
#' # Finally, can commence the optimization loop
#' par <- rb0
#' for (iter in 1:3) {
#'   res <- mize_step(opt, par, rosenbrock_fg)
#'   par <- res$par
#'   opt <- res$opt
#' }
mize_init <- function(
  opt,
  par,
  fg,
  max_iter = Inf,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  abs_tol = NULL,
  rel_tol = abs_tol,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
) {
  opt <- register_hooks(opt)
  opt$iter <- 0
  opt <- life_cycle_hook("opt", "init", opt, par, fg, opt$iter)
  if (is.null(opt$convergence)) {
    opt$convergence <- list()
  }
  if (is.null(opt$convergence$max_iter)) {
    opt$convergence$max_iter <- max_iter
  }
  if (is.null(opt$convergence$max_fn)) {
    opt$convergence$max_fn <- max_fn
  }
  if (is.null(opt$convergence$max_gr)) {
    opt$convergence$max_gr <- max_gr
  }
  if (is.null(opt$convergence$max_fg)) {
    opt$convergence$max_fg <- max_fg
  }
  if (is.null(opt$convergence$abs_tol)) {
    opt$convergence$abs_tol <- abs_tol
  }
  if (is.null(opt$convergence$rel_tol)) {
    opt$convergence$rel_tol <- rel_tol
  }
  if (is.null(opt$convergence$grad_tol)) {
    opt$convergence$grad_tol <- grad_tol
  }
  if (is.null(opt$convergence$ginf_tol)) {
    opt$convergence$ginf_tol <- ginf_tol
  }
  if (is.null(opt$convergence$step_tol)) {
    opt$convergence$step_tol <- step_tol
  }
  opt <- opt_clear_cache(opt)
  opt$is_initialized <- TRUE
  opt
}
