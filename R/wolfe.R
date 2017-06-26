# Functions for line searches

# p62 of Nocedal & Wright defines a "loose" line search as c1 = 1.e-4, c2 = 0.9
# But note that CG and SD methods are not considered suitable for loose line
# search because of the search directions are not well-scaled. c2 = 0.1 is
# suggested for CG on p34. With the Strong Wolfe conditions, reducing c2 makes
# the line search stricter (i.e. forces it closer to a minimum).

# More-Thuente ------------------------------------------------------------
more_thuente_ls <- function(c1 = c2 / 2, c2 = 0.1,
                            max_alpha_mult = Inf,
                            min_step_size = .Machine$double.eps,
                            initializer = "s",
                            initial_step_length = 1,
                            try_newton_step = FALSE,
                            stop_at_min = TRUE,
                            max_fn = Inf,
                            max_gr = Inf,
                            max_fg = Inf,
                            approx_armijo = FALSE,
                            strong_curvature = TRUE,
                            debug = FALSE) {
  if (!is_in_range(c1, 0, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c1 must be between 0 and 1")
  }
  if (!is.null(c2) && !is_in_range(c2, c1, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c2 must be between c1 and 1")
  }
  max_ls_fn <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(more_thuente(c1 = c1, c2 = c2,
                           max_fn = max_ls_fn,
                           strong_curvature = strong_curvature,
                           approx_armijo = approx_armijo),
              name = "more-thuente",
              max_alpha_mult = max_alpha_mult,
              min_step_size = min_step_size, stop_at_min = stop_at_min,
              initializer = initializer,
              initial_step_length = initial_step_length,
              try_newton_step = try_newton_step,
              debug = debug)
}


# Rasmussen ---------------------------------------------------------------

rasmussen_ls <- function(c1 = c2 / 2, c2 = 0.1, int = 0.1, ext = 3.0,
                         max_alpha_mult = Inf,
                         min_step_size = .Machine$double.eps,
                         initializer = "s",
                         initial_step_length = 1,
                         try_newton_step = FALSE,
                         stop_at_min = TRUE, eps = .Machine$double.eps,
                         max_fn = Inf,
                         max_gr = Inf,
                         max_fg = Inf,
                         strong_curvature = TRUE,
                         approx_armijo = FALSE,
                         debug = FALSE) {
  if (!is_in_range(c1, 0, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c1 must be between 0 and 1")
  }
  if (!is.null(c2) && !is_in_range(c2, c1, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c2 must be between c1 and 1")
  }

  max_ls_fn <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(rasmussen(c1 = c1, c2 = c2, int = int, ext = ext,
                        max_fn = max_ls_fn,
                        strong_curvature = strong_curvature,
                        approx_armijo = approx_armijo),
              name = "rasmussen",
              max_alpha_mult = max_alpha_mult,
              min_step_size = min_step_size, stop_at_min = stop_at_min,
              initializer = initializer,
              initial_step_length = initial_step_length,
              try_newton_step = try_newton_step,
              eps = eps,
              debug = debug)
}


# Schmidt (minfunc) -------------------------------------------------------

schmidt_ls <- function(c1 = c2 / 2, c2 = 0.1,
                         max_alpha_mult = Inf,
                         min_step_size = .Machine$double.eps,
                         initializer = "s",
                         initial_step_length = "schmidt",
                         try_newton_step = FALSE,
                         stop_at_min = TRUE, eps = .Machine$double.eps,
                         max_fn = Inf,
                         max_gr = Inf,
                         max_fg = Inf,
                         strong_curvature = TRUE,
                         approx_armijo = FALSE,
                         debug = FALSE) {
  if (!is_in_range(c1, 0, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c1 must be between 0 and 1")
  }
  if (!is.null(c2) && !is_in_range(c2, c1, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c2 must be between c1 and 1")
  }

  max_ls_fn <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(schmidt(c1 = c1, c2 = c2, max_fn = max_ls_fn,
                      strong_curvature = strong_curvature,
                      approx_armijo = approx_armijo),
              name = "schmidt",
              max_alpha_mult = max_alpha_mult,
              min_step_size = min_step_size, stop_at_min = stop_at_min,
              initializer = initializer,
              initial_step_length = initial_step_length,
              try_newton_step = try_newton_step,
              eps = eps,
              debug = debug)
}


schmidt_armijo_ls <- function(c1 = 0.005,
                       max_alpha_mult = Inf,
                       min_step_size = .Machine$double.eps,
                       initializer = "s",
                       initial_step_length = "schmidt",
                       try_newton_step = FALSE,
                       step_down = NULL,
                       stop_at_min = TRUE, eps = .Machine$double.eps,
                       max_fn = Inf,
                       max_gr = Inf,
                       max_fg = Inf,
                       debug = FALSE) {
  if (!is_in_range(c1, 0, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c1 must be between 0 and 1")
  }

  max_ls_fn <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(schmidt_armijo_backtrack(c1 = c1, max_fn = max_ls_fn,
                                       step_down = step_down),
              name = "schmidt_armijo",
              max_alpha_mult = max_alpha_mult,
              min_step_size = min_step_size, stop_at_min = stop_at_min,
              initializer = initializer,
              initial_step_length = initial_step_length,
              try_newton_step = try_newton_step,
              eps = eps,
              debug = debug)
}


# Hager-Zhang -------------------------------------------------------------

hager_zhang_ls <- function(c1 = c2 / 2, c2 = 0.1,
                           max_alpha_mult = Inf,
                           min_step_size = .Machine$double.eps,
                           initializer = "hz",
                           initial_step_length = "hz",
                           try_newton_step = FALSE,
                           stop_at_min = TRUE, eps = .Machine$double.eps,
                           max_fn = Inf,
                           max_gr = Inf,
                           max_fg = Inf,
                           strong_curvature = FALSE,
                           approx_armijo = TRUE,
                           debug = FALSE) {
  if (!is_in_range(c1, 0, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c1 must be between 0 and 1")
  }
  if (!is.null(c2) && !is_in_range(c2, c1, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c2 must be between c1 and 1")
  }

  max_ls_fn <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(hager_zhang(c1 = c1, c2 = c2, max_fn = max_ls_fn,
                          strong_curvature = strong_curvature,
                          approx_armijo = approx_armijo),
              name = "hager-zhang",
              max_alpha_mult = max_alpha_mult,
              min_step_size = min_step_size, stop_at_min = stop_at_min,
              initializer = initializer,
              initial_step_length = initial_step_length,
              try_newton_step = try_newton_step,
              eps = eps,
              debug = debug)

}

# Line Search -------------------------------------------------------------

line_search <- function(ls_fn,
                        name,
                        initializer = "slope ratio",
                        try_newton_step = FALSE,
                        initial_step_length = 1,
                        max_alpha_mult = Inf,
                        min_step_size = .Machine$double.eps,
                        stop_at_min = TRUE,
                        debug = FALSE,
                        eps = .Machine$double.eps) {

  initializer <- match.arg(tolower(initializer),
                           c("slope ratio", "quadratic", "hz", "hager-zhang"))
  if (initializer == "hager-zhang") {
    initializer <- "hz"
  }

  if (!is.numeric(initial_step_length)) {
    initializer0 <- match.arg(tolower(initial_step_length),
                                     c("rasmussen", "scipy", "schmidt",
                                       "hz", "hager-zhang"))
    if (initializer0 == "hager-zhang") {
      initializer0 <- "hz"
    }
  }
  else {
    initializer0 <- initial_step_length
  }

  make_step_size(list(
    name = name,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Initializing Wolfe line search for ", stage$type)

      if (!is_first_stage(opt, stage)) {
        # Requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn on eager updating
        #message("Wolfe: setting stage updating to eager")
        opt$eager_update <- TRUE
      }

      sub_stage$value <- NULL
      sub_stage$alpha0 <- NULL
      sub_stage$d0 <- NULL
      sub_stage$f0 <- NULL
      sub_stage$df <- NULL

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {

      pm <- stage$direction$value
      if (norm2(pm) < .Machine$double.eps) {
        sub_stage$value <- 0
        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, opt$cache$fn_curr, iter)
        }
        return(list(opt = opt, sub_stage = sub_stage))
      }

      if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
#        message(sub_stage$name, ": fetching fn_curr from cache ", formatC(opt$cache$fn_curr))
        f0 <- opt$cache$fn_curr
      }
      else {
        opt <- calc_fn(opt, par, fg$fn)
        f0 <- opt$fn
      }

      #message("gr = ", vec_formatC(opt$cache$gr_curr), " pm = ", vec_formatC(pm))
      step0 <- list(
        alpha = 0,
        f = f0,
        df = opt$cache$gr_curr,
        d = dot(opt$cache$gr_curr, pm)
      )

      alpha_prev <- sub_stage$value

      phi_alpha <- make_phi_alpha(par, fg, pm,
                                  calc_gradient_default = TRUE, debug = debug)


      # described on p59 of Nocedal and Wright
      if (initializer == "slope ratio" && !is.null(sub_stage$d0)) {
        sub_stage$value <- step_next_slope_ratio(alpha_prev, sub_stage$d0,
                                            step0, eps, max_alpha_mult)
      }
      else if (initializer == "quadratic" && !is.null(sub_stage$f0)) {
        # quadratic interpolation
        sub_stage$value <- step_next_quad_interp(sub_stage$f0, step0,
                                            try_newton_step = try_newton_step)
      }
      else if (initializer == "hz" && !is.null(alpha_prev)) {
        step_next_res <- step_next_hz(phi_alpha, alpha_prev, step0)
        sub_stage$value <- step_next_res$alpha
        opt$counts$fn <- opt$counts$fn + step_next_res$fn
      }

      if (is.null(sub_stage$value) || sub_stage$value <= 0) {
        sub_stage$value <- guess_alpha0(initializer0,
                                        par,
                                        step0$f,
                                        step0$df,
                                        step0$d,
                                        try_newton_step)
      }

      sub_stage$alpha0 <- sub_stage$value
      sub_stage$d0 <- step0$d
      sub_stage$f0 <- step0$f

      max_fn <- Inf
      max_gr <- Inf
      max_fg <- Inf
      if (!is.null(opt$convergence$max_fn) && is.finite(opt$convergence$max_fn)) {
        max_fn <- opt$convergence$max_fn - opt$counts$fn
      }
      if (!is.null(opt$convergence$max_gr) && is.finite(opt$convergence$max_gr)) {
        max_gr <- opt$convergence$max_gr - opt$counts$gr
      }
      if (!is.null(opt$convergence$max_fg) && is.finite(opt$convergence$max_fg)) {
        max_fg <- opt$convergence$max_fg - (opt$counts$fn + opt$counts$gr)
      }

      if (max_fn <= 0 || max_gr <= 0 || max_fg <= 0) {
        sub_stage$value <- 0
        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, step0$f, iter)
          sub_stage$df <- step0$df
        }
      }
      else {
        ls_result <- ls_fn(phi_alpha, step0, sub_stage$value,
                           total_max_fn = max_fn, total_max_gr = max_gr,
                           total_max_fg = max_fg, pm = pm)
        sub_stage$value <- ls_result$step$alpha
        opt$counts$fn <- opt$counts$fn + ls_result$nfn
        opt$counts$gr <- opt$counts$gr + ls_result$ngr

        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, ls_result$step$f, iter)
          if (is.null(ls_result$step$df)) {
            sub_stage$df <- rep(sub_stage$eps, length(par))
          }
          else {
            sub_stage$df <- ls_result$step$df
          }
        }
      }

      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {
      if (opt$ok && is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
        opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
      }

      if (opt$ok && is_single_stage(opt)) {
        opt <- set_gr_curr(opt, sub_stage$df, iter + 1)
      }

      list(opt = opt)
    },
    depends = c("gradient")
  ))
}

make_phi_alpha <- function(par, fg, pm,
                           calc_gradient_default = FALSE, debug = FALSE) {
  # LS functions are responsible for updating fn and gr count
  function(alpha, calc_gradient = calc_gradient_default) {
    y_alpha <- par + (alpha * pm)

    if (!is.null(fg$fg) && calc_gradient) {
      fg_res <- fg$fg(y_alpha)
      f <- fg_res$fn
      g <- fg_res$gr

      step <- list(
        alpha = alpha,
        f = f,
        df = g,
        d = dot(g, pm)
      )
    }

    else {
      f <- fg$fn(y_alpha)
      step <- list(
        alpha = alpha,
        f = f
      )

      if (calc_gradient) {
        step$df <- fg$gr(y_alpha)
        step$d <- dot(step$df, pm)
      }
    }

    step$par <- y_alpha

    if (debug) {
      message(format_list(step))
    }
    step
  }
}

# Ensure Valid Step Size
#
# Given an initial step size, if either the function value or the directional
# derivative is non-finite (NaN or infinite), reduce the step size until
# finite values are found.
#
# @param phi Line function.
# @param alpha Initial step size.
# @param min_alpha Minimum step size.
# @param max_fn Maximum number of function evaluations allowed.
# @return List containing:
# \itemize{
#   \item step Valid step size or the last step size evaluated, or NULL if
#     max_fn == 0.
#   \item nfn Number of function evaluations.
#   \item ok If a valid step was found
# }
find_finite <- function(phi, alpha, min_alpha = 0, max_fn = 20) {
  nfn <- 0
  ok <- FALSE
  step <- NULL
  while (nfn < max_fn && alpha >= min_alpha) {
    step <- phi(alpha)
    nfn <- nfn + 1
    if (step_is_finite(step)) {
      ok <- TRUE
      break
    }
    alpha <- (min_alpha + alpha) / 2
  }
  list(step = step, nfn = nfn, ok = ok)
}

step_is_finite <- function(step) {
  is.finite(step$f) && is.finite(step$df)
}


# Initial Step Length -----------------------------------------------------

# Set the initial step length. If initial_step_length is a numeric scalar,
# then use that as-is. Otherwise, use one of several variations based around
# the only thing we know (the directional derivative)
guess_alpha0 <- function(guess_type, x0, f0, gr0, d0,
                           try_newton_step = FALSE) {
  if (is.numeric(guess_type)) {
    return(guess_type)
  }

  s <- switch(guess_type,
    rasmussen = step0_rasmussen(d0),
    scipy = step0_scipy(gr0, d0),
    schmidt = step0_schmidt(gr0),
    hz = step0_hz(x0, f0, gr0, psi0 = 0.01)
  )

  if (try_newton_step) {
    s <- min(1, 1.01 * s)
  }
  s
}

# From minimize.m
step0_rasmussen <- function(d0) {
  1 / (1 - d0)
}

# found in _minimize_bfgs in optimize.py with this comment:
# # Sets the initial step guess to dx ~ 1
# actually sets f_old to f0 + 0.5 * ||g||2 then uses f_old in the quadratic
# update formula. If you do the algebra,  you get -||g||2 / d
# Assuming steepest descent for step0, this can be simplified further to
# 1 / sqrt(-d0), but may as well not assume that
step0_scipy <- function(gr0, d0) {
  -norm2(gr0) / d0
}

# Mark Schmidt's minFunc.m uses reciprocal of the one-norm
step0_schmidt <- function(gr0) {
  1 / norm1(gr0)
}

# I0 in the 'initial' routine in CG_DESCENT paper
step0_hz <- function(x0, f0, gr0, psi0 = 0.01) {
  alpha <- 1
  if (is.null(x0)) {
    return(alpha)
  }
  ginf_norm <- norm_inf(gr0)
  if (ginf_norm != 0) {
    xinf_norm <- norm_inf(x0)
    if (xinf_norm != 0) {
      alpha <- psi0 * (xinf_norm / ginf_norm)
    }
    else if (is_finite_numeric(f0) && f0 != 0) {
      g2_norm2 <- sqnorm2(gr0)
      if (is_finite_numeric(g2_norm2) && g2_norm2 != 0) {
        alpha <- psi0 * (abs(f0) / ginf_norm ^ 2)
      }
    }
  }
  alpha
}

# Next Step Length --------------------------------------------------------

# described on p59 of Nocedal and Wright
# slope ratio method
step_next_slope_ratio <- function(alpha_prev, d0_prev, step0, eps,
                                  max_alpha_mult) {
  # NB the p vector must be a descent direction or the directional
  # derivative will be positive => a negative initial step size!
  slope_ratio <- d0_prev / (step0$d + eps)
  s <- alpha_prev * min(max_alpha_mult, slope_ratio)
  max(s, eps)
}

# quadratic interpolation
step_next_quad_interp <- function(f0_prev, step0, try_newton_step = FALSE) {
  s <- 2  * (step0$f - f0_prev) / step0$d
  if (try_newton_step) {
    s <- min(1, 1.01 * s)
  }
  s
  max(.Machine$double.eps, s)
}

# steps I1-2 in the routine 'initial' of the CG_DESCENT paper
step_next_hz <- function(phi, alpha_prev, step0, psi1 = 0.1, psi2 = 2) {
  if (alpha_prev < .Machine$double.eps) {
    return(list(alpha = .Machine$double.eps, fn = 0))
  }

  # I2: use if QuadStep fails
  alpha <- alpha_prev * psi2
  nfn <- 0
  # I1: QuadStep
  # If function is reduced at the initial guess, carry out a quadratic
  # interpolation. If the resulting quadratic is strongly convex, use the
  # minimizer of the quadratic. Otherwise, use I2.
  step_psi <- phi(alpha_prev * psi1, calc_gradient = FALSE)
  nfn <- 1
  if (step_psi$f <= step0$f) {
    alpha_q <- quadratic_interpolate_step(step0, step_psi)
    # A 1D quadratic Ax^2 + Bx + C is strongly convex if A > 0. Second clause in
    # if statement is A expressed in terms of the minimizer (this is easy to
    # derive by looking at Nocedal & Wright 2nd Edition, equations 3.57 and
    # 3.58)
    if (alpha_q > 0 && -step0$d / (2 * alpha_q) > 0) {
      alpha <- alpha_q
    }
  }
  alpha <- max(.Machine$double.eps, alpha)
  list(alpha = alpha, fn = nfn)
}

# Line Search Checks -------------------------------------------------------

# Armijo Rule (or Sufficient Decrease Condition)
#
# Line search test.
#
# The sufficient decrease condition is met if the line search step length yields
# a decrease in the function value that is sufficiently large (relative to the
# size of the step).
#
# This test prevents large step sizes that, while representing a function value
# decrease, don't reduce it by very much, which could indicate that the
# function minimum has been stepped over and you're now going back up the slope.
# Also, this condition can always be met by taking a sufficiently small step,
# so line searches involving this condition can always terminate. The downside
# is that you can end up taking very small steps, so it's usual to combine this
# condition with one that encourages larger step sizes.
#
# @param f0 Value of function at starting point of line search.
# @param d0 Directional derivative at starting point of line search.
# @param alpha the step length.
# @param fa Value of function at alpha.
# @param c1 the sufficient decrease constant. Should take a value between 0 and
#   1.
# @return \code{TRUE} if the step \code{alpha} represents a sufficient decrease.
armijo_ok <- function(f0, d0, alpha, fa, c1) {
  fa <= f0 + c1 * alpha * d0
}

# Armijo Rule (or Sufficient Decrease Condition)
#
# Line search test.
#
# The sufficient decrease condition is met if the line search step length yields
# a decrease in the function value that is sufficiently large (relative to the
# size of the step).
#
# This test prevents large step sizes that, while representing a function value
# decrease, don't reduce it by very much, which could indicate that the
# function minimum has been stepped over and you're now going back up the slope.
# Also, this condition can always be met by taking a sufficiently small step,
# so line searches involving this condition can always terminate. The downside
# is that you can end up taking very small steps, so it's usual to combine this
# condition with one that encourages larger step sizes.
#
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param c1 the sufficient decrease constant. Should take a value between 0 and
#   1.
# @return \code{TRUE} if the step represents a sufficient decrease.
armijo_ok_step <- function(step0, step, c1) {
  armijo_ok(step0$f, step0$d, step$alpha, step$f, c1)
}


# Curvature Condition
#
# Line search test.
#
# Ensures that the directional derivative of the line search direction at a
# candidate step size is greater than a specified fraction of the slope of the
# line at the starting point of the search. This condition is used to stop step
# sizes being too small.
#
# In combination with the sufficient decrease condition \code{\link{armijo_ok}}
# these conditions make up the Wolfe conditions.
#
# @param d0 Directional derivative at starting point.
# @param da Directional derivative at step alpha.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
curvature_ok <- function(d0, da, c2) {
  da > c2 * d0
}

# Curvature Condition
#
# Line search test.
#
# Ensures that the directional derivative of the line search direction at a
# candidate step size is greater than a specified fraction of the slope of the
# line at the starting point of the search. This condition is used to stop step
# sizes being too small.
#
# In combination with the sufficient decrease condition \code{\link{armijo_ok}}
# these conditions make up the Wolfe conditions.
#
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
curvature_ok_step <- function(step0, step, c2) {
  curvature_ok(step0$d, step$d, c2)
}

# Strong Curvature Condition
#
# Line search test.
#
# Ensures that the value of the directional derivative of the line search
# direction at a candidate step size is equal to or greater than a specified
# fraction of the slope of the line at the starting point of the search, while
# having the same direction. This condition is used to make the step size lie
# close to a stationary point. Unlike the normal curvature condition, a step
# size where the sign of the gradient changed (e.g. the minimum had been
# skipped) would not be acceptable for the strong curvature condition.
#
# In combination with the sufficient decrease condition \code{\link{armijo_ok}}
# these conditions make up the Strong Wolfe conditions.
#
# @param d0 Directional derivative at starting point.
# @param da Directrional derivative at step alpha.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the strong curvature condition is met.
strong_curvature_ok <- function(d0, da, c2) {
  abs(da) <= -c2 * d0
}

# Strong Curvature Condition
#
# Line search test.
#
# Ensures that the value of the directional derivative of the line search
# direction at a candidate step size is equal to or greater than a specified
# fraction of the slope of the line at the starting point of the search, while
# having the same direction. This condition is used to make the step size lie
# close to a stationary point. Unlike the normal curvature condition, a step
# size where the sign of the gradient changed (e.g. the minimum had been
# skipped) would not be acceptable for the strong curvature condition.
#
# In combination with the sufficient decrease condition \code{\link{armijo_ok}}
# these conditions make up the Strong Wolfe conditions.
#
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
strong_curvature_ok_step <- function(step0, step, c2) {
  strong_curvature_ok(step0$d, step$d, c2)
}

# Are the Strong Wolfe Conditions Met?
#
# Step size check.
#
# Returns true if the Strong Wolfe conditions are met, consisting of the
# sufficient decrease conditions and the strong curvature condition.
#
# @param f0 Function value at starting point.
# @param d0 Directional derivative value at starting point.
# @param alpha Step length.
# @param fa Function value at alpha.
# @param da Directional derivative at alpha.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @return TRUE if the Strong Wolfe condition is met by the step size.
strong_wolfe_ok <- function(f0, d0, alpha, fa, da, c1, c2) {
  armijo_ok(f0, d0, alpha, fa, c1) &&
    strong_curvature_ok(d0, da, c2)
}

# Are the Strong Wolfe Conditions Met for the Given Step?
#
# Line search test.
#
# Returns true if the candidate step size meets the Strong Wolfe conditions,
# consisting of the sufficient decrease conditions and the strong curvature
# condition.
#
# @param step0 Line search values at starting point of line search.
# @param step Line search value at candiate step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @return TRUE if the Strong Wolfe condition is met by the step size.
strong_wolfe_ok_step <- function(step0, step, c1, c2) {
  armijo_ok_step(step0, step, c1) && strong_curvature_ok_step(step0, step, c2)
}

# Are the Wolfe Conditions Met for the Given Step?
wolfe_ok_step <- function(step0, step, c1, c2) {
  armijo_ok_step(step0, step, c1) && curvature_ok_step(step0, step, c2)
}

# Create a Wolfe Conditions check function allowing for either approximate or
# exact Armijo condition and weak or strong curvature condition
make_wolfe_ok_step_fn <- function(approx_armijo = FALSE,
                                  strong_curvature = TRUE, eps = 1e-6) {

  approx_armijo_ok_fn <- make_approx_armijo_ok_step(eps)

  function(step0, step, c1, c2) {
    if (approx_armijo) {
      ok <- approx_armijo_ok_fn(step0, step, c1)
    }
    else {
      ok <- armijo_ok_step(step0, step, c1)
    }

    if (ok) {
      if (strong_curvature) {
        ok <- strong_curvature_ok_step(step0, step, c2)
      }
      else {
        ok <- curvature_ok_step(step0, step, c2)
      }
    }
    ok
  }
}

# Create Approximation Armijo check function:
# Checks approximate Armijo conditions only if exact Armijo check fails and
# if function value is sufficiently close to step0 value according to eps
make_approx_armijo_ok_step <- function(eps) {
  function(step0, step, c1) {
    eps_k <- eps * abs(step0$f)

    if (armijo_ok_step(step0, step, c1)) {
      return(TRUE)
    }
    (step$f <= step0$f + eps_k) && approx_armijo_ok_step(step0, step, c1)
  }
}

# Is the approximate Armijo condition met?
#
# Suggested by Hager and Zhang (2005) as part of the Approximate Wolfe
# Conditions. The second of these conditions is identical to the (weak)
# curvature condition.
#
# The first condition applies the armijo condition to a quadratic approximation
# to the function, which allows for higher precision in finding the minimizer.
#
# It is suggested that the approximate version of the Armijo condition be used
# when fa is 'close' to f0, e.g. fa <= f0 + eps * |f0| where eps = 1e-6
#
# c1 should be < 0.5
approx_armijo_ok <- function(d0, da, c1) {
  (2 * c1 - 1) * d0 >= da
}

# Is the approximate Armijo condition met for the given step?
approx_armijo_ok_step <- function(step0, step, c1) {
  approx_armijo_ok(step0$d, step$d, c1)
}


# Bracket -----------------------------------------------------------------

bracket_is_finite <- function(bracket) {
  all(is.finite(c(bracket_props(bracket, c('f', 'd')))))
}

# extracts all the properties (e.g. 'f', 'df', 'd' or 'alpha') from all members
# of the bracket. Works if there are one or two bracket members. Can get
# multiple properties at once, by providing an array of the properties,
# e.g. bracket_props(bracket, c('f', 'd'))
bracket_props <- function(bracket, prop) {
  unlist(sapply(bracket, `[`, prop))
}

bracket_width <- function(bracket) {
  bracket_range <- bracket_props(bracket, 'alpha')
  abs(bracket_range[2] - bracket_range[1])
}

bracket_min_alpha <- function(bracket) {
  min(bracket_props(bracket, 'alpha'))
}

best_bracket_step <- function(bracket) {
  LOpos <- which.min(bracket_props(bracket, 'f'))
  bracket[[LOpos]]
}


is_in_bracket <- function(bracket, alpha) {
  is_in_range(alpha, bracket[[1]]$alpha, bracket[[2]]$alpha)
}

format_bracket <- function(bracket) {
  paste0("[", formatC(bracket[[1]]$alpha), ", ", formatC(bracket[[2]]$alpha),
         "]")
}
