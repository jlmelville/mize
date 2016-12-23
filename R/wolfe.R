# Functions for line searches

# p62 of Nocedal & Wright defines a "loose" line search as c1 = 1.e-4, c2 = 0.9
# But note that CG and SD methods are not considered suitable for loose line
# search because of the search directions are not well-scaled. c2 = 0.1 is
# suggested for CG on p34. With the Strong Wolfe conditions, reducing c2 makes
# the line search stricter (i.e. forces it closer to a minimum).

# More-Thuente ------------------------------------------------------------
more_thuente_ls <- function(c1 = c2 / 2, c2 = 0.1,
                            max_alpha_mult = 10,
                            min_step_size = .Machine$double.eps,
                            initializer = "s",
                            initial_step_length = 1,
                            try_newton_step = FALSE,
                            stop_at_min = TRUE,
                            max_fn = Inf,
                            max_gr = Inf,
                            max_fg = Inf,
                            debug = FALSE) {

  line_search(more_thuente(c1 = c1, c2 = c2,
                           max_fn = min(max_fn, max_gr, max_fg)),
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
                         max_alpha_mult = 10,
                         min_step_size = .Machine$double.eps,
                         initializer = "s",
                         initial_step_length = 1,
                         try_newton_step = FALSE,
                         stop_at_min = TRUE, eps = .Machine$double.eps,
                         max_fn = Inf,
                         max_gr = Inf,
                         max_fg = Inf,
                         debug = FALSE) {

  line_search(rasmussen(c1 = c1, c2 = c2, int = int, ext = ext,
                        max_fn = min(max_fn, max_gr, max_fg)),
              name = "rasmussen",
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
                        initializer = c("slope ratio", "quadratic"),
                        try_newton_step = FALSE,
                        initial_step_length = 1,
                        max_alpha_mult = 10,
                        min_step_size = .Machine$double.eps,
                        stop_at_min = TRUE,
                        debug = FALSE,
                        eps = .Machine$double.eps) {

  initializer <- match.arg(initializer)

  make_step_size(list(
    name = name,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Initializing Wolfe line search for ", stage$type)

      if (!is_first_stage(opt, stage)) {
        # Requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn eager updating
        #message("Wolfe: setting stage updating to eager")
        opt$eager_update <- TRUE
      }
      list(opt = opt)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {

      pm <- stage$direction$value
      if (norm2(pm) < sqrt(sub_stage$eps)) {
        sub_stage$value <- 0
        return(list(sub_stage = sub_stage))
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

      old_step_length <- sub_stage$value

      phi_alpha <- make_phi_alpha(par, fg, pm,
                                  calc_gradient_default = TRUE, debug = debug)


      # described on p59 of Nocedal and Wright
      if (initializer == "slope ratio" && !is.null(sub_stage$d0)) {
        sub_stage$value <- step_slope_ratio(old_step_length, sub_stage$d0,
                                            step0, eps, max_alpha_mult)
      }
      else if (initializer == "quadratic" && !is.null(sub_stage$f0)) {
        # quadratic interpolation
        sub_stage$value <- step_quad_interp(sub_stage$f0, step0,
                                            try_newton_step = try_newton_step)
      }

      if (is.null(sub_stage$value) || sub_stage$value <= 0) {
        sub_stage$value <- make_step_zero(initial_step_length, step0$d,
                                          try_newton_step)
      }

      sub_stage$alpha0 <- sub_stage$value
      ls_result <- ls_fn(phi_alpha, step0, sub_stage$value)
      sub_stage$d0 <- step0$d
      sub_stage$f0 <- step0$f
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

# Set the initial step length. If initial_step_length is a numeric scalar,
# then use that as-is. Otherwise, use one of several variations based around
# the only thing we know (the directional derivative)
make_step_zero <- function(initial_step_length, d0,
                           try_newton_step = FALSE) {
  if (is.numeric(initial_step_length)) {
    return(initial_step_length)
  }

  if (initial_step_length == "r") { # Rasmussen default from minimize.m
    s <- 1 / (1 - d0)
  }
  else if (initial_step_length == "s") { # scipy
    # found in _minimize_bfgs in optimize.py with this comment:
    # # Sets the initial step guess to dx ~ 1
    # actually sets f_old to f0 + 0.5 * ||g||2 then uses f_old in the quadratic
    # update formula. If you do the algebra, you get 1 / sqrt(-d)
    # (2 norm of g is sqrt(d) when starting with steepest descent)
    s <- 1 / sqrt(-d0)
  }
  else if (initial_step_length == "m") {
    # Mark Schmidt's minFunc.m uses reciprocal of the one-norm
    s <- 1 / sum(abs(d0))
  }
  else {
    stop("Unknown initial step method '", initial_step_length, "'")
  }

  if (try_newton_step) {
    s <- min(1, 1.01 * s)
  }
  s
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

    if (debug) {
      message(format_list(step))
    }
    step
  }
}

# described on p59 of Nocedal and Wright
# slope ratio method
step_slope_ratio <- function(old_step_length, d0, step0, eps, max_alpha_mult) {
  # NB the p vector must be a descent direction or the directional
  # derivative will be positive => a negative initial step size!
  slope_ratio <- d0 / (step0$d + eps)
  s <- old_step_length * min(max_alpha_mult, slope_ratio)
  max(s, eps)
}

# quadratic interpolation
step_quad_interp <- function(f0, step0, try_newton_step = FALSE) {
  s <- 2  * (step0$f - f0) / step0$d
  if (try_newton_step) {
    s <- min(1, 1.01 * s)
  }
  s
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

# Are the Strong Wolfe Conditions Met?
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


