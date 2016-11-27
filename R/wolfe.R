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
                             stop_at_min = TRUE, debug = FALSE) {
  if (!requireNamespace("rconjgrad",
                        quietly = TRUE,
                        warn.conflicts = FALSE)) {
    stop("Using More-Thuente line search requires 'rconjgrad' package")
  }
  rcg_line_search(rconjgrad::more_thuente(c1 = c1, c2 = c2),
                  name = "more-thuente",
                   max_alpha_mult = max_alpha_mult,
                   min_step_size = min_step_size, stop_at_min = stop_at_min,
                   initializer = initializer,
                   initial_step_length = initial_step_length,
                   debug = debug)
}


# Rasmussen ---------------------------------------------------------------

rasmussen_ls <- function(c1 = c2 / 2, c2 = 0.1, int = 0.1, ext = 3.0,
                          max_alpha_mult = 10,
                          min_step_size = .Machine$double.eps,
                          initializer = "s",
                          initial_step_length = 1,
                          stop_at_min = TRUE, debug = FALSE) {
  if (!requireNamespace("rconjgrad",
                        quietly = TRUE,
                        warn.conflicts = FALSE)) {
    stop("Using Rasmussen line search requires 'rconjgrad' package")
  }
  rcg_line_search(rconjgrad::rasmussen(c1 = c1, c2 = c2, int = int, ext = ext),
                  name = "rasmussen",
                  max_alpha_mult = max_alpha_mult,
                   min_step_size = min_step_size, stop_at_min = stop_at_min,
                   initializer = initializer,
                   initial_step_length = initial_step_length,
                   debug = debug)
}

rcg_line_search <- function(ls_fn,
                            name,
                            initializer = "s",
                            initial_step_length = 1,
                            max_alpha_mult = 10,
                            min_step_size = .Machine$double.eps,
                            stop_at_min = TRUE,
                            debug = FALSE,
                            eps = .Machine$double.eps) {
  make_step_size(list(
    name = name,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      #message("Initializing Wolfe line search for ", stage$type)

      if (!is_first_stage(opt, stage)) {
        # Requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn eager updating
        #message("Wolfe: setting stage updating to eager")
        opt$eager_update <- TRUE
      }
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {

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
        opt <- calc_fn(opt, par, fn)
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
      sub_stage$value <- initial_step_length

      phi_alpha <- make_phi_alpha(par, fn, gr, pm,
                                  calc_gradient_default = TRUE, debug = debug)


      # described on p59 of Nocedal and Wright
      if (initializer == "s" && !is.null(sub_stage$d0)) {
        # slope ratio method
        # NB the p vector must be a descent direction or the directional
        # derivative will be positive => a negative initial step size!
        slope_ratio <- sub_stage$d0 / (step0$d + eps)
        s <- old_step_length * min(max_alpha_mult, slope_ratio)
        sub_stage$value <- max(s, eps)
      }
      else if (initializer == "q" && !is.null(sub_stage$f0)) {
        # quadratic interpolation

        s <- 2  * (step0$f - sub_stage$f0) / step0$d
        #message("step0$f = ", step0$f, " sub_stage$f0 = ", sub_stage$f0,
        #        " step0$d = ", step0$d, " s = ", s)
        sub_stage$value <- min(1, 1.01 * s)
      }

      sub_stage$alpha0 <- sub_stage$value
      ls_result <- ls_fn(phi_alpha, step0, sub_stage$value)
      sub_stage$d0 <- step0$d
      sub_stage$f0 <- step0$f
      sub_stage$value <- ls_result$step$alpha

      opt$counts$fn <- opt$counts$fn + ls_result$nfn
      opt$counts$gr <- opt$counts$gr + ls_result$ngr

      if (is_last_stage(opt, stage)) {
#        message("wolfe: setting fn_step for iter ", iter, " fn_step = "
#                , formatC(ls_result$step$f))
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
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {

      if (is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
#        message("wolfe: setting next f_old from f for iter ", iter, " "
#                , formatC(opt$cache$fn_new)
#                )
        opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
      }

      if (is_single_stage(opt)) {
        opt <- set_gr_curr(opt, sub_stage$df, iter + 1)
      }

      list(opt = opt)
    },
    depends = c("gradient")
  ))
}

make_phi_alpha <- function(par, fn, gr, pm,
                            calc_gradient_default = FALSE, debug = FALSE) {
  function(alpha, calc_gradient = calc_gradient_default) {
    y_alpha <- par + (alpha * pm)

    f <- fn(y_alpha)
    step <- list(
      alpha = alpha,
      f = f
    )

    if (calc_gradient) {
      step$df <- gr(y_alpha)
      step$d <- dot(step$df, pm)
      if (debug) {
        message("alpha = ", formatC(alpha)
                , " f = ", formatC(f)
                , " d = ", formatC(step$d)
                , " df = ", format_vec(step$df)
                , " pm = ", format_vec(pm)
        )
      }
    }
    else if (debug) {
      message("alpha = ", formatC(alpha), " f = ", formatC(f))
    }
    step
  }
}
