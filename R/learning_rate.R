# Adaptive algorithms, normally used in the neural network/deep learning
# communities, and normally associated with stochastic gradient descent.

# Implementation-wise treated as step size methods, but they actually modify
# both the direction and step size simultaneously to allow each parameter to
# be updated at a different rate. The direction method should always be
# steepest descent.

# Delta-Bar-Delta ----------------------------------------------------------
# A modification of Jacobs' Delta Bar Delta method is used to optimize
# the objective function in t-SNE.
#
# The heuristic is to look at successive directions of the parameter: if the
# direction is the same as the previous iteration, the minimum has yet to
# be encountered, so increase the step in that direction. If the direction
# has changed, then the minimum has been skipped, so reduce the step size.
#
# The t-SNE version differs from the version in the paper by using the update
# vector stored for use in momentum rather than storing a separate vector of
# exponentially weighted gradients.
#
# Default arguments are similar to bold driver, in that the learning rate values
# are multiplied, but for the authentic DBD experience (and also as used in
# t-SNE), you can specify kappa_fun to be `+` to add kappa to the learning rate
# when increasing the learning rate.
#
# Idea behind combining momentum with adaptive step size: slide 25 of
# http://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf
# notes that
# "Use the agreement in sign between the current gradient for a weight and the
# velocity for that weight (Jacobs, 1989)."
#
# kappa amount to increase the learning rate by.
# kappa_fun - operator to apply to kappa and the current learning rate when
#  increasing the learning rate. Set it to `+` to add to the current learning
#  rate rather than scaling up proportionally to the current value.
delta_bar_delta <- function(kappa = 1.1, kappa_fun = `*`,
                            phi = 0.5, epsilon = 1,
                            min_eps = 0,
                            theta = 0.1,
                            use_momentum = FALSE) {

  if (kappa <= 0) {
    stop("kappa must be positive")
  }
  if (!is_in_range(phi, 0, 1)) {
    stop("phi must be between 0 and 1")
  }
  if (!is_in_range(theta, 0, 1)) {
    stop("theta must be between 0 and 1")
  }

  make_step_size(list(
    name = "delta_bar_delta",
    kappa = kappa,
    kappa_fun = kappa_fun,
    phi = phi,
    min_eps = min_eps,
    theta = theta,
    epsilon = epsilon,
    use_momentum = use_momentum,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$delta_bar_old <- rep(0, length(par))
      sub_stage$gamma_old <- rep(1, length(par))
      sub_stage$value <- rep(sub_stage$init_eps, length(par))
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      delta <- opt$cache$gr_curr

      if (!is.numeric(sub_stage$epsilon)) {
        d0 <- dot(delta, stage$direction$value)
        sub_stage$epsilon <- guess_alpha0(sub_stage$epsilon,
                                          x0 = NULL,
                                          f0 = NULL,
                                          gr0 = delta,
                                          d0 = d0,
                                          try_newton_step = FALSE)
      }

      if (use_momentum && !is.null(opt$cache$update_old)) {
        # previous update includes -eps*grad_old, so reverse sign
        delta_bar_old <- -opt$cache$update_old
      }
      else {
        delta_bar_old <- sub_stage$delta_bar_old
      }
      # technically delta_bar_delta = delta_bar * delta
      # but only its sign matters, so just compare signs of delta_bar and delta
      # Force step size increase on first stage to be like the t-SNE
      # implementation
      if (all(delta_bar_old == 0)) {
        delta_bar_delta <- TRUE
      }
      else {
        delta_bar_delta <- sign(delta_bar_old) == sign(delta)
      }
      kappa <- sub_stage$kappa
      phi <- sub_stage$phi
      gamma_old <- sub_stage$gamma_old

      # signs of delta_bar and delta are the same, increase step size
      # if they're not, decrease.
      gamma <-
        (kappa_fun(gamma_old,kappa)) * abs(delta_bar_delta) +
        (gamma_old * phi) * abs(!delta_bar_delta)

      sub_stage$value <- clamp(sub_stage$epsilon * gamma,
                               min_val = sub_stage$min_eps)
      if (!use_momentum || is.null(opt$cache$update_old)) {
        theta <- sub_stage$theta
        sub_stage$delta_bar_old <- ((1 - theta) * delta) + (theta * delta_bar_old)
      }
      sub_stage$gamma_old <- gamma
      list(opt = opt, sub_stage = sub_stage)
    },
    depends = c("gradient")
  ))
}
