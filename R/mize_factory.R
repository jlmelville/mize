mize_is_numeric_scalar <- function(value) {
  is.numeric(value) &&
    length(value) == 1 &&
    is.null(dim(value)) &&
    !is.na(value)
}

mize_validate_count <- function(
  value,
  argument,
  minimum = 0,
  allow_inf = FALSE
) {
  valid_inf <- allow_inf && mize_is_numeric_scalar(value) && value == Inf
  valid_integer <- mize_is_numeric_scalar(value) &&
    is.finite(value) &&
    value >= minimum &&
    value == floor(value)

  if (!valid_inf && !valid_integer) {
    if (allow_inf) {
      stop(
        argument,
        " must be a non-negative integer or Inf",
        call. = FALSE
      )
    }
    qualifier <- if (minimum == 0) "a finite non-negative" else "a positive"
    stop(argument, " must be ", qualifier, " integer", call. = FALSE)
  }

  invisible(value)
}

mize_validate_tolerance <- function(value, argument) {
  if (
    !is.null(value) &&
      (!mize_is_numeric_scalar(value) || !is.finite(value) || value < 0)
  ) {
    stop(
      argument,
      " must be NULL or a finite non-negative numeric scalar",
      call. = FALSE
    )
  }

  invisible(value)
}

mize_validate_convergence_controls <- function(convergence) {
  mize_validate_count(convergence$max_iter, "max_iter", allow_inf = TRUE)
  for (argument in c("max_fn", "max_gr", "max_fg")) {
    mize_validate_count(convergence[[argument]], argument, allow_inf = TRUE)
  }
  for (argument in c(
    "abs_tol",
    "rel_tol",
    "grad_tol",
    "ginf_tol",
    "step_tol"
  )) {
    mize_validate_tolerance(convergence[[argument]], argument)
  }

  invisible(convergence)
}

mize_validate_range <- function(
  value,
  argument,
  left,
  right,
  left_open = FALSE,
  right_open = FALSE
) {
  valid <- mize_is_numeric_scalar(value) && is.finite(value)
  if (valid) {
    valid_left <- if (left_open) value > left else value >= left
    valid_right <- if (right_open) value < right else value <= right
    valid <- valid_left && valid_right
  }
  if (!valid) {
    stop(
      argument,
      " must be a finite numeric scalar in the required range",
      call. = FALSE
    )
  }

  invisible(value)
}

mize_validate_positive_numeric <- function(value, argument, allow_inf = FALSE) {
  valid <- mize_is_numeric_scalar(value) &&
    value > 0 &&
    (is.finite(value) || (allow_inf && value == Inf))
  if (!valid) {
    qualifier <- if (allow_inf) {
      "a positive numeric scalar"
    } else {
      "a positive finite numeric scalar"
    }
    stop(argument, " must be ", qualifier, call. = FALSE)
  }

  invisible(value)
}

mize_validate_initial_par <- function(par) {
  if (!is.numeric(par) || !is.null(dim(par)) || length(par) == 0) {
    stop("par must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(par))) {
    stop("par must contain only finite values", call. = FALSE)
  }

  invisible(par)
}

#' Create an Optimizer
#'
#' Factory function for creating a (possibly uninitialized) optimizer.
#'
#' If the function to be optimized and starting point are not present at
#' creation time, then the optimizer should be initialized using
#' [mize_init()] before being used with [mize_step()].
#'
#' See the documentation to [mize()] for an explanation of all the
#' parameters.
#'
#' Details of the `fg` list containing the function to be optimized and its
#' gradient can be found in the 'Details' section of [mize()]. It is
#' optional for this function, but if it is passed to this function, along with
#' the vector of initial values, `par`, the optimizer will be returned
#' already initialized for this function. Otherwise, [mize_init()]
#' must be called before optimization begins.
#'
#' Additionally, optional convergence parameters may also be passed here, for
#' use with [check_mize_convergence()]. They are optional here if you
#' plan to call [mize_init()] later, or if you want to do your own
#' convergence checking.
#'
#' @param method Optimization method. See 'Details' of [mize()].
#' @param norm_direction If `TRUE`, then the steepest descent direction is
#'   normalized to unit length. Useful for adaptive step size methods where the
#'   previous step size is used to initialize the next iteration.
#' @param scale_hess if `TRUE`, the approximation to the inverse Hessian is
#'   scaled according to the method described by Nocedal and Wright
#'   (approximating an eigenvalue). Applies only to the methods `BFGS`
#'   (where the scaling is applied only during the first step) and `L-BFGS`
#'   (where the scaling is applied during every iteration). Ignored otherwise.
#' @param memory The number of updates to store if using the `L-BFGS`
#'   method. Ignored otherwise. Must be a positive integer.
#' @param cg_update Type of update to use for the `"CG"` method. For
#'   details see the "CG" subsection of the "Optimization Methods" section.
#'   Ignored if `method` is not `"CG"`.
#' @param preconditioner Type of preconditioner to use in Truncated Newton.
#'   Leave blank or set to  `"L-BFGS"` to use a limited memory BFGS
#'   preconditioner. Use the `"memory"` parameter to control the number of
#'   updates to store. Applies only if `method = "TN"`, or `"CG"`,
#'   ignored otherwise.
#' @param tn_init Type of initialization to use in inner loop of Truncated
#'   Newton. Use `0` to use the zero vector (the usual TN initialization),
#'   or `"previous"` to use the final result from the previous iteration,
#'   as suggested by Martens (2010). Applies only if `method = "TN"`,
#'   ignored otherwise.
#' @param tn_exit Type of exit criterion to use when terminating the inner CG
#'   loop of Truncated Newton method. Either `"curvature"` to use the
#'   standard negative curvature test, or `"strong"` to use the modified
#'   "strong" curvature test in TNPACK (Xie and Schlick, 1999). Applies only
#'   if `method = "TN"`, ignored otherwise.
#' @param nest_q NAG momentum parameter. Must take a value between 0 (largest
#'   momentum) and 1 (zero momentum, giving steepest descent). Ignored unless
#'   the `method` is `"NAG"` and `nest_convex_approx` is `FALSE`.
#' @param nest_convex_approx If `TRUE`, then use an approximation due to
#'   Sutskever for calculating the momentum parameter in the NAG method. Only
#'   applies if `method` is `"NAG"`.
#' @param nest_burn_in Number of iterations to wait before using a non-zero
#'   momentum. Only applies if using the `"NAG"` method or setting the
#'   `momentum_type` to "Nesterov".
#' @param step_up Value by which to increase the step size for the `"bold"`
#'   step size method or the `"DBD"` method.
#' @param step_up_fun Operator to use when combining the current step size with
#'   `step_up`. Can be one of `"*"` (to multiply the current step size
#'   with `step_up`) or `"+"` (to add).
#' @param step_down Multiplier to reduce the step size by if using the
#'   `"DBD"` method or the `"bold"`. Can also be used with the
#'   `"back"` line search method, but is optional. Should be a positive
#'   value less than 1.
#' @param dbd_weight Numeric weighting value used by the `"DBD"` method only,
#'   and only if no momentum scheme is provided. Must be between 0 and 1,
#'   inclusive.
#' @param line_search Type of line search to use. See 'Details' of
#'   [mize()].
#' @param c1 Sufficient decrease parameter for Wolfe-type line searches. Should
#'   be a value between 0 and 1.
#' @param c2 Sufficient curvature parameter for line search for Wolfe-type line
#'   searches. Should be a value between `c1` and 1.
#' @param step0 Initial value for the line search on the first step. See
#'   'Details' of [mize()].
#' @param step_next_init For Wolfe-type line searches only, how to initialize
#'   the line search on iterations after the first. See 'Details' of
#'   [mize()].
#' @param try_newton_step For Wolfe-type line searches only, try the line step
#'   value of 1 as the initial step size whenever `step_next_init` suggests
#'   a step size > 1. Defaults to `TRUE` for quasi-Newton methods such as
#'   BFGS and L-BFGS, `FALSE` otherwise.
#' @param ls_max_fn Maximum number of function evaluations allowed during a line
#'   search.
#' @param ls_max_gr Maximum number of gradient evaluations allowed during a line
#'   search.
#' @param ls_max_fg Maximum number of function or gradient evaluations allowed
#'   during a line search.
#' @param ls_max_alpha Maximum value of alpha allowed during line search. Only
#'   applies for `line_search = "more-thuente"`.
#' @param ls_max_alpha_mult The maximum value that can be attained by the ratio
#'   of the initial guess for alpha for the current line search, to the final
#'   value of alpha of the previous line search. Used to stop line searches
#'   diverging due to very large initial guesses. Only applies for Wolfe-type
#'   line searches.
#' @param ls_safe_cubic (Optional). If `TRUE`, check that cubic
#'   interpolation in the Wolfe line search does not produce too small a value.
#'   Only applies for `line_search = "more-thuente"`.
#' @param strong_curvature (Optional). If `TRUE` use the strong
#'   curvature condition in Wolfe line search. See the 'Line Search' section of
#'   [mize()] for details.
#' @param approx_armijo (Optional). If `TRUE` use the approximate Armijo
#'   condition in Wolfe line search. See the 'Line Search' section of
#'   [mize()] for details.
#' @param mom_type Momentum type, either `"classical"` or
#'   `"nesterov"`.
#' @param mom_schedule Momentum schedule. See 'Details' of [mize()].
#' @param mom_init Initial momentum value.
#' @param mom_final Final momentum value.
#' @param mom_switch_iter For `mom_schedule` `"switch"` only, the
#'   iteration when `mom_init` is changed to `mom_final`.
#' @param mom_linear_weight If `TRUE`, the gradient contribution to the
#'   update is weighted using momentum contribution.
#' @param use_init_mom If `TRUE`, then the momentum coefficient on the
#'   first iteration is non-zero. Otherwise, it's zero. Only applies if using a
#'   momentum schedule.
#' @param restart Momentum restart type. Can be one of "fn" or "gr". See
#'   'Details' of [mize()].
#' @param restart_wait Number of iterations to wait between restarts. Ignored if
#'   `restart` is `NULL`.
#' @param par (Optional) Initial values for the function to be optimized over.
#' @param fg (Optional). Function and gradient list. See 'Details' of
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
#' @export
#' @examples
#' # Function to optimize and starting point
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
#' # Create an optimizer and initialize it for use with the Rosenbrock function
#' opt <- make_mize(method = "L-BFGS", par = rb0, fg = rosenbrock_fg)
#'
#' # Create optimizer without initialization
#' opt <- make_mize(method = "L-BFGS")
#'
#' # Need to call mize_init separately:
#' opt <- mize_init(opt, rb0, rosenbrock_fg)
make_mize <- function(
  method = "L-BFGS",
  norm_direction = FALSE,
  # BFGS
  scale_hess = TRUE,
  memory = 5,
  # CG
  cg_update = "PR+",
  preconditioner = "",
  # TN
  tn_init = 0,
  tn_exit = "curvature",
  # NAG
  nest_q = 0,
  nest_convex_approx = FALSE,
  nest_burn_in = 0,
  # DBD
  step_up = 1.1,
  step_up_fun = c("*", "+"),
  step_down = NULL,
  dbd_weight = 0.1,
  # Line Search
  line_search = "More-Thuente",
  c1 = 1e-4,
  c2 = NULL,
  step0 = NULL,
  step_next_init = NULL,
  try_newton_step = NULL,
  ls_max_fn = 20,
  ls_max_gr = Inf,
  ls_max_fg = Inf,
  ls_max_alpha_mult = Inf,
  ls_max_alpha = Inf,
  ls_safe_cubic = FALSE,
  strong_curvature = NULL,
  approx_armijo = NULL,
  # Momentum
  mom_type = NULL,
  mom_schedule = NULL,
  mom_init = NULL,
  mom_final = NULL,
  mom_switch_iter = NULL,
  mom_linear_weight = FALSE,
  use_init_mom = FALSE,
  restart = NULL,
  restart_wait = 10,
  par = NULL,
  fg = NULL,
  max_iter = 100,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  abs_tol = NULL,
  rel_tol = abs_tol,
  grad_tol = NULL,
  ginf_tol = NULL,
  step_tol = NULL
) {
  # Gradient Descent Direction configuration
  dir_type <- NULL
  method <- match.arg(
    tolower(method),
    c(
      "sd",
      "newton",
      "phess",
      "cg",
      "bfgs",
      "sr1",
      "l-bfgs",
      "nag",
      "momentum",
      "dbd",
      "tn"
    )
  )

  convergence <- list(
    max_iter = max_iter,
    max_fn = max_fn,
    max_gr = max_gr,
    max_fg = max_fg,
    abs_tol = abs_tol,
    rel_tol = rel_tol,
    grad_tol = grad_tol,
    ginf_tol = ginf_tol,
    step_tol = step_tol
  )
  mize_validate_convergence_controls(convergence)

  if (method %in% c("cg", "tn")) {
    if (
      !is.character(preconditioner) ||
        length(preconditioner) != 1 ||
        is.na(preconditioner)
    ) {
      stop(
        "preconditioner must be '' or 'l-bfgs' for method ",
        toupper(method),
        call. = FALSE
      )
    }
    preconditioner <- tolower(preconditioner)
    if (!(preconditioner %in% c("", "l-bfgs"))) {
      stop(
        "preconditioner must be '' or 'l-bfgs' for method ",
        toupper(method),
        call. = FALSE
      )
    }
  }
  if (
    method == "l-bfgs" ||
      (method %in% c("cg", "tn") && preconditioner == "l-bfgs")
  ) {
    mize_validate_count(memory, "memory", minimum = 1)
  }

  switch(
    method,
    sd = {
      dir_type <- sd_direction(normalize = norm_direction)
    },
    newton = {
      dir_type <- newton_direction()
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    },
    phess = {
      dir_type <- partial_hessian_direction()
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    },
    cg = {
      cg_update <- match.arg(
        tolower(cg_update),
        c(
          "fr",
          "cd",
          "dy",
          "hs",
          "hs+",
          "pr",
          "pr+",
          "ls",
          "hz",
          "hz+",
          "prfr"
        )
      )
      cg_update_fn <- switch(
        cg_update,
        fr = fr_update,
        cd = cd_update,
        dy = dy_update,
        hs = hs_update,
        "hs+" = hs_plus_update,
        pr = pr_update,
        "pr+" = pr_plus_update,
        ls = ls_update,
        hz = hz_update,
        "hz+" = hz_plus_update,
        prfr = prfr_update
      )
      dir_type <- cg_direction(
        cg_update = cg_update_fn,
        preconditioner = preconditioner,
        memory = memory
      )
    },
    bfgs = {
      dir_type <- bfgs_direction(scale_inverse = scale_hess)
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    },
    sr1 = {
      dir_type <- sr1_direction(scale_inverse = scale_hess)
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    },
    "l-bfgs" = {
      dir_type <- lbfgs_direction(memory = memory, scale_inverse = scale_hess)
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    },
    nag = {
      dir_type <- sd_direction(normalize = norm_direction)
    },
    momentum = {
      dir_type <- sd_direction(normalize = norm_direction)
    },
    dbd = {
      dir_type <- sd_direction(normalize = norm_direction)
    },
    tn = {
      if (
        mize_is_numeric_scalar(tn_init) &&
          is.finite(tn_init) &&
          tn_init == 0
      ) {
        tn_init <- 0
      } else if (
        is.character(tn_init) &&
          length(tn_init) == 1 &&
          !is.na(tn_init) &&
          tolower(tn_init) %in% c("previous", "prev")
      ) {
        tn_init <- "previous"
      } else {
        stop(
          "tn_init must be numeric 0, 'previous', or 'prev' for method TN",
          call. = FALSE
        )
      }
      tn_exit <- match.arg(
        tolower(tn_exit),
        c("curvature", "strong")
      )

      dir_type <- tn_direction(
        init = tn_init,
        exit_criterion = tn_exit,
        preconditioner = preconditioner,
        memory = memory
      )
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    },
    stop("Unknown method: '", method, "'")
  )

  # If it's not already been turned on, turn off the Newton step option
  if (is.null(try_newton_step)) {
    try_newton_step <- FALSE
  }

  # Line Search configuration
  step_type <- NULL
  if (method == "dbd") {
    mize_validate_positive_numeric(step_up, "step_up")
    step_up_fun <- match.arg(step_up_fun)
    if (!is.null(step_down)) {
      mize_validate_range(step_down, "step_down", 0, 1)
    }
    mize_validate_range(dbd_weight, "dbd_weight", 0, 1)

    if (is.character(step0) || is.numeric(step0)) {
      eps_init <- step0
    } else {
      eps_init <- "rasmussen"
    }
    if (step_up_fun == "*") {
      step_up_fun <- `*`
    } else if (step_up_fun == "+") {
      step_up_fun <- `+`
    } else {
      stop("Unknown delta-bar-delta step_up function '", step_up_fun, "'")
    }
    if (is.null(step_down)) {
      step_down <- 0.5
    }
    step_type <- delta_bar_delta(
      epsilon = eps_init,
      kappa = step_up,
      kappa_fun = step_up_fun,
      phi = step_down,
      theta = dbd_weight,
      use_momentum = !is.null(mom_schedule)
    )
  } else {
    line_search <- match.arg(
      tolower(line_search),
      c(
        "more-thuente",
        "mt",
        "rasmussen",
        "bold driver",
        "backtracking",
        "constant",
        "schmidt",
        "minfunc",
        "armijo",
        "hager-zhang",
        "hz"
      )
    )
    if (line_search == "hager-zhang") {
      line_search <- "hz"
    }
    if (line_search == "more-thuente") {
      line_search <- "mt"
    }
    if (line_search == "minfunc") {
      line_search <- "schmidt"
    }
    if (line_search == "armijo") {
      line_search <- "backtracking"
    }

    if (method %in% c("newton", "phess", "bfgs", "l-bfgs", "tn")) {
      if (is.null(c2)) {
        c2 <- 0.9
      }
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    } else {
      if (is.null(c2)) {
        c2 <- 0.1
      }
      if (is.null(try_newton_step)) {
        try_newton_step <- FALSE
      }
    }

    wolfe_searches <- c("mt", "rasmussen", "schmidt", "hz")
    gradient_searches <- c(wolfe_searches, "backtracking")

    if (line_search %in% gradient_searches) {
      mize_validate_range(c1, "c1", 0, 1, left_open = TRUE, right_open = TRUE)
      if (line_search %in% wolfe_searches && !is.null(c2)) {
        mize_validate_range(
          c2,
          "c2",
          c1,
          1,
          left_open = TRUE,
          right_open = TRUE
        )
      }
      mize_validate_count(ls_max_fn, "ls_max_fn", allow_inf = TRUE)
      mize_validate_count(ls_max_gr, "ls_max_gr", allow_inf = TRUE)
      mize_validate_count(ls_max_fg, "ls_max_fg", allow_inf = TRUE)
      mize_validate_positive_numeric(
        ls_max_alpha_mult,
        "ls_max_alpha_mult",
        allow_inf = TRUE
      )
      if (is.numeric(step_next_init)) {
        mize_validate_positive_numeric(step_next_init, "step_next_init")
      }
    } else if (line_search == "bold driver") {
      mize_validate_positive_numeric(step_up, "step_up")
      if (!is.null(step_down)) {
        mize_validate_range(step_down, "step_down", 0, 1)
      }
      mize_validate_count(ls_max_fn, "ls_max_fn", allow_inf = TRUE)
    } else {
      if (
        !mize_is_numeric_scalar(step0) ||
          !is.finite(step0)
      ) {
        stop(
          "step0 must be a finite numeric scalar when line_search = 'constant'",
          call. = FALSE
        )
      }
    }
    if (line_search == "mt") {
      mize_validate_positive_numeric(
        ls_max_alpha,
        "ls_max_alpha",
        allow_inf = TRUE
      )
    }
    if (line_search == "backtracking" && !is.null(step_down)) {
      mize_validate_range(step_down, "step_down", 0, 1)
    }

    if (line_search == "bold driver") {
      if (is.null(step_down)) {
        step_down <- 0.5
      }
    }

    # Set Wolfe line search termination defaults
    # Most Wolfe Line Searches use the standard Strong Wolfe conditions
    if (
      line_search %in%
        c(
          "more-thuente",
          "mt",
          "rasmussen",
          "schmidt",
          "minfunc"
        )
    ) {
      if (is.null(strong_curvature)) {
        strong_curvature <- TRUE
      }
      if (is.null(approx_armijo)) {
        approx_armijo <- FALSE
      }
    }

    # Hager-Zhang uses weak Wolfe condtions with an approximation to the
    # Armijo condition. Also use the step initialization methods used in
    # CG_DESCENT by default
    if (line_search == "hz") {
      if (is.null(strong_curvature)) {
        strong_curvature <- FALSE
      }
      if (is.null(approx_armijo)) {
        approx_armijo <- TRUE
      }
      if (is.null(step_next_init)) {
        step_next_init <- "hz"
      }
      if (is.null(step0)) {
        step0 <- "hz"
      }
    } else if (line_search %in% gradient_searches) {
      if (is.null(step0)) {
        step0 <- "rasmussen"
      }
      if (is.null(step_next_init)) {
        step_next_init <- "quad"
      }
    }

    if (line_search %in% gradient_searches && !is.numeric(step_next_init)) {
      step_next_init <- tolower(step_next_init)
    }
    step_type <- switch(
      line_search,
      mt = more_thuente_ls(
        c1 = c1,
        c2 = c2,
        initializer = step_next_init,
        initial_step_length = step0,
        try_newton_step = try_newton_step,
        max_fn = ls_max_fn,
        max_gr = ls_max_gr,
        max_fg = ls_max_fg,
        max_alpha = ls_max_alpha,
        max_alpha_mult = ls_max_alpha_mult,
        strong_curvature = strong_curvature,
        approx_armijo = approx_armijo,
        safeguard_cubic = ls_safe_cubic
      ),
      rasmussen = rasmussen_ls(
        c1 = c1,
        c2 = c2,
        initializer = step_next_init,
        initial_step_length = step0,
        try_newton_step = try_newton_step,
        max_fn = ls_max_fn,
        max_gr = ls_max_gr,
        max_fg = ls_max_fg,
        max_alpha_mult = ls_max_alpha_mult,
        strong_curvature = strong_curvature,
        approx_armijo = approx_armijo
      ),
      "bold driver" = bold_driver(
        inc_mult = step_up,
        dec_mult = step_down,
        max_fn = ls_max_fn
      ),
      constant = constant_step_size(value = step0),
      schmidt = schmidt_ls(
        c1 = c1,
        c2 = c2,
        initializer = step_next_init,
        initial_step_length = step0,
        try_newton_step = try_newton_step,
        max_fn = ls_max_fn,
        max_gr = ls_max_gr,
        max_fg = ls_max_fg,
        max_alpha_mult = ls_max_alpha_mult,
        strong_curvature = strong_curvature,
        approx_armijo = approx_armijo
      ),
      backtracking = schmidt_armijo_ls(
        c1 = c1,
        initializer = step_next_init,
        initial_step_length = step0,
        try_newton_step = try_newton_step,
        step_down = step_down,
        max_fn = ls_max_fn,
        max_gr = ls_max_gr,
        max_fg = ls_max_fg,
        max_alpha_mult = ls_max_alpha_mult
      ),
      hz = hager_zhang_ls(
        c1 = c1,
        c2 = c2,
        initializer = step_next_init,
        initial_step_length = step0,
        try_newton_step = try_newton_step,
        max_fn = ls_max_fn,
        max_gr = ls_max_gr,
        max_fg = ls_max_fg,
        max_alpha_mult = ls_max_alpha_mult,
        strong_curvature = strong_curvature,
        approx_armijo = approx_armijo
      )
    )
  }

  # Create Gradient Descent stage
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = dir_type,
        step_size = step_type
      )
    )
  )

  # Momentum Configuration
  if (is.null(mom_type)) {
    mom_type <- "classical"
  }
  mom_type <- match.arg(tolower(mom_type), c("classical", "nesterov"))

  mom_direction <- momentum_direction()

  if (method == "nag") {
    # Nesterov Accelerated Gradient
    mom_type <- "classical"
    if (is.null(mom_schedule)) {
      mom_schedule <- "nsconvex"
    }
    mom_direction <- nesterov_momentum_direction()
  } else if (method == "momentum") {
    # Default momentum values
    if (mom_type == "nesterov") {
      mom_direction <- nesterov_momentum_direction()
    }
    if (is.null(mom_schedule)) {
      mom_schedule <- 0.9
    }
  }

  # Momentum configuration
  if (!is.null(mom_schedule)) {
    if (is.numeric(mom_schedule)) {
      mom_step <- make_momentum_step(
        mu_fn = make_constant(value = mom_schedule),
        use_init_mom = use_init_mom
      )
    } else if (is.function(mom_schedule)) {
      mom_step <- make_momentum_step(
        mu_fn = make_user_momentum_schedule(mom_schedule),
        use_init_mom = use_init_mom
      )
    } else {
      mom_schedule <- match.arg(
        tolower(mom_schedule),
        c("ramp", "switch", "nsconvex")
      )
      if (mom_schedule == "switch") {
        mize_validate_count(mom_switch_iter, "mom_switch_iter")
      }
      if (mom_schedule == "nsconvex") {
        mize_validate_count(nest_burn_in, "nest_burn_in")
        if (!nest_convex_approx) {
          mize_validate_range(nest_q, "nest_q", 0, 1)
        }
      }

      mom_step <- switch(
        mom_schedule,
        ramp = make_momentum_step(
          make_ramp(
            init_value = mom_init,
            final_value = mom_final,
            wait = ifelse(use_init_mom, 0, 1)
          ),
          use_init_mom = use_init_mom
        ),
        "switch" = make_momentum_step(
          make_switch(
            init_value = mom_init,
            final_value = mom_final,
            switch_iter = mom_switch_iter
          ),
          use_init_mom = use_init_mom
        ),
        nsconvex = nesterov_step(
          burn_in = nest_burn_in,
          q = if (nest_convex_approx) 0 else nest_q,
          use_approx = nest_convex_approx,
          use_init_mu = use_init_mom
        )
      )
    }

    mom_stage <- momentum_stage(
      direction = mom_direction,
      step_size = mom_step
    )

    opt <- append_stage(opt, mom_stage)

    if (mom_linear_weight) {
      opt <- append_stage(opt, momentum_correction_stage())
    }
  }

  # Adaptive Restart
  if (!is.null(restart)) {
    restart <- match.arg(tolower(restart), c("none", "fn", "gr", "speed"))
    if (restart != "none") {
      mize_validate_count(restart_wait, "restart_wait", minimum = 1)
      opt <- adaptive_restart(opt, restart, wait = restart_wait)
    }
  }

  opt$convergence <- convergence

  # Initialize for specific dataset if par and fg are provided
  if (!is.null(par) && !is.null(fg)) {
    opt <- mize_init(opt, par, fg)
  }

  opt
}
