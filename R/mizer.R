#' Numerical Optimization
#'
#' Numerical optimization including conjugate gradient,
#' Broyden-Fletcher-Goldfarb-Shanno (BFGS), and the limited memory BFGS.
#'
#' The function to be optimized should be passed as a list to the \code{fg}
#' parameter. This should consist of:
#' \itemize{
#' \item{\code{fn}}. The function to be optimized. Takes a vector of parameters
#'   and returns a scalar.
#' \item{\code{gr}}. The gradient of the function. Takes a vector of parameters
#' and returns a vector with the same length as the input parameter vector.
#' \item{\code{fg}}. Optional function which calculates the function and
#' gradient in thesame routine. Takes a vector of parameters and returns a list
#' containing the function result as \code{fn} and the gradient result as
#' \code{gr}.
#' }
#'
#' The \code{fg} function is optional, but for some methods (e.g. line search
#' methods based on the Wolfe criteria), both the function and gradient values
#' are needed for the same parameter value. Calculating them in the same
#' function can save time if there is a lot of shared work.
#'
#' The \code{method} specifies the optimization method:
#'
#' \itemize{
#' \item \code{"SD"} is plain steepest descent. Not very effective on its own,
#' but can be combined with various momentum approaches.
#' \item \code{"BFGS"} is the Broyden-Fletcher-Goldfarb-Shanno quasi-Newton
#' method. This stores an approximation to the inverse of the Hessian of the
#' function being minimized, which requires storage proportional to the
#' square of the length of \code{par}, so is unsuitable for large problems.
#' \item \code{"L-BFGS"} is the Limited memory Broyden-Fletcher-Goldfarb-Shanno
#' quasi-Newton method. This does not store the inverse Hessian approximation
#' directly and so can scale to larger-sized problems than \code{"BFGS"}. The
#' amount of memory used can be controlled with the \code{memory} parameter.
#' \item \code{"CG"} is the conjugate gradient method. The \code{cg_update}
#' parameter allows for different methods for choosing the next direction:
#'   \itemize{
#'     \item \code{"FR"} The method of Fletcher and Reeves.
#'     \item \code{"PR"} The method of Polak and Ribiere.
#'     \item \code{"PR+"} The method of Polak and Ribiere with a restart to
#'     steepest descent if conjugacy is lost. The default.
#'     \item \code{"HS"} The method of Hestenes and Stiefel.
#'     \item \code{"DY"} The method of Dai and Yuan.
#'   }
#' \item \code{"NAG"} is the Nesterov Accelerated Gradient method. The exact
#' form of the momentum update in this method can be controlled with the
#' following parameters:
#'   \itemize{
#'   \item{\code{nest_q}} Strong convexity parameter. Must take a value
#'   between 0 (strongly convex) and 1 (zero momentum). Ignored if
#'   \code{nest_convex_approx} is \code{TRUE}.
#'   \item{\code{nest_convex_approx}} If \code{TRUE}, then use an approximation
#'   due to Sutskever for calculating the momentum parameter.
#'   \item{\code{use_nest_mu_zero}} If \code{TRUE}, then the momentum on
#'   iteration zero is set to 0.4. Otherwise, it's zero. Ignored if
#'   \code{nest_convex_approx} is \code{FALSE}.
#'   \item{\code{nest_burn_in}} Number of iterations to wait before using a
#'   non-zero momentum.
#'   }
#' \item \code{"DBD"} is the Delta-Bar-Delta method of Jacobs.
#' \item \code{"MOM"} is steepest descent with momentum. See below for momentum
#' options.
#' }
#'
#' For more details on gradient-based optimization in general, and the BFGS,
#' L-BFGS and CG methods, see Nocedal and Wright.
#'
#' The parameter \code{line_search} determines the line search to be carried
#' out:
#'
#' \itemize{
#'   \item If a numeric scalar is provided, then a constant value will be used
#'   for the line search. Note that this value will be multiplied by the
#'   magnitude of the direction vector used in the gradient descent method.
#'   For method \code{"SD"} only, setting the \code{norm_direction} parameter to
#'   \code{TRUE} will scale the direction vector so it has unit length.
#'   \item \code{"RAS"} carries out a line search using the strong Wolfe
#'   conditions and the method of Rasmussen.
#'   \item \code{"MT"} carries out a line search using the strong Wolfe
#'   conditions and the method of More-Thuente.
#'   \item \code{"BOLD"} carries out a back tracking line search until a
#'   reduction in the function value is found.
#' }
#'
#' If using one of the methods: \code{"BFGS"}, \code{"L-BFGS"}, \code{"CG"} or
#' \code{"NAG"}, one of the Wolfe line searches, \code{"RAS"} or \code{"MT"}
#' should be used, otherwise very poor performance is likely to be encountered.
#' The following parameters can be used to control the line search:
#'
#'  \itemize{
#'    \item{\code{c1}} The sufficient decrease condition. Normally left at its
#'    default value of 1e-4.
#'    \item{\code{c2}} The sufficient curvature condition. Defaults to 0.9 if
#'    using the methods \code{"BFGS"} and \code{"L-BFGS"}, and to 0.1 for
#'    every other method, more or less in line with the recommendations given
#'    by Nocedal and Wright. The smaller the value of \code{c2}, the stricter
#'    the line search, but it should not be set to smaller than \code{c1}.
#'    \item{\code{step0}} Initial value for the line search on the first step.
#'    If a positive numeric value is passed as an argument, that value is used
#'    as-is. Otherwise, by passing a character as an argument, a guess is made
#'    based on the gradient at the starting point:
#'    \itemize{
#'      \item{\code{"r"}} As used by Rasmussen in \code{minimize.m}:
#'      \deqn{\frac{1}{1+\left|g\right|^2}}{1 / 1 + (|g|^2)}
#'      \item{\code{"s"}} As used in scipy's \code{optimize.py}
#'      \deqn{\frac{1}{\left|g\right|}}{1 / |g|}
#'      \item{\code{"m"}} As used by Schmidt in \code{minFunc.m}
#'      (the reciprocal of the l1 norm of g)
#'      \deqn{\frac{1}{\left|g\right|_1}}{1 / |g|1}
#'    }
#'    \item{\code{ls_initializer}} How to initialize subsequent line searches
#'    after the first, using results from the previous line search,
#'    based on two suggestions mentioned by Nocedal and Wright:
#'    \itemize{
#'      \item{\code{"r"}} Slope ratio method.
#'      \item{\code{"q"}} Quadratic interpolation method.
#'    }
#'    \item{\code{try_newton_step}} For quasi-Newton methods (\code{"BFGS"} and
#'    \code{"L-BFGS"}), setting this to \code{TRUE} will try the "natural" step
#'    size of 1, whenever the \code{ls_initializer} method suggests an initial
#'    step size larger than that. On by default for BFGS and L-BFGS, off for
#'    everything else.
#' }
#'
#' If the \code{"DBD"} is used for the optimization \code{"method"}, then the
#' \code{line_search} parameter is ignored, because this method controls both
#' the direction of the search and the step size simultaneously. The following
#' parameters can be used to control the step size:
#'
#' \itemize{
#'   \item{\code{kappa}} The amount by which to increase the step size in a
#'   direction where the current step size is deemed to be too short. This
#'   should be a positive scalar.
#'   \item{\code{phi}} The amount by which to decrease the step size in a
#'   direction where the currents step size is deemed to be too long. This
#'   should be a positive scalar smaller than 1.
#'   \item{\code{kappa_fun}} How to increase the step size: either the method of
#'   Jacobs (addition of \code{kappa}) or Janet and co-workers (multiplication
#'   by \code{kappa}). Note that the step size decrease \code{phi} is always
#'   a multiplication.
#' }
#'
#' The \code{"BOLD"} line search also uses the \code{kappa} and \code{phi}
#' parameters with similar meanings to their use with the \code{"DBD"} method:
#' the backtracking portion reduces the step size by a factor of \code{phi}.
#' Once a satisfactory step size has been found, the line search for the
#' next iteration is initialized by multiplying the previously found step size
#' by \code{kappa}.
#'
#' For \code{method} \code{"MOM"}, momentum schemes can be accessed through the
#' momentum arguments:
#'
#' \itemize{
#' \item{\code{mom_type}} Momentum type, either \code{"classical"} or
#'  \code{"nesterov"}. Using "Nesterov" applies the momentum step before the
#'   gradient descent as suggested by Sutskever, emulating the behavior of the
#'   Nesterov Accelerated Gradient method.
#' \item{\code{mom_schedule}} How the momentum changes over the course of the
#'   optimization:
#'   \itemize{
#'   \item{If a numerical scalar is provided, a constant momentum will be
#'     applied throughout.}
#'   \item{\code{"nesterov"}} Use the momentum schedule from the Nesterov
#'   Accelerated Gradient method. Parameters which control the NAG momentum
#'   can also be used in combination with this option.
#'   \item{\code{"switch"}} Switch from one momentum value (specified via
#'   \code{mom_init}) to another (\code{mom_final}) at a
#'   a specified iteration (\code{mom_switch_iter}).
#'   \item{\code{"ramp"}} Linearly increase from one momentum value
#'   (\code{mom_init}) to another (\code{mom_final}) over the specified
#'   period (\code{max_iter}).
#'   }
#' }
#'
#' The \code{restart} parameter provides a way to restart the momentum if the
#' optimization appears to be not be making progress, using the method of
#' O'Donoghue and Candes. There are two strategies:
#' \itemize{
#'   \item{\code{"fn"}} A restart is applied if the function does not decrease
#'   on consecutive iterations.
#'   \item{\code{"gr"}} A restart is applied if the direction of the
#'   optimization is not a descent direction.
#' }
#'
#' The effect of the restart is to "forget" any previous momentum update vector,
#' and, for those momentum schemes that change with iteration number, to
#' effectively reset the iteration number back to zero. If the \code{mom_type}
#' is \code{"nesterov"}, the gradient-based restart is not available.
#'
#' If \code{method} type \code{"MOM"} is specified with no other values, the
#' momentum scheme will default to a constant value of \code{0.9}, with a
#' function-based restart.
#'
#' @param par Initial values for the function to be optimized over.
#' @param fg Function and gradient list. See 'Details'.
#' @param method Optimization method. See 'Details'.
#' @param norm_direction If \code{TRUE}, then the steepest descent direction
#' is normalized to unit length. Useful for adaptive step size methods where
#' the previous step size is used to initialize the next iteration.
#' @param scale_hess if \code{TRUE}, the approximation to the inverse Hessian
#' is scaled according to the method described by Nocedal and Wright
#' (approximating an eigenvalue). Applies only to the methods \code{BFGS}
#' (where the scaling is applied only during the first step) and \code{L-BFGS}
#' (where the scaling is applied during every iteration). Ignored otherwise.
#' @param memory The number of updates to store if using the \code{L-BFGS}
#' method. Ignored otherwise. Must be a positive integer.
#' @param cg_update Type of update to use for the \code{CG} method. Can be
#' one of \code{"FR"} (Fletcher-Reeves), \code{"PR"} (Polak-Ribiere),
#' \code{"PR+"} (Polak-Ribiere with a reset to steepest descent), \code{"HS"}
#' (Hestenes-Stiefel), or \code{"DY"} (Dai-Yuan). Ignored if \code{method} is
#' not \code{"CG"}.
#' @param nest_q Strong convexity parameter for the NAG
#' momentum term. Must take a value between 0 (strongly convex) and 1
#' (zero momentum). Only applies using the NAG method or a momentum method with
#' Nesterov momentum schedule. Also does nothing if \code{nest_convex_approx}
#' is \code{TRUE}.
#' @param nest_convex_approx If \code{TRUE}, then use an approximation due to
#' Sutskever for calculating the momentum parameter in the NAG method. Only
#' applies using the NAG method or a momentum method with Nesterov momentum
#' schedule.
#' @param nest_burn_in Number of iterations to wait before using a non-zero
#' momentum. Only applies using the NAG method or a momentum method with
#' Nesterov momentum schedule.
#' @param use_nest_mu_zero If \code{TRUE}, then the momentum on iteration zero
#' is set to 0.4. Otherwise, it's zero. Only applies using the NAG method or a
#' momentum method with Nesterov momentum schedule.
#' @param kappa Value by which to increase the step size for the \code{"bold"}
#' step size method or the \code{"DBD"} method.
#' @param kappa_fun Operator to use when combining the current step size with
#' \code{kappa}. Can be one of \code{"*"} (to multiply the current step size
#' with \code{kappa}) or \code{"+"} (to add).
#' @param phi Multiplier to reduce the step size by if using the \code{"DBD"}
#' method or the \code{"bold"} or \code{"back"} line search method. Should be
#' a positive value less than 1.
#' @param theta Weighting parameter used by the \code{"DBD"} method only, and
#' only if no momentum scheme is provided. Must be an integer between 0 and 1.
#' @param line_search Type of line search to use. See 'Details'.
#' @param c1 Sufficient decrease parameter for Wolfe-type line searches. Should
#' be a value between 0 and 1.
#' @param c2 Sufficient curvature parameter for line search for Wolfe-type line
#' searches. Should be a value between \code{c1} and 1.
#' @param step0 Initial value for the line search on the first step. See
#' 'Details'.
#' @param ls_initializer For Wolfe-type line searches only, how to initialize
#' the line search on iterations after the first. See 'Details'.
#' @param try_newton_step For Wolfe-type line searches only, try the
#' line step value of 1 as the initial step size whenever \code{ls_initializer}
#' suggests a step size > 1. Defaults to \code{TRUE} for quasi-Newton methods
#' such as BFGS and L-BFGS, \code{FALSE} otherwise.
#' @param mom_type Momentum type, either \code{"classical"} or
#' \code{"nesterov"}. See 'Details'.
#' @param mom_schedule Momentum schedule. See 'Details'.
#' @param mom_init Initial momentum value.
#' @param mom_final Final momentum value.
#' @param mom_switch_iter For \code{mom_schedule} \code{"switch"} only, the
#' iteration when \code{mom_init} is changed to \code{mom_final}.
#' @param mom_linear_weight If \code{TRUE}, the gradient contribution to the
#' update is weighted using momentum contribution.
#' @param restart Momentum restart type. Can be one of "fn" or "gr". See
#' 'Details'. Ignored if no momentum scheme is being used.
#' @param max_iter Maximum number of iterations to optimize for. Defaults to
#' 100.
#' @param max_fn Maximum number of function evaluations.
#' @param max_gr Maximum number of gradient evaluations.
#' @param max_fg Maximum number of function or gradient evaluations.
#' @param abs_tol Absolute tolerance for comparing two function evaluations.
#' @param rel_tol Relative tolerance for comparing two function evaluations.
#' @param grad_tol Absolute tolerance for the length (l2-norm) of the gradient
#' vector.
#' @param verbose If \code{TRUE}, log information about the progress of the
#' optimization to the console.
#' @param store_progress If \code{TRUE} store information about the progress
#' of the optimization in a data frame, and include it as part of the return
#' value.
#' @return A list with components:
#'\itemize{
#'  \item{\code{par}} Optimized parameters.
#'  \item{\code{nf}} Ttotal number of function evaluations carried out .
#'  \item{\code{ng}} Running total number of gradient evaluations carried out since
#'    iteration 1.
#'  \item{\code{f}} Value of the function, evaluated at the returned
#'    value of \code{par}.
#'  \item{\code{g2n}} The length (l2-norm) of the gradient vector, evaluated
#'    at the returned value of \code{par}.
#'  \item{\code{iter}} The number of iterations the optimization was carried
#'    out for.
#'  \item{\code{terminate}} List containing items: \code{what}, indicating what
#'  convergence criterion was met, and \code{val} specifying the value at
#'  convergence.
#'  \item{\code{progress}} Optional data frame containing information on the
#'  value of the function, gradient, momentum, and step sizes evaluated at each
#'  iteration. Only present if \code{store_progress} is set to \code{TRUE}.
#'}
#' @references
#' Jacobs, R. A. (1988).
#' Increased rates of convergence through learning rate adaptation.
#' \emph{Neural networks}, \emph{1}(4), 295-307.
#'
#' Janet, J. A., Scoggins, S. M., Schultz, S. M., Snyder, W. E., White, M. W.,
#' & Sutton, J. C. (1998, May).
#' Shocking: An approach to stabilize backprop training with greedy adaptive
#' learning rates.
#' In \emph{1998 IEEE International Joint Conference on Neural Networks Proceedings.}
#' (Vol. 3, pp. 2218-2223). IEEE.
#'
#' Nocedal, J., & Wright, S. (2006).
#' Numerical optimization.
#' Springer Science & Business Media.
#'
#' O'Donoghue, B., & Candes, E. (2013).
#' Adaptive restart for accelerated gradient schemes.
#' \emph{Foundations of computational mathematics}, \emph{15}(3), 715-732.
#'
#' Sutskever, I., Martens, J., Dahl, G., & Hinton, G. (2013).
#' On the importance of initialization and momentum in deep learning.
#' In \emph{Proceedings of the 30th international conference on machine learning (ICML-13)}
#' (pp. 1139-1147).
#' @examples
#' # Function to optimize and starting point defined after creating optimizer
#' rosenbrock_fg <- list(
#'   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
#'   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                          200 *        (x[2] - x[1] * x[1])) })
#' rb0 <- c(-1.2, 1)
#'
#' # Minimize using L-BFGS
#' res <- mizer(rb0, rosenbrock_fg, method = "L-BFGS")
#'
#' # Conjugate gradient with Fletcher-Reeves update, tight Wolfe line search
#' res <- mizer(rb0, rosenbrock_fg, method = "CG", cg_update = "FR", c2 = 0.1)
#'
#' # Steepest decent with constant momentum = 0.9
#' res <- mizer(rb0, rosenbrock_fg, method = "SD", mom_type = "classical",
#'              mom_schedule = 0.9)
#'
#' # Steepest descent with constant momentum in the Nesterov style as described
#' # by Sutskever and co-workers
#' res <- mizer(rb0, rosenbrock_fg, method = "SD", mom_type = "nesterov",
#'              mom_schedule = 0.9)
#'
#' # Nesterov momentum with adaptive restart comparing function values
#' res <- mizer(rb0, rosenbrock_fg, method = "SD", mom_type = "nesterov",
#'              mom_schedule = 0.9, restart = "fn")
#' @export
mizer <- function(par, fg,
                  method = "L-BFGS",
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
                  c2 = NULL,
                  step0 = NULL,
                  ls_initializer = NULL,
                  try_newton_step = NULL,
                  # Momentum
                  mom_type = NULL,
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
                    try_newton_step = try_newton_step,
                    mom_type = mom_type,
                    mom_schedule = mom_schedule,
                    mom_init = mom_init,
                    mom_final = mom_final,
                    mom_switch_iter = mom_switch_iter,
                    mom_linear_weight = mom_linear_weight,
                    max_iter = max_iter,
                    restart = restart)

  res <- opt_loop(opt, par, fg,
          max_iter = max_iter,
          max_fn = max_fn, max_gr = max_gr, max_fg = max_fg,
          abs_tol = abs_tol, rel_tol = rel_tol, grad_tol = grad_tol,
          store_progress = store_progress,
          verbose = verbose)

  Filter(Negate(is.null),
         res[c("f", "g2n", "nf", "ng", "par", "iter", "terminate")])
}

#' Create an Optimizer
#'
#' Factory function for creating a (possibly uninitialized) optimizer.
#'
#' If the function to be optimized and starting point are not present at
#' creation time, then the optimizer should be initialized using
#' \code{\link{mizer_init}} before being used with \code{\link{mizer_step}}.
#'
#' See the documentation to \code{\link{mizer}} for an explanation of all
#' the parameters.
#'
#' Details of the \code{fg} list containing the function to be optimized and
#' its gradient can be found in the 'Details' section of \code{\link{mizer}}.
#' It is optional for this function, but if it is passed to this function,
#' along with the vector of initial values, \code{par}, the optimizer will be
#' returned already initialized for this function. Otherwise,
#' \code{\link{mizer_init}} must be called before optimization begins.
#'
#' @param method Optimization method. See 'Details'.
#' @param norm_direction If \code{TRUE}, then the steepest descent direction
#' is normalized to unit length. Useful for adaptive step size methods where
#' the previous step size is used to initialize the next iteration.
#' @param scale_hess if \code{TRUE}, the approximation to the inverse Hessian
#' is scaled according to the method described by Nocedal and Wright
#' (approximating an eigenvalue). Applies only to the methods \code{BFGS}
#' (where the scaling is applied only during the first step) and \code{L-BFGS}
#' (where the scaling is applied during every iteration). Ignored otherwise.
#' @param memory The number of updates to store if using the \code{L-BFGS}
#' method. Ignored otherwise. Must be a positive integer.
#' @param cg_update Type of update to use for the \code{CG} method. Can be
#' one of \code{"FR"} (Fletcher-Reeves), \code{"PR"} (Polak-Ribiere),
#' \code{"PR+"} (Polak-Ribiere with a reset to steepest descent), \code{"HS"}
#' (Hestenes-Stiefel), or \code{"DY"} (Dai-Yuan). Ignored if \code{method} is
#' not \code{"CG"}.
#' @param nest_q Strong convexity parameter for the \code{"NAG"} method's
#' momentum term. Must take a value between 0 (strongly convex) and 1 (results
#' in steepest descent).Ignored unless the \code{method} is \code{"NAG"} and
#' \code{nest_convex_approx} is \code{FALSE}.
#' @param nest_convex_approx If \code{TRUE}, then use an approximation due to
#' Sutskever for calculating the momentum parameter in the NAG method. Only
#' applies if \code{method} is \code{"NAG"}.
#' @param nest_burn_in Number of iterations to wait before using a non-zero
#' momentum. Only applies if using the \code{"NAG"} method or setting the
#' \code{momentum_type} to "Nesterov".
#' @param use_nest_mu_zero If \code{TRUE}, then the momentum on iteration zero
#' is set to 0.4. Otherwise, it's zero. Applies only if
#' \code{nest_convex_approx} is \code{TRUE}.
#' @param kappa Value by which to increase the step size for the \code{"bold"}
#' step size method or the \code{"DBD"} method.
#' @param kappa_fun Operator to use when combining the current step size with
#' \code{kappa}. Can be one of \code{"*"} (to multiply the current step size
#' with \code{kappa}) or \code{"+"} (to add).
#' @param phi Multiplier to reduce the step size by if using the \code{"DBD"}
#' method or the \code{"bold"} or \code{"back"} line search method. Should be
#' a positive value less than 1.
#' @param theta Weighting parameter used by the \code{"DBD"} method only, and
#' only if no momentum scheme is provided. Must be an integer between 0 and 1.
#' @param line_search Type of line search to use. See 'Details'.
#' @param c1 Sufficient decrease parameter for Wolfe-type line searches. Should
#' be a value between 0 and 1.
#' @param c2 Sufficient curvature parameter for line search for Wolfe-type line
#' searches. Should be a value between \code{c1} and 1.
#' @param step0 Initial value for the line search on the first step. See
#' 'Details'.
#' @param ls_initializer For Wolfe-type line searches only, how to initialize
#' the line search on iterations after the first. See 'Details'.
#' @param try_newton_step For Wolfe-type line searches only, try the
#' line step value of 1 as the initial step size whenever \code{ls_initializer}
#' suggests a step size > 1. Defaults to \code{TRUE} for quasi-Newton methods
#' such as BFGS and L-BFGS, \code{FALSE} otherwise.
#' @param mom_type Momentum type, either \code{"classical"} or
#' \code{"nesterov"}.
#' @param mom_schedule Momentum schedule. See 'Details'.
#' @param mom_init Initial momentum value.
#' @param mom_final Final momentum value.
#' @param mom_switch_iter For \code{mom_schedule} \code{"switch"} only, the
#' iteration when \code{mom_init} is changed to \code{mom_final}.
#' @param mom_linear_weight If \code{TRUE}, the gradient contribution to the
#' update is weighted using momentum contribution.
#' @param max_iter Maximum number of iterations the optimization will be carried
#' out over. Used only if \code{mom_schedule} is set to \code{"ramp"}.
#' @param restart Momentum restart type. Can be one of "fn" or "gr". See
#' 'Details'.
#' @param par Initial values for the function to be optimized over. Optional.
#' @param fg Function and gradient list. See 'Details'. Optional.
#' @export
#' @examples
#' # Function to optimize and starting point
#' rosenbrock_fg <- list(
#'   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
#'   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                          200 *        (x[2] - x[1] * x[1])) })
#' rb0 <- c(-1.2, 1)
#'
#' # Create an optimizer and initialize it for use with the Rosenbrock function
#' opt <- make_mizer(method = "l-bfgs", par = rb0, fg = rosenbrock_fg)
#'
#' # Create optimizer without initialization
#' opt <- make_mizer(method = "l-bfgs")
#'
#' # Need to call mizer_init separately:
#' opt <- mizer_init(opt, rb0, rosenbrock_fg)
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
                       c1 = 1e-4, c2 = NULL,
                       step0 = NULL,
                       ls_initializer = NULL,
                       try_newton_step = NULL,
                       # Momentum
                       mom_type = NULL,
                       mom_schedule = NULL,
                       mom_init = NULL,
                       mom_final = NULL,
                       mom_switch_iter = NULL,
                       mom_linear_weight = FALSE,
                       max_iter = NULL,
                       restart = NULL,
                       par = NULL,
                       fg = NULL) {
  dir_type <- NULL
  method <- toupper(method)
  if (method == "SD") {
    dir_type <- sd_direction(normalize = norm_direction)
  }
  else if (method == "NEWTON") {
    dir_type <- newton_direction()
    if (is.null(try_newton_step)) {
      try_newton_step <- TRUE
    }
  }
  else if (method == "PHESS") {
    dir_type <- partial_hessian_direction()
    if (is.null(try_newton_step)) {
      try_newton_step <- TRUE
    }
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
    if (is.null(try_newton_step)) {
      try_newton_step <- TRUE
    }
  }
  else if (method == "L-BFGS") {
    dir_type <- lbfgs_direction(memory = memory, scale_inverse = scale_hess)
    if (is.null(try_newton_step)) {
      try_newton_step <- TRUE
    }
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

  # If it's not already been turned on, turn off the Newton step option
  if (is.null(try_newton_step)) {
    try_newton_step <- FALSE
  }

  step_type <- NULL
  line_search <- toupper(line_search)
  if (method == "DBD") {
    if (is.character(step0) || is.numeric(step0)) {
      eps_init <- step0
    }
    else {
      eps_init <- "r"
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
    if (method %in% c("NEWTON", "PHESS", "BFGS", "L-BFGS")) {
      if (is.null(c2)) {
        c2 <- 0.9
      }
      if (is.null(step0)) {
        step0 <- 1
      }
      if (is.null(ls_initializer)) {
        ls_initializer <- "q"
      }
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    }
    else {
      if (is.null(c2)) {
        c2 <- 0.1
      }
      if (is.null(step0)) {
        step0 <- "r"
      }
      if (is.null(ls_initializer)) {
        ls_initializer <- "s"
      }
      if (is.null(try_newton_step)) {
        try_newton_step <- FALSE
      }
    }

    if (line_search == "MT") {
      step_type <- more_thuente_ls(c1 = c1, c2 = c2,
                                   initializer = tolower(ls_initializer),
                                   initial_step_length = step0,
                                   try_newton_step = try_newton_step)
    }
    else if (line_search == "RAS") {
      step_type <- rasmussen_ls(c1 = c1, c2 = c2,
                                initializer = tolower(ls_initializer),
                                initial_step_length = step0,
                                try_newton_step = try_newton_step)
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
        step_size = step_type)))

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

  if (method == "MOM") {
    if (is.null(mom_type)) {
      mom_type <- "classical"
    }
    if (is.null(mom_schedule)) {
      mom_schedule <- 0.9
    }
    if (is.null(restart)) {
      restart <- "fn"
    }
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
#' Details of the \code{fg} list can be found in the 'Details' section of
#' \code{\link{mizer}}.
#'
#' @param opt Optimizer, created by \code{\link{make_mizer}}.
#' @param par Vector of initial values for the function to be optimized over.
#' @param fg Function and gradient list. See the documentaion of
#' \code{\link{mizer}}.
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
#' @examples
#' rosenbrock_fg <- list(
#'   fn = function(x) {
#'     100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2
#'   },
#'   gr = function(x) {
#'     c(
#'      -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'       200 *        (x[2] - x[1] * x[1]))
#'  })
#'  rb0 <- c(-1.2, 1)
#'
#'  opt <- make_mizer(method = "sd", line_search = "const", step0 = 0.0001,
#'                    par = rb0, fg = rosenbrock_fg)
#'  par <- rb0
#'  for (iter in 1:3) {
#'    res <- mizer_step(opt, par, rosenbrock_fg, iter)
#'    par <- res$par
#'    opt <- res$opt
#'  }
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
#' Should be called after creating an optimizer with \code{\link{make_mizer}}
#' and before beginning any optimization with \code{\link{mizer_step}}. Note
#' that if \code{fg} and \code{par} are available at the time
#' \code{\link{mizer_step}} is called, they can be passed to that function
#' and initialization will be carried out automatically, avoiding the need to
#' call \code{mizer_init}.
#'
#' Details of the \code{fg} list can be found in the 'Details' section of
#' \code{\link{mizer}}.
#'
#' @param opt Optimizer, created by \code{\link{make_mizer}}.
#' @param par Vector of initial values for the function to be optimized over.
#' @param fg Function and gradient list. See the documentaion of
#' \code{\link{mizer}}.
#' @return Initialized optimizer.
#' @export
#' @examples
#'
#' # Create an optimizer
#' opt <- make_mizer(method = "l-bfgs")
#'
#' # Function to optimize and starting point defined after creating optimizer
#' rosenbrock_fg <- list(
#'   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
#'   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                          200 *        (x[2] - x[1] * x[1])) })
#' rb0 <- c(-1.2, 1)
#'
#' # Initialize with function and starting point before commencing optimization
#' opt <- mizer_init(opt, rb0, rosebrock_fg)
#'
#' # Finally, can commence the optimization loop
#' par <- rb0
#' for (iter in 1:3) {
#'   res <- mizer_step(opt, par, rosenbrock_fg, iter)
#'   par <- res$par
#'   opt <- res$opt
#' }
#'
mizer_init <- function(opt, par, fg) {
  opt <- register_hooks(opt)
  opt <- life_cycle_hook("opt", "init", opt, par, fg, 0)
  opt
}
