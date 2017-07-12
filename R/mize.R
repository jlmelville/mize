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
#' gradient in the same routine. Takes a vector of parameters and returns a list
#' containing the function result as \code{fn} and the gradient result as
#' \code{gr}.
#' }
#'
#' The \code{fg} function is optional, but for some methods (e.g. line search
#' methods based on the Wolfe criteria), both the function and gradient values
#' are needed for the same parameter value. Calculating them in the same
#' function can save time if there is a lot of shared work.
#'
#' @section Optimization Methods:
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
#'     \item \code{"HZ"} The method of Hager and Zhang.
#'     \item \code{"HZ+"} The method of Hager and Zhang with restart, as used
#'     in CG_DESCENT.
#'   }
#' The \code{"PR+"} and \code{"HZ+"} are likely to be most robust in practice.
#' Other updates are available more for curiosity purposes.
#' \item \code{"NAG"} is the Nesterov Accelerated Gradient method. The exact
#' form of the momentum update in this method can be controlled with the
#' following parameters:
#'   \itemize{
#'   \item{\code{nest_q}} Strong convexity parameter. Must take a value
#'   between 0 (strongly convex) and 1 (zero momentum). Ignored if
#'   \code{nest_convex_approx} is \code{TRUE}.
#'   \item{\code{nest_convex_approx}} If \code{TRUE}, then use an approximation
#'   due to Sutskever for calculating the momentum parameter.
#'   \item{\code{nest_burn_in}} Number of iterations to wait before using a
#'   non-zero momentum.
#'   }
#' \item \code{"DBD"} is the Delta-Bar-Delta method of Jacobs.
#' \item \code{"Momentum"} is steepest descent with momentum. See below for
#' momentum options.
#' }
#'
#' For more details on gradient-based optimization in general, and the BFGS,
#' L-BFGS and CG methods, see Nocedal and Wright.
#'
#' @section Line Search:
#' The parameter \code{line_search} determines the line search to be carried
#' out:
#'
#' \itemize{
#'   \item \code{"Rasmussen"} carries out a line search using the strong Wolfe
#'   conditions as implemented by Carl Edward Rasmussen's minimize.m routines.
#'   \item \code{"More-Thuente"} carries out a line search using the strong Wolfe
#'   conditions and the method of More-Thuente. Can be abbreviated to
#'   \code{"MT"}.
#'   \item \code{"Schmidt"} carries out a line search using the strong Wolfe
#'   conditions as implemented in Mark Schmidt's minFunc routines.
#'   \item \code{"Backtracking"} carries out a back tracking line search using
#'   the sufficient decrease (Armijo) condition. By default, cubic interpolation
#'   is used to find an acceptable step size. A constant step size reduction
#'   can be used by specifying a value for \code{step_down} between 0 and 1
#'   (e.g. step size will be halved if \code{step_down} is set to \code{0.5}).
#'   \item \code{"Bold Driver"} carries out a back tracking line search until a
#'   reduction in the function value is found.
#'   \item \code{"Constant"} uses a constant line search, the value of which
#'   should be provided with \code{step0}. Note that this value will be
#'   multiplied by the magnitude of the direction vector used in the gradient
#'   descent method. For method \code{"SD"} only, setting the
#'   \code{norm_direction} parameter to \code{TRUE} will scale the direction
#'   vector so it has unit length.
#' }
#'
#' If using one of the methods: \code{"BFGS"}, \code{"L-BFGS"}, \code{"CG"} or
#' \code{"NAG"}, one of the Wolfe line searches: \code{"Rasmussen"} or
#' \code{"More-Thuente"}, \code{"Schmidt"} or \code{"Hager-Zhang"} should be
#' used, otherwise very poor performance is likely to be encountered. The
#' following parameters can be used to control the line search:
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
#'    as-is. Otherwise, by passing a string as an argument, a guess is made
#'    based on values of the gradient, function or parameters, at the starting
#'    point:
#'    \itemize{
#'      \item{\code{"rasmussen"}} As used by Rasmussen in \code{minimize.m}:
#'      \deqn{\frac{1}{1+\left|g\right|^2}}{1 / 1 + (|g|^2)}
#'      \item{\code{"scipy"}} As used in scipy's \code{optimize.py}
#'      \deqn{\frac{1}{\left|g\right|}}{1 / |g|}
#'      \item{\code{"schmidt"}} As used by Schmidt in \code{minFunc.m}
#'      (the reciprocal of the l1 norm of g)
#'      \deqn{\frac{1}{\left|g\right|_1}}{1 / |g|1}
#'      \item{\code{"hz"}} The method suggested by Hager and Zhang (2006) for
#'      the CG_DESCENT software.
#'    }
#'    These arguments can be abbreviated.
#'    \item{\code{step_next_init}} How to initialize subsequent line searches
#'    after the first, using results from the previous line search:
#'    \itemize{
#'      \item{\code{"slope ratio"}} Slope ratio method.
#'      \item{\code{"quadratic"}} Quadratic interpolation method.
#'      \item{\code{"hz"}} The QuadStep method of Hager and Zhang (2006) for
#'      the CG_DESCENT software.
#'    }
#'    These arguments can be abbreviated. Details on the first two methods
#'    are provided by Nocedal and Wright.
#'    \item{\code{try_newton_step}} For quasi-Newton methods (\code{"BFGS"} and
#'    \code{"L-BFGS"}), setting this to \code{TRUE} will try the "natural" step
#'    size of 1, whenever the \code{step_next_init} method suggests an initial
#'    step size larger than that. On by default for BFGS and L-BFGS, off for
#'    everything else.
#'    \item{\code{strong_curvature}} If \code{TRUE}, then the strong curvature
#'    condition will be used to check termination in Wolfe line search methods.
#'    If \code{FALSE}, then the standard curvature condition will be used. The
#'    default is \code{NULL} which lets the different Wolfe line searches choose
#'    whichever is their default behavior. This option is ignored if not using
#'    a Wolfe line search method.
#'    \item{\code{approx_armijo}} If \code{TRUE}, then the approximate Armijo
#'    sufficient decrease condition (Hager and Zhang, 2005) will be used to
#'    check termination in Wolfe line search methods. If \code{FALSE}, then the
#'    exact curvature condition will be used. The default is \code{NULL} which
#'    lets the different Wolfe line searches choose whichever is their default
#'    behavior. This option is ignored if not using a Wolfe line search method.
#' }
#'
#' For the Wolfe line searches, the methods of \code{"Rasmussen"},
#' \code{"Schmidt"} and \code{"More-Thuente"} default to using the strong
#' curvature condition and the exact Armijo condition to terminate the line
#' search (i.e. Strong Wolfe conditions). The default step size initialization
#' methods use the Rasmussen method for the first iteration and quadratic
#' interpolation for subsequent iterations.
#'
#' The \code{"Hager-Zhang"} Wolfe line search method defaults to the standard
#' curvature condition and the approximate Armijo condition (i.e. approximate
#' Wolfe conditions). The default step size initialization methods are those
#' used by Hager and Zhang (2006) in the description of CG_DESCENT.
#'
#' If the \code{"DBD"} is used for the optimization \code{"method"}, then the
#' \code{line_search} parameter is ignored, because this method controls both
#' the direction of the search and the step size simultaneously. The following
#' parameters can be used to control the step size:
#'
#' \itemize{
#'   \item{\code{step_up}} The amount by which to increase the step size in a
#'   direction where the current step size is deemed to be too short. This
#'   should be a positive scalar.
#'   \item{\code{step_down}} The amount by which to decrease the step size in a
#'   direction where the currents step size is deemed to be too long. This
#'   should be a positive scalar smaller than 1. Default is 0.5.
#'   \item{\code{step_up_fun}} How to increase the step size: either the method of
#'   Jacobs (addition of \code{step_up}) or Janet and co-workers (multiplication
#'   by \code{step_up}). Note that the step size decrease \code{step_down} is always
#'   a multiplication.
#' }
#'
#' The \code{"bold driver"} line search also uses the \code{step_up} and
#' \code{step_down} parameters with similar meanings to their use with the
#' \code{"DBD"} method: the backtracking portion reduces the step size by a
#' factor of \code{step_down}. Once a satisfactory step size has been found, the
#' line search for the next iteration is initialized by multiplying the
#' previously found step size by \code{step_up}.
#'
#' @section Momentum:
#' For \code{method} \code{"Momentum"}, momentum schemes can be accessed
#' through the momentum arguments:
#'
#' \itemize{
#' \item{\code{mom_type}} Momentum type, either \code{"classical"} or
#'   \code{"nesterov"} (case insensitive, can be abbreviated). Using
#'   \code{"nesterov"} applies the momentum step before the
#'   gradient descent as suggested by Sutskever, emulating the behavior of the
#'   Nesterov Accelerated Gradient method.
#' \item{\code{mom_schedule}} How the momentum changes over the course of the
#'   optimization:
#'   \itemize{
#'   \item{If a numerical scalar is provided, a constant momentum will be
#'     applied throughout.}
#'   \item{\code{"nsconvex"}} Use the momentum schedule from the Nesterov
#'   Accelerated Gradient method suggested for non-strongly convex functions.
#'   Parameters which control the NAG momentum
#'   can also be used in combination with this option.
#'   \item{\code{"switch"}} Switch from one momentum value (specified via
#'   \code{mom_init}) to another (\code{mom_final}) at a
#'   a specified iteration (\code{mom_switch_iter}).
#'   \item{\code{"ramp"}} Linearly increase from one momentum value
#'   (\code{mom_init}) to another (\code{mom_final}).
#'   \item{If a function is provided, this will be invoked to provide a momentum
#'   value. It must take one argument (the current iteration number) and return
#'   a scalar.}
#'   }
#'   String arguments are case insensitive and can be abbreviated.
#' }
#'
#' The \code{restart} parameter provides a way to restart the momentum if the
#' optimization appears to be not be making progress, inspired by the method
#' of O'Donoghue and Candes (2013) and Su and co-workers (2014). There are three
#' strategies:
#' \itemize{
#'   \item{\code{"fn"}} A restart is applied if the function does not decrease
#'   on consecutive iterations.
#'   \item{\code{"gr"}} A restart is applied if the direction of the
#'   optimization is not a descent direction.
#'   \item{\code{"speed"}} A restart is applied if the update vector is not
#'   longer (as measured by Euclidean 2-norm) in consecutive iterations.
#' }
#'
#' The effect of the restart is to "forget" any previous momentum update vector,
#' and, for those momentum schemes that change with iteration number, to
#' effectively reset the iteration number back to zero. If the \code{mom_type}
#' is \code{"nesterov"}, the gradient-based restart is not available. The
#' \code{restart_wait} parameter controls how many iterations to wait after a
#' restart, before allowing another restart. Must be a positive integer. Default
#' is 10, as used by Su and co-workers (2014). Setting this too low could
#' cause premature convergence. These methods were developed specifically
#' for the NAG method, but can be employed with any momentum type and schedule.
#'
#' If \code{method} type \code{"momentum"} is specified with no other values,
#' the momentum scheme will default to a constant value of \code{0.9}.
#'
#' @section Convergence:
#'
#' There are several ways for the optimization to terminate. The type of
#' termination is communicated by a two-item list \code{terminate} in the return
#' value, consisting of \code{what}, a short string describing what caused the
#' termination, and \code{val}, the value of the termination criterion that
#' caused termination.
#'
#' The following parameters control various stopping criteria:
#'
#' \itemize{
#'   \item{\code{max_iter}} Maximum number of iterations to calculate. Reaching
#'   this limit is indicated by \code{terminate$what} being \code{"max_iter"}.
#'   \item{\code{max_fn}} Maximum number of function evaluations allowed.
#'   Indicated by \code{terminate$what} being \code{"max_fn"}.
#'   \item{\code{max_gr}} Maximum number of gradient evaluations allowed.
#'   Indicated by \code{terminate$what} being \code{"max_gr"}.
#'   \item{\code{max_fg}} Maximum number of gradient evaluations allowed.
#'   Indicated by \code{terminate$what} being \code{"max_fg"}.
#'   \item{\code{abs_tol}} Absolute tolerance of the function value. If the
#'   absolute value of the function falls below this threshold,
#'   \code{terminate$what} will be \code{"abs_tol"}. Will only be triggered if
#'   the objective function has a minimum value of zero.
#'   \item{\code{rel_tol}} Relative tolerance of the function value, comparing
#'   consecutive function evaluation results. Indicated by \code{terminate$what}
#'   being \code{"rel_tol"}.
#'   \item{\code{grad_tol}} Absolute tolerance of the l2 (Euclidean) norm of
#'   the gradient. Indicated by \code{terminate$what} being \code{"grad_tol"}.
#'   Note that the gradient norm is not a very reliable stopping criterion
#'   (see Nocedal and co-workers 2002), but is quite commonly used, so this
#'   might be useful for comparison with results from other optimizers.
#'   \item{\code{ginf_tol}} Absolute tolerance of the infinity norm (maximum
#'   absolute component) of the gradient. Indicated by \code{terminate$what}
#'   being \code{"ginf_tol"}.
#'   \item{\code{step_tol}} Absolute tolerance of the step size, i.e. the
#'   Euclidean distance between values of \code{par} fell below the specified
#'   value. Indicated by \code{terminate$what} being \code{"step_tol"}.
#'   For those optimization methods which allow for abandoning the result of an
#'   iteration and restarting using the previous iteration's value of
#'   \code{par} an iteration, \code{step_tol} will not be triggered.
#' }
#'
#' Convergence is checked between specific iterations. How often is determined
#' by the \code{check_conv_every} parameter, which specifies the number of
#' iterations between each check. By default, this is set for every iteration.
#'
#' Be aware that if \code{abs_tol} or \code{rel_tol} are non-\code{NULL}, this
#' requires the function to have been evaluated at the current position at the
#' end of each iteration. If the function at that position has not been
#' calculated, it will be calculated and will contribute to the total reported
#' in the \code{counts} list in the return value. The calculated function value
#' is cached for use by the optimizer in the next iteration, so if the optimizer
#' would have needed to calculate the function anyway (e.g. use of the strong
#' Wolfe line search methods), there is no significant cost accrued by
#' calculating it earlier for convergence calculations. However, for methods
#' that don't use the function value at that location, this could represent a
#' lot of extra function evaluations. On the other hand, not checking
#' convergence could result in a lot of extra unnecessary iterations.
#' Similarly, if \code{grad_tol} or \code{ginf_tol} is non-\code{NULL}, then
#' the gradient will be calculated if needed.
#'
#' If extra function or gradient evaluations is an issue, set
#' \code{check_conv_every} to a higher value, but be aware that this can cause
#' convergence limits to be exceeded by a greater amount.
#'
#' Note also that if the \code{verbose} parameter is \code{TRUE}, then a summary
#' of the results so far will be logged to the console whenever a convergence
#' check is carried out. If the \code{store_progress} parameter is \code{TRUE},
#' then the same information will be returned as a data frame in the return
#' value. For a long optimization this could be a lot of data, so by default it
#' is not stored.
#'
#' Other ways for the optimization to terminate is if an iteration generates a
#' non-finite (i.e. \code{Inf} or \code{NaN}) gradient or function value.
#' Some, but not all, line-searches will try to recover from the latter, by
#' reducing the step size, but a non-finite gradient calculation during the
#' gradient descent portion of optimization is considered catastrophic by mize,
#' and it will give up. Termination under non-finite gradient or function
#' conditions will result in \code{terminate$what} being \code{"gr_inf"} or
#' \code{"fn_inf"} respectively. Unlike the convergence criteria, the
#' optimization will detect these error conditions and terminate even if a
#' convergence check would not be carried out for this iteration.
#'
#' The value of \code{par} in the return value should be the parameters which
#' correspond to the lowest value of the function that has been calculated
#' during the optimization. As discussed above however, determining which set
#' of parameters requires a function evaluation at the end of each iteration,
#' which only happens if either the optimization method calculates it as part
#' of its own operation or if a convergence check is being carried out during
#' this iteration. Therefore, if your method does not carry out function
#' evaluations and \code{check_conv_every} is set to be so large that no
#' convergence calculation is carried out before \code{max_iter} is reached,
#' then the returned value of \code{par} is the last value encountered.
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
#' @param step_up Value by which to increase the step size for the \code{"bold"}
#' step size method or the \code{"DBD"} method.
#' @param step_up_fun Operator to use when combining the current step size with
#' \code{step_up}. Can be one of \code{"*"} (to multiply the current step size
#' with \code{step_up}) or \code{"+"} (to add).
#' @param step_down Multiplier to reduce the step size by if using the \code{"DBD"}
#' method or the \code{"bold"} line search method. Should be a positive value
#' less than 1. Also optional for use with the \code{"back"} line search method.
#' @param dbd_weight Weighting parameter used by the \code{"DBD"} method only, and
#' only if no momentum scheme is provided. Must be an integer between 0 and 1.
#' @param line_search Type of line search to use. See 'Details'.
#' @param c1 Sufficient decrease parameter for Wolfe-type line searches. Should
#' be a value between 0 and 1.
#' @param c2 Sufficient curvature parameter for line search for Wolfe-type line
#' searches. Should be a value between \code{c1} and 1.
#' @param step0 Initial value for the line search on the first step. See
#' 'Details'.
#' @param step_next_init For Wolfe-type line searches only, how to initialize
#' the line search on iterations after the first. See 'Details'.
#' @param try_newton_step For Wolfe-type line searches only, try the
#' line step value of 1 as the initial step size whenever \code{step_next_init}
#' suggests a step size > 1. Defaults to \code{TRUE} for quasi-Newton methods
#' such as BFGS and L-BFGS, \code{FALSE} otherwise.
#' @param ls_max_fn Maximum number of function evaluations allowed during a
#' line search.
#' @param ls_max_gr Maximum number of gradient evaluations allowed during a
#' line search.
#' @param ls_max_fg Maximum number of function or gradient evaluations allowed
#' during a line search.
#' @param ls_max_alpha_mult Maximum multiplier for alpha between iterations.
#' Only applies for Wolfe-type line searches and if \code{step_next_init} is
#' set to \code{"slope"}
#' @param strong_curvature (Optional). If \code{TRUE} use the strong
#' curvature condition in Wolfe line search. See the 'Line Search' section
#' for details.
#' @param approx_armijo (Optional). If \code{TRUE} use the approximate Armijo
#' condition in Wolfe line search. See the 'Line Search' section for details.
#' @param mom_type Momentum type, either \code{"classical"} or
#' \code{"nesterov"}. See 'Details'.
#' @param mom_schedule Momentum schedule. See 'Details'.
#' @param mom_init Initial momentum value.
#' @param mom_final Final momentum value.
#' @param mom_switch_iter For \code{mom_schedule} \code{"switch"} only, the
#' iteration when \code{mom_init} is changed to \code{mom_final}.
#' @param use_init_mom If \code{TRUE}, then the momentum coefficient on
#' the first iteration is non-zero. Otherwise, it's zero. Only applies if
#' using a momentum schedule.
#' @param mom_linear_weight If \code{TRUE}, the gradient contribution to the
#' update is weighted using momentum contribution.
#' @param restart Momentum restart type. Can be one of "fn", "gr" or "speed".
#' See Details'. Ignored if no momentum scheme is being used.
#' @param restart_wait Number of iterations to wait between restarts. Ignored
#' if \code{restart} is \code{NULL}.
#' @param max_iter Maximum number of iterations to optimize for. Defaults to
#' 100. See the 'Convergence' section for details.
#' @param max_fn Maximum number of function evaluations. See the 'Convergence'
#' section for details.
#' @param max_gr Maximum number of gradient evaluations. See the 'Convergence'
#' section for details.
#' @param max_fg Maximum number of function or gradient evaluations. See the
#' 'Convergence' section for details.
#' @param abs_tol Absolute tolerance for comparing two function evaluations.
#' See the 'Convergence' section for details.
#' @param rel_tol Relative tolerance for comparing two function evaluations.
#' See the 'Convergence' section for details.
#' @param grad_tol Absolute tolerance for the length (l2-norm) of the gradient
#' vector. See the 'Convergence' section for details.
#' @param ginf_tol Absolute tolerance for the infinity norm (maximum absolute
#' component) of the gradient vector. See the 'Convergence' section for details.
#' @param step_tol Absolute tolerance for the size of the parameter update.
#' See the 'Convergence' section for details.
#' @param check_conv_every Positive integer indicating how often to check
#' convergence. Default is 1, i.e. every iteration. See the 'Convergence'
#' section for details.
#' @param log_every Positive integer indicating how often to log convergence
#' results to the console. Ignored if \code{verbose} is \code{FALSE}.
#' If not an integer multiple of \code{check_conv_every}, it will be set to
#' \code{check_conv_every}.
#' @param verbose If \code{TRUE}, log information about the progress of the
#' optimization to the console.
#' @param store_progress If \code{TRUE} store information about the progress
#' of the optimization in a data frame, and include it as part of the return
#' value.
#' @return A list with components:
#'\itemize{
#'  \item{\code{par}} Optimized parameters. Normally, this is the best set of
#'  parameters seen during optimization, i.e. the set that produced the minimum
#'  function value. This requires that convergence checking with is carried out,
#'  including function evaluation where necessary. See the 'Convergence'
#'  section for details.
#'  \item{\code{nf}} Total number of function evaluations carried out. This
#'  includes any extra evaluations required for convergence calculations. Also,
#'  a function evaluation may be required to calculate the value of \code{f}
#'  returned in this list (see below). Additionally, if the \code{verbose}
#'  parameter is \code{TRUE}, then function and gradient information for the
#'  initial value of \code{par} will be logged to the console. These values
#'  are cached for subsequent use by the optimizer.
#'  \item{\code{ng}} Total number of gradient evaluations carried out. This
#'  includes any extra evaluations required for convergence calculations using
#'  \code{grad_tol}. As with \code{nf}, additional gradient calculations beyond
#'  what you're expecting may have been needed for logging, convergence and
#'  calculating the value of \code{g2} or \code{ginf} (see below).
#'  \item{\code{f}} Value of the function, evaluated at the returned
#'  value of \code{par}.
#'  \item{\code{g2}} Optional: the length (Euclidean or l2-norm) of the
#'  gradient vector, evaluated at the returned value of \code{par}. Calculated
#'  only if \code{grad_tol} is non-null.
#'  \item{\code{ginf}} Optional: the infinity norm (maximum absolute component)
#'  of the gradient vector, evaluated at the returned value of \code{par}.
#'  Calculated only if \code{ginf_tol} is non-null.
#'  \item{\code{iter}} The number of iterations the optimization was carried
#'  out for.
#'  \item{\code{terminate}} List containing items: \code{what}, indicating what
#'  convergence criterion was met, and \code{val} specifying the value at
#'  convergence. See the 'Convergence' section for more details.
#'  \item{\code{progress}} Optional data frame containing information on the
#'  value of the function, gradient, momentum, and step sizes evaluated at each
#'  iteration where convergence is checked. Only present if
#'  \code{store_progress} is set to \code{TRUE}. Could get quite large if the
#'  optimization is long and the convergence is checked regularly.
#'}
#' @references
#'
#' Hager, W. W., & Zhang, H. (2005).
#' A new conjugate gradient method with guaranteed descent and an efficient line search.
#' \emph{SIAM Journal on Optimization}, \emph{16}(1), 170-192.
#'
#' Hager, W. W., & Zhang, H. (2006).
#' Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{32}(1), 113-137.
#'
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
#' More', J. J., & Thuente, D. J. (1994).
#' Line search algorithms with guaranteed sufficient decrease.
#' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{20}(3), 286-307.
#'
#' Nocedal, J., Sartenaer, A., & Zhu, C. (2002).
#' On the behavior of the gradient norm in the steepest descent method.
#' \emph{Computational Optimization and Applications}, \emph{22}(1), 5-35.
#'
#' Nocedal, J., & Wright, S. (2006).
#' Numerical optimization.
#' Springer Science & Business Media.
#'
#' O'Donoghue, B., & Candes, E. (2013).
#' Adaptive restart for accelerated gradient schemes.
#' \emph{Foundations of computational mathematics}, \emph{15}(3), 715-732.
#'
#' Schmidt, M. (2005).
#' minFunc: unconstrained differentiable multivariate optimization in Matlab.
#' \url{http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html}
#'
#' Su, W., Boyd, S., & Candes, E. (2014).
#' A differential equation for modeling Nesterov's accelerated gradient method: theory and insights.
#' In \emph{Advances in Neural Information Processing Systems} (pp. 2510-2518).
#'
#' Sutskever, I. (2013).
#' \emph{Training recurrent neural networks}
#' (Doctoral dissertation, University of Toronto).
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
#' res <- mize(rb0, rosenbrock_fg, method = "L-BFGS")
#'
#' # Conjugate gradient with Fletcher-Reeves update, tight Wolfe line search
#' res <- mize(rb0, rosenbrock_fg, method = "CG", cg_update = "FR", c2 = 0.1)
#'
#' # Steepest decent with constant momentum = 0.9
#' res <- mize(rb0, rosenbrock_fg, method = "MOM", mom_schedule = 0.9)
#'
#' # Steepest descent with constant momentum in the Nesterov style as described
#' # in papers by Sutskever and Bengio
#' res <- mize(rb0, rosenbrock_fg, method = "MOM", mom_type = "nesterov",
#'              mom_schedule = 0.9)
#'
#' # Nesterov momentum with adaptive restart comparing function values
#' res <- mize(rb0, rosenbrock_fg, method = "MOM", mom_type = "nesterov",
#'              mom_schedule = 0.9, restart = "fn")
#' @export
mize <- function(par, fg,
                 method = "L-BFGS",
                 norm_direction = FALSE,
                 # L-BFGS
                 memory = 5,
                 scale_hess = TRUE,
                 # CG
                 cg_update = "PR+",
                 # NAG
                 nest_q = 0, # 1 - SD,
                 nest_convex_approx = FALSE,
                 nest_burn_in = 0,
                 # DBD
                 step_up = 1.1,
                 step_up_fun = "*",
                 step_down = NULL,
                 dbd_weight = 0.1,
                 # Line Search configuration
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
                 # Adaptive Restart
                 restart = NULL,
                 restart_wait = 10,
                 # Termination criterion
                 max_iter = 100,
                 max_fn = Inf,
                 max_gr = Inf,
                 max_fg = Inf,
                 abs_tol = sqrt(.Machine$double.eps),
                 rel_tol = abs_tol,
                 grad_tol = NULL,
                 ginf_tol = NULL,
                 step_tol = sqrt(.Machine$double.eps),
                 check_conv_every = 1,
                 log_every = check_conv_every,
                 verbose = FALSE,
                 store_progress = FALSE) {

  opt <- make_mize(method = method,
                   norm_direction = norm_direction,
                   scale_hess = scale_hess,
                   memory = memory,
                   cg_update = cg_update,
                   nest_q = nest_q, nest_convex_approx = nest_convex_approx,
                   nest_burn_in = nest_burn_in,
                   use_init_mom = use_init_mom,
                   step_up = step_up,
                   step_up_fun = step_up_fun,
                   step_down = step_down,
                   dbd_weight = dbd_weight,
                   line_search = line_search, step0 = step0, c1 = c1, c2 = c2,
                   step_next_init = step_next_init,
                   try_newton_step = try_newton_step,
                   ls_max_fn = ls_max_fn, ls_max_gr = ls_max_gr,
                   ls_max_fg = ls_max_fg,
                   ls_max_alpha_mult = ls_max_alpha_mult,
                   strong_curvature = strong_curvature,
                   approx_armijo = approx_armijo,
                   mom_type = mom_type,
                   mom_schedule = mom_schedule,
                   mom_init = mom_init,
                   mom_final = mom_final,
                   mom_switch_iter = mom_switch_iter,
                   mom_linear_weight = mom_linear_weight,
                   max_iter = max_iter,
                   restart = restart,
                   restart_wait = restart_wait)
  if (max_iter < 0) {
    stop("max_iter must be non-negative")
  }
  if (max_fn < 0) {
    stop("max_fn must be non-negative")
  }
  if (max_gr < 0) {
    stop("max_gr must be non-negative")
  }
  if (max_fg < 0) {
    stop("max_fg must be non-negative")
  }
  if (store_progress && is.null(check_conv_every)) {
    stop("check_conv_every must be non-NULL if store_progress is TRUE")
  }

  res <- opt_loop(opt, par, fg,
          max_iter = max_iter,
          max_fn = max_fn, max_gr = max_gr, max_fg = max_fg,
          abs_tol = abs_tol, rel_tol = rel_tol,
          grad_tol = grad_tol, ginf_tol = ginf_tol,
          step_tol = step_tol,
          check_conv_every = check_conv_every,
          log_every = log_every,
          store_progress = store_progress,
          verbose = verbose)

  Filter(Negate(is.null),
         res[c("f", "g2n", "ginfn", "nf", "ng", "par", "iter", "terminate",
               "progress")])
}

#' Create an Optimizer
#'
#' Factory function for creating a (possibly uninitialized) optimizer.
#'
#' If the function to be optimized and starting point are not present at
#' creation time, then the optimizer should be initialized using
#' \code{\link{mize_init}} before being used with \code{\link{mize_step}}.
#'
#' See the documentation to \code{\link{mize}} for an explanation of all the
#' parameters.
#'
#' Details of the \code{fg} list containing the function to be optimized and its
#' gradient can be found in the 'Details' section of \code{\link{mize}}. It is
#' optional for this function, but if it is passed to this function, along with
#' the vector of initial values, \code{par}, the optimizer will be returned
#' already initialized for this function. Otherwise, \code{\link{mize_init}}
#' must be called before optimization begins.
#'
#' Additionally, optional convergence parameters may also be passed here, for
#' use with \code{\link{check_mize_convergence}}. They are optional here if you
#' plan to call \code{\link{mize_init}} later, or if you want to do your own
#' convergence checking.
#'
#' @param method Optimization method. See 'Details' of \code{\link{mize}}.
#' @param norm_direction If \code{TRUE}, then the steepest descent direction is
#'   normalized to unit length. Useful for adaptive step size methods where the
#'   previous step size is used to initialize the next iteration.
#' @param scale_hess if \code{TRUE}, the approximation to the inverse Hessian is
#'   scaled according to the method described by Nocedal and Wright
#'   (approximating an eigenvalue). Applies only to the methods \code{BFGS}
#'   (where the scaling is applied only during the first step) and \code{L-BFGS}
#'   (where the scaling is applied during every iteration). Ignored otherwise.
#' @param memory The number of updates to store if using the \code{L-BFGS}
#'   method. Ignored otherwise. Must be a positive integer.
#' @param cg_update Type of update to use for the \code{CG} method. Can be one
#'   of \code{"FR"} (Fletcher-Reeves), \code{"PR"} (Polak-Ribiere), \code{"PR+"}
#'   (Polak-Ribiere with a reset to steepest descent), \code{"HS"}
#'   (Hestenes-Stiefel), or \code{"DY"} (Dai-Yuan). Ignored if \code{method} is
#'   not \code{"CG"}.
#' @param nest_q Strong convexity parameter for the \code{"NAG"} method's
#'   momentum term. Must take a value between 0 (strongly convex) and 1 (results
#'   in steepest descent).Ignored unless the \code{method} is \code{"NAG"} and
#'   \code{nest_convex_approx} is \code{FALSE}.
#' @param nest_convex_approx If \code{TRUE}, then use an approximation due to
#'   Sutskever for calculating the momentum parameter in the NAG method. Only
#'   applies if \code{method} is \code{"NAG"}.
#' @param nest_burn_in Number of iterations to wait before using a non-zero
#'   momentum. Only applies if using the \code{"NAG"} method or setting the
#'   \code{momentum_type} to "Nesterov".
#' @param step_up Value by which to increase the step size for the \code{"bold"}
#'   step size method or the \code{"DBD"} method.
#' @param step_up_fun Operator to use when combining the current step size with
#'   \code{step_up}. Can be one of \code{"*"} (to multiply the current step size
#'   with \code{step_up}) or \code{"+"} (to add).
#' @param step_down Multiplier to reduce the step size by if using the
#'   \code{"DBD"} method or the \code{"bold"}. Can also be used with the
#'   \code{"back"} line search method, but is optional. Should be a positive
#'   value less than 1.
#' @param dbd_weight Weighting parameter used by the \code{"DBD"} method only,
#'   and only if no momentum scheme is provided. Must be an integer between 0
#'   and 1.
#' @param line_search Type of line search to use. See 'Details' of
#'   \code{\link{mize}}.
#' @param c1 Sufficient decrease parameter for Wolfe-type line searches. Should
#'   be a value between 0 and 1.
#' @param c2 Sufficient curvature parameter for line search for Wolfe-type line
#'   searches. Should be a value between \code{c1} and 1.
#' @param step0 Initial value for the line search on the first step. See
#'   'Details' of \code{\link{mize}}.
#' @param step_next_init For Wolfe-type line searches only, how to initialize
#'   the line search on iterations after the first. See 'Details' of
#'   \code{\link{mize}}.
#' @param try_newton_step For Wolfe-type line searches only, try the line step
#'   value of 1 as the initial step size whenever \code{step_next_init} suggests
#'   a step size > 1. Defaults to \code{TRUE} for quasi-Newton methods such as
#'   BFGS and L-BFGS, \code{FALSE} otherwise.
#' @param ls_max_fn Maximum number of function evaluations allowed during a line
#'   search.
#' @param ls_max_gr Maximum number of gradient evaluations allowed during a line
#'   search.
#' @param ls_max_fg Maximum number of function or gradient evaluations allowed
#'   during a line search.
#' @param ls_max_alpha_mult Maximum multiplier for alpha between iterations.
#'   Only applies for Wolfe-type line searches and if \code{step_next_init} is
#'   set to \code{"slope"}
#' @param strong_curvature (Optional). If \code{TRUE} use the strong
#'   curvature condition in Wolfe line search. See the 'Line Search' section of
#'   \code{\link{mize}} for details.
#' @param approx_armijo (Optional). If \code{TRUE} use the approximate Armijo
#'   condition in Wolfe line search. See the 'Line Search' section of
#'   \code{\link{mize}} for details.
#' @param mom_type Momentum type, either \code{"classical"} or
#'   \code{"nesterov"}.
#' @param mom_schedule Momentum schedule. See 'Details' of \code{\link{mize}}.
#' @param mom_init Initial momentum value.
#' @param mom_final Final momentum value.
#' @param mom_switch_iter For \code{mom_schedule} \code{"switch"} only, the
#'   iteration when \code{mom_init} is changed to \code{mom_final}.
#' @param mom_linear_weight If \code{TRUE}, the gradient contribution to the
#'   update is weighted using momentum contribution.
#' @param use_init_mom If \code{TRUE}, then the momentum coefficient on the
#'   first iteration is non-zero. Otherwise, it's zero. Only applies if using a
#'   momentum schedule.
#' @param restart Momentum restart type. Can be one of "fn" or "gr". See
#'   'Details' of \code{\link{mize}}.
#' @param restart_wait Number of iterations to wait between restarts. Ignored if
#'   \code{restart} is \code{NULL}.
#' @param par (Optional) Initial values for the function to be optimized over.
#' @param fg (Optional). Function and gradient list. See 'Details' of
#'   \code{\link{mize}}.
#' @param max_iter (Optional). Maximum number of iterations. See the
#'   'Convergence' section of \code{\link{mize}} for details.
#' @param max_fn (Optional). Maximum number of function evaluations. See the
#'   'Convergence' section of \code{\link{mize}} for details.
#' @param max_gr (Optional). Maximum number of gradient evaluations. See the
#'   'Convergence' section of \code{\link{mize}} for details.
#' @param max_fg (Optional). Maximum number of function or gradient evaluations.
#'   See the 'Convergence' section of \code{\link{mize}} for details.
#' @param abs_tol (Optional). Absolute tolerance for comparing two function
#'   evaluations. See the 'Convergence' section of \code{\link{mize}} for
#'   details.
#' @param rel_tol (Optional). Relative tolerance for comparing two function
#'   evaluations. See the 'Convergence' section of \code{\link{mize}} for
#'   details.
#' @param grad_tol (Optional). Absolute tolerance for the length (l2-norm) of
#'   the gradient vector. See the 'Convergence' section of \code{\link{mize}}
#'   for details.
#' @param ginf_tol (Optional). Absolute tolerance for the infinity norm (maximum
#'   absolute component) of the gradient vector. See the 'Convergence' section
#'   of \code{\link{mize}} for details.
#' @param step_tol (Optional). Absolute tolerance for the size of the parameter
#'   update. See the 'Convergence' section of \code{\link{mize}} for details.
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
#' opt <- make_mize(method = "L-BFGS", par = rb0, fg = rosenbrock_fg)
#'
#' # Create optimizer without initialization
#' opt <- make_mize(method = "L-BFGS")
#'
#' # Need to call mizer_init separately:
#' opt <- mize_init(opt, rb0, rosenbrock_fg)
make_mize <- function(method = "L-BFGS",
                      norm_direction = FALSE,
                      # BFGS
                      scale_hess = TRUE,
                      memory = 5,
                      # CG
                      cg_update = "PR+",
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
                      c1 = 1e-4, c2 = NULL,
                      step0 = NULL,
                      step_next_init = NULL,
                      try_newton_step = NULL,
                      ls_max_fn = 20,
                      ls_max_gr = Inf,
                      ls_max_fg = Inf,
                      ls_max_alpha_mult = Inf,
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
                      max_fn = Inf, max_gr = Inf, max_fg = Inf,
                      abs_tol = NULL,
                      rel_tol = abs_tol, grad_tol = NULL, ginf_tol = NULL,
                      step_tol = NULL) {

  if (memory < 1) {
    stop("memory must be > 0")
  }
  if (!is_in_range(nest_q, 0, 1)) {
    stop("nest_q must be between 0 and 1")
  }
  if (nest_burn_in < 0) {
    stop("nest_burn_in must be non-negative")
  }
  if (step_up <= 0) {
    stop("step_up must be positive")
  }
  step_up_fun <- match.arg(step_up_fun)
  if (!is.null(step_down) && !is_in_range(step_down, 0, 1)) {
    stop("step_down must be between 0 and 1")
  }
  if (!is_in_range(dbd_weight, 0, 1)) {
    stop("dbd_weight must be between 0 and 1")
  }
  if (!is_in_range(c1, 0, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c1 must be between 0 and 1")
  }
  if (!is.null(c2) && !is_in_range(c2, c1, 1, lopen = FALSE, ropen = FALSE)) {
    stop("c2 must be between c1 and 1")
  }
  if (ls_max_fn < 0) {
    stop("ls_max_fn must be non-negative")
  }
  if (ls_max_gr < 0) {
    stop("ls_max_gr must be non-negative")
  }
  if (ls_max_fg < 0) {
    stop("ls_max_fg must be non-negative")
  }
  if (ls_max_alpha_mult <= 0) {
    stop("ls_max_alpha_mult must be positive")
  }
  if (restart_wait < 1) {
    stop("restart_wait must be a positive integer")
  }

  # Gradient Descent Direction configuration
  dir_type <- NULL
  method <- match.arg(tolower(method), c("sd", "newton", "phess", "cg", "bfgs",
                                "l-bfgs", "nag", "momentum", "dbd"))
  switch(method,
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
      cg_update <- match.arg(tolower(cg_update),
                             c("fr", "cd", "dy",
                               "hs", "hs+", "pr", "pr+", "ls", "hz", "hz+"))
      cg_update_fn <- switch(cg_update,
        fr = fr_update,
        cd = cd_update,
        dy = dy_update,
        hs = hs_update,
        "hs+" = hs_plus_update,
        pr = pr_update,
        "pr+" = pr_plus_update,
        ls = ls_update,
        hz = hz_update,
        "hz+" = hz_plus_update
      )
      dir_type <- cg_direction(cg_update = cg_update_fn)
    },
    bfgs = {
      dir_type <- bfgs_direction(scale_inverse = scale_hess)
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
    stop("Unknown method: '", method, "'")
  )

  # If it's not already been turned on, turn off the Newton step option
  if (is.null(try_newton_step)) {
    try_newton_step <- FALSE
  }

  # Line Search configuration
  step_type <- NULL
  line_search <- tolower(line_search)
  if (method == "dbd") {
    if (is.character(step0) || is.numeric(step0)) {
      eps_init <- step0
    }
    else {
      eps_init <- "rasmussen"
    }
    if (step_up_fun == "*") {
      step_up_fun <- `*`
    }
    else if (step_up_fun == "+") {
      step_up_fun <- `+`
    }
    else {
      stop("Unknown delta-bar-delta step_up function '", step_up_fun, "'")
    }
    if (is.null(step_down)) {
      step_down <- 0.5
    }
    step_type <- delta_bar_delta(epsilon = eps_init,
                                 kappa = step_up, kappa_fun = step_up_fun,
                                 phi = step_down, theta = dbd_weight,
                                 use_momentum = is.null(mom_schedule))
  }
  else {
    if (method %in% c("newton", "phess", "bfgs", "l-bfgs")) {
      if (is.null(c2)) {
        c2 <- 0.9
      }
      if (is.null(try_newton_step)) {
        try_newton_step <- TRUE
      }
    }
    else {
      if (is.null(c2)) {
        c2 <- 0.1
      }
      if (is.null(try_newton_step)) {
        try_newton_step <- FALSE
      }
    }

    line_search <- match.arg(tolower(line_search),
                             c("more-thuente", "mt", "rasmussen",
                               "bold driver",
                               "backtracking", "constant",
                               "schmidt", "minfunc", "armijo",
                               "hager-zhang", "hz"))
    if (line_search == "hager-zhang") {
      line_search <- "hz"
    }
    if (line_search == "more-thuente") {
      line_search <- "mt"
    }
    if (line_search == "minfunc") {
      line_search <- "schmidt"
    }

    if (line_search == "bold driver") {
      if (is.null(step_down)) {
        step_down <- 0.5
      }
    }

    # Set Wolfe line search termination defaults
    # Most Wolfe Line Searches use the standard Strong Wolfe conditions
    if (line_search %in% c("more-thuente", "mt", "rasmussen", "schmidt",
                           "minfunc")) {
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
    }
    else {
      if (is.null(step0)) {
        step0 <- "rasmussen"
      }
      if (is.null(step_next_init)) {
        step_next_init <- "quad"
      }
    }

    step_type <- switch(line_search,
      mt = more_thuente_ls(c1 = c1, c2 = c2,
                           initializer = tolower(step_next_init),
                           initial_step_length = step0,
                           try_newton_step = try_newton_step,
                           max_fn = ls_max_fn,
                           max_gr = ls_max_gr,
                           max_fg = ls_max_fg,
                           max_alpha_mult = ls_max_alpha_mult,
                           strong_curvature = strong_curvature,
                           approx_armijo = approx_armijo),
      rasmussen = rasmussen_ls(c1 = c1, c2 = c2,
                              initializer = tolower(step_next_init),
                              initial_step_length = step0,
                              try_newton_step = try_newton_step,
                              max_fn = ls_max_fn,
                              max_gr = ls_max_gr,
                              max_fg = ls_max_fg,
                              max_alpha_mult = ls_max_alpha_mult,
                              strong_curvature = strong_curvature,
                              approx_armijo = approx_armijo),
      "bold driver" = bold_driver(inc_mult = step_up, dec_mult = step_down,
                                  max_fn = ls_max_fn),
      constant = constant_step_size(value = step0),
      schmidt = schmidt_ls(c1 = c1, c2 = c2,
                           initializer = tolower(step_next_init),
                           initial_step_length = step0,
                           try_newton_step = try_newton_step,
                           max_fn = ls_max_fn,
                           max_gr = ls_max_gr,
                           max_fg = ls_max_fg,
                           max_alpha_mult = ls_max_alpha_mult,
                           strong_curvature = strong_curvature,
                           approx_armijo = approx_armijo),
      backtracking = schmidt_armijo_ls(c1 = c1,
                          initializer = tolower(step_next_init),
                          initial_step_length = step0,
                          try_newton_step = try_newton_step,
                          step_down = step_down,
                          max_fn = ls_max_fn,
                          max_gr = ls_max_gr,
                          max_fg = ls_max_fg,
                          max_alpha_mult = ls_max_alpha_mult),
      hz =  hager_zhang_ls(c1 = c1, c2 = c2,
                           initializer = tolower(step_next_init),
                           initial_step_length = step0,
                           try_newton_step = try_newton_step,
                           max_fn = ls_max_fn,
                           max_gr = ls_max_gr,
                           max_fg = ls_max_fg,
                           max_alpha_mult = ls_max_alpha_mult,
                           strong_curvature = strong_curvature,
                           approx_armijo = approx_armijo)
    )
  }

  # Create Gradient Descent stage
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = dir_type,
        step_size = step_type)))

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
  }
  else if (method == "momentum") {
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
        use_init_mom = use_init_mom)
    }
    else if (is.function(mom_schedule)) {
      mom_step <- make_momentum_step(mu_fn = mom_schedule,
                                     use_init_mom = use_init_mom)
    }
    else {
      mom_schedule <- match.arg(tolower(mom_schedule),
                                c("ramp", "switch", "nsconvex"))

      mom_step <- switch(mom_schedule,
        ramp = make_momentum_step(
          make_ramp(init_value = mom_init,
                    final_value = mom_final,
                    wait = ifelse(use_init_mom, 0, 1)),
          use_init_mom = use_init_mom),
        "switch" = make_momentum_step(
          make_switch(
            init_value = mom_init,
            final_value = mom_final,
            switch_iter = mom_switch_iter),
          use_init_mom = use_init_mom),
        nsconvex = nesterov_step(burn_in = nest_burn_in, q = nest_q,
                                 use_approx = nest_convex_approx,
                                 use_init_mu = use_init_mom)
        )
    }

    mom_stage <- momentum_stage(
      direction = mom_direction,
      step_size = mom_step)

    opt <- append_stage(opt, mom_stage)

    if (mom_linear_weight) {
      opt <- append_stage(opt, momentum_correction_stage())
    }
  }

  # Adaptive Restart
  if (!is.null(restart)) {
    restart <- match.arg(tolower(restart), c("none", "fn", "gr", "speed"))
    if (restart != "none") {
      opt <- adaptive_restart(opt, restart, wait = restart_wait)
    }
  }

  # Initialize for specific dataset if par and fg are provided
  if (!is.null(par) && !is.null(fg)) {
    opt <- mize_init(opt, par, fg, max_iter = max_iter,
                     max_fn = max_fn, max_gr = max_gr, max_fg = max_fg,
                     abs_tol = abs_tol, rel_tol = rel_tol,
                     grad_tol = grad_tol, ginf_tol = ginf_tol,
                     step_tol = step_tol)
  }

  opt
}

#'One Step of Optimization
#'
#'Performs one iteration of optimization using a specified optimizer.
#'
#'This function returns both the (hopefully) optimized vector of parameters, and
#'an updated version of the optimizer itself. This is intended to be used when
#'you want more control over the optimization process compared to the more black
#'box approach of the \code{\link{mize}} function. In return for having to
#'manually call this function every time you want the next iteration of
#'optimization, you gain the ability to do your own checks for convergence,
#'logging and so on, as well as take other action between iterations, e.g.
#'visualization.
#'
#'Normally calling this function should return a more optimized vector of
#'parameters than the input, or at  least leave the parameters unchanged if no
#'improvement was found, although this is determined by how the optimizer was
#'configured by \code{\link{make_mize}}. It is very possible to create an
#'optimizer that can cause a solution to diverge. It is the responsibility of
#'the caller to check that the result of the optimization step has actually
#'reduced the value returned from function being optimized.
#'
#'Details of the \code{fg} list can be found in the 'Details' section of
#'\code{\link{mize}}.
#'
#'@param opt Optimizer, created by \code{\link{make_mize}}.
#'@param par Vector of initial values for the function to be optimized over.
#'@param fg Function and gradient list. See the documentation of
#'  \code{\link{mize}}.
#'@return Result of the current optimization step, a list with components:
#'  \itemize{
#'
#'  \item{\code{opt}}. Updated version of the optimizer passed to the \code{opt}
#'  argument Should be passed as the \code{opt} argument in the next iteration.
#'
#'  \item{\code{par}}. Updated version of the parameters passed to the
#'  \code{par} argument. Should be passed as the \code{par} argument in the next
#'  iteration.
#'
#'  \item{\code{nf}}. Running total number of function evaluations carried out
#'  since iteration 1.
#'
#'  \item{\code{ng}}. Running total number of gradient evaluations carried out
#'  since iteration 1.
#'
#'  \item{\code{f}}. Optional. The new value of the function, evaluated at the
#'  returned value of \code{par}. Only present if calculated as part of the
#'  optimization step (e.g. during a line search calculation).
#'
#'  \item{\code{g}}. Optional. The gradient vector, evaluated at the returned
#'  value of \code{par}. Only present if the gradient was calculated as part of
#'  the optimization step (e.g. during a line search calculation.)}
#'
#'@seealso \code{\link{make_mize}} to create a value to pass to \code{opt},
#'  \code{\link{mize_init}} to initialize \code{opt} before passing it to this
#'  function for the first time. \code{\link{mize}} creates an optimizer and
#'  carries out a full optimization with it.
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
#'  opt <- make_mize(method = "SD", line_search = "const", step0 = 0.0001,
#'                   par = rb0, fg = rosenbrock_fg)
#'  par <- rb0
#'  for (iter in 1:3) {
#'    res <- mize_step(opt, par, rosenbrock_fg)
#'    par <- res$par
#'    opt <- res$opt
#'  }
#'@export
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
    }
    else {
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
    opt <- life_cycle_hook("validation", "before", opt, par, fg, iter,
                           par0, step_result)
    opt <- life_cycle_hook("validation", "during", opt, par, fg, iter,
                           par0, step_result)
  }
  # If the this solution was vetoed or the catastrophe happened,
  # roll back to the previous one.
  if (!is.null(opt$terminate) || !opt$ok) {
    par <- par0
  }

  if (is.null(opt$terminate)) {
    opt <- life_cycle_hook("step", "after", opt, par, fg, iter, par0,
                         step_result)
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
#' Should be called after creating an optimizer with \code{\link{make_mize}} and
#' before beginning any optimization with \code{\link{mize_step}}. Note that if
#' \code{fg} and \code{par} are available at the time \code{\link{mize_step}} is
#' called, they can be passed to that function and initialization will be
#' carried out automatically, avoiding the need to call \code{mize_init}.
#'
#' Optional convergence parameters may also be passed here, for use with
#' \code{\link{check_mize_convergence}}. They are optional if you do your own
#' convergence checking.
#'
#' Details of the \code{fg} list can be found in the 'Details' section of
#' \code{\link{mize}}.
#'
#' @param opt Optimizer, created by \code{\link{make_mize}}.
#' @param par Vector of initial values for the function to be optimized over.
#' @param fg Function and gradient list. See the documentation of
#'   \code{\link{mize}}.
#' @param max_iter (Optional). Maximum number of iterations. See the
#'   'Convergence' section of \code{\link{mize}} for details.
#' @param max_fn (Optional). Maximum number of function evaluations. See the
#'   'Convergence' section of \code{\link{mize}} for details.
#' @param max_gr (Optional). Maximum number of gradient evaluations. See the
#'   'Convergence' section of \code{\link{mize}} for details.
#' @param max_fg (Optional). Maximum number of function or gradient evaluations.
#'   See the 'Convergence' section of \code{\link{mize}} for details.
#' @param abs_tol (Optional). Absolute tolerance for comparing two function
#'   evaluations. See the 'Convergence' section of \code{\link{mize}} for
#'   details.
#' @param rel_tol (Optional). Relative tolerance for comparing two function
#'   evaluations. See the 'Convergence' section of \code{\link{mize}} for
#'   details.
#' @param grad_tol (Optional). Absolute tolerance for the length (l2-norm) of
#'   the gradient vector. See the 'Convergence' section of \code{\link{mize}}
#'   for details.
#' @param ginf_tol (Optional). Absolute tolerance for the infinity norm (maximum
#'   absolute component) of the gradient vector. See the 'Convergence' section
#'   of \code{\link{mize}} for details.
#' @param step_tol (Optional). Absolute tolerance for the size of the parameter
#'   update. See the 'Convergence' section of \code{\link{mize}} for details.
#' @return Initialized optimizer.
#' @export
#' @examples
#'
#' # Create an optimizer
#' opt <- make_mize(method = "L-BFGS")
#'
#' # Function to optimize and starting point defined after creating optimizer
#' rosenbrock_fg <- list(
#'   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
#'   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                          200 *        (x[2] - x[1] * x[1])) })
#' rb0 <- c(-1.2, 1)
#'
#' # Initialize with function and starting point before commencing optimization
#' opt <- mize_init(opt, rb0, rosebrock_fg)
#'
#' # Finally, can commence the optimization loop
#' par <- rb0
#' for (iter in 1:3) {
#'   res <- mize_step(opt, par, rosenbrock_fg)
#'   par <- res$par
#'   opt <- res$opt
#' }
#'
mize_init <- function(opt, par, fg,
                      max_iter = Inf,
                      max_fn = Inf, max_gr = Inf, max_fg = Inf,
                      abs_tol = NULL,
                      rel_tol = abs_tol, grad_tol = NULL, ginf_tol = NULL,
                      step_tol = NULL) {
  opt <- register_hooks(opt)
  opt$iter <- 0
  opt <- life_cycle_hook("opt", "init", opt, par, fg, opt$iter)
  opt$convergence <- list(
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
  opt$is_initialized <- TRUE
  opt
}

#'Mize Step Summary
#'
#'Produces a result summary for an optimization iteration. Information such as
#'function value, gradient norm and step size may be returned.
#'
#'By default, convergence tolerance parameters will be used to determine what
#'function and gradient data is returned. The function value will be returned if
#'it was already calculated and cached in the optimization iteration. Otherwise,
#'it will be calculated only if a non-null absolute or relative tolerance value
#'was asked for. A gradient norm will be returned only if a non-null gradient
#'tolerance was specified, even if the gradient is available.
#'
#'Note that if a function tolerance was specified, but was not calculated for
#'the relevant value of \code{par}, they will be calculated here and the
#'calculation does contribute to the total function count (and will be cached
#'for potential use in the next iteration). The same applies for gradient
#'tolerances and gradient calculation. Function and gradient calculation can
#'also be forced here by setting the \code{calc_fn} and \code{calc_gr}
#'(respectively) parameters to \code{TRUE}.
#'
#'@param opt Optimizer to generate summary for, from return value of
#'  \code{\link{mize_step}}.
#'@param par Vector of parameters at the end of the iteration, from return value
#'  of \code{\link{mize_step}}.
#'@param fg Function and gradient list. See the documentation of
#'  \code{\link{mize}}.
#'@param par_old (Optional). Vector of parameters at the end of the previous
#'  iteration. Used to calculate step size.
#'@param calc_fn (Optional). If \code{TRUE}, force calculation of function if
#'  not already cached in \code{opt}, even if it would not be needed for
#'  convergence checking.
#'@return A list with the following items: \itemize{
#'
#'  \item \code{opt} Optimizer with updated state (e.g. function and gradient
#'  counts).
#'
#'  \item \code{iter} Iteration number.
#'
#'  \item \code{f} Function value at \code{par}.
#'
#'  \item \code{g2n} 2-norm of the gradient at \code{par}.
#'
#'  \item \code{ginfn} Infinity-norm of the gradient at \code{par}.
#'
#'  \item \code{nf} Number of function evaluations so far.
#'
#'  \item \code{ng} Number of gradient evaluations so far.
#'
#'  \item \code{step} Size of the step between \code{par_old} and \code{par},
#'  if \code{par_old} is provided.
#'
#'  \item \code{alpha} Step length of the gradient descent part of the step.
#'
#'  \item \code{mu} Momentum coefficient for this iteration}
#'@export
#'@examples
#' rb_fg <- list(
#'   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
#'   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                          200 *        (x[2] - x[1] * x[1])) })
#' rb0 <- c(-1.2, 1)
#'
#' opt <- make_mize(method = "BFGS", par = rb0, fg = rb_fg, max_iter = 30)
#' mize_res <- mize_step(opt = opt, par = rb0, fg = rb_fg)
#' # Get info about first step, use rb0 to compare new par with initial value
#' step_info <- mize_step_summary(mize_res$opt, mize_res$par, rb_fg, rb0)
mize_step_summary <- function(opt, par, fg, par_old = NULL, calc_fn = NULL) {

  iter <- opt$iter
  # An internal flag useful for unit tests: if FALSE, doesn't count any
  # fn/gr calculations towards their counts. Can still get back fn/gr values
  # without confusing the issue of the expected number of fn/gr evaluations
  if (!is.null(opt$count_res_fg)) {
    count_fg <- opt$count_res_fg
  }
  else {
    count_fg <- TRUE
  }

  # Whether and what convergence info to calculate if fn/gr calculation not
  # explicitly asked for
  if (is.null(calc_fn)) {
    calc_fn <- is_finite_numeric(opt$convergence$abs_tol) ||
      is_finite_numeric(opt$convergence$rel_tol)
  }

  gr_norms <- c()
  if (is_finite_numeric(opt$convergence$grad_tol)) {
    gr_norms <- c(gr_norms, 2)
  }
  if (is_finite_numeric(opt$convergence$ginf_tol)) {
    gr_norms <- c(gr_norms, Inf)
  }
  calc_gr <- length(gr_norms) > 0

  f <- NULL
  if (calc_fn || has_fn_curr(opt, iter + 1)) {
    if (!has_fn_curr(opt, iter + 1)) {
      f <- fg$fn(par)
      if (count_fg) {
        opt <- set_fn_curr(opt, f, iter + 1)
        opt$counts$fn <- opt$counts$fn + 1
      }
    }
    else {
      f <- opt$cache$fn_curr
    }
  }

  g2n <- NULL
  ginfn <- NULL
  if (calc_gr || has_gr_curr(opt, iter + 1)) {
    if (!has_gr_curr(opt, iter + 1)) {
      g <- fg$gr(par)
      if (grad_is_first_stage(opt) && count_fg) {
        opt <- set_gr_curr(opt, g, iter + 1)
        opt$counts$gr <- opt$counts$gr + 1
      }
    }
    else {
      g <- opt$cache$gr_curr
    }
    if (2 %in% gr_norms) {
      g2n <- norm2(g)
    }
    if (Inf %in% gr_norms) {
      ginfn <- norm_inf(g)
    }
  }

  if (!is.null(par_old)) {
    step_size <- norm2(par - par_old)
  }
  else {
    step_size <- 0
  }

  alpha <- 0
  if (!is.null(opt$stages[["gradient_descent"]]$step_size$value)) {
    alpha <- norm2(opt$stages[["gradient_descent"]]$step_size$value)
    if (is.null(alpha)) {
      alpha <- 0
    }
  }

  res <- list(
    opt = opt,
    f = f,
    g2n = g2n,
    ginfn = ginfn,
    nf = opt$counts$fn,
    ng = opt$counts$gr,
    step = step_size,
    alpha = alpha,
    iter = iter
  )

  if ("momentum" %in% names(opt$stages)) {
    res$mu <- opt$stages[["momentum"]]$step_size$value
    if (is.null(res$mu)) {
      res$mu <- 0
    }
  }

  Filter(Negate(is.null), res)
}


#' Check Optimization Convergence
#'
#' Updates the optimizer with information about convergence or termination,
#' signaling if the optimization process should stop.
#'
#' On returning from this function, the updated value of \code{opt} will
#' contain: \itemize{
#'
#' \item A boolean value \code{is_terminated} which is \code{TRUE} if
#' termination has been indicated, and \code{FALSE} otherwise.
#'
#' \item A list \code{terminate} if \code{is_terminated} is \code{TRUE}. This
#' contains two items: \code{what}, a short string describing what caused the
#' termination, and \code{val}, the value of the termination criterion that
#' caused termination. This list will not be present if \code{is_terminated} is
#' \code{FALSE}.}
#'
#' Convergence criteria are only checked here. To set these criteria, use
#' \code{\link{make_mize}} or \code{\link{mize_init}}.
#'
#' @param mize_step_info Step info for this iteration, created by
#'   \code{\link{mize_step_summary}}
#' @return \code{opt} updated with convergence and termination data. See
#'   'Details'.
#' @export
#'@examples
#' rb_fg <- list(
#'   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
#'   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                          200 *        (x[2] - x[1] * x[1])) })
#' rb0 <- c(-1.2, 1)
#'
#' opt <- make_mize(method = "BFGS", par = rb0, fg = rb_fg, max_iter = 30)
#' mize_res <- mize_step(opt = opt, par = rb0, fg = rb_fg)
#' step_info <- mize_step_summary(mize_res$opt, mize_res$par, rb_fg, rb0)
#' # check convergence by looking at opt$is_terminated
#' opt <- check_mize_convergence(step_info)
check_mize_convergence <- function(mize_step_info) {

  opt <- mize_step_info$opt

  convergence <- opt$convergence

  terminate <- check_counts(opt, convergence$max_fn, convergence$max_gr,
                            convergence$max_fg)
  if (!is.null(terminate)) {
    opt$terminate <- terminate
    opt$is_terminated <- TRUE
    return(opt)
  }

  terminate <- check_step_conv(opt, opt$iter, mize_step_info$step,
                               convergence$step_tol)
  if (!is.null(terminate)) {
    opt$terminate <- terminate
    opt$is_terminated <- TRUE
    return(opt)
  }

  terminate <- check_gr_conv(opt, convergence$grad_tol, convergence$ginf_tol)
  if (!is.null(terminate)) {
    opt$terminate <- terminate
    opt$is_terminated <- TRUE
    return(opt)
  }

  if (!is.null(opt$cache$fn_curr)) {
    fn_new <- opt$cache$fn_curr
    fn_old <- convergence$fn_new
    convergence$fn_new <- fn_new
    opt$convergence <- convergence

    terminate <- check_fn_conv(opt, opt$iter, fn_old, fn_new,
                               convergence$abs_tol, convergence$rel_tol)
    if (!is.null(terminate)) {
      opt$is_terminated <- TRUE
      opt$terminate <- terminate
      return(opt)
    }
  }

  if (opt$iter == convergence$max_iter) {
    opt$is_terminated <- TRUE
    opt$terminate <- list(what = "max_iter", val = convergence$max_iter)
  }

  opt
}
