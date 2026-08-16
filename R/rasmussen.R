# Rasmussen Line Search
#
# Line search algorithm originally written by Carl Edward Rasmussen in his
# conjugate gradient routine. It consists of two main parts:
#
# 1. Using cubic extrapolation from an initial starting guess for the step
#    size until either the sufficient decrease condition is not met or the
#    curvature condition is met.
# 2. Interpolation (quadratic or cubic) between that point and the start
#    point of the search until either a step size is found which meets the
#    Strong Wolfe conditions or the maximum number of allowed function
#    evaluations is reached.
#
# The extrapolation and interpolation steps are bounded at each stage to ensure
# they don't represent too large or small a change to the step size.
#
# @seealso This implementation was informed by Carl Edward Rasmussen's
#  `minimize.m` routine in the Matlab
#  [GPML](https://www.gaussianprocess.org/gpml/code/matlab/doc/) package.
rasmussen_core <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  direction,
  method_policy
) {
  run_bracket_zoom(
    evaluator = evaluator,
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    condition_policy = condition_policy,
    direction = direction,
    policy = method_policy
  )
}
