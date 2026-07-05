# Various Gradient-based optimization routines: steepest descent, conjugate
# gradient, BFGS etc.

# Gradient Direction -----------------------------------------------------------

# Creates a direction sub stage
make_direction <- function(sub_stage) {
  make_sub_stage(sub_stage, "direction")
}

# Steepest Descent
#
# normalize - If TRUE, then the returned direction vector is normalized to unit
# length. This can be useful for some adaptive line search methods, so that the
# total step length is purely determined by the line search value, rather than
# the product of the line search value and the magnitude of the direction.
sd_direction <- function(normalize = FALSE) {
  make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      opt$cache$gr_curr <- NULL
      list(opt = opt)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$value <- -opt$cache$gr_curr

      if (sub_stage$normalize) {
        sub_stage$value <- normalize(sub_stage$value)
      }

      list(sub_stage = sub_stage)
    },
    normalize = normalize
  ))
}
