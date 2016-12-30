# Adaptive Restart --------------------------------------------------------

# Adds Adaptive Restart for optimizers which are using a momentum scheme.
#
# Candes and O'Donoghue suggested a restart scheme for Nesterov Accelerated
# Gradient schemes to avoid oscillatory behavior. It effectively restarts the
# momentum part of the optimization and can hence be applied to any optimization
# that uses momentum.
#
# There are two ways to check if a restart is needed: comparing the direction
# of the optimization with the gradient to see if the direction is a descent
# direction (gradient-based validation), or by comparing function evaluation
# before and after the step (function-based validation). Normally, the gradient
# based validation is cheaper, because the gradient has to be calculated
# anyway, but if using the momentum-version of NAG, this isn't available,
# because the gradient is calculated at the wrong position. If using a
# line search that calculates the function, then you probably get the function
# based validation for free.
#
# opt An optimizer
# validation_type - one of "fn" or "gr" for function- or gradient-based
# validation, respectively.
# wait - number of iterations after a restart to wait before validating again.
adaptive_restart <- function(opt, validation_type, wait = 1) {
  stage_names <- names(opt$stages)
  if (!("momentum" %in% stage_names)) {
    stop("Adaptive restart is only applicable to optimizers which use momentum")
  }
  if (stage_names[[2]] != "momentum" && validation_type == "gr") {
    stop("Can't use gradient-based adaptive restart with Nesterov momentum")
  }
  opt$restart_wait <- wait
  append_depends(opt, "momentum", "direction",
                 c('adaptive_restart', paste0('validate_', validation_type)))
}

# Function-based adaptive restart
adaptive_restart_fn <- function(opt) {
  adaptive_restart(opt, "fn")
}

# Gradient-based adaptive restart
adaptive_restart_gr <- function(opt) {
  adaptive_restart(opt, "gr")
}


# Replace the usual momentum after step event with one which restarts
# the momentum if validation failed
require_adaptive_restart <- function(opt, par, fg, iter, par0, update) {
  if (!opt$ok && can_restart(opt, iter)) {
    opt <- life_cycle_hook("momentum", "init", opt, par, fg, iter)
    opt$restart_at <- iter
  }
  else {
    opt$cache$update_old <- update
  }
  opt
}
attr(require_adaptive_restart, 'event') <- 'after step'
# Should have the same name as normal update old: we want to replace that hook
attr(require_adaptive_restart, 'name') <- 'update_old'
attr(require_adaptive_restart, 'depends') <- 'update_old_init'

# Add a depend function to one of opt, a stage or sub stage
append_depends <- function(opt, stage_type = NULL, sub_stage_type = NULL,
                           new_depends) {
  if (!is.null(sub_stage_type)) {
    if (is.null(stage_type)) {
      stop("Must provide stage for sub_stage '", sub_stage_type, "'")
    }
    stage <- opt$stages[[stage_type]]
    if (is.null(stage)) {
      stop("No stage '", stage_type, "' exists for this optimizer")
    }
    depends <- stage[[sub_stage_type]]$depends
  }
  else if (!is.null(stage)) {
    stage <- opt$stages[[stage_type]]
    if (is.null(stage)) {
      stop("No stage '", stage_type, "' exists for this optimizer")
    }
    depends <- stage$depends
  }
  else {
    depends <- opt$depends
  }
  if (is.null(depends)) {
    depends <- c()
  }

  depends <- c(depends, new_depends)

  if (!is.null(sub_stage_type)) {
    opt$stages[[stage_type]][[sub_stage_type]]$depends <- depends
  }
  else if (!is.null(stage)) {
    opt$stages[[stage_type]]$depends <- depends
  }
  else {
    opt$depends <- depends
  }

  opt
}

# True if we aren't currently waiting between restarts
can_restart <- function(opt, iter) {
  is.null(opt$restart_at) || iter - opt$restart_at > opt$restart_wait
}

