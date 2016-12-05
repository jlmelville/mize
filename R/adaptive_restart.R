# Adaptive Restart --------------------------------------------------------

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

adaptive_restart_fn <- function(opt) {
  adaptive_restart(opt, "fn")
}

adaptive_restart_gr <- function(opt) {
  adaptive_restart(opt, "gr")
}

require_adaptive_restart <- function(opt, par, fn, gr, iter, par0, update) {
  if (!opt$ok
      && (is.null(opt$restart_at)
          || iter - opt$restart_at > opt$restart_wait)) {
    opt <- life_cycle_hook("momentum", "init", opt, par, fn, gr, iter)
    #message("adaptive restart: restarting momentum")
    opt$restart_at <- iter
  }
  else {
    opt$cache$update_old <- update
    #    message("adaptive restart: continuing")
  }
  opt
}
attr(require_adaptive_restart, 'event') <- 'after step'
# Should have the same name as normal update old: we want to replace that hook
attr(require_adaptive_restart, 'name') <- 'update_old'
attr(require_adaptive_restart, 'depends') <- 'update_old_init'
