# Adaptive Restart --------------------------------------------------------

append_depends <- function(opt, stage = NULL, sub_stage = NULL, new_depends) {
  if (!is.null(sub_stage)) {
    if (is.null(stage)) {
      stop("Must provide stage for sub_stage '", sub_stage, "'")
    }
    depends <- opt$stages[[stage]][[sub_stage]]$depends
  }
  else if (!is.null(stage)) {
    depends <- opt$stages[[stage]]$depends
  }
  else {
    depends <- opt$depends
  }
  if (is.null(depends)) {
    depends <- c()
  }

  depends <- c(depends, new_depends)

  if (!is.null(sub_stage)) {
    opt$stages[[stage]][[sub_stage]]$depends <- depends
  }
  else if (!is.null(stage)) {
    opt$stages[[stage]]$depends <- depends
  }
  else {
    opt$depends <- depends
  }

  opt
}

adaptive_restart <- function(opt, validation_type) {
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
  if (!opt$ok) {
    opt <- life_cycle_hook("momentum", "init", opt, par, fn, gr, iter)
#    message("adaptive restart: restarting momentum")
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
