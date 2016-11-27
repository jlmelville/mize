# Life Cycle Constructors ------------------------------------------------------

# [before|during|after] [init|step|gradient_descent|momentum|validation] [direction|step_size]

# initialization
# post-initialization
# step
# stages
# validation

# initialization
# post-initialization

# before
# during
# after

# finalization

# step
#   stages
#     stage
#     gradient_descent
#       direction
#       step_size
#     momentum
#       direction
#       step_size
#   validation


make_life_cycle_hook <- function(life_cycle_name) {
  function(opt, par, fn, gr, iter, ...) {
    # Optimizer hook
    if (!is.null(opt[[life_cycle_name]])) {
      opt <- opt[[life_cycle_name]](opt, par, fn, gr, iter, ...)
    }

    # Stage hook
    for (i in 1:length(opt$stages)) {
      stage <- opt$stages[[i]]
      if (!is.null(stage[[life_cycle_name]])) {
        stage_res <- stage[[life_cycle_name]](opt, stage, par, fn, gr, iter,
                                              ...)
        if (!is.null(stage_res$opt)) {
          opt <- stage_res$opt
        }
        if (!is.null(stage_res$stage)) {
          stage <- stage_res$stage
        }
      }

      # Sub stage hooks
      for (sub_stage_name in c("direction", "step_size")) {
        if (!is.null(stage[[sub_stage_name]][[life_cycle_name]])) {
          #message("Calling ", stage$type, "#", sub_stage_name, "#", life_cycle_name)
          sub_stage_res <-
            stage[[sub_stage_name]][[life_cycle_name]](opt, stage, par, fn, gr,
                                                       iter, ...)
          if (!is.null(sub_stage_res$opt)) {
            opt <- sub_stage_res$opt
          }
          if (!is.null(sub_stage_res$stage)) {
            stage <- sub_stage_res$stage
          }
        }
      }

      opt$stages[[i]] <- stage
    }
    opt
  }
}

# Triggers around each stage
make_life_cycle_stage_hook <- function(life_cycle_name) {
  function(opt, stage, par, fn, gr, iter, ...) {
    # Optimizer hook
    # Defining a stage hook for the opt means this function will trigger
    # for every stage
    if (!is.null(opt[[life_cycle_name]])) {
      opt <- opt[[life_cycle_name]](opt, stage, par, fn, gr, iter, ...)
    }

    # Stage hook
    if (!is.null(stage[[life_cycle_name]])) {
      stage_res <- stage[[life_cycle_name]](opt, stage, par, fn, gr, iter,
                                            ...)
      if (!is.null(stage_res$opt)) {
        opt <- stage_res$opt
      }
      if (!is.null(stage_res$stage)) {
        stage <- stage_res$stage
      }
    }

    # Sub stage hooks
    for (sub_stage_name in c("direction", "step_size")) {
      if (!is.null(stage[[sub_stage_name]][[life_cycle_name]])) {
        sub_stage_res <-
          stage[[sub_stage_name]][[life_cycle_name]](opt, stage, par, fn, gr,
                                                     iter, ...)
        if (!is.null(sub_stage_res$opt)) {
          opt <- sub_stage_res$opt
        }
        if (!is.null(sub_stage_res$stage)) {
          stage <- sub_stage_res$stage
        }
      }
    }
    list(opt = opt, stage = stage)
  }
}

life_cycle_before_stage <- make_life_cycle_stage_hook("before_stage")

add_stage_validate <- function(stage, validate, sub_stage = NULL) {
  if (!is.null(sub_stage)) {
    phase <- stage[[sub_stage]]
  }
  else {
    phase <- stage
  }

  if (!is.null(phase$validate)) {
    old_validate <- phase$validate
    new_validate <- function(opt, par, fn, gr, iter, par0, update) {
      old_res <- old_validate(opt, stage, par, fn, gr, iter, par0, update)
      if (!is.null(old_res$opt)) {
        opt <- old_res$opt
      }
      if (!is.null(old_res$stage)) {
        stage <- old_res$stage
      }
      validate(opt, stage, par, fn, gr, iter, par0, update)
    }
  }
  else {
    new_validate <- validate
  }

  if (!is.null(sub_stage)) {
    stage[[sub_stage]]$validate <- new_validate
  }
  else {
    stage$validate <- new_validate
  }
  stage
}

# If sub_stage is NULL, then the function is added to the stage directly
# Otherwise should be one of "gradient_descent" or "momentum"
add_stage_after_step <- function(stage, after_step, sub_stage = NULL) {

  if (!is.null(sub_stage)) {
    phase <- stage[[sub_stage]]
  }
  else {
    phase <- stage
  }

  if (!is.null(phase$after_step)) {
    old_after_step <- phase$after_step

    new_after_step <- function(opt, stage, par, fn, gr, iter, par0, update) {
      old_res <- old_after_step(opt, stage, par, fn, gr, iter, par0, update)
      if (!is.null(old_res$opt)) {
        opt <- old_res$opt
      }
      if (!is.null(old_res$stage)) {
        stage <- old_res$stage
      }
      after_step(opt, stage, par, fn, gr, iter, par0, update)
    }
  }
  else {
    new_after_step <- after_step
  }

  if (!is.null(sub_stage)) {
    stage[[sub_stage]]$after_step <- new_after_step
  }
  else {
    stage$after_step <- new_after_step
  }
  stage
}


# phase:
#   step
#     stage: gradient_descent
#       sub_stage: direction
#       sub_stage: step_size
#     stage: momentum
#     stage? validate
# event:
#   after
#   before
#   during?

# in init function
# create a nested list
# life_cycle:
#   phase
#     event
#       hook_name => function

# Calls all hooks registered with the phase firing this event
life_cycle_hook <- function(phase, advice_type, opt, par, fn, gr, iter, ...) {
#  message("life cycle advice_type: '", advice_type, "' phase: '", phase, "' triggered")
  hooks <- opt$hooks
  if (is.null(hooks)) {
    return(opt)
  }
  hooks <- hooks[[phase]]
  if (is.null(hooks)) {
    return(opt)
  }
  hooks <- hooks[[advice_type]]
  if (is.null(hooks)) {
    return(opt)
  }

  for (name in names(hooks)) {
    hook <- hooks[[name]]
    #message("life_cycle_hook: ", name)
    opt <- hook(opt, par, fn, gr, iter, ...)
  }

  opt
}

# An event is something like "before step"
# Borrowing terminology from Aspect Oriented Programming:
# The first part is the "advice type": before, after.
# The second part is the "join point": the life cycle function that is going
#  to fire.
# Event can consist of the just the join point, e.g. "init" and the advice
#   type is then implicitly assumed to be "during".
register_hook <- function(opt, hook,
                          stage_type = NULL,
                          sub_stage_type = NULL) {
  name <- attr(hook, "name")
  if (is.null(name)) {
    stop("No 'name' attribute for hook function")
  }
  if (!is.null(sub_stage_type)) {
    # Deambiguate the specific sub stage function by adding stage type
    # e.g. could choose bold driver for both gradient descent and momentum
    name <- paste0(stage_type, " ", name)
  }

  event <- attr(hook, "event")
  if (is.null(event)) {
    stop("hook function ", name, " no 'event' attribute")
  }

  event_tok <- strsplit(event, "\\s+")[[1]]
  if (length(event_tok) == 2) {
    advice_type <- event_tok[1]
    join_point <- event_tok[2]
  }
  else {
    advice_type <- "during"
    join_point <- event_tok[1]
  }

  if (!is.null(sub_stage_type)) {
    if (is.null(stage_type)) {
      stop("sub stage type '", sub_stage_type, "' but stage type is NULL")
    }
    if (join_point == "direction" || join_point == "step_size") {
      join_point <- paste0(stage_type, " ", join_point)
    #message("wrapping sub stage hook '", name, "' for event '", join_point,
    #        "' from phase: '", stage_type, " ", sub_stage_type, "'")
    }
    hook <- wrap_sub_stage_hook(hook, stage_type, sub_stage_type)
  }
  else if (!is.null(stage_type)) {
    #message("wrapping stage hook '", name, "' for event '", event,
    #        "' from phase: '", stage_type, "'")
    hook <- wrap_stage_hook(hook, stage_type)
  }

  # store the hook
  if (is.null(opt$hooks[[join_point]])) {
    opt$hooks[[join_point]] <- list()
  }
  join_point_hooks <- opt$hooks[[join_point]]

  if (is.null(join_point_hooks[[advice_type]])) {
    join_point_hooks[[advice_type]] <- list()
  }
  advice <- join_point_hooks[[advice_type]]

#message("registering name: '", name, "' advice_type: '", advice_type,
#          "' join_point: '", join_point, "'")
  advice[[name]] <- hook
  join_point_hooks[[advice_type]] <- advice
  opt$hooks[[join_point]] <- join_point_hooks

  depends <- attr(hook, "depends")
  if (!is.null(depends)) {
    depends <- strsplit(depends, "\\s+")[[1]]
    for (depend in depends) {
      # eventually recursively calls this function
      opt <- depend_to_hook(opt, depend)
    }
  }

  opt
}

wrap_stage_hook <- function(stage_hook, stage_type) {
  callback <- stage_hook
  function(opt, par, fn, gr, iter, ...) {
    stage <- opt$stages[[stage_type]]
    res <- callback(opt, stage, par, fn, gr, iter, ...)

    if (!is.null(res$opt)) {
      opt <- res$opt
    }
    if (!is.null(res$stage)) {
      opt$stages[[stage_type]] <- res$stage
    }
    opt
  }
}

wrap_sub_stage_hook <- function(sub_stage_hook, stage_type, sub_stage_type) {
  callback <- sub_stage_hook
  function(opt, par, fn, gr, iter, ...) {
    stage <- opt$stages[[stage_type]]
    sub_stage <- stage[[sub_stage_type]]
    res <- callback(opt, stage, sub_stage, par, fn, gr, iter, ...)

    if (!is.null(res$opt)) {
      opt <- res$opt
    }

    if (!is.null(res$stage)) {
      stage <- res$stage
    }

    if (!is.null(res$sub_stage)) {
      sub_stage <- res$sub_stage
    }

    stage[[sub_stage_type]] <- sub_stage
    opt$stages[[stage_type]] <- stage

    opt
  }
}

# registers all hook functions:
#   functions in the phase (opt, stage, sub_stage) object itself with attribute
#     'name'. These are intended to be private functions for housekeeping of
#     the phase itself and any calculated values would not be shared. These
#     stage callbacks will be passed the stage, and sub_stage callbacks the
#     stage and the sub_stage so local data can be stored in these objects
#   functions in the depends list. These functions are NOT passed the stage or
#     sub_stage objects, because they are intended to be shared functionality
#     between arbitrary parts of the optimizer tree (e.g. storing the gradient)
register_hooks <- function(opt) {
  # Optimizer hook
  opt$hooks <- list()

  for (name in names(opt)) {
    if (!is.null(attr(opt[[name]], 'event'))) {
      message("registering ", attr(opt[[name]], 'event'))
      opt <- register_hook(opt, opt[[name]])
    }
  }
  if (!is.null(opt$depends)) {
    opt <- depends_to_hooks(opt, opt)
  }


  # Stage hook
  for (i in 1:length(opt$stages)) {
    stage <- opt$stages[[i]]
    #message("Looking for stage hooks for stage ", i, " '", stage$type, "'")
    for (stage_prop_name in names(stage)) {
      if (!is.null(attr(stage[[stage_prop_name]], 'event'))) {
        #message("Found hook '", stage_prop_name, "'")
        opt <- register_hook(opt, stage[[stage_prop_name]], stage$type)
      }
    }

    if (!is.null(stage$depends)) {
      opt <- depends_to_hooks(opt, stage)
    }

    # Sub stage hooks

    for (sub_stage_type in c("direction", "step_size")) {
      #message("Looking for sub stage hooks for stage ", i, " '", stage$type, " ", sub_stage_type, "'")

      sub_stage <- stage[[sub_stage_type]]
      for (sub_stage_prop_name in names(sub_stage)) {
        if (!is.null(attr(sub_stage[[sub_stage_prop_name]], 'event'))) {
          #message("Found hook '", sub_stage_prop_name, "'")

          opt <- register_hook(opt, sub_stage[[sub_stage_prop_name]],
                               stage$type, sub_stage_type)
        }
      }

      if (!is.null(sub_stage$depends)) {
        opt <- depends_to_hooks(opt, sub_stage)
      }
    }
  }
  opt
}


depends_to_hooks <- function(opt, phase, stage_type = NULL,
                             sub_stage_type = NULL) {
  if (is.null(phase$depends)) {
    return(opt)
  }

  for (name in (phase$depends)) {
    opt <- depend_to_hook(opt, name, stage_type, sub_stage_type)
  }

  opt
}

depend_to_hook <- function(opt, depend, stage_type = NULL,
                           sub_stage_type = NULL) {
  f_name <- paste0("require_", depend)
  f <- get0(f_name)
  if (!is.null(f)) {
   # message("registering: ", f_name)
    opt <- register_hook(opt, f, stage_type, sub_stage_type)
  }

  opt
}
