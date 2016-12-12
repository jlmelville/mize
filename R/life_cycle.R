# Life Cycle ------------------------------------------------------
# Various internals for registering functions that should fire during certain
# points of the optimization. You don't want to look to closely at any of this.

# Calls all hooks registered with the phase firing this event
life_cycle_hook <- function(phase, advice_type, opt, par, fg, iter, ...) {
  handler <- life_cycle_handler(phase, advice_type, opt)
  if (is.null(handler)) {
    opt <- default_handler(phase, advice_type, opt, par, fg, iter, ...)
  }
  else {
    opt <- handler(opt, par, fg, iter, ...)
  }
  opt
}

# Look for a handler function that will deal with this event
life_cycle_handler <- function(phase, advice_type, opt) {
  handlers <- opt$handlers
  if (is.null(handlers)) {
    return(NULL)
  }
  handlers <- handlers[[phase]]
  if (is.null(handlers)) {
    return(NULL)
  }
  handlers[[advice_type]]
}

# A default function to handle an event by just iterating over ever function
# that was registered and invoke them in turn
default_handler <- function(phase, advice_type, opt, par, fg, iter, ...) {
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
    opt <- hook(opt, par, fg, iter, ...)
  }
  opt
}

# registers all hook functions:
register_hooks <- function(opt) {
  # Optimizer hook

  for (name in names(opt)) {
    if (!is.null(attr(opt[[name]], 'event'))) {
      opt <- register_hook(opt, opt[[name]])
    }
  }
  if (!is.null(opt$depends)) {
    opt <- depends_to_hooks(opt, opt)
  }

  # Stage hook
  for (i in 1:length(opt$stages)) {
    stage <- opt$stages[[i]]
    for (stage_prop_name in names(stage)) {
      stage_prop <- stage[[stage_prop_name]]
      if (!is.null(attr(stage_prop, 'event'))) {
        opt <- register_hook(opt, stage_prop, stage$type)
      }
    }

    if (!is.null(stage$depends)) {
      opt <- depends_to_hooks(opt, stage)
    }

    # Sub stage hooks
    for (sub_stage_type in c("direction", "step_size")) {
      sub_stage <- stage[[sub_stage_type]]
      for (sub_stage_prop_name in names(sub_stage)) {
        sub_stage_prop <- sub_stage[[sub_stage_prop_name]]
        if (!is.null(attr(sub_stage_prop, 'event'))) {
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
    stop("No 'name' attribute for function")
  }
  if (!is.null(sub_stage_type) && name != "handler") {
    # Deambiguate the specific sub stage function by adding stage type
    # e.g. could choose bold driver for both gradient descent and momentum
    name <- paste0(stage_type, " ", name)
  }

  event <- attr(hook, "event")
  if (is.null(event)) {
    stop("hook function ", name, " no 'event' attribute")
  }

  event_tok <- strsplit(event, "\\s+")[[1]]

  advice_type <- event_tok[1]
  join_point <- event_tok[2]
  if (join_point == "stage"
      || join_point == "gradient_descent"
      || join_point == "momentum"
      && is.null(stage_type)) {
    stage_type <- join_point
  }

  # Functions defined outside of the stage/sub stage constructors can
  # define events completely e.g. "init momentum step_size"
  # but functions defined inside a sub stage don't know the stage
  # that they are defind for and only define e.g. "init step_size"
  # For the free functions, temporarily redefine the join point to just be e.g.
  # "step_size"
  if (length(event_tok) == 3) {
    stage_type <- event_tok[2]
    sub_stage_type <- event_tok[3]
    join_point <- sub_stage_type
  }

  if (!is.null(sub_stage_type)) {
    if (is.null(stage_type)) {
      stop("sub stage type '", sub_stage_type, "' but stage type is NULL")
    }
    if (join_point == "direction" || join_point == "step_size") {
      join_point <- paste0(stage_type, " ", join_point)
    }
    hook <- wrap_sub_stage_hook(hook, stage_type, sub_stage_type)
  }
  else if (!is.null(stage_type)) {
    hook <- wrap_stage_hook(hook, stage_type)
  }

  if (name == "handler") {
    opt <- store_handler(opt, join_point, advice_type, hook)
  }
  else {
    # store the hook
    opt <- store_hook(opt, join_point, advice_type, name, hook)
  }

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

# Puts hook in the correct sub list
store_hook <- function(opt, join_point, advice_type, name, hook) {
  if (is.null(opt$hooks[[join_point]])) {
    opt$hooks[[join_point]] <- list()
  }
  join_point_hooks <- opt$hooks[[join_point]]

  if (is.null(join_point_hooks[[advice_type]])) {
    join_point_hooks[[advice_type]] <- list()
  }
  advice <- join_point_hooks[[advice_type]]

  advice[[name]] <- hook
  join_point_hooks[[advice_type]] <- advice
  opt$hooks[[join_point]] <- join_point_hooks

  opt
}

# Puts the handler in the correct sub list
store_handler <- function(opt, join_point, advice_type, handler) {
  if (is.null(opt$handlers[[join_point]])) {
    opt$handlers[[join_point]] <- list()
  }
  join_point_handlers <- opt$handlers[[join_point]]

  if (is.null(join_point_handlers[[advice_type]])) {
    join_point_handlers[[advice_type]] <- list()
  }

  join_point_handlers[[advice_type]] <- handler
  opt$handlers[[join_point]] <- join_point_handlers

  opt
}

# Wraps a hook that should be fired for a specific stage
# stage_type can be "gradient_descent", "momentum" etc., but also "stage"
# if it should be fired for every stage.
wrap_stage_hook <- function(stage_hook, stage_type) {
  callback <- stage_hook
  function(opt, par, fg, iter, ...) {
    if (stage_type == "stage") {
      stage <- opt$stages[[opt$stage_i]]
    }
    else {
      stage <- opt$stages[[stage_type]]
    }

    res <- callback(opt, stage, par, fg, iter, ...)

    if (!is.null(res$opt)) {
      opt <- res$opt
    }
    if (!is.null(res$stage)) {
      if (stage_type == "stage") {
        opt$stages[[opt$stage_i]] <- res$stage
      }
      else {
        opt$stages[[stage_type]] <- res$stage
      }
    }
    opt
  }
}

# Wraps a hook that should be fired for a specific sub stage.
# stage_type can be "gradient_descent", "momentum" etc., but also "stage"
# if it should be fired for every stage.
# sub_stage should be one of "direction" or "step_size"
wrap_sub_stage_hook <- function(sub_stage_hook, stage_type, sub_stage_type) {
  callback <- sub_stage_hook
  function(opt, par, fg, iter, ...) {
    if (stage_type == "stage") {
      stage <- opt$stages[[opt$stage_i]]
    }
    else {
      stage <- opt$stages[[stage_type]]
    }

    sub_stage <- stage[[sub_stage_type]]
    res <- callback(opt, stage, sub_stage, par, fg, iter, ...)
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

    if (stage_type == "stage") {
      opt$stages[[opt$stage_i]] <- stage
    }
    else {
      opt$stages[[stage_type]] <- stage
    }

    opt
  }
}

# Convert all functions named in the depends vector of a phase into a hook
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

# Convert a specific named function (in depend) into a hook
depend_to_hook <- function(opt, depend, stage_type = NULL,
                           sub_stage_type = NULL) {
  f_name <- paste0("require_", depend)
  f <- get0(f_name)
  if (!is.null(f)) {
    opt <- register_hook(opt, f, stage_type, sub_stage_type)
  }

  opt
}

# Lists all functions and the phases/events they should fire for.
list_hooks <- function(opt) {
  message("handlers")
  if (!is.null(opt$handlers)) {
    handlers <- opt$handlers
    for (phase in names(handlers)) {
      phandlers <- handlers[[phase]]
      for (advice in names(phandlers)) {
        message(advice, " ", phase)
      }
    }
  }

  message("hooks")
  if (!is.null(opt$hooks)) {
    hooks <- opt$hooks
    for (phase in names(hooks)) {
      phooks <- hooks[[phase]]
      for (advice in names(phooks)) {
        aphooks <- phooks[[advice]]
        for (name in names(aphooks)) {
          message(advice, " ", phase, ": ", name)
        }
      }
    }
  }
}
