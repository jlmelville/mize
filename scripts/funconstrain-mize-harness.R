# Sourceable helpers for adapting optional funconstrain benchmark cases.

resolve_funconstrain_start <- function(problem, n, name) {
  if (!is.function(problem$x0)) {
    par <- problem$x0
    return(list(
      par = par,
      source = "fixed",
      resolution = "fixed",
      requested_dimension = n,
      requested_dimension_accepted = NA,
      actual_dimension = length(par)
    ))
  }

  x0_args <- names(formals(problem$x0))
  if ("n" %in% x0_args) {
    requested <- tryCatch(problem$x0(n = n), error = identity)
    if (!inherits(requested, "error")) {
      return(list(
        par = requested,
        source = "function",
        resolution = "requested_dimension",
        requested_dimension = n,
        requested_dimension_accepted = TRUE,
        actual_dimension = length(requested)
      ))
    }

    fallback <- tryCatch(problem$x0(), error = identity)
    if (inherits(fallback, "error")) {
      stop(requested)
    }
    message(
      "funconstrain problem ",
      name,
      " rejected n = ",
      n,
      "; using its default dimension ",
      length(fallback),
      "."
    )
    return(list(
      par = fallback,
      source = "function",
      resolution = "default_after_rejected_request",
      requested_dimension = n,
      requested_dimension_accepted = FALSE,
      actual_dimension = length(fallback)
    ))
  }

  par <- problem$x0()
  list(
    par = par,
    source = "function",
    resolution = "function_default",
    requested_dimension = n,
    requested_dimension_accepted = NA,
    actual_dimension = length(par)
  )
}

funconstrain_x0 <- function(problem, n, name) {
  resolve_funconstrain_start(problem, n, name)$par
}

funconstrain_factory_arguments <- function(maker, supplied_arguments) {
  if (!is.list(supplied_arguments)) {
    stop("factory arguments must be supplied as a list", call. = FALSE)
  }
  if (
    length(supplied_arguments) > 0L &&
      (is.null(names(supplied_arguments)) ||
        any(!nzchar(names(supplied_arguments))))
  ) {
    stop("factory arguments must be named", call. = FALSE)
  }

  factory_formals <- formals(maker)
  unknown <- setdiff(names(supplied_arguments), names(factory_formals))
  if (length(unknown) > 0L) {
    stop(
      "unknown factory argument(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  effective_arguments <- list()
  evaluation_environment <- new.env(parent = environment(maker))
  for (argument_name in names(factory_formals)) {
    if (argument_name %in% names(supplied_arguments)) {
      value <- supplied_arguments[[argument_name]]
    } else {
      default <- factory_formals[[argument_name]]
      if (identical(default, quote(expr = ))) {
        next
      }
      value <- eval(default, envir = evaluation_environment)
    }
    effective_arguments[argument_name] <- list(value)
    assign(argument_name, value, envir = evaluation_environment)
  }
  effective_arguments
}

funconstrain_provenance <- function(source_commit = NULL) {
  description <- utils::packageDescription("funconstrain")
  installed_commit <- NULL
  for (field in c("RemoteSha", "GithubSHA1")) {
    candidate <- description[[field]]
    if (
      !is.null(candidate) &&
        length(candidate) == 1L &&
        !is.na(candidate) &&
        nzchar(candidate)
    ) {
      installed_commit <- candidate
      break
    }
  }

  provided_commit <- !is.null(source_commit) &&
    length(source_commit) == 1L &&
    !is.na(source_commit) &&
    nzchar(source_commit)
  if (provided_commit) {
    commit <- source_commit
    commit_source <- "provided"
  } else if (!is.null(installed_commit)) {
    commit <- installed_commit
    commit_source <- "installed_package_metadata"
  } else {
    commit <- NA_character_
    commit_source <- "unavailable"
  }

  list(
    package = "funconstrain",
    version = as.character(description[["Version"]]),
    source_commit = commit,
    source_commit_source = commit_source
  )
}

funconstrain_reference_metadata <- function(name, problem, dimension, maker) {
  fmin <- problem[["fmin"]]
  xmin <- problem[["xmin"]]
  fmin_available <- !is.null(fmin) && length(fmin) > 0L && !all(is.na(fmin))
  xmin_available <- !is.null(xmin) && length(xmin) > 0L
  xmin_usable <- xmin_available && is.numeric(xmin) && !anyNA(xmin)

  fmin_applicable <- if (fmin_available) NA else FALSE
  xmin_applicable <- if (xmin_available) NA else FALSE
  fmin_basis <- if (fmin_available) "not_encoded" else "unavailable"
  xmin_basis <- if (xmin_available) "not_encoded" else "unavailable"

  if (xmin_available && !xmin_usable) {
    xmin_applicable <- FALSE
    xmin_basis <- "no_single_minimizer"
  } else if (xmin_usable && length(xmin) != dimension) {
    xmin_applicable <- FALSE
    xmin_basis <- "dimension_mismatch"
  }

  fixed_configuration <- !is.function(problem$x0) &&
    length(formals(maker)) == 0L
  if (fixed_configuration) {
    if (fmin_available) {
      fmin_applicable <- TRUE
      fmin_basis <- "sole_fixed_configuration"
    }
    if (xmin_usable && length(xmin) == dimension) {
      xmin_applicable <- TRUE
      xmin_basis <- "sole_fixed_configuration"
    }
  } else if (name %in% c("ex_rosen", "var_dim")) {
    if (fmin_available) {
      fmin_applicable <- TRUE
      fmin_basis <- "documented_all_valid_dimensions"
    }
    reference_dimension <- if (name == "ex_rosen") 8L else 6L
    if (xmin_usable && length(xmin) == dimension) {
      xmin_applicable <- dimension == reference_dimension
      xmin_basis <- paste0("documented_dimension_", reference_dimension)
    }
  } else if (name == "chebyquad") {
    if (fmin_available) {
      fmin_applicable <- dimension == 8L
      fmin_basis <- "documented_dimension_8"
    }
    if (xmin_usable && length(xmin) == dimension) {
      xmin_applicable <- dimension == 8L
      xmin_basis <- "documented_dimension_8"
    }
  }

  list(
    fmin = fmin,
    xmin = xmin,
    fmin_applicable = fmin_applicable,
    xmin_applicable = xmin_applicable,
    fmin_basis = fmin_basis,
    xmin_basis = xmin_basis
  )
}

validate_funconstrain_case <- function(case) {
  case_label <- paste0("funconstrain case ", case$name)
  fail <- function(...) {
    stop(case_label, " ", ..., call. = FALSE)
  }

  if (!is.list(case$fg)) {
    fail("must provide callbacks in a list")
  }

  required_callbacks <- c("fn", "gr", "hs")
  missing_callbacks <- setdiff(required_callbacks, names(case$fg))
  if (length(missing_callbacks) > 0L) {
    fail(
      "is missing required callback(s): ",
      paste0("fg$", missing_callbacks, collapse = ", ")
    )
  }

  invalid_callbacks <- required_callbacks[
    !vapply(case$fg[required_callbacks], is.function, logical(1))
  ]
  if (length(invalid_callbacks) > 0L) {
    fail(
      "requires function callback(s): ",
      paste0("fg$", invalid_callbacks, collapse = ", ")
    )
  }

  par <- case$par
  if (!is.numeric(par) || !is.null(dim(par)) || length(par) == 0L) {
    fail("starting parameters must be a non-empty numeric vector")
  }
  dimension <- length(par)

  evaluate_callback <- function(callback_name) {
    tryCatch(
      case$fg[[callback_name]](par),
      error = function(err) {
        fail(
          "callback fg$",
          callback_name,
          " failed at the resolved start: ",
          conditionMessage(err)
        )
      }
    )
  }

  objective <- evaluate_callback("fn")
  valid_objective <- is.numeric(objective) &&
    is.null(dim(objective)) &&
    length(objective) == 1L &&
    is.finite(objective)
  if (!valid_objective) {
    fail(
      "callback fg$fn must return a finite numeric scalar at the resolved start"
    )
  }

  gradient <- evaluate_callback("gr")
  valid_gradient <- is.numeric(gradient) &&
    is.null(dim(gradient)) &&
    length(gradient) == dimension
  if (!valid_gradient) {
    fail(
      "callback fg$gr must return a numeric vector of length ",
      dimension,
      " at the resolved start"
    )
  }

  hessian <- evaluate_callback("hs")
  if (!is.matrix(hessian) || !is.numeric(hessian)) {
    fail("callback fg$hs must return a numeric matrix at the resolved start")
  }
  if (!identical(dim(hessian), c(dimension, dimension))) {
    fail(
      "callback fg$hs must return a ",
      dimension,
      " x ",
      dimension,
      " matrix at the resolved start"
    )
  }
  if (!all(is.finite(hessian))) {
    fail("callback fg$hs must return only finite values at the resolved start")
  }

  hessian_scale <- max(1, max(abs(hessian)))
  relative_asymmetry <- max(abs(hessian - t(hessian))) / hessian_scale
  symmetry_tolerance <- sqrt(.Machine$double.eps)
  if (relative_asymmetry > symmetry_tolerance) {
    fail(
      "callback fg$hs must return an approximately symmetric matrix at the resolved start"
    )
  }

  invisible(case)
}

counted_fg <- function(fg) {
  counts <- new.env(parent = emptyenv())
  counts$n_fn_call <- 0L
  counts$n_gr_call <- 0L
  counts$n_fg_call <- 0L
  counts$n_hs_call <- 0L

  wrap_callback <- function(callback_name, count_name) {
    callback <- fg[[callback_name]]
    if (!is.function(callback)) {
      return(callback)
    }

    function(...) {
      counts[[count_name]] <- counts[[count_name]] + 1L
      callback(...)
    }
  }

  counted <- fg
  counted$fn <- wrap_callback("fn", "n_fn_call")
  counted$gr <- wrap_callback("gr", "n_gr_call")
  counted$fg <- wrap_callback("fg", "n_fg_call")
  counted$hs <- wrap_callback("hs", "n_hs_call")

  list(
    fg = counted,
    counts = function() {
      list(
        n_fn_call = counts$n_fn_call,
        n_gr_call = counts$n_gr_call,
        n_fg_call = counts$n_fg_call,
        n_hs_call = counts$n_hs_call
      )
    }
  )
}

funconstrain_problem_case <- function(
  name,
  n,
  factory_arguments = list(),
  source_commit = NULL
) {
  maker <- getExportedValue("funconstrain", name)
  effective_arguments <- funconstrain_factory_arguments(
    maker,
    factory_arguments
  )
  problem <- do.call(maker, factory_arguments)
  required <- c("fn", "gr", "x0")
  if (!all(required %in% names(problem))) {
    stop("funconstrain problem ", name, " lacks fn, gr, or x0")
  }

  start <- resolve_funconstrain_start(problem, n, name)
  fg <- list(fn = problem$fn, gr = problem$gr)
  if (is.function(problem$fg)) {
    fg$fg <- problem$fg
  }
  if (is.function(problem$he)) {
    fg$hs <- problem$he
  }

  case <- list(
    name = paste0("funconstrain-", name),
    source = "funconstrain",
    par = start$par,
    fg = fg,
    start = start,
    factory = list(
      name = name,
      arguments = effective_arguments,
      supplied_arguments = factory_arguments
    ),
    provenance = funconstrain_provenance(source_commit),
    reference = funconstrain_reference_metadata(
      name = name,
      problem = problem,
      dimension = start$actual_dimension,
      maker = maker
    )
  )
  validate_funconstrain_case(case)
  case
}

optional_funconstrain_cases <- function(
  problem_names,
  n,
  factory_arguments_by_problem = list(),
  source_commit = NULL
) {
  if (!requireNamespace("funconstrain", quietly = TRUE)) {
    message("funconstrain is not installed; skipping funconstrain cases.")
    return(list())
  }

  exports <- getNamespaceExports("funconstrain")
  cases <- list()
  for (name in problem_names) {
    if (!name %in% exports) {
      message("funconstrain does not export ", name, "; skipping.")
      next
    }
    problem_arguments <- factory_arguments_by_problem[[name]]
    if (is.null(problem_arguments)) {
      problem_arguments <- list()
    }
    case <- tryCatch(
      funconstrain_problem_case(
        name,
        n = n,
        factory_arguments = problem_arguments,
        source_commit = source_commit
      ),
      error = function(err) {
        message(
          "Could not adapt funconstrain problem ",
          name,
          ": ",
          conditionMessage(err)
        )
        NULL
      }
    )
    if (!is.null(case)) {
      cases[[case$name]] <- case
    }
  }
  cases
}
