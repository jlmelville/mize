# Sourceable helpers for adapting optional funconstrain benchmark cases.

funconstrain_x0 <- function(problem, n, name) {
  if (!is.function(problem$x0)) {
    return(problem$x0)
  }
  x0_args <- names(formals(problem$x0))
  if ("n" %in% x0_args) {
    requested <- tryCatch(problem$x0(n = n), error = identity)
    if (!inherits(requested, "error")) {
      return(requested)
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
    return(fallback)
  }
  problem$x0()
}

funconstrain_problem_case <- function(name, n) {
  maker <- getExportedValue("funconstrain", name)
  problem <- maker()
  required <- c("fn", "gr", "x0")
  if (!all(required %in% names(problem))) {
    stop("funconstrain problem ", name, " lacks fn, gr, or x0")
  }

  fg <- list(fn = problem$fn, gr = problem$gr)
  if (is.function(problem$fg)) {
    fg$fg <- problem$fg
  }

  list(
    name = paste0("funconstrain-", name),
    source = "funconstrain",
    par = funconstrain_x0(problem, n, name),
    fg = fg
  )
}

optional_funconstrain_cases <- function(problem_names, n) {
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
    case <- tryCatch(
      funconstrain_problem_case(name, n = n),
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
