#!/usr/bin/env Rscript

usage <- function(status = 0L) {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/benchmark-optimizers.R [options]",
      "",
      "Options:",
      "  --cases CASES              Comma-separated cases:",
      "                             rosenbrock,spd-quadratic,mmds-eurodist.",
      "  --max-iter N               Maximum iterations per optimizer run.",
      "  --reps N                   Repetitions per optimizer/configuration.",
      "  --seed N                   RNG seed used to build reproducible cases.",
      "  --spd-n N                  Dimension for the SPD quadratic case.",
      "  --line-search VALUES       Comma-separated mize line searches.",
      "  --step0 VALUES             Comma-separated mize first-step settings.",
      "                             Use default for each line search default.",
      "  --funconstrain-cases NAMES Comma-separated funconstrain problem names.",
      "  --no-funconstrain          Skip optional funconstrain cases.",
      "  --out PATH                 Write CSV results to PATH instead of stdout.",
      "  --help                     Show this help.",
      "",
      "Example:",
      "  Rscript scripts/benchmark-optimizers.R --cases rosenbrock \\",
      "    --max-iter 20 --line-search More-Thuente --step0 default",
      sep = "\n"
    )
  )
  quit(save = "no", status = status)
}

split_arg <- function(value) {
  parts <- strsplit(value, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

parse_int <- function(value, name) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  parsed
}

parse_args <- function(args) {
  config <- list(
    cases = c("rosenbrock", "spd-quadratic", "mmds-eurodist"),
    max_iter = 100L,
    reps = 1L,
    seed = 42L,
    spd_n = 50L,
    line_searches = c("More-Thuente", "Hager-Zhang"),
    step0 = c("default", "1"),
    include_funconstrain = TRUE,
    funconstrain_cases = c("rosen", "chebyquad"),
    out = NULL
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    read_value <- function() {
      if (i == length(args)) {
        stop(arg, " requires a value", call. = FALSE)
      }
      args[[i + 1L]]
    }

    if (arg %in% c("--help", "-h")) {
      usage()
    } else if (arg %in% c("--cases", "--case")) {
      config$cases <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--max-iter") {
      config$max_iter <- parse_int(read_value(), "--max-iter")
      i <- i + 1L
    } else if (arg == "--reps") {
      config$reps <- parse_int(read_value(), "--reps")
      i <- i + 1L
    } else if (arg == "--seed") {
      config$seed <- parse_int(read_value(), "--seed")
      i <- i + 1L
    } else if (arg == "--spd-n") {
      config$spd_n <- parse_int(read_value(), "--spd-n")
      i <- i + 1L
    } else if (arg == "--line-search") {
      config$line_searches <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--step0") {
      config$step0 <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--funconstrain-cases") {
      config$funconstrain_cases <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--no-funconstrain") {
      config$include_funconstrain <- FALSE
    } else if (arg == "--out") {
      config$out <- read_value()
      i <- i + 1L
    } else {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  config
}

load_local_mize <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
    return("pkgload::load_all(.)")
  }
  if (requireNamespace("mize", quietly = TRUE)) {
    suppressPackageStartupMessages(library(mize))
    return("installed mize")
  }
  stop(
    "Install pkgload or install mize before running this harness.",
    call. = FALSE
  )
}

norm2 <- function(x) {
  sqrt(sum(x * x))
}

rosenbrock_case <- function() {
  fg <- list(
    fn = function(x) {
      100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
    },
    gr = function(x) {
      c(
        -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
        200 * (x[2] - x[1] * x[1])
      )
    },
    fg = function(x) {
      x1 <- x[1]
      x2 <- x[2]
      a <- x2 - x1 * x1
      b <- 1 - x1
      list(
        fn = 100 * a * a + b * b,
        gr = c(-400 * x1 * a - 2 * b, 200 * a)
      )
    }
  )

  list(
    name = "rosenbrock",
    source = "built-in",
    par = c(-1.2, 1),
    fg = fg
  )
}

rposdef <- function(n, ev = stats::runif(n, 0, 10)) {
  z <- matrix(stats::rnorm(n * n), ncol = n)
  decomp <- qr(z)
  q <- qr.Q(decomp)
  r <- qr.R(decomp)
  ph <- diag(r) / abs(diag(r))
  o <- q %*% diag(ph)
  t(o) %*% diag(ev) %*% o
}

quad_matrix_fg <- function(a, b) {
  minimizer <- solve(a, b)
  f_min <- drop(0.5 * t(minimizer) %*% a %*% minimizer - t(b) %*% minimizer)

  list(
    fn = function(par) {
      drop(0.5 * t(par) %*% a %*% par - t(b) %*% par) - f_min
    },
    gr = function(par) {
      as.vector(0.5 * (t(a) %*% par + a %*% par) - b)
    },
    fg = function(par) {
      list(
        fn = drop(0.5 * t(par) %*% a %*% par - t(b) %*% par) - f_min,
        gr = as.vector(0.5 * (t(a) %*% par + a %*% par) - b)
      )
    },
    hs = function(par) {
      a
    },
    minimizer = minimizer
  )
}

spd_quadratic_case <- function(n = 50L) {
  eig <- stats::runif(n = n, min = 0.001, max = 1)
  a <- rposdef(n = n, ev = eig)
  b <- stats::rnorm(n = n, mean = 0, sd = sqrt(25))

  list(
    name = paste0("spd-quadratic-n", n),
    source = "built-in",
    par = rep(0, n),
    fg = quad_matrix_fg(a, b)
  )
}

make_mmds_fg <- function(distmat) {
  r <- as.matrix(distmat)
  cost_fun <- function(ref, d) {
    sum((ref - d)^2) * 0.5
  }
  cost_grad <- function(ref, d, y) {
    k <- (ref - d) / (d + 1.e-10)
    g <- matrix(nrow = nrow(y), ncol = ncol(y))

    for (i in seq_len(nrow(y))) {
      dyij <- sweep(-y, 2, -y[i, ])
      g[i, ] <- colSums(dyij * k[, i])
    }

    as.vector(t(g)) * -2
  }

  list(
    fn = function(par) {
      y <- matrix(par, ncol = 2, byrow = TRUE)
      d <- as.matrix(stats::dist(y))
      cost_fun(r, d)
    },
    gr = function(par) {
      y <- matrix(par, ncol = 2, byrow = TRUE)
      d <- as.matrix(stats::dist(y))
      cost_grad(r, d, y)
    },
    fg = function(par) {
      y <- matrix(par, ncol = 2, byrow = TRUE)
      d <- as.matrix(stats::dist(y))
      list(
        fn = cost_fun(r, d),
        gr = cost_grad(r, d, y)
      )
    }
  )
}

mmds_case <- function() {
  list(
    name = "mmds-eurodist",
    source = "built-in",
    par = stats::rnorm(attr(datasets::eurodist, "Size") * 2),
    fg = make_mmds_fg(datasets::eurodist)
  )
}

funconstrain_x0 <- function(problem, n) {
  if (!is.function(problem$x0)) {
    return(problem$x0)
  }
  x0_args <- names(formals(problem$x0))
  if ("n" %in% x0_args) {
    return(problem$x0(n = n))
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
    par = funconstrain_x0(problem, n),
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

benchmark_cases <- function(config) {
  allowed <- c("rosenbrock", "spd-quadratic", "mmds-eurodist")
  unknown <- setdiff(config$cases, allowed)
  if (length(unknown) > 0L) {
    stop("Unknown built-in case(s): ", paste(unknown, collapse = ", "))
  }

  set.seed(config$seed)
  cases <- list()
  if ("rosenbrock" %in% config$cases) {
    cases$rosenbrock <- rosenbrock_case()
  }
  if ("spd-quadratic" %in% config$cases) {
    case <- spd_quadratic_case(n = config$spd_n)
    cases[[case$name]] <- case
  }
  if ("mmds-eurodist" %in% config$cases) {
    cases$`mmds-eurodist` <- mmds_case()
  }
  if (config$include_funconstrain) {
    cases <- c(
      cases,
      optional_funconstrain_cases(
        problem_names = config$funconstrain_cases,
        n = config$spd_n
      )
    )
  }

  cases
}

counted_fg <- function(fg) {
  counts <- new.env(parent = emptyenv())
  counts$nf <- 0L
  counts$ng <- 0L

  list(
    fn = function(par) {
      counts$nf <- counts$nf + 1L
      fg$fn(par)
    },
    gr = function(par) {
      counts$ng <- counts$ng + 1L
      fg$gr(par)
    },
    counts = function() {
      list(nf = counts$nf, ng = counts$ng)
    }
  )
}

run_timed <- function(thunk) {
  warnings <- character()
  result <- NULL
  elapsed <- system.time({
    result <- tryCatch(
      withCallingHandlers(
        thunk(),
        warning = function(warn) {
          warnings <<- c(warnings, conditionMessage(warn))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(err) err
    )
  })[["elapsed"]]

  list(
    result = result,
    elapsed = unname(elapsed),
    warnings = unique(warnings)
  )
}

safe_final_metrics <- function(fg, par) {
  if (inherits(par, "error") || is.null(par)) {
    return(list(final_f = NA_real_, grad_norm = NA_real_))
  }

  tryCatch(
    {
      final_f <- fg$fn(par)
      grad_norm <- norm2(fg$gr(par))
      list(final_f = final_f, grad_norm = grad_norm)
    },
    error = function(err) {
      list(final_f = NA_real_, grad_norm = NA_real_)
    }
  )
}

row_result <- function(
  case,
  rep,
  optimizer,
  method,
  line_search,
  step0,
  max_iter,
  final_f,
  grad_norm,
  nf,
  ng,
  elapsed,
  iter,
  termination,
  failure_mode,
  warnings,
  error
) {
  data.frame(
    case = case$name,
    source = case$source,
    rep = rep,
    optimizer = optimizer,
    method = method,
    line_search = line_search,
    step0 = step0,
    max_iter = max_iter,
    final_f = final_f,
    grad_norm = grad_norm,
    nf = nf,
    ng = ng,
    elapsed_sec = elapsed,
    iter = iter,
    termination = termination,
    failure_mode = failure_mode,
    warnings = paste(warnings, collapse = " | "),
    error = error,
    stringsAsFactors = FALSE
  )
}

optim_failure <- function(res) {
  if (inherits(res, "error")) {
    return("error")
  }
  if (isTRUE(res$convergence == 0)) {
    return("ok")
  }
  paste0("optim_convergence_", res$convergence)
}

run_optim_case <- function(case, rep, method, max_iter) {
  counted <- counted_fg(case$fg)
  timed <- run_timed(function() {
    stats::optim(
      par = case$par,
      fn = counted$fn,
      gr = counted$gr,
      method = method,
      control = list(maxit = max_iter)
    )
  })
  counts <- counted$counts()
  res <- timed$result
  error <- if (inherits(res, "error")) conditionMessage(res) else ""
  par <- if (inherits(res, "error")) NULL else res$par
  metrics <- safe_final_metrics(case$fg, par)
  termination <- if (inherits(res, "error")) {
    "error"
  } else {
    paste0("convergence=", res$convergence)
  }

  row_result(
    case = case,
    rep = rep,
    optimizer = "stats::optim",
    method = method,
    line_search = NA_character_,
    step0 = NA_character_,
    max_iter = max_iter,
    final_f = metrics$final_f,
    grad_norm = metrics$grad_norm,
    nf = counts$nf,
    ng = counts$ng,
    elapsed = timed$elapsed,
    iter = NA_integer_,
    termination = termination,
    failure_mode = optim_failure(res),
    warnings = timed$warnings,
    error = error
  )
}

mize_step0_value <- function(step0) {
  if (identical(tolower(step0), "default")) {
    return(NULL)
  }
  numeric_value <- suppressWarnings(as.numeric(step0))
  if (!is.na(numeric_value)) {
    return(numeric_value)
  }
  step0
}

mize_failure <- function(res) {
  if (inherits(res, "error")) {
    return("error")
  }

  what <- res$terminate$what
  if (what %in% c("abs_tol", "rel_tol", "grad_tol", "ginf_tol", "step_tol")) {
    return("ok")
  }
  if (what %in% c("max_iter", "max_fn", "max_gr", "max_fg")) {
    return("budget")
  }
  if (what %in% c("fn_inf", "gr_inf")) {
    return("nonfinite")
  }
  what
}

run_mize_case <- function(
  case,
  rep,
  method,
  cg_update,
  line_search,
  step0,
  max_iter
) {
  step0_value <- mize_step0_value(step0)
  timed <- run_timed(function() {
    mize(
      par = case$par,
      fg = case$fg,
      method = method,
      cg_update = cg_update,
      line_search = line_search,
      step0 = step0_value,
      max_iter = max_iter
    )
  })
  res <- timed$result
  error <- if (inherits(res, "error")) conditionMessage(res) else ""
  par <- if (inherits(res, "error")) NULL else res$par
  metrics <- safe_final_metrics(case$fg, par)
  termination <- if (inherits(res, "error")) {
    "error"
  } else {
    paste0(res$terminate$what, "=", res$terminate$val)
  }
  method_label <- if (identical(method, "CG")) {
    paste0(method, ":", cg_update)
  } else {
    method
  }

  row_result(
    case = case,
    rep = rep,
    optimizer = "mize",
    method = method_label,
    line_search = line_search,
    step0 = step0,
    max_iter = max_iter,
    final_f = metrics$final_f,
    grad_norm = metrics$grad_norm,
    nf = if (inherits(res, "error")) NA_integer_ else res$nf,
    ng = if (inherits(res, "error")) NA_integer_ else res$ng,
    elapsed = timed$elapsed,
    iter = if (inherits(res, "error")) NA_integer_ else res$iter,
    termination = termination,
    failure_mode = mize_failure(res),
    warnings = timed$warnings,
    error = error
  )
}

run_case_grid <- function(case, config) {
  optim_methods <- c("BFGS", "CG", "L-BFGS-B")
  mize_methods <- list(
    list(method = "L-BFGS", cg_update = "PR+"),
    list(method = "BFGS", cg_update = "PR+"),
    list(method = "CG", cg_update = "PR+"),
    list(method = "CG", cg_update = "HZ+"),
    list(method = "TN", cg_update = "PR+")
  )

  rows <- list()
  for (rep in seq_len(config$reps)) {
    for (method in optim_methods) {
      rows[[length(rows) + 1L]] <- run_optim_case(
        case = case,
        rep = rep,
        method = method,
        max_iter = config$max_iter
      )
    }
    for (method_config in mize_methods) {
      for (line_search in config$line_searches) {
        for (step0 in config$step0) {
          rows[[length(rows) + 1L]] <- run_mize_case(
            case = case,
            rep = rep,
            method = method_config$method,
            cg_update = method_config$cg_update,
            line_search = line_search,
            step0 = step0,
            max_iter = config$max_iter
          )
        }
      }
    }
  }

  do.call(rbind, rows)
}

run_benchmarks <- function(config) {
  cases <- benchmark_cases(config)
  if (length(cases) == 0L) {
    stop("No benchmark cases selected", call. = FALSE)
  }

  rows <- lapply(cases, run_case_grid, config = config)
  do.call(rbind, rows)
}

write_results <- function(results, out) {
  if (is.null(out)) {
    utils::write.csv(results, file = stdout(), row.names = FALSE, na = "")
    return(invisible(NULL))
  }

  utils::write.csv(results, file = out, row.names = FALSE, na = "")
  message("Wrote benchmark results to ", out)
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  loader <- load_local_mize()
  message("Loaded mize via ", loader)
  results <- run_benchmarks(config)
  write_results(results, config$out)
}

main()
