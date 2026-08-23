#!/usr/bin/env Rscript

source(file.path("scripts", "funconstrain-mize-harness.R"))

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
      "  --warmup N                 Unreported warmups per configuration.",
      "  --seed N                   RNG seed used to build reproducible cases.",
      "  --spd-n N                  Dimension for the SPD quadratic case.",
      "  --methods VALUES           Comma-separated optimizer method profiles.",
      "                             stats-bfgs,stats-cg,stats-l-bfgs-b,",
      "                             mize-l-bfgs,mize-bfgs,mize-cg-pr+,",
      "                             mize-cg-hz+,mize-tn,mize-newton,",
      "                             mize-sr1,mize-tn-l-bfgs,",
      "                             mize-newton-safeguarded, and",
      "                             mize-cg-pr+-l-bfgs.",
      "                             The default retains the legacy grid.",
      "  --line-search VALUES       Comma-separated mize line searches.",
      "  --step0 VALUES             Comma-separated mize first-step settings.",
      "                             Use default for each line search default.",
      "  --callbacks VALUES         Comma-separated mize callback modes:",
      "                             combined,separate (default: combined).",
      "  --funconstrain-cases NAMES Comma-separated funconstrain problem names.",
      "  --funconstrain-commit SHA  Optional source commit provenance.",
      "  --no-funconstrain          Skip optional funconstrain cases.",
      "  --out PATH                 Write CSV results to PATH instead of stdout.",
      "  --progress-out PATH        Write opt-in per-iteration progress CSV.",
      "  --case-out PATH            Write opt-in case/target manifest CSV.",
      "  --summary-out PATH         Write opt-in resource summary CSV.",
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

benchmark_method_profiles <- function() {
  list(
    `stats-bfgs` = list(optimizer = "stats::optim", method = "BFGS"),
    `stats-cg` = list(optimizer = "stats::optim", method = "CG"),
    `stats-l-bfgs-b` = list(
      optimizer = "stats::optim",
      method = "L-BFGS-B"
    ),
    `mize-l-bfgs` = list(
      optimizer = "mize",
      method = "L-BFGS",
      cg_update = "PR+"
    ),
    `mize-bfgs` = list(
      optimizer = "mize",
      method = "BFGS",
      cg_update = "PR+"
    ),
    `mize-cg-pr+` = list(
      optimizer = "mize",
      method = "CG",
      cg_update = "PR+"
    ),
    `mize-cg-hz+` = list(
      optimizer = "mize",
      method = "CG",
      cg_update = "HZ+"
    ),
    `mize-tn` = list(
      optimizer = "mize",
      method = "TN",
      cg_update = "PR+",
      preconditioner = ""
    ),
    `mize-newton` = list(
      optimizer = "mize",
      method = "NEWTON",
      cg_update = "PR+"
    ),
    `mize-sr1` = list(
      optimizer = "mize",
      method = "SR1",
      cg_update = "PR+"
    ),
    `mize-tn-l-bfgs` = list(
      optimizer = "mize",
      method = "TN",
      cg_update = "PR+",
      preconditioner = "L-BFGS"
    ),
    `mize-newton-safeguarded` = list(
      optimizer = "mize",
      method = "NEWTON",
      cg_update = "PR+",
      safeguarded_newton = TRUE
    ),
    `mize-cg-pr+-l-bfgs` = list(
      optimizer = "mize",
      method = "CG",
      cg_update = "PR+",
      preconditioner = "L-BFGS"
    )
  )
}

benchmark_default_method_profiles <- function() {
  c(
    "stats-bfgs",
    "stats-cg",
    "stats-l-bfgs-b",
    "mize-l-bfgs",
    "mize-bfgs",
    "mize-cg-pr+",
    "mize-cg-hz+",
    "mize-tn"
  )
}

parse_int <- function(value, name) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  parsed
}

parse_nonnegative_int <- function(value, name) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 0L) {
    stop(name, " must be a non-negative integer", call. = FALSE)
  }
  parsed
}

parse_args <- function(args) {
  config <- list(
    cases = c("rosenbrock", "spd-quadratic", "mmds-eurodist"),
    max_iter = 100L,
    reps = 1L,
    warmup = 0L,
    seed = 42L,
    spd_n = 50L,
    methods = benchmark_default_method_profiles(),
    line_searches = c("More-Thuente", "Hager-Zhang"),
    step0 = c("default", "1"),
    callback_modes = "combined",
    include_funconstrain = TRUE,
    funconstrain_cases = c("rosen", "chebyquad"),
    funconstrain_commit = NULL,
    out = NULL,
    progress_out = NULL,
    case_out = NULL,
    summary_out = NULL
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
    } else if (arg == "--warmup") {
      config$warmup <- parse_nonnegative_int(read_value(), "--warmup")
      i <- i + 1L
    } else if (arg == "--seed") {
      config$seed <- parse_int(read_value(), "--seed")
      i <- i + 1L
    } else if (arg == "--spd-n") {
      config$spd_n <- parse_int(read_value(), "--spd-n")
      i <- i + 1L
    } else if (arg == "--methods") {
      config$methods <- tolower(split_arg(read_value()))
      i <- i + 1L
    } else if (arg == "--line-search") {
      config$line_searches <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--step0") {
      config$step0 <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--callbacks") {
      config$callback_modes <- tolower(split_arg(read_value()))
      i <- i + 1L
    } else if (arg == "--funconstrain-cases") {
      config$funconstrain_cases <- split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--funconstrain-commit") {
      config$funconstrain_commit <- read_value()
      i <- i + 1L
    } else if (arg == "--no-funconstrain") {
      config$include_funconstrain <- FALSE
    } else if (arg == "--out") {
      config$out <- read_value()
      i <- i + 1L
    } else if (arg == "--progress-out") {
      config$progress_out <- read_value()
      i <- i + 1L
    } else if (arg == "--case-out") {
      config$case_out <- read_value()
      i <- i + 1L
    } else if (arg == "--summary-out") {
      config$summary_out <- read_value()
      i <- i + 1L
    } else {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  if (length(config$methods) == 0L) {
    stop("--methods must select at least one method profile", call. = FALSE)
  }
  unknown_methods <- setdiff(
    config$methods,
    names(benchmark_method_profiles())
  )
  if (length(unknown_methods) > 0L) {
    stop(
      "Unknown method profile(s): ",
      paste(unknown_methods, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_callback_modes <- setdiff(
    config$callback_modes,
    c("combined", "separate")
  )
  if (length(unknown_callback_modes) > 0L) {
    stop(
      "Unknown callback mode(s): ",
      paste(unknown_callback_modes, collapse = ", "),
      call. = FALSE
    )
  }

  benchmark_preflight_output_paths(list(
    `--out` = config$out,
    `--progress-out` = config$progress_out,
    `--case-out` = config$case_out,
    `--summary-out` = config$summary_out
  ))

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

  fg <- quad_matrix_fg(a, b)
  list(
    name = paste0("spd-quadratic-n", n),
    source = "built-in",
    par = rep(0, n),
    fg = fg,
    start = list(
      resolution = "fixed",
      requested_dimension = n,
      requested_dimension_accepted = TRUE,
      actual_dimension = n
    ),
    reference = list(
      fmin = 0,
      xmin = fg$minimizer,
      fmin_applicable = TRUE,
      xmin_applicable = TRUE,
      fmin_basis = "exact_spd_minimum",
      xmin_basis = "exact_spd_minimizer"
    )
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
        n = config$spd_n,
        source_commit = config$funconstrain_commit
      )
    )
  }

  lapply(cases, function(case) {
    case$dimension <- length(case$par)
    case$initial_f <- case$fg$fn(case$par)
    case$initial_grad_norm <- norm2(case$fg$gr(case$par))
    case
  })
}

benchmark_scalar <- function(value, default = NA) {
  if (is.null(value) || length(value) != 1L || is.na(value)) {
    return(default)
  }
  value
}

benchmark_applicability <- function(value) {
  if (is.null(value) || length(value) != 1L || !is.logical(value)) {
    return(NA)
  }
  value
}

benchmark_quality_targets <- function(case) {
  reference <- case$reference
  if (is.null(reference)) {
    reference <- list()
  }
  fmin <- benchmark_scalar(reference$fmin, NA_real_)
  fmin_applicable <- identical(reference$fmin_applicable, TRUE) &&
    is.numeric(fmin) &&
    is.finite(fmin)

  list(
    objective_reduction_target = case$initial_f -
      0.99 * max(abs(case$initial_f), 1e-12),
    gradient_reduction_target = max(
      1e-8,
      1e-6 * case$initial_grad_norm
    ),
    reference_gap_target = if (fmin_applicable) {
      fmin + 1e-8 * max(1, abs(case$initial_f - fmin))
    } else {
      NA_real_
    },
    fmin = fmin,
    fmin_applicable = fmin_applicable
  )
}

benchmark_case_manifest <- function(cases) {
  rows <- lapply(cases, function(case) {
    start <- case$start
    if (is.null(start)) {
      start <- list(
        resolution = "fixed",
        requested_dimension = case$dimension,
        requested_dimension_accepted = NA,
        actual_dimension = case$dimension
      )
    }
    reference <- case$reference
    if (is.null(reference)) {
      reference <- list()
    }
    provenance <- case$provenance
    if (is.null(provenance)) {
      provenance <- list()
    }
    targets <- benchmark_quality_targets(case)

    data.frame(
      case = case$name,
      source = case$source,
      dimension = case$dimension,
      start_resolution = start$resolution,
      requested_dimension = start$requested_dimension,
      requested_dimension_accepted = start$requested_dimension_accepted,
      actual_dimension = start$actual_dimension,
      initial_f = case$initial_f,
      initial_grad_norm = case$initial_grad_norm,
      objective_reduction_target = targets$objective_reduction_target,
      gradient_reduction_target = targets$gradient_reduction_target,
      fmin = targets$fmin,
      fmin_applicable = benchmark_applicability(
        reference$fmin_applicable
      ),
      fmin_target_applicable = targets$fmin_applicable,
      fmin_basis = benchmark_scalar(
        reference$fmin_basis,
        "unavailable"
      ),
      reference_gap_target = targets$reference_gap_target,
      xmin_applicable = benchmark_applicability(
        reference$xmin_applicable
      ),
      xmin_basis = benchmark_scalar(
        reference$xmin_basis,
        "unavailable"
      ),
      package = benchmark_scalar(provenance$package, NA_character_),
      package_version = benchmark_scalar(
        provenance$version,
        NA_character_
      ),
      source_commit = benchmark_scalar(
        provenance$source_commit,
        NA_character_
      ),
      source_commit_source = benchmark_scalar(
        provenance$source_commit_source,
        NA_character_
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
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
  callback_mode,
  max_iter,
  final_f,
  grad_norm,
  nf,
  ng,
  n_fn_call,
  n_gr_call,
  n_fg_call,
  n_hs_call,
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
    callback_mode = callback_mode,
    dimension = case$dimension,
    initial_f = case$initial_f,
    max_iter = max_iter,
    final_f = final_f,
    grad_norm = grad_norm,
    nf = nf,
    ng = ng,
    n_fn_call = n_fn_call,
    n_gr_call = n_gr_call,
    n_fg_call = n_fg_call,
    n_hs_call = n_hs_call,
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
      fn = counted$fg$fn,
      gr = counted$fg$gr,
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
    callback_mode = "separate",
    max_iter = max_iter,
    final_f = metrics$final_f,
    grad_norm = metrics$grad_norm,
    nf = counts$n_fn_call,
    ng = counts$n_gr_call,
    n_fn_call = counts$n_fn_call,
    n_gr_call = counts$n_gr_call,
    n_fg_call = counts$n_fg_call,
    n_hs_call = counts$n_hs_call,
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

mize_fg_for_callback_mode <- function(fg, callback_mode) {
  if (callback_mode == "combined") {
    if (!is.function(fg$fg)) {
      stop("Combined callback mode requires an fg$fg function", call. = FALSE)
    }
    return(fg)
  }

  fg$fg <- NULL
  fg
}

benchmark_trace_callbacks <- function(fg) {
  trace <- new.env(parent = emptyenv())
  trace$objectives <- list()
  trace$gradients <- list()

  record_objective <- function(par, value) {
    if (is.numeric(value)) {
      trace$objectives[[length(trace$objectives) + 1L]] <- list(
        par = as.numeric(par),
        value = value
      )
    }
  }
  record_gradient <- function(par, value) {
    if (is.numeric(value)) {
      trace$gradients[[length(trace$gradients) + 1L]] <- list(
        par = as.numeric(par),
        value = value
      )
    }
  }
  wrap_callback <- function(callback_name) {
    callback <- fg[[callback_name]]
    if (!is.function(callback)) {
      return(callback)
    }

    function(...) {
      arguments <- list(...)
      par <- arguments[[1L]]
      value <- callback(...)
      if (callback_name == "fn") {
        record_objective(par, value)
      } else if (callback_name == "gr") {
        record_gradient(par, value)
      } else if (callback_name == "fg" && is.list(value)) {
        record_objective(par, value[["fn"]])
        record_gradient(par, value[["gr"]])
      }
      value
    }
  }

  traced <- fg
  traced$fn <- wrap_callback("fn")
  traced$gr <- wrap_callback("gr")
  traced$fg <- wrap_callback("fg")
  list(
    fg = traced,
    callbacks = function() {
      list(objectives = trace$objectives, gradients = trace$gradients)
    }
  )
}

benchmark_attach_gradient_trace <- function(
  progress,
  callbacks,
  initial_par,
  initial_grad_norm
) {
  if (!is.data.frame(progress) || nrow(progress) == 0L) {
    return(progress)
  }

  progress$g2n <- NA_real_
  selected_par <- as.numeric(initial_par)
  selected_grad_norm <- initial_grad_norm
  for (index in seq_len(nrow(progress))) {
    step <- progress$step[[index]]
    if (index > 1L && is.finite(step) && step > 0) {
      objective_limit <- min(
        as.integer(progress$nf[[index]]),
        length(callbacks$objectives)
      )
      objective_candidates <- if (objective_limit > 0L) {
        seq_len(objective_limit)
      } else {
        integer()
      }
      objective_candidates <- objective_candidates[vapply(
        callbacks$objectives[objective_candidates],
        function(record) {
          value_matches <- isTRUE(record$value == progress$f[[index]])
          candidate_step <- norm2(record$par - selected_par)
          step_tolerance <- 100 *
            .Machine$double.eps *
            max(1, abs(step), abs(candidate_step))
          value_matches &&
            is.finite(candidate_step) &&
            abs(candidate_step - step) <= step_tolerance
        },
        logical(1)
      )]
      if (length(objective_candidates) == 0L) {
        stop(
          "objective trace does not identify the selected point",
          call. = FALSE
        )
      }
      selected_par <- callbacks$objectives[[tail(
        objective_candidates,
        1L
      )]]$par

      gradient_limit <- min(
        as.integer(progress$ng[[index]]),
        length(callbacks$gradients)
      )
      gradient_candidates <- if (gradient_limit > 0L) {
        seq_len(gradient_limit)
      } else {
        integer()
      }
      gradient_candidates <- gradient_candidates[vapply(
        callbacks$gradients[gradient_candidates],
        function(record) identical(record$par, selected_par),
        logical(1)
      )]
      if (length(gradient_candidates) == 0L) {
        stop(
          "gradient trace does not identify the selected point",
          call. = FALSE
        )
      }
      selected_gradient <- callbacks$gradients[[tail(
        gradient_candidates,
        1L
      )]]$value
      selected_grad_norm <- norm2(selected_gradient)
    }
    progress$g2n[[index]] <- selected_grad_norm
  }
  progress
}

run_safeguarded_newton <- function(
  par,
  fg,
  cg_update,
  preconditioner,
  memory,
  line_search,
  step0,
  max_iter,
  store_progress
) {
  namespace <- asNamespace("mize")
  newton_constructor <- get(
    "newton_direction",
    envir = namespace,
    inherits = FALSE
  )
  optimizer_loop <- get("opt_loop", envir = namespace, inherits = FALSE)
  abs_tol <- sqrt(.Machine$double.eps)
  opt <- make_mize(
    method = "NEWTON",
    cg_update = cg_update,
    preconditioner = preconditioner,
    memory = memory,
    line_search = line_search,
    step0 = step0,
    ls_max_fn = 20,
    ls_max_gr = Inf,
    ls_max_fg = Inf,
    ls_max_alpha_mult = Inf,
    ls_max_alpha = Inf,
    ls_safe_cubic = FALSE,
    max_iter = max_iter,
    max_fn = Inf,
    max_gr = Inf,
    max_fg = Inf,
    abs_tol = abs_tol,
    rel_tol = abs_tol,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = abs_tol
  )
  opt$stages$gradient_descent$direction <- newton_constructor(
    try_safe_chol = TRUE
  )
  optimizer_loop(
    opt = opt,
    par = par,
    fg = fg,
    max_iter = max_iter,
    max_fn = Inf,
    max_gr = Inf,
    max_fg = Inf,
    abs_tol = abs_tol,
    rel_tol = abs_tol,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = abs_tol,
    check_conv_every = 1,
    log_every = 1,
    store_progress = store_progress
  )
}

run_mize_optimizer <- function(
  par,
  fg,
  method,
  cg_update,
  preconditioner,
  memory,
  safeguarded_newton,
  line_search,
  step0,
  max_iter,
  store_progress
) {
  if (safeguarded_newton) {
    return(run_safeguarded_newton(
      par = par,
      fg = fg,
      cg_update = cg_update,
      preconditioner = preconditioner,
      memory = memory,
      line_search = line_search,
      step0 = step0,
      max_iter = max_iter,
      store_progress = store_progress
    ))
  }

  mize(
    par = par,
    fg = fg,
    method = method,
    cg_update = cg_update,
    preconditioner = preconditioner,
    memory = memory,
    line_search = line_search,
    step0 = step0,
    ls_max_fn = 20,
    ls_max_gr = Inf,
    ls_max_fg = Inf,
    ls_max_alpha_mult = Inf,
    ls_max_alpha = Inf,
    ls_safe_cubic = FALSE,
    max_iter = max_iter,
    max_fn = Inf,
    max_gr = Inf,
    max_fg = Inf,
    check_conv_every = 1,
    log_every = 1,
    store_progress = store_progress
  )
}

mize_method_label <- function(
  method,
  cg_update,
  preconditioner,
  safeguarded_newton
) {
  if (safeguarded_newton) {
    return("NEWTON:safeguarded")
  }
  label <- if (identical(method, "CG")) {
    paste0(method, ":", cg_update)
  } else {
    method
  }
  if (identical(tolower(preconditioner), "l-bfgs")) {
    label <- paste0(label, ":L-BFGS-preconditioned")
  }
  label
}

benchmark_progress_rows <- function(
  progress,
  case,
  rep,
  profile,
  method,
  line_search,
  step0,
  callback_mode,
  max_iter
) {
  if (!is.data.frame(progress) || nrow(progress) == 0L) {
    return(NULL)
  }
  progress_iter <- suppressWarnings(as.integer(rownames(progress)))
  metadata <- data.frame(
    case = base::rep(case$name, nrow(progress)),
    source = base::rep(case$source, nrow(progress)),
    rep = base::rep(rep, nrow(progress)),
    profile = base::rep(profile, nrow(progress)),
    optimizer = base::rep("mize", nrow(progress)),
    method = base::rep(method, nrow(progress)),
    line_search = base::rep(line_search, nrow(progress)),
    step0 = base::rep(step0, nrow(progress)),
    callback_mode = base::rep(callback_mode, nrow(progress)),
    dimension = base::rep(case$dimension, nrow(progress)),
    max_iter = base::rep(max_iter, nrow(progress)),
    progress_iter = progress_iter,
    stringsAsFactors = FALSE
  )
  rownames(progress) <- NULL
  data.frame(metadata, progress, stringsAsFactors = FALSE)
}

run_mize_case <- function(
  case,
  rep,
  profile,
  method,
  cg_update,
  preconditioner,
  memory,
  safeguarded_newton,
  line_search,
  step0,
  callback_mode,
  max_iter,
  store_progress = FALSE
) {
  step0_value <- mize_step0_value(step0)
  run_config <- function(fg, capture_progress) {
    run_mize_optimizer(
      par = case$par,
      fg = fg,
      method = method,
      cg_update = cg_update,
      preconditioner = preconditioner,
      memory = memory,
      safeguarded_newton = safeguarded_newton,
      line_search = line_search,
      step0 = step0_value,
      max_iter = max_iter,
      store_progress = capture_progress
    )
  }

  counted <- counted_fg(case$fg)
  callback_fg <- mize_fg_for_callback_mode(counted$fg, callback_mode)
  timed <- run_timed(function() {
    run_config(callback_fg, capture_progress = FALSE)
  })
  counts <- counted$counts()
  res <- timed$result
  error <- if (inherits(res, "error")) conditionMessage(res) else ""
  par <- if (inherits(res, "error")) NULL else res$par
  metrics <- safe_final_metrics(case$fg, par)
  termination <- if (inherits(res, "error")) {
    "error"
  } else {
    paste0(res$terminate$what, "=", res$terminate$val)
  }
  method_label <- mize_method_label(
    method = method,
    cg_update = cg_update,
    preconditioner = preconditioner,
    safeguarded_newton = safeguarded_newton
  )

  row <- row_result(
    case = case,
    rep = rep,
    optimizer = "mize",
    method = method_label,
    line_search = line_search,
    step0 = step0,
    callback_mode = callback_mode,
    max_iter = max_iter,
    final_f = metrics$final_f,
    grad_norm = metrics$grad_norm,
    nf = if (inherits(res, "error")) NA_integer_ else res$nf,
    ng = if (inherits(res, "error")) NA_integer_ else res$ng,
    n_fn_call = counts$n_fn_call,
    n_gr_call = counts$n_gr_call,
    n_fg_call = counts$n_fg_call,
    n_hs_call = counts$n_hs_call,
    elapsed = timed$elapsed,
    iter = if (inherits(res, "error")) NA_integer_ else res$iter,
    termination = termination,
    failure_mode = mize_failure(res),
    warnings = timed$warnings,
    error = error
  )
  progress <- NULL
  if (store_progress && !inherits(res, "error")) {
    replay_counted <- counted_fg(case$fg)
    replay_traced <- benchmark_trace_callbacks(replay_counted$fg)
    replay_fg <- mize_fg_for_callback_mode(replay_traced$fg, callback_mode)
    replay_timed <- run_timed(function() {
      run_config(replay_fg, capture_progress = TRUE)
    })
    replay <- replay_timed$result
    replay_projection <- if (inherits(replay, "error")) {
      replay
    } else {
      replay[c("par", "nf", "ng", "iter", "terminate")]
    }
    primary_projection <- res[c("par", "nf", "ng", "iter", "terminate")]
    replay_matches <- !inherits(replay, "error") &&
      isTRUE(all.equal(
        primary_projection,
        replay_projection,
        tolerance = 0,
        check.attributes = TRUE
      )) &&
      identical(counts, replay_counted$counts())
    if (!replay_matches) {
      stop(
        "unmeasured progress replay did not match the timed run",
        call. = FALSE
      )
    }
    progress <- benchmark_attach_gradient_trace(
      replay$progress,
      replay_traced$callbacks(),
      case$par,
      case$initial_grad_norm
    )
  }
  attr(row, "benchmark_progress") <- benchmark_progress_rows(
    progress = progress,
    case = case,
    rep = rep,
    profile = profile,
    method = method_label,
    line_search = line_search,
    step0 = step0,
    callback_mode = callback_mode,
    max_iter = max_iter
  )
  row
}

benchmark_profile_value <- function(profile, name, default) {
  value <- profile[[name]]
  if (is.null(value)) default else value
}

benchmark_rbind_fill <- function(values) {
  if (length(values) == 0L) {
    return(data.frame())
  }
  column_names <- Reduce(union, lapply(values, names))
  prototypes <- lapply(column_names, function(name) {
    source <- values[[which(vapply(
      values,
      function(value) name %in% names(value),
      logical(1)
    ))[[1L]]]][[name]]
    if (is.character(source)) {
      NA_character_
    } else if (is.integer(source)) {
      NA_integer_
    } else if (is.logical(source)) {
      NA
    } else {
      NA_real_
    }
  })
  names(prototypes) <- column_names

  aligned <- lapply(values, function(value) {
    for (name in setdiff(column_names, names(value))) {
      value[[name]] <- rep(prototypes[[name]], nrow(value))
    }
    value[column_names]
  })
  do.call(rbind, aligned)
}

run_case_grid <- function(case, config) {
  method_profiles <- benchmark_method_profiles()[config$methods]
  store_progress <- !is.null(config$progress_out) ||
    !is.null(config$summary_out)

  rows <- list()
  resource_rows <- list()
  progress_rows <- list()
  for (rep in seq_len(config$reps)) {
    for (profile in names(method_profiles)) {
      method_config <- method_profiles[[profile]]
      if (identical(method_config$optimizer, "stats::optim")) {
        if (rep == 1L && config$warmup > 0L) {
          for (warmup in seq_len(config$warmup)) {
            run_optim_case(
              case = case,
              rep = -warmup,
              method = method_config$method,
              max_iter = config$max_iter
            )
          }
        }
        row <- run_optim_case(
          case = case,
          rep = rep,
          method = method_config$method,
          max_iter = config$max_iter
        )
        rows[[length(rows) + 1L]] <- row
        resource_rows[[length(resource_rows) + 1L]] <- data.frame(
          profile = profile,
          row,
          stringsAsFactors = FALSE
        )
        next
      }

      for (line_search in config$line_searches) {
        for (step0 in config$step0) {
          for (callback_mode in config$callback_modes) {
            run_profile <- function(run_rep, capture_progress) {
              run_mize_case(
                case = case,
                rep = run_rep,
                profile = profile,
                method = method_config$method,
                cg_update = benchmark_profile_value(
                  method_config,
                  "cg_update",
                  "PR+"
                ),
                preconditioner = benchmark_profile_value(
                  method_config,
                  "preconditioner",
                  ""
                ),
                memory = benchmark_profile_value(
                  method_config,
                  "memory",
                  5L
                ),
                safeguarded_newton = benchmark_profile_value(
                  method_config,
                  "safeguarded_newton",
                  FALSE
                ),
                line_search = line_search,
                step0 = step0,
                callback_mode = callback_mode,
                max_iter = config$max_iter,
                store_progress = capture_progress
              )
            }
            if (rep == 1L && config$warmup > 0L) {
              for (warmup in seq_len(config$warmup)) {
                run_profile(-warmup, capture_progress = FALSE)
              }
            }
            row <- run_profile(rep, capture_progress = store_progress)
            progress <- attr(row, "benchmark_progress")
            attr(row, "benchmark_progress") <- NULL
            rows[[length(rows) + 1L]] <- row
            resource_rows[[length(resource_rows) + 1L]] <- data.frame(
              profile = profile,
              row,
              stringsAsFactors = FALSE
            )
            if (!is.null(progress)) {
              progress_rows[[length(progress_rows) + 1L]] <- progress
            }
          }
        }
      }
    }
  }

  list(
    results = do.call(rbind, rows),
    resource_results = do.call(rbind, resource_rows),
    progress = if (length(progress_rows) == 0L) {
      data.frame()
    } else {
      benchmark_rbind_fill(progress_rows)
    }
  )
}

run_benchmarks <- function(config) {
  cases <- benchmark_cases(config)
  if (length(cases) == 0L) {
    stop("No benchmark cases selected", call. = FALSE)
  }

  grids <- lapply(cases, run_case_grid, config = config)
  progress <- lapply(grids, function(grid) grid$progress)
  progress <- Filter(function(value) nrow(value) > 0L, progress)
  list(
    results = do.call(rbind, lapply(grids, function(grid) grid$results)),
    resource_results = do.call(
      rbind,
      lapply(grids, function(grid) grid$resource_results)
    ),
    progress = if (length(progress) == 0L) {
      data.frame()
    } else {
      benchmark_rbind_fill(progress)
    },
    cases = benchmark_case_manifest(cases)
  )
}

benchmark_first_target <- function(progress, value_name, target) {
  unavailable <- list(
    applicable = FALSE,
    hit = FALSE,
    iter = NA_integer_,
    nf = NA_integer_,
    ng = NA_integer_
  )
  if (
    !is.numeric(target) ||
      length(target) != 1L ||
      !is.finite(target)
  ) {
    return(unavailable)
  }
  unavailable$applicable <- TRUE
  if (
    !is.data.frame(progress) ||
      nrow(progress) == 0L ||
      !value_name %in% names(progress)
  ) {
    return(unavailable)
  }
  values <- progress[[value_name]]
  hits <- which(is.finite(values) & values <= target)
  if (length(hits) == 0L) {
    return(unavailable)
  }

  index <- hits[[1L]]
  list(
    applicable = TRUE,
    hit = TRUE,
    iter = as.integer(progress$progress_iter[[index]]),
    nf = as.integer(progress$nf[[index]]),
    ng = as.integer(progress$ng[[index]])
  )
}

benchmark_frequency <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    return("")
  }
  frequencies <- table(values)
  paste(
    paste(names(frequencies), as.integer(frequencies), sep = "="),
    collapse = ";"
  )
}

benchmark_progress_count <- function(progress, name, value) {
  if (!name %in% names(progress)) {
    return(0L)
  }
  sum(!is.na(progress[[name]]) & progress[[name]] == value)
}

benchmark_resource_summary_row <- function(row, progress, case_manifest) {
  case_row <- case_manifest[
    case_manifest$case == row$case[[1L]],
    ,
    drop = FALSE
  ]
  if (nrow(case_row) != 1L) {
    stop("resource summary requires one matching case-manifest row")
  }

  objective_target <- benchmark_first_target(
    progress,
    "f",
    case_row$objective_reduction_target[[1L]]
  )
  gradient_target <- benchmark_first_target(
    progress,
    "g2n",
    case_row$gradient_reduction_target[[1L]]
  )
  reference_target <- benchmark_first_target(
    progress,
    "f",
    case_row$reference_gap_target[[1L]]
  )
  nonstationary_no_step <- if (
    nrow(progress) > 0L &&
      all(c("ls_outcome", "g2n") %in% names(progress))
  ) {
    !is.na(progress$ls_outcome) &
      progress$ls_outcome == "no_step" &
      is.finite(progress$g2n) &
      progress$g2n > case_row$gradient_reduction_target[[1L]]
  } else {
    logical()
  }
  non_descent_escape <- if (
    nrow(progress) > 0L &&
      all(c("ls_reason", "ls_outcome", "alpha") %in% names(progress))
  ) {
    !is.na(progress$ls_reason) &
      progress$ls_reason == "non_descent_direction" &
      ((is.finite(progress$alpha) & progress$alpha > 0) |
        (!is.na(progress$ls_outcome) & progress$ls_outcome != "no_step"))
  } else {
    logical()
  }
  termination_what <- sub("=.*$", "", row$termination[[1L]])
  convergence_termination <- termination_what %in%
    c(
      "abs_tol",
      "rel_tol",
      "grad_tol",
      "ginf_tol",
      "step_tol"
    )
  objective_tolerance <- 100 *
    .Machine$double.eps *
    max(1, abs(case_row$initial_f[[1L]]))

  data.frame(
    row,
    objective_reduction_target = case_row$objective_reduction_target[[1L]],
    objective_reduction_hit = objective_target$hit,
    objective_reduction_iter = objective_target$iter,
    objective_reduction_nf = objective_target$nf,
    objective_reduction_ng = objective_target$ng,
    gradient_reduction_target = case_row$gradient_reduction_target[[1L]],
    gradient_reduction_hit = gradient_target$hit,
    gradient_reduction_iter = gradient_target$iter,
    gradient_reduction_nf = gradient_target$nf,
    gradient_reduction_ng = gradient_target$ng,
    reference_gap_target_applicable = reference_target$applicable,
    reference_gap_target = case_row$reference_gap_target[[1L]],
    reference_gap_hit = reference_target$hit,
    reference_gap_iter = reference_target$iter,
    reference_gap_nf = reference_target$nf,
    reference_gap_ng = reference_target$ng,
    ls_outcome_wolfe = benchmark_progress_count(
      progress,
      "ls_outcome",
      "wolfe"
    ),
    ls_outcome_armijo = benchmark_progress_count(
      progress,
      "ls_outcome",
      "armijo"
    ),
    ls_outcome_improving_fallback = benchmark_progress_count(
      progress,
      "ls_outcome",
      "improving_fallback"
    ),
    ls_outcome_no_step = benchmark_progress_count(
      progress,
      "ls_outcome",
      "no_step"
    ),
    ls_reason_frequency = if ("ls_reason" %in% names(progress)) {
      benchmark_frequency(progress$ls_reason)
    } else {
      ""
    },
    direction_reason_frequency = if ("direction_reason" %in% names(progress)) {
      benchmark_frequency(progress$direction_reason)
    } else {
      ""
    },
    nonstationary_no_step_count = sum(nonstationary_no_step),
    false_convergence = convergence_termination &&
      any(nonstationary_no_step),
    objective_worse_than_start = is.finite(row$final_f[[1L]]) &&
      row$final_f[[1L]] > case_row$initial_f[[1L]] + objective_tolerance,
    non_descent_escape_count = sum(non_descent_escape),
    final_metric_nonfinite = !is.finite(row$final_f[[1L]]) ||
      !is.finite(row$grad_norm[[1L]]),
    callback_or_optimizer_error = nzchar(row$error[[1L]]) ||
      identical(row$failure_mode[[1L]], "nonfinite"),
    stringsAsFactors = FALSE
  )
}

benchmark_resource_summary <- function(
  resource_results,
  progress,
  case_manifest
) {
  same_value <- function(values, target) {
    if (length(target) != 1L || is.na(target)) {
      return(is.na(values))
    }
    !is.na(values) & values == target
  }
  rows <- lapply(seq_len(nrow(resource_results)), function(index) {
    row <- resource_results[index, , drop = FALSE]
    run_progress <- if (nrow(progress) == 0L) {
      data.frame()
    } else {
      progress[
        same_value(progress$case, row$case[[1L]]) &
          same_value(progress$rep, row$rep[[1L]]) &
          same_value(progress$profile, row$profile[[1L]]) &
          same_value(progress$line_search, row$line_search[[1L]]) &
          same_value(progress$step0, row$step0[[1L]]) &
          same_value(progress$callback_mode, row$callback_mode[[1L]]),
        ,
        drop = FALSE
      ]
    }
    benchmark_resource_summary_row(row, run_progress, case_manifest)
  })
  do.call(rbind, rows)
}

benchmark_resolved_output_path <- function(path, option) {
  if (
    !is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(path)
  ) {
    stop(option, " must be a non-empty path", call. = FALSE)
  }

  expanded <- path.expand(path)
  link_target <- Sys.readlink(expanded)
  if (!is.na(link_target) && nzchar(link_target)) {
    stop(option, " must not be a symbolic link", call. = FALSE)
  }
  if (file.exists(expanded)) {
    return(normalizePath(expanded, winslash = "/", mustWork = TRUE))
  }

  resolved_parent <- normalizePath(
    dirname(expanded),
    winslash = "/",
    mustWork = FALSE
  )
  file.path(resolved_parent, basename(expanded))
}

benchmark_existing_file_identity <- function(path, option) {
  if (!file.exists(path)) {
    stop("existing-file identity requires an existing path", call. = FALSE)
  }

  stat_command <- Sys.which("stat")
  system_name <- Sys.info()[["sysname"]]
  stat_args <- if (
    nzchar(stat_command) &&
      system_name %in% c("Darwin", "FreeBSD", "OpenBSD", "NetBSD")
  ) {
    c("-f", "%d:%i", shQuote(path))
  } else if (nzchar(stat_command) && identical(system_name, "Linux")) {
    c("-Lc", "%d:%i", shQuote(path))
  } else {
    NULL
  }
  if (!is.null(stat_args)) {
    stat_output <- suppressWarnings(system2(
      stat_command,
      args = stat_args,
      stdout = TRUE,
      stderr = TRUE
    ))
    stat_status <- attr(stat_output, "status")
    if (is.null(stat_status)) {
      stat_status <- 0L
    }
    if (
      identical(stat_status, 0L) &&
        length(stat_output) == 1L &&
        grepl("^[0-9]+:[0-9]+$", stat_output)
    ) {
      return(stat_output)
    }
  }

  if (requireNamespace("fs", quietly = TRUE)) {
    info <- fs::file_info(path)
    device <- as.character(info$device_id[[1L]])
    inode <- as.character(info$inode[[1L]])
    if (!is.na(device) && !is.na(inode)) {
      return(paste(device, inode, sep = ":"))
    }
  }

  stop(
    "could not establish existing-file identity for ",
    option,
    "; install the fs package or provide unused output paths",
    call. = FALSE
  )
}

benchmark_directory_case_sensitive <- function(directory) {
  if (!dir.exists(directory)) {
    stop("output parent directory does not exist: ", directory, call. = FALSE)
  }

  probe <- tempfile(pattern = ".mize-case-sensitivity-", tmpdir = directory)
  case_variant <- file.path(directory, toupper(basename(probe)))
  if (file.exists(case_variant)) {
    stop(
      "could not choose an unused output case-sensitivity probe",
      call. = FALSE
    )
  }
  if (!file.create(probe, showWarnings = FALSE)) {
    stop(
      "could not probe output-directory case sensitivity: ",
      directory,
      call. = FALSE
    )
  }
  on.exit(unlink(probe, force = TRUE), add = TRUE)

  if (!file.exists(case_variant)) {
    return(TRUE)
  }
  probe_identity <- benchmark_existing_file_identity(
    probe,
    "case-sensitivity probe"
  )
  variant_identity <- benchmark_existing_file_identity(
    case_variant,
    "case-sensitivity probe"
  )
  !identical(probe_identity, variant_identity)
}

benchmark_case_aliased_paths <- function(resolved) {
  case_aliased <- rep(FALSE, length(resolved))
  parents <- dirname(resolved)
  for (parent in unique(parents)) {
    indices <- which(parents == parent)
    folded <- tolower(basename(resolved[indices]))
    candidates <- duplicated(folded) |
      duplicated(folded, fromLast = TRUE)
    if (
      any(candidates) &&
        !benchmark_directory_case_sensitive(parent)
    ) {
      case_aliased[indices[candidates]] <- TRUE
    }
  }
  case_aliased
}

benchmark_preflight_output_paths <- function(paths) {
  paths <- paths[!vapply(paths, is.null, logical(1))]
  if (length(paths) == 0L) {
    return(invisible(character()))
  }
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("output paths must be named by their command-line options")
  }

  resolved <- Map(
    benchmark_resolved_output_path,
    path = paths,
    option = names(paths)
  )
  resolved <- vapply(resolved, identity, character(1))
  path_aliased <- duplicated(resolved) |
    duplicated(resolved, fromLast = TRUE)
  if (any(path_aliased)) {
    stop(
      "active output paths must resolve to different destinations: ",
      paste(names(resolved)[path_aliased], collapse = ", "),
      call. = FALSE
    )
  }

  case_aliased <- benchmark_case_aliased_paths(resolved)
  if (any(case_aliased)) {
    stop(
      "active output paths must not be case-equivalent aliases: ",
      paste(names(resolved)[case_aliased], collapse = ", "),
      call. = FALSE
    )
  }

  existing <- which(file.exists(resolved))
  if (length(existing) > 1L) {
    identities <- vapply(
      existing,
      function(index) {
        benchmark_existing_file_identity(
          resolved[[index]],
          names(resolved)[[index]]
        )
      },
      character(1)
    )
    identity_aliased <- duplicated(identities) |
      duplicated(identities, fromLast = TRUE)
    if (any(identity_aliased)) {
      stop(
        "active output paths must identify different existing files: ",
        paste(
          names(resolved)[existing[identity_aliased]],
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }

  invisible(resolved)
}

write_results <- function(results, out) {
  if (is.null(out)) {
    utils::write.csv(results, file = stdout(), row.names = FALSE, na = "")
    return(invisible(NULL))
  }

  utils::write.csv(results, file = out, row.names = FALSE, na = "")
  message("Wrote benchmark results to ", out)
}

write_benchmark_artifact <- function(results, out, label) {
  if (is.null(out)) {
    return(invisible(NULL))
  }
  utils::write.csv(results, file = out, row.names = FALSE, na = "")
  message("Wrote ", label, " to ", out)
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  loader <- load_local_mize()
  message("Loaded mize via ", loader)
  benchmark <- run_benchmarks(config)
  write_results(benchmark$results, config$out)
  write_benchmark_artifact(
    benchmark$progress,
    config$progress_out,
    "optimizer progress"
  )
  write_benchmark_artifact(
    benchmark$cases,
    config$case_out,
    "benchmark case manifest"
  )
  if (!is.null(config$summary_out)) {
    summary <- benchmark_resource_summary(
      benchmark$resource_results,
      benchmark$progress,
      benchmark$cases
    )
    write_benchmark_artifact(
      summary,
      config$summary_out,
      "resource summary"
    )
  }
}

if (sys.nframe() == 0L) {
  main()
}
