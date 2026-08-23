#!/usr/bin/env Rscript

source(file.path("scripts", "benchmark-optimizers.R"))

packet5a_usage <- function(status = 0L) {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/benchmark-hager-zhang-rescue.R [options]",
      "",
      "Options:",
      "  --scope SCOPE             tranche or all35.",
      "  --funconstrain-commit SHA Required funconstrain source commit.",
      "  --out PATH                Measured run CSV.",
      "  --progress-out PATH       Unmeasured progress replay CSV.",
      "  --case-out PATH           Case and target manifest CSV.",
      "  --summary-out PATH        Per-run resource summary CSV.",
      "  --initializer-out PATH    Later-HZ initializer evidence CSV.",
      "  --manifest-out PATH       Provenance and artifact manifest CSV.",
      "  --help                    Show this help.",
      "",
      "The scientific grid is fixed by the Packet 5A contract.",
      sep = "\n"
    )
  )
  quit(save = "no", status = status)
}

packet5a_problem_names <- function(scope) {
  if (identical(scope, "tranche")) {
    return(c("brown_bs", "gulf", "meyer"))
  }
  c(
    "bard",
    "beale",
    "biggs_exp6",
    "box_3d",
    "brown_al",
    "brown_bs",
    "brown_den",
    "broyden_band",
    "broyden_tri",
    "chebyquad",
    "disc_bv",
    "disc_ie",
    "ex_powell",
    "ex_rosen",
    "freud_roth",
    "gauss",
    "gulf",
    "helical",
    "jenn_samp",
    "kow_osb",
    "linfun_fr",
    "linfun_r1",
    "linfun_r1z",
    "meyer",
    "osborne_1",
    "osborne_2",
    "penalty_1",
    "penalty_2",
    "powell_bs",
    "powell_s",
    "rosen",
    "trigon",
    "var_dim",
    "watson",
    "wood"
  )
}

packet5a_profiles <- function() {
  c(
    "mize-newton",
    "mize-newton-safeguarded",
    "mize-tn",
    "mize-tn-l-bfgs",
    "mize-bfgs",
    "mize-l-bfgs",
    "mize-sr1",
    "mize-cg-pr+",
    "mize-cg-hz+",
    "mize-cg-pr+-l-bfgs"
  )
}

packet5a_policies <- function() {
  list(
    `current-hz-20` = list(
      initializer = "hz",
      ls_max_fn = 20L,
      scale_rescue = FALSE
    ),
    `rescue-hz-20` = list(
      initializer = "hz",
      ls_max_fn = 20L,
      scale_rescue = TRUE
    ),
    `quadratic-20` = list(
      initializer = "quadratic",
      ls_max_fn = 20L,
      scale_rescue = FALSE
    ),
    `slope-ratio-20` = list(
      initializer = "slope ratio",
      ls_max_fn = 20L,
      scale_rescue = FALSE
    ),
    `current-hz-100` = list(
      initializer = "hz",
      ls_max_fn = 100L,
      scale_rescue = FALSE
    )
  )
}

packet5a_parse_args <- function(args) {
  config <- list(
    scope = NULL,
    funconstrain_commit = NULL,
    out = NULL,
    progress_out = NULL,
    case_out = NULL,
    summary_out = NULL,
    initializer_out = NULL,
    manifest_out = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      packet5a_usage()
    }
    if (i == length(args)) {
      stop(arg, " requires a value", call. = FALSE)
    }
    value <- args[[i + 1L]]
    if (arg == "--scope") {
      config$scope <- tolower(value)
    } else if (arg == "--funconstrain-commit") {
      config$funconstrain_commit <- value
    } else if (arg == "--out") {
      config$out <- value
    } else if (arg == "--progress-out") {
      config$progress_out <- value
    } else if (arg == "--case-out") {
      config$case_out <- value
    } else if (arg == "--summary-out") {
      config$summary_out <- value
    } else if (arg == "--initializer-out") {
      config$initializer_out <- value
    } else if (arg == "--manifest-out") {
      config$manifest_out <- value
    } else {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    i <- i + 2L
  }

  if (!config$scope %in% c("tranche", "all35")) {
    stop("--scope must be tranche or all35", call. = FALSE)
  }
  if (
    is.null(config$funconstrain_commit) ||
      !nzchar(config$funconstrain_commit)
  ) {
    stop("--funconstrain-commit is required", call. = FALSE)
  }
  output_paths <- list(
    `--out` = config$out,
    `--progress-out` = config$progress_out,
    `--case-out` = config$case_out,
    `--summary-out` = config$summary_out,
    `--initializer-out` = config$initializer_out,
    `--manifest-out` = config$manifest_out
  )
  if (any(vapply(output_paths, is.null, logical(1)))) {
    stop("all six Packet 5A output paths are required", call. = FALSE)
  }
  benchmark_preflight_output_paths(output_paths)
  config
}

packet5a_case_config <- function(config) {
  list(
    cases = "spd-quadratic",
    seed = 20260823L,
    spd_n = 20L,
    include_funconstrain = TRUE,
    funconstrain_cases = packet5a_problem_names(config$scope),
    funconstrain_commit = config$funconstrain_commit
  )
}

packet5a_method_c2 <- function(method) {
  if (tolower(method) %in% c("newton", "phess", "bfgs", "l-bfgs", "tn")) {
    0.9
  } else {
    0.1
  }
}

packet5a_try_newton_step <- function(method) {
  !identical(tolower(method), "cg")
}

packet5a_initializer_recorder <- function() {
  records <- list()
  list(
    observe = function(record) {
      records[[length(records) + 1L]] <<- record
      invisible(NULL)
    },
    values = function() records
  )
}

packet5a_run_optimizer <- function(
  par,
  fg,
  method,
  cg_update,
  preconditioner,
  memory,
  safeguarded_newton,
  policy,
  store_progress,
  initialization_observer = NULL
) {
  namespace <- asNamespace("mize")
  optimizer_loop <- get("opt_loop", envir = namespace, inherits = FALSE)
  abs_tol <- sqrt(.Machine$double.eps)
  opt <- make_mize(
    method = method,
    cg_update = cg_update,
    preconditioner = preconditioner,
    memory = memory,
    line_search = "Hager-Zhang",
    step0 = "hz",
    step_next_init = policy$initializer,
    ls_max_fn = policy$ls_max_fn,
    ls_max_gr = Inf,
    ls_max_fg = Inf,
    ls_max_alpha_mult = Inf,
    ls_max_alpha = Inf,
    ls_safe_cubic = FALSE,
    max_iter = 100L,
    max_fn = Inf,
    max_gr = Inf,
    max_fg = Inf,
    abs_tol = abs_tol,
    rel_tol = abs_tol,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = abs_tol
  )
  if (safeguarded_newton) {
    newton_constructor <- get(
      "newton_direction",
      envir = namespace,
      inherits = FALSE
    )
    opt$stages$gradient_descent$direction <- newton_constructor(
      try_safe_chol = TRUE
    )
  }
  if (policy$scale_rescue || is.function(initialization_observer)) {
    hager_zhang_constructor <- get(
      "hager_zhang_ls",
      envir = namespace,
      inherits = FALSE
    )
    opt$stages$gradient_descent$step_size <- hager_zhang_constructor(
      c1 = 1e-4,
      c2 = packet5a_method_c2(method),
      initializer = policy$initializer,
      initial_step_length = "hz",
      try_newton_step = packet5a_try_newton_step(method),
      max_fn = policy$ls_max_fn,
      max_gr = Inf,
      max_fg = Inf,
      max_alpha_mult = Inf,
      strong_curvature = FALSE,
      approx_armijo = TRUE,
      scale_rescue = policy$scale_rescue,
      initialization_observer = initialization_observer
    )
  }

  optimizer_loop(
    opt = opt,
    par = par,
    fg = fg,
    max_iter = 100L,
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

packet5a_initializer_rows <- function(
  records,
  case,
  rep,
  profile,
  method,
  policy,
  callback_mode
) {
  if (length(records) == 0L) {
    return(NULL)
  }
  rows <- lapply(records, function(record) {
    data.frame(
      case = case$name,
      source = case$source,
      rep = rep,
      profile = profile,
      method = method,
      policy = policy,
      callback_mode = callback_mode,
      dimension = case$dimension,
      iter = as.integer(record$iter),
      previous_alpha = record$previous_alpha,
      previous_slope = record$previous_slope,
      current_slope = record$current_slope,
      unrescued_alpha = record$unrescued_alpha,
      scale_estimate = record$scale_estimate,
      trial_evaluation_budget = record$trial_evaluation_budget,
      required_contractions = record$required_contractions,
      available_contractions = record$available_contractions,
      selected_initial_alpha = record$selected_initial_alpha,
      rescue_applied = record$rescue_applied,
      observer_provider_callbacks = 0L,
      probe_function_evaluations = record$probe_function_evaluations,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

packet5a_run_one <- function(
  case,
  rep,
  profile,
  method_config,
  policy_name,
  policy,
  callback_mode,
  capture_progress
) {
  method <- method_config$method
  cg_update <- benchmark_profile_value(method_config, "cg_update", "PR+")
  preconditioner <- benchmark_profile_value(
    method_config,
    "preconditioner",
    ""
  )
  memory <- benchmark_profile_value(method_config, "memory", 5L)
  safeguarded_newton <- benchmark_profile_value(
    method_config,
    "safeguarded_newton",
    FALSE
  )
  run_config <- function(fg, store_progress, observer = NULL) {
    packet5a_run_optimizer(
      par = case$par,
      fg = fg,
      method = method,
      cg_update = cg_update,
      preconditioner = preconditioner,
      memory = memory,
      safeguarded_newton = safeguarded_newton,
      policy = policy,
      store_progress = store_progress,
      initialization_observer = observer
    )
  }

  counted <- counted_fg(case$fg)
  callback_fg <- mize_fg_for_callback_mode(counted$fg, callback_mode)
  timed <- run_timed(function() run_config(callback_fg, FALSE))
  counts <- counted$counts()
  result <- timed$result
  error <- if (inherits(result, "error")) conditionMessage(result) else ""
  par <- if (inherits(result, "error")) NULL else result$par
  metrics <- safe_final_metrics(case$fg, par)
  termination <- if (inherits(result, "error")) {
    "error"
  } else {
    paste0(result$terminate$what, "=", result$terminate$val)
  }
  method_label <- mize_method_label(
    method,
    cg_update,
    preconditioner,
    safeguarded_newton
  )
  row <- row_result(
    case = case,
    rep = rep,
    optimizer = "mize",
    method = method_label,
    line_search = "Hager-Zhang",
    step0 = "hz",
    callback_mode = callback_mode,
    max_iter = 100L,
    final_f = metrics$final_f,
    grad_norm = metrics$grad_norm,
    nf = if (inherits(result, "error")) NA_integer_ else result$nf,
    ng = if (inherits(result, "error")) NA_integer_ else result$ng,
    n_fn_call = counts$n_fn_call,
    n_gr_call = counts$n_gr_call,
    n_fg_call = counts$n_fg_call,
    n_hs_call = counts$n_hs_call,
    elapsed = timed$elapsed,
    iter = if (inherits(result, "error")) NA_integer_ else result$iter,
    termination = termination,
    failure_mode = mize_failure(result),
    warnings = timed$warnings,
    error = error
  )
  row <- data.frame(
    profile = profile,
    policy = policy_name,
    later_initializer = policy$initializer,
    ls_max_fn = policy$ls_max_fn,
    scale_rescue = policy$scale_rescue,
    row,
    stringsAsFactors = FALSE
  )

  progress <- NULL
  initializer_rows <- NULL
  if (capture_progress && !inherits(result, "error")) {
    replay_counted <- counted_fg(case$fg)
    replay_traced <- benchmark_trace_callbacks(replay_counted$fg)
    replay_fg <- mize_fg_for_callback_mode(replay_traced$fg, callback_mode)
    recorder <- if (policy$initializer == "hz") {
      packet5a_initializer_recorder()
    } else {
      NULL
    }
    replay_timed <- run_timed(function() {
      run_config(
        replay_fg,
        TRUE,
        if (is.null(recorder)) NULL else recorder$observe
      )
    })
    replay <- replay_timed$result
    replay_projection <- if (inherits(replay, "error")) {
      replay
    } else {
      replay[c("par", "nf", "ng", "iter", "terminate")]
    }
    primary_projection <- result[c("par", "nf", "ng", "iter", "terminate")]
    replay_matches <- !inherits(replay, "error") &&
      length(replay_timed$warnings) == 0L &&
      isTRUE(all.equal(
        primary_projection,
        replay_projection,
        tolerance = 0,
        check.attributes = TRUE
      )) &&
      identical(counts, replay_counted$counts())
    if (!replay_matches) {
      stop(
        "unmeasured Packet 5A replay did not match the timed run",
        call. = FALSE
      )
    }
    replay_progress <- benchmark_attach_gradient_trace(
      replay$progress,
      replay_traced$callbacks(),
      case$par,
      case$initial_grad_norm
    )
    progress <- benchmark_progress_rows(
      progress = replay_progress,
      case = case,
      rep = rep,
      profile = profile,
      method = method_label,
      line_search = "Hager-Zhang",
      step0 = "hz",
      callback_mode = callback_mode,
      max_iter = 100L
    )
    progress <- data.frame(
      policy = policy_name,
      later_initializer = policy$initializer,
      ls_max_fn = policy$ls_max_fn,
      scale_rescue = policy$scale_rescue,
      progress,
      stringsAsFactors = FALSE
    )
    if (!is.null(recorder)) {
      initializer_rows <- packet5a_initializer_rows(
        records = recorder$values(),
        case = case,
        rep = rep,
        profile = profile,
        method = method_label,
        policy = policy_name,
        callback_mode = callback_mode
      )
    }
  }

  list(row = row, progress = progress, initializers = initializer_rows)
}

packet5a_rbind <- function(values) {
  values <- Filter(Negate(is.null), values)
  benchmark_rbind_fill(values)
}

packet5a_run_grid <- function(cases) {
  profiles <- benchmark_method_profiles()[packet5a_profiles()]
  policies <- packet5a_policies()
  callback_modes <- c("combined", "separate")
  rows <- list()
  progress_rows <- list()
  initializer_rows <- list()
  for (case in cases) {
    for (profile in names(profiles)) {
      for (policy_name in names(policies)) {
        for (callback_mode in callback_modes) {
          for (rep in seq_len(3L)) {
            if (rep == 1L) {
              packet5a_run_one(
                case,
                rep = -1L,
                profile,
                profiles[[profile]],
                policy_name,
                policies[[policy_name]],
                callback_mode,
                capture_progress = FALSE
              )
            }
            run <- packet5a_run_one(
              case,
              rep,
              profile,
              profiles[[profile]],
              policy_name,
              policies[[policy_name]],
              callback_mode,
              capture_progress = TRUE
            )
            rows[[length(rows) + 1L]] <- run$row
            progress_rows[[length(progress_rows) + 1L]] <- run$progress
            initializer_rows[[length(initializer_rows) + 1L]] <-
              run$initializers
          }
        }
      }
    }
  }
  list(
    runs = packet5a_rbind(rows),
    progress = packet5a_rbind(progress_rows),
    initializers = packet5a_rbind(initializer_rows)
  )
}

packet5a_select_run_progress <- function(progress, row) {
  progress[
    progress$case == row$case &
      progress$rep == row$rep &
      progress$profile == row$profile &
      progress$policy == row$policy &
      progress$callback_mode == row$callback_mode,
    ,
    drop = FALSE
  ]
}

packet5a_summary <- function(runs, progress, case_manifest) {
  rows <- lapply(seq_len(nrow(runs)), function(index) {
    row <- runs[index, , drop = FALSE]
    selected <- packet5a_select_run_progress(progress, row)
    summary <- benchmark_resource_summary_row(row, selected, case_manifest)
    summary
  })
  do.call(rbind, rows)
}

packet5a_git_value <- function(args) {
  value <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(value, "status")) || length(value) != 1L) {
    return(NA_character_)
  }
  value
}

packet5a_sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", output[[1L]])
}

packet5a_manifest <- function(config, artifacts) {
  description <- utils::packageDescription("mize")
  funconstrain <- utils::packageDescription("funconstrain")
  session <- sessionInfo()
  rows <- lapply(names(artifacts), function(name) {
    path <- artifacts[[name]]
    data.frame(
      scope = config$scope,
      artifact = name,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      data_rows = max(0L, length(readLines(path, warn = FALSE)) - 1L),
      sha256 = packet5a_sha256(path),
      mize_version = as.character(description[["Version"]]),
      mize_head = packet5a_git_value(c("rev-parse", "HEAD")),
      funconstrain_version = as.character(funconstrain[["Version"]]),
      funconstrain_commit = config$funconstrain_commit,
      r_version = R.version.string,
      platform = R.version$platform,
      running = session$running,
      blas = session$BLAS,
      lapack = session$LAPACK,
      seed = 20260823L,
      requested_dimension = 20L,
      max_iter = 100L,
      repetitions = 3L,
      warmups = 1L,
      callback_modes = "combined;separate",
      profiles = paste(packet5a_profiles(), collapse = ";"),
      policies = paste(names(packet5a_policies()), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

packet5a_main <- function() {
  config <- packet5a_parse_args(commandArgs(trailingOnly = TRUE))
  loader <- load_local_mize()
  message("Loaded mize via ", loader)
  cases <- benchmark_cases(packet5a_case_config(config))
  expected_case_count <- if (config$scope == "tranche") 4L else 36L
  if (length(cases) != expected_case_count) {
    stop("Packet 5A case resolution is incomplete", call. = FALSE)
  }
  manifest <- benchmark_case_manifest(cases)
  result <- packet5a_run_grid(cases)
  summary <- packet5a_summary(result$runs, result$progress, manifest)

  write_benchmark_artifact(result$runs, config$out, "Packet 5A runs")
  write_benchmark_artifact(
    result$progress,
    config$progress_out,
    "Packet 5A progress"
  )
  write_benchmark_artifact(manifest, config$case_out, "Packet 5A cases")
  write_benchmark_artifact(
    summary,
    config$summary_out,
    "Packet 5A resource summary"
  )
  write_benchmark_artifact(
    result$initializers,
    config$initializer_out,
    "Packet 5A initializer evidence"
  )
  artifact_manifest <- packet5a_manifest(
    config,
    c(
      runs = config$out,
      progress = config$progress_out,
      cases = config$case_out,
      summary = config$summary_out,
      initializers = config$initializer_out
    )
  )
  write_benchmark_artifact(
    artifact_manifest,
    config$manifest_out,
    "Packet 5A manifest"
  )
}

if (sys.nframe() == 0L) {
  packet5a_main()
}
